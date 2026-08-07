#!/usr/bin/env python
"""
omero_nuclei_flow.py
====================

High-content ImageXpress -> OMERO Screen analysis pipeline.

Pipeline
--------
1.  Connect to OMERO, walk Screen -> Plate -> Well -> WellSample (field) -> Image.
2.  Segment nuclei on the DAPI channel with StarDist (default) or Cellpose.
    StarDist's star-convex polygons + NMS are the strongest option for *dense,
    touching* nuclei, which is why it is the default here.
3.  Measure per-nucleus intensity (mean / median / integrated) in *every*
    channel, plus morphology, with optional per-image background subtraction.
4.  Gate each marker channel (everything that is not DAPI) into +/- using an
    automatic threshold (Otsu / GMM / quantile) or negative-control wells.
5.  Classify every cell as negative / single / double / triple / ... positive
    based on how many marker channels it is positive for. The number of
    classes adapts to the number of channels present in each screen.
6.  Export flow-cytometry-style plots: asinh-scaled histograms with gates,
    biaxial density scatters with quadrant statistics, population stacked bars,
    combination bar charts and plate heatmaps.
7.  Optionally push masks back to OMERO as ROIs and attach the CSV/plots.

Quick start
-----------
    python omero_nuclei_flow.py \
        --host omero.myinstitute.org --user jdoe \
        --screen-id 1234 \
        --outdir ./results

    # re-gate / re-plot without touching the server
    python omero_nuclei_flow.py --from-csv results/Plate1/per_nucleus_measurements.csv.gz \
        --outdir ./results_regated --threshold-method gmm

Author: generated scaffold - review thresholds and QC filters before trusting
        the numbers on a new assay.
"""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import math
import os
import re
import sys
from dataclasses import dataclass, asdict, field
from itertools import combinations
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from scipy import ndimage as ndi
from skimage.filters import threshold_otsu
from skimage.measure import regionprops_table
from skimage.segmentation import clear_border, find_boundaries

LOG = logging.getLogger("omero_nuclei_flow")

DAPI_PATTERNS = ("dapi", "hoechst", "405", "dna", "nuc", "w1")
POSITIVITY_NAMES = {
    0: "negative",
    1: "single",
    2: "double",
    3: "triple",
    4: "quadruple",
    5: "quintuple",
}

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #


@dataclass
class Config:
    # connection
    host: str = "localhost"
    port: int = 4064
    user: Optional[str] = None
    password: Optional[str] = None
    secure: bool = True
    group: str = "-1"

    # data selection
    screen_id: Optional[int] = None
    plate_ids: List[int] = field(default_factory=list)
    wells: Optional[List[str]] = None
    max_fields: Optional[int] = None

    # imaging
    dapi_channel: Optional[str] = None       # index or name fragment
    z_mode: str = "mid"                      # mid | mip
    timepoint: int = 0
    server_tile_height: int = 0              # 0 = whole plane; else strip rows

    # segmentation
    model: str = "stardist"                  # stardist | cellpose
    stardist_model: str = "2D_versatile_fluo"
    cellpose_model: str = "nuclei"
    prob_thresh: Optional[float] = None
    nms_thresh: Optional[float] = None
    diameter: Optional[float] = None
    n_tiles: Optional[int] = None
    gpu: bool = True
    gpu_memory_limit: Optional[int] = None   # MB cap for the GPU
    gc_every: int = 50                       # garbage-collect every N images
    reload_model_every: int = 0              # rebuild the model every N images

    # QC filters
    min_area: int = 30
    max_area: int = 100000
    exclude_border: bool = True
    min_solidity: float = 0.0

    # dead / dying cell removal (bright, condensed DAPI)
    filter_dead: bool = False
    dead_method: str = "mad"                 # mad | quantile | otsu | manual
    dead_stat: str = "mean"                  # mean | median | integrated
    dead_k: float = 3.0                      # robust SDs above the median
    dead_quantile: float = 0.99
    dead_threshold: Optional[float] = None
    dead_scope: str = "plate"                # screen | plate | well | image
    dead_require_small: bool = False         # bright AND condensed
    dead_area_quantile: float = 0.25
    dead_action: str = "remove"              # remove | flag

    # measurement
    background: str = "median"               # none | median
    intensity_stat: str = "mean"             # mean | median | integrated

    # foci detection (e.g. gamma-H2AX, 53BP1)
    foci_channels: List[str] = field(default_factory=list)
    pixel_size_um: Optional[float] = None    # override; else read from OMERO
    foci_diameter_um: float = 0.8            # expected physical focus diameter
    foci_auto_scale: bool = True             # derive px params from pixel size
    foci_upsample: int = 0                   # 0 = auto (2-3x for small foci)
    foci_bg_clip: bool = True                # sigma-clip the in-nucleus background
    save_foci_crops: int = 0                 # N zoomed QC montages per plate
    foci_method: str = "log"                 # log | tophat | hmax
    foci_tophat_radius: float = 5.0
    foci_min_sigma: float = 1.0
    foci_max_sigma: float = 4.0
    foci_num_sigma: int = 5
    foci_log_threshold: float = 0.02
    foci_k: float = 3.0                      # robust SDs over in-nucleus background
    foci_abs_threshold: Optional[float] = None
    foci_h: float = 0.0                      # h-maxima depth (0 = auto)
    foci_min_area: int = 3
    foci_max_area: int = 400
    foci_min_count: int = 3                  # >= N foci => foci-positive
    foci_as_marker: bool = False             # use foci count, not intensity, for +/-

    # OMERO key-value pairs (map annotations)
    fetch_kv: bool = True
    kv_levels: List[str] = field(default_factory=lambda: ["screen", "plate", "well"])
    kv_prefix: str = "kv_"
    group_by: Optional[str] = None

    # tables
    per_nucleus_format: str = "csv.gz"       # csv | csv.gz | parquet | xlsx | none
    save_foci_table: bool = False
    export_explorer: bool = False            # write the interactive HTML explorer
    explorer_no_embed: bool = False          # don't inline the data

    # gating
    threshold_method: str = "otsu"           # otsu | gmm | quantile | control | manual
    threshold_scope: str = "plate"           # screen | plate | well | image
    quantile: float = 0.99
    control_wells: List[str] = field(default_factory=list)
    control_sd: float = 3.0
    manual_thresholds: Dict[str, float] = field(default_factory=dict)
    marker_channels: Optional[List[str]] = None
    subset_analyses: bool = False            # per-panel re-gating
    subset_sizes: List[int] = field(default_factory=lambda: [1, 2])
    subset_panels: List[str] = field(default_factory=list)

    # plotting
    cofactor: Optional[float] = None
    max_points: int = 50000
    per_well_plots: bool = False
    save_overlays: int = 0                   # number of QC overlays per plate
    overlay_line_width: int = 2              # outline thickness in px

    # output
    outdir: Path = Path("./results")
    upload_rois: bool = False
    attach_results: bool = False
    from_csv: Optional[Path] = None


# --------------------------------------------------------------------------- #
# OMERO access
# --------------------------------------------------------------------------- #


def connect(cfg: Config):
    """Open a BlitzGateway connection (omero-py imported lazily)."""
    from omero.gateway import BlitzGateway

    password = cfg.password or os.environ.get("OMERO_PASSWORD")
    if not password:
        import getpass
        password = getpass.getpass(f"OMERO password for {cfg.user}@{cfg.host}: ")

    conn = BlitzGateway(cfg.user, password, host=cfg.host, port=cfg.port,
                        secure=cfg.secure)
    if not conn.connect():
        raise RuntimeError(f"Could not connect to OMERO at {cfg.host}:{cfg.port}")
    conn.c.enableKeepAlive(60)
    conn.SERVICE_OPTS.setOmeroGroup(cfg.group)
    LOG.info("Connected to OMERO as %s", conn.getUser().getName())
    return conn


def map_annotations(obj) -> Dict[str, str]:
    """All key-value pairs (MapAnnotations) attached to one OMERO object.

    OMERO allows repeated keys within and across annotations, so duplicates are
    joined with '; ' rather than silently overwritten.
    """
    from omero.gateway import MapAnnotationWrapper

    out: Dict[str, str] = {}
    try:
        anns = list(obj.listAnnotations())
    except Exception as exc:                          # noqa: BLE001
        LOG.debug("Could not list annotations on %s: %s", obj, exc)
        return out

    for ann in anns:
        if not isinstance(ann, MapAnnotationWrapper):
            continue
        for pair in ann.getValue() or []:
            try:
                k, v = pair[0], pair[1]
            except (TypeError, IndexError):
                continue
            k, v = str(k).strip(), str(v).strip()
            if not k:
                continue
            if k in out:
                if v and v not in out[k].split("; "):
                    out[k] = f"{out[k]}; {v}"
            else:
                out[k] = v
    return out


def kv_column_name(key: str, prefix: str) -> Optional[str]:
    """Turn an arbitrary OMERO key into a safe, predictable column name."""
    s = re.sub(r"[^0-9A-Za-z]+", "_", str(key)).strip("_")
    return f"{prefix}{s}" if s else None


def kv_to_columns(kv: Dict[str, str], prefix: str) -> Dict[str, str]:
    cols: Dict[str, str] = {}
    for k, v in kv.items():
        col = kv_column_name(k, prefix)
        if col is None:
            continue
        if col in cols and cols[col] != v:
            LOG.debug("Key-value collision on column %s - keeping first", col)
            continue
        cols[col] = v
    return cols


def coerce_numeric_kv(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    """Convert key-value columns that are entirely numeric into real numbers,
    so downstream code can sort and plot concentrations without parsing."""
    for c in [c for c in df.columns if c.startswith(prefix)]:
        conv = pd.to_numeric(df[c], errors="coerce")
        if df[c].notna().any() and conv.notna().sum() == df[c].notna().sum():
            df[c] = conv
    return df


def iter_plates(conn, cfg: Config):
    """Yield PlateWrapper objects from a screen id or explicit plate ids."""
    if cfg.screen_id is not None:
        screen = conn.getObject("Screen", cfg.screen_id)
        if screen is None:
            raise RuntimeError(f"Screen {cfg.screen_id} not found (check group/permissions)")
        LOG.info("Screen %s: %s", cfg.screen_id, screen.getName())
        for plate in screen.listChildren():
            yield plate
    for pid in cfg.plate_ids:
        plate = conn.getObject("Plate", pid)
        if plate is None:
            LOG.warning("Plate %s not found - skipping", pid)
            continue
        yield plate


def well_label(well) -> str:
    """A1-style label, robust across omero-py versions."""
    try:
        return well.getWellPos()
    except Exception:
        return f"{chr(ord('A') + well.row)}{well.column + 1}"


def iter_images(plate, cfg: Config):
    """Yield (well_label, row, col, field_index, ImageWrapper, WellWrapper)."""
    wanted = {w.strip().upper() for w in cfg.wells} if cfg.wells else None
    for well in plate.listChildren():
        label = well_label(well)
        if wanted and label.upper() not in wanted:
            continue
        n_fields = well.countWellSample()
        if cfg.max_fields:
            n_fields = min(n_fields, cfg.max_fields)
        for idx in range(n_fields):
            ws = well.getWellSample(idx)
            if ws is None:
                continue
            img = ws.getImage()
            if img is None:
                continue
            yield label, well.row, well.column, idx, img, well


def channel_names(image) -> List[str]:
    """Channel labels WITHOUT initialising the server-side rendering engine.

    Plain image.getChannels() prepares the rendering engine, which makes the
    server load pixel data just to compute display defaults. On Bio-Formats
    backed plates (ImageXpress/CellWorx read in place) that is a whole-plane
    Java heap allocation per call, for information we already have in the
    logical-channel metadata. noRE=True skips all of it.
    """
    try:
        channels = image.getChannels(noRE=True)
    except TypeError:                     # very old omero-py without the kwarg
        channels = image.getChannels()
    names = []
    for i, ch in enumerate(channels or []):
        name = None
        try:
            name = ch.getLabel() or ch.getName()
        except Exception:                 # noqa: BLE001
            pass
        names.append(str(name).strip() if name else f"C{i + 1}")
    return names


def _is_server_oom(exc: BaseException) -> bool:
    """Server-side Java heap exhaustion, reported as an OMERO InternalException."""
    text = f"{type(exc).__name__}: {exc}"
    return ("OutOfMemoryError" in text or "Java heap space" in text
            or "InternalException" in text)


def resolve_dapi(names: Sequence[str], requested: Optional[str]) -> int:
    """Return the index of the nuclear channel."""
    if requested is not None:
        if re.fullmatch(r"\d+", str(requested)):
            return int(requested)
        for i, n in enumerate(names):
            if str(requested).lower() in n.lower():
                return i
        raise ValueError(f"Channel matching '{requested}' not found in {names}")
    for i, n in enumerate(names):
        if any(p in n.lower() for p in DAPI_PATTERNS):
            return i
    LOG.warning("No DAPI-like channel name in %s - falling back to channel 0", names)
    return 0


def read_stack(image, cfg: Config, tile_height: int = 0) -> np.ndarray:
    """Return a (H, W, C) float32 array for the chosen z/t.

    tile_height > 0 reads each plane in horizontal strips via getTile() instead
    of one getPlane() call. The server allocates a buffer per request, so on
    Bio-Formats backed data (where every read re-opens the original file) strips
    are the difference between a 16 MB allocation and a 200 KB one.
    """
    pix = image.getPrimaryPixels()
    size_c, size_z = image.getSizeC(), image.getSizeZ()
    size_x, size_y = image.getSizeX(), image.getSizeY()
    t = min(cfg.timepoint, image.getSizeT() - 1)

    if cfg.z_mode == "mip" and size_z > 1:
        zs = list(range(size_z))
    else:
        zs = [size_z // 2]

    def _plane(z: int, c: int) -> np.ndarray:
        if not tile_height:
            return np.asarray(pix.getPlane(z, c, t), dtype=np.float32)
        strips = []
        for y0 in range(0, size_y, tile_height):
            h = min(tile_height, size_y - y0)
            strips.append(np.asarray(pix.getTile(z, c, t, tile=(0, y0, size_x, h)),
                                     dtype=np.float32))
        return np.vstack(strips)

    if not tile_height:
        zct = [(z, c, t) for c in range(size_c) for z in zs]
        planes = [np.asarray(p, dtype=np.float32) for p in pix.getPlanes(zct)]
    else:
        planes = [_plane(z, c) for c in range(size_c) for z in zs]

    out = np.empty((planes[0].shape[0], planes[0].shape[1], size_c), dtype=np.float32)
    k = 0
    for c in range(size_c):
        block = np.stack(planes[k:k + len(zs)], axis=0)
        k += len(zs)
        out[..., c] = block.max(axis=0) if len(zs) > 1 else block[0]
    return out


# --------------------------------------------------------------------------- #
# Segmentation
# --------------------------------------------------------------------------- #


def normalize_percentile(img: np.ndarray, lo: float = 1.0,
                         hi: float = 99.8) -> np.ndarray:
    p_lo, p_hi = np.percentile(img, [lo, hi])
    if p_hi <= p_lo:
        p_hi = p_lo + 1.0
    return np.clip((img - p_lo) / (p_hi - p_lo), 0, 1).astype(np.float32)


def _is_oom(exc: BaseException) -> bool:
    """TensorFlow and PyTorch report GPU exhaustion under several names."""
    text = f"{type(exc).__name__}: {exc}".lower()
    return any(s in text for s in (
        "dst tensor is not initialized",     # TF host->device copy failed
        "oom when allocating",
        "resourceexhausted",
        "out of memory",
        "cuda_error_out_of_memory",
        "cudnn_status_alloc_failed",
    ))


class Segmenter:
    """Thin wrapper so StarDist and Cellpose share one call signature."""

    def __init__(self, cfg: Config):
        self.cfg = cfg
        self.kind = cfg.model.lower()
        self._model = None
        self._calls = 0
        self._oom_events = 0

    # ---- GPU memory management -------------------------------------------

    def _configure_gpu(self):
        """Enable incremental GPU allocation before the first op runs.

        By default TensorFlow reserves essentially all VRAM up front, then
        fragments it over thousands of predict() calls until a copy fails with
        'Dst tensor is not initialized'. Memory growth makes it allocate as
        needed, which both survives fragmentation better and leaves room for
        anything else using the card.
        """
        if not self.cfg.gpu or self.kind != "stardist":
            return
        try:
            import tensorflow as tf
            gpus = tf.config.list_physical_devices("GPU")
            for g in gpus:
                try:
                    tf.config.experimental.set_memory_growth(g, True)
                except RuntimeError:
                    pass          # already initialised; harmless
            if gpus and self.cfg.gpu_memory_limit:
                tf.config.set_logical_device_configuration(
                    gpus[0], [tf.config.LogicalDeviceConfiguration(
                        memory_limit=int(self.cfg.gpu_memory_limit))])
                LOG.info("GPU memory capped at %d MB", self.cfg.gpu_memory_limit)
            if gpus:
                LOG.info("GPU memory growth enabled on %d device(s)", len(gpus))
        except Exception as exc:                       # noqa: BLE001
            LOG.debug("Could not configure GPU memory: %s", exc)

    def _free_memory(self, reload_model: bool = False):
        """Release what we can after an OOM, optionally rebuilding the model."""
        import gc
        if reload_model and self.kind == "stardist":
            try:
                from tensorflow.keras import backend as K
                self._model = None
                K.clear_session()
            except Exception:                          # noqa: BLE001
                self._model = None
        gc.collect()
        if self.kind == "cellpose":
            try:
                import torch
                torch.cuda.empty_cache()
            except Exception:                          # noqa: BLE001
                pass

    def _load(self):
        if self._model is not None:
            return
        if self.kind == "stardist":
            # TensorFlow >=2.16 ships Keras 3, which StarDist does not support.
            # If the user installed tf-keras (the legacy Keras 2 package), point
            # tf.keras at it. Must happen before TensorFlow is first imported;
            # the import below is the first one, since it is lazy. Harmless on
            # TF <2.16, where the variable is ignored.
            import importlib.util
            if importlib.util.find_spec("tf_keras") is not None:
                if os.environ.setdefault("TF_USE_LEGACY_KERAS", "1") == "1":
                    LOG.info("tf-keras found -> TF_USE_LEGACY_KERAS=1 "
                             "(StarDist needs Keras 2)")
            try:
                from stardist.models import StarDist2D
            except ModuleNotFoundError as exc:
                if "pkg_resources" in str(exc):
                    raise RuntimeError(
                        "StarDist needs pkg_resources, which setuptools 82 "
                        "removed (Feb 2026). Fix it with:\n"
                        "    pip install \"setuptools<81\"\n"
                        "The pinned environment.yml already does this; run "
                        "--check-env to verify the rest of the stack."
                    ) from exc
                raise RuntimeError(
                    f"Could not import StarDist ({exc}). Install it with "
                    "'pip install stardist tensorflow==2.10.1', or use "
                    "--model cellpose instead."
                ) from exc
            self._configure_gpu()
            LOG.info("Loading StarDist model '%s'", self.cfg.stardist_model)
            self._model = StarDist2D.from_pretrained(self.cfg.stardist_model)
        elif self.kind == "cellpose":
            try:
                from cellpose import models as cp_models
            except ModuleNotFoundError as exc:
                raise RuntimeError(
                    f"Could not import Cellpose ({exc}). Install it with "
                    "'pip install cellpose', or use --model stardist."
                ) from exc
            LOG.info("Loading Cellpose model '%s'", self.cfg.cellpose_model)
            try:  # Cellpose >= 4 (cpsam generalist)
                self._model = cp_models.CellposeModel(gpu=self.cfg.gpu)
                self._cp_api = 4
            except TypeError:
                self._model = cp_models.Cellpose(gpu=self.cfg.gpu,
                                                 model_type=self.cfg.cellpose_model)
                self._cp_api = 3
        else:
            raise ValueError(f"Unknown model '{self.cfg.model}'")

    def _segment(self, dapi: np.ndarray, tiles: Optional[int] = None,
                 force_cpu: bool = False) -> np.ndarray:
        if self.kind == "stardist":
            kw = {}
            if self.cfg.prob_thresh is not None:
                kw["prob_thresh"] = self.cfg.prob_thresh
            if self.cfg.nms_thresh is not None:
                kw["nms_thresh"] = self.cfg.nms_thresh
            n = tiles if tiles is not None else self.cfg.n_tiles
            if n:
                kw["n_tiles"] = (int(n), int(n))
            img = normalize_percentile(dapi)
            if force_cpu:
                import tensorflow as tf
                with tf.device("/CPU:0"):
                    labels, _ = self._model.predict_instances(img, **kw)
            else:
                labels, _ = self._model.predict_instances(img, **kw)
            return labels.astype(np.int32)

        # cellpose
        diam = self.cfg.diameter or 0
        res = self._model.eval(normalize_percentile(dapi), diameter=diam or None,
                               channels=[0, 0])
        return np.asarray(res[0], dtype=np.int32)

    def __call__(self, dapi: np.ndarray) -> np.ndarray:
        self._load()
        self._calls += 1

        # Periodic housekeeping: TensorFlow/Keras accumulate memory across many
        # predict() calls, so a run that starts fine can fail thousands of
        # images later. Cheap insurance against a mid-screen crash.
        if self.cfg.gc_every and self._calls % self.cfg.gc_every == 0:
            self._free_memory()
        if (self.cfg.reload_model_every
                and self._calls % self.cfg.reload_model_every == 0):
            LOG.info("Reloading segmentation model after %d images "
                     "(--reload-model-every)", self._calls)
            self._free_memory(reload_model=True)
            self._load()

        try:
            return self._segment(dapi)
        except Exception as exc:                       # noqa: BLE001
            if not _is_oom(exc):
                raise
            self._oom_events += 1
            base = self.cfg.n_tiles or 1
            more = max(base * 2, 2)
            LOG.warning("GPU out of memory on this field (event %d). Freeing "
                        "memory and retrying with %dx%d tiling.",
                        self._oom_events, more, more)
            self._free_memory()          # cheap first: collect, keep the model
            try:
                return self._segment(dapi, tiles=more)
            except Exception as exc2:                  # noqa: BLE001
                if not _is_oom(exc2):
                    raise
                LOG.warning("Still out of GPU memory - segmenting this field on "
                            "the CPU (slow, but the run continues). Consider "
                            "--reload-model-every 200 or --gpu-memory-limit.")
                self._free_memory(reload_model=True)
                self._load()
                return self._segment(dapi, tiles=more, force_cpu=True)


def filter_labels(labels: np.ndarray, cfg: Config) -> np.ndarray:
    """Drop border objects and objects outside the accepted size/shape range."""
    if cfg.exclude_border:
        labels = clear_border(labels)
    ids, counts = np.unique(labels, return_counts=True)
    keep = ids[(ids != 0) & (counts >= cfg.min_area) & (counts <= cfg.max_area)]
    if len(keep) == len(ids) - 1:
        return labels
    mask = np.isin(labels, keep)
    return np.where(mask, labels, 0).astype(np.int32)


# --------------------------------------------------------------------------- #
# Measurement
# --------------------------------------------------------------------------- #


def measure(labels: np.ndarray, stack: np.ndarray, names: Sequence[str],
            cfg: Config) -> pd.DataFrame:
    """Per-nucleus intensities in every channel + morphology."""
    ids = np.unique(labels)
    ids = ids[ids != 0]
    if ids.size == 0:
        return pd.DataFrame()

    morph = pd.DataFrame(regionprops_table(
        labels,
        properties=("label", "area", "centroid", "eccentricity",
                    "solidity", "perimeter"),
    ))
    morph = morph.rename(columns={"centroid-0": "y", "centroid-1": "x"})
    morph = morph.set_index("label").loc[ids].reset_index()

    background = {}
    bg_mask = labels == 0
    for c, name in enumerate(names):
        plane = stack[..., c]
        bg = float(np.median(plane[bg_mask])) if (cfg.background == "median"
                                                  and bg_mask.any()) else 0.0
        background[name] = bg
        corrected = plane - bg
        morph[f"mean_{name}"] = ndi.mean(corrected, labels, ids)
        morph[f"median_{name}"] = ndi.median(corrected, labels, ids)
        morph[f"integrated_{name}"] = ndi.sum_labels(corrected, labels, ids)
        morph[f"bg_{name}"] = bg

    if cfg.min_solidity > 0:
        morph = morph[morph["solidity"] >= cfg.min_solidity]
    return morph


# --------------------------------------------------------------------------- #
# Foci detection (gamma-H2AX / 53BP1 style punctate signals)
# --------------------------------------------------------------------------- #


def resolve_channels(names: Sequence[str], requested: Sequence[str]) -> List[int]:
    """Map channel indices or name fragments onto channel positions."""
    out = []
    for r in requested:
        r = str(r).strip()
        if not r:
            continue
        if re.fullmatch(r"\d+", r):
            out.append(int(r))
            continue
        hit = [i for i, n in enumerate(names) if r.lower() in n.lower()]
        if not hit:
            raise ValueError(f"Foci channel '{r}' not found in {list(names)}")
        out.append(hit[0])
    return out


def _label_lut(labels: np.ndarray, ids: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Broadcast one value per label back onto the pixel grid."""
    lut = np.zeros(int(labels.max()) + 1, dtype=np.float64)
    lut[ids] = values
    return lut[labels]


@dataclass
class FociParams:
    """Detection parameters in pixels, derived from physical size where possible."""
    min_sigma: float
    max_sigma: float
    num_sigma: int
    tophat_radius: int
    min_area: int
    max_area: int
    upsample: int
    seed_merge_radius: int
    log_threshold: float
    k: float
    method: str
    abs_threshold: Optional[float]
    h: float
    bg_clip: bool
    source: str = "manual"

    def describe(self) -> str:
        return (f"sigma {self.min_sigma:.2f}-{self.max_sigma:.2f} px, "
                f"tophat r={self.tophat_radius}, area {self.min_area}-{self.max_area} px, "
                f"upsample x{self.upsample} [{self.source}]")


def derive_foci_params(cfg: Config, pixel_um: Optional[float]) -> FociParams:
    """Scale the detector to the actual focus size in pixels.

    At 20x an ~0.8 um gamma-H2AX focus spans barely 1-3 pixels, where defaults
    tuned for 60x images miss most of them: the Laplacian-of-Gaussian scales are
    far too coarse, the top-hat structuring element is too large to flatten the
    background around such small objects, and the minimum-area filter discards
    genuine single-pixel foci. Deriving all three from the pixel size fixes the
    whole family of problems at once.
    """
    manual = FociParams(
        min_sigma=cfg.foci_min_sigma, max_sigma=cfg.foci_max_sigma,
        num_sigma=cfg.foci_num_sigma,
        tophat_radius=int(max(1, round(cfg.foci_tophat_radius))),
        min_area=cfg.foci_min_area, max_area=cfg.foci_max_area,
        upsample=max(1, cfg.foci_upsample),
        seed_merge_radius=max(1, int(round(cfg.foci_tophat_radius / 4))),
        log_threshold=cfg.foci_log_threshold,
        k=cfg.foci_k, method=cfg.foci_method,
        abs_threshold=cfg.foci_abs_threshold, h=cfg.foci_h,
        bg_clip=cfg.foci_bg_clip, source="manual")

    px = cfg.pixel_size_um or pixel_um
    if not cfg.foci_auto_scale or not px or px <= 0:
        return manual

    r_px = 0.5 * cfg.foci_diameter_um / px          # focus radius in pixels
    area = math.pi * r_px ** 2
    sigma = r_px / math.sqrt(2.0)                   # LoG peak response scale

    upsample = cfg.foci_upsample
    if upsample <= 0:
        # Work at a scale where LoG is numerically stable (sigma >~ 2 px).
        # Below that, discretisation of the second derivative dominates and
        # small foci are simply not detected reliably at any threshold.
        upsample = int(np.clip(math.ceil(2.0 / max(sigma, 1e-6)), 1, 3))

    s = sigma * upsample                            # scale in the working grid
    min_sigma = max(0.5, 0.6 * s)
    return FociParams(
        min_sigma=min_sigma,
        max_sigma=max(min_sigma + 0.4, 1.8 * s),
        num_sigma=cfg.foci_num_sigma,
        tophat_radius=int(max(3, round(4 * r_px))),
        min_area=int(max(1, round(0.4 * area))),
        max_area=int(max(9, round(8 * area))),
        upsample=upsample,
        seed_merge_radius=int(round(0.75 * r_px)),
        log_threshold=cfg.foci_log_threshold, k=cfg.foci_k,
        method=cfg.foci_method, abs_threshold=cfg.foci_abs_threshold,
        h=cfg.foci_h, bg_clip=cfg.foci_bg_clip,
        source=f"auto from {px:.3f} um/px, focus {cfg.foci_diameter_um} um")


def _clipped_bg(th: np.ndarray, labels: np.ndarray, ids: np.ndarray,
                k: float, iters: int = 2) -> np.ndarray:
    """Per-label median + k*robust-SD, sigma-clipped so that the foci themselves
    don't inflate the background they are being compared against.

    This matters most in heavily damaged nuclei: at high foci burden the naive
    median sits inside the foci distribution and the threshold runs away
    upward, so the most interesting cells are the ones that get undercounted."""
    med = np.asarray(ndi.median(th, labels, ids), dtype=float)
    mad = np.asarray(ndi.median(np.abs(th - _label_lut(labels, ids, med)),
                                labels, ids), dtype=float)
    for _ in range(max(0, iters)):
        hi = _label_lut(labels, ids, med + 3.0 * 1.4826 * np.maximum(mad, 1e-9))
        kept = np.where(th <= hi, labels, 0)
        new_med = np.asarray(ndi.median(th, kept, ids), dtype=float)
        new_mad = np.asarray(ndi.median(np.abs(th - _label_lut(labels, ids, new_med)),
                                        kept, ids), dtype=float)
        ok = np.isfinite(new_med) & np.isfinite(new_mad) & (new_mad > 0)
        med = np.where(ok, new_med, med)
        mad = np.where(ok, new_mad, mad)
    return med + k * 1.4826 * mad


def remove_small(mask: np.ndarray, min_size: int) -> np.ndarray:
    """Drop connected components below min_size (version-proof)."""
    lbl, n = ndi.label(mask)
    if n == 0:
        return mask
    sizes = np.bincount(lbl.ravel())
    return mask & (sizes[lbl] >= min_size)


def detect_foci(plane: np.ndarray, nuc_labels: np.ndarray,
                cfg: Config, params: Optional[FociParams] = None
                ) -> Tuple[np.ndarray, pd.DataFrame]:
    """Find punctate foci *inside* nuclei.

    Returns (foci_label_image, per-focus DataFrame with a 'nucleus_label' column).

    Stages: white top-hat to flatten the diffuse nuclear signal; a threshold
    computed per nucleus (foci brightness varies enormously between cells);
    optional 2-3x upsampling so that under-sampled foci land in the range where
    Laplacian-of-Gaussian scale selection is numerically stable; then seeded
    watershed to split touching foci.
    """
    from skimage.morphology import white_tophat, disk, h_maxima
    from skimage.segmentation import watershed

    if params is None:
        params = derive_foci_params(cfg, None)

    empty = pd.DataFrame(columns=["nucleus_label", "focus_id", "area", "y", "x",
                                  "mean_intensity", "max_intensity",
                                  "integrated_intensity"])
    nuc_mask = nuc_labels > 0
    ids = np.unique(nuc_labels)
    ids = ids[ids != 0]
    if ids.size == 0:
        return np.zeros_like(nuc_labels), empty

    th = white_tophat(plane, disk(params.tophat_radius))

    # ---- adaptive threshold, per nucleus ----------------------------------
    if params.abs_threshold is not None:
        thr_map = np.full(plane.shape, float(params.abs_threshold))
    elif params.bg_clip:
        thr_map = _label_lut(nuc_labels, ids,
                             _clipped_bg(th, nuc_labels, ids, params.k))
    else:
        med = np.asarray(ndi.median(th, nuc_labels, ids), dtype=float)
        dev = np.abs(th - _label_lut(nuc_labels, ids, med))
        mad = np.asarray(ndi.median(dev, nuc_labels, ids), dtype=float)
        thr_map = _label_lut(nuc_labels, ids, med + params.k * 1.4826 * mad)

    mask = (th > thr_map) & nuc_mask
    if not mask.any():
        return np.zeros_like(nuc_labels), empty
    mask = remove_small(mask, max(1, params.min_area))

    # ---- seeds -------------------------------------------------------------
    method = params.method.lower()
    if method == "tophat":
        foci_labels, _ = ndi.label(mask)
    else:
        if method == "log":
            from skimage.feature import blob_log
            # Robust normalisation: a percentile-based scale is dominated by the
            # foci themselves in heavily damaged nuclei, which shrinks the LoG
            # response below a fixed threshold and silently loses counts.
            vals = th[nuc_mask]
            med_g = float(np.median(vals))
            mad_g = float(np.median(np.abs(vals - med_g)))
            scale = max(med_g + 5.0 * 1.4826 * mad_g, 1e-6)

            src = th
            if params.upsample > 1:
                src = ndi.zoom(th, params.upsample, order=1)
            norm = np.clip(src / scale, 0, 1)

            blobs = blob_log(norm, min_sigma=params.min_sigma,
                             max_sigma=params.max_sigma,
                             num_sigma=int(params.num_sigma),
                             threshold=params.log_threshold, overlap=0.5)
            if len(blobs):
                coords = blobs[:, :2] / params.upsample
                coords = np.clip(np.round(coords).astype(int), 0,
                                 np.array(th.shape) - 1)
                inside = nuc_mask[coords[:, 0], coords[:, 1]]
                coords = coords[inside]
            else:
                coords = np.empty((0, 2), int)

            # Merge seeds that are closer together than one focus radius:
            # multi-scale LoG routinely fires twice on the same small focus,
            # and each extra seed would otherwise split it into two "foci".
            seeds = np.zeros(plane.shape, dtype=bool)
            if len(coords):
                seeds[coords[:, 0], coords[:, 1]] = True
                if params.seed_merge_radius >= 1:
                    seeds = ndi.binary_dilation(
                        seeds, structure=disk(params.seed_merge_radius))
            markers, _ = ndi.label(seeds)
        else:  # hmax
            h = params.h or float(np.percentile(th[nuc_mask], 99) -
                                  np.percentile(th[nuc_mask], 50)) / 2.0
            markers, _ = ndi.label(h_maxima(th, max(h, 1e-6)) & nuc_mask)

        markers = np.where(mask, markers, 0)      # seeds must be real signal
        if markers.max() == 0:
            foci_labels, _ = ndi.label(mask)
        else:
            foci_labels = watershed(-th, markers=markers, mask=mask)

    # ---- size filter -------------------------------------------------------
    fids, counts = np.unique(foci_labels, return_counts=True)
    keep = fids[(fids != 0) & (counts >= params.min_area) &
                (counts <= params.max_area)]
    if keep.size == 0:
        return np.zeros_like(nuc_labels), empty
    foci_labels = np.where(np.isin(foci_labels, keep), foci_labels, 0).astype(np.int32)

    # ---- assign each focus to the nucleus holding most of its pixels -------
    sel = foci_labels > 0
    fl, nl = foci_labels[sel].astype(np.int64), nuc_labels[sel].astype(np.int64)
    span = int(nuc_labels.max()) + 1
    pair, cnt = np.unique(fl * span + nl, return_counts=True)
    order = np.lexsort((-cnt, pair // span))
    pair, cnt = pair[order], cnt[order]
    first = np.unique(pair // span, return_index=True)[1]
    focus_ids = (pair // span)[first]
    nucleus_of = (pair % span)[first]

    # ---- per-focus measurements -------------------------------------------
    props = pd.DataFrame(regionprops_table(
        foci_labels, intensity_image=plane,
        properties=("label", "area", "centroid", "intensity_mean", "intensity_max"),
    )).rename(columns={"label": "focus_id", "centroid-0": "y", "centroid-1": "x",
                       "intensity_mean": "mean_intensity",
                       "intensity_max": "max_intensity"})
    props["integrated_intensity"] = props["mean_intensity"] * props["area"]
    props = props.merge(pd.DataFrame({"focus_id": focus_ids,
                                      "nucleus_label": nucleus_of}),
                        on="focus_id", how="inner")
    props = props[props["nucleus_label"] != 0]
    return foci_labels, props


def summarise_foci(props: pd.DataFrame, nuclei: pd.DataFrame,
                   chan: str) -> pd.DataFrame:
    """Collapse the per-focus table into per-nucleus columns."""
    cols = {
        f"foci_n_{chan}": ("focus_id", "count"),
        f"foci_area_total_{chan}": ("area", "sum"),
        f"foci_area_mean_{chan}": ("area", "mean"),
        f"foci_mean_intensity_{chan}": ("mean_intensity", "mean"),
        f"foci_max_intensity_{chan}": ("max_intensity", "max"),
        f"foci_integrated_intensity_{chan}": ("integrated_intensity", "sum"),
    }
    if props.empty:
        agg = pd.DataFrame(columns=["label"] + list(cols))
    else:
        agg = (props.groupby("nucleus_label")
                    .agg(**{k: pd.NamedAgg(column=c, aggfunc=f)
                            for k, (c, f) in cols.items()})
                    .reset_index().rename(columns={"nucleus_label": "label"}))

    out = nuclei.merge(agg, on="label", how="left")
    out[f"foci_n_{chan}"] = out[f"foci_n_{chan}"].fillna(0).astype(int)
    for c in list(cols)[1:]:
        out[c] = out[c].fillna(0.0)
    if "area" in out.columns:
        out[f"foci_area_fraction_{chan}"] = (out[f"foci_area_total_{chan}"] /
                                             out["area"].replace(0, np.nan))
    return out


def foci_channel_columns(df: pd.DataFrame) -> List[str]:
    return [c[len("foci_n_"):] for c in df.columns if c.startswith("foci_n_")]


# --------------------------------------------------------------------------- #
# Gating
# --------------------------------------------------------------------------- #


def _gmm_threshold(values: np.ndarray) -> float:
    """Split a log-intensity distribution with a 2-component GMM."""
    try:
        from sklearn.mixture import GaussianMixture
    except ImportError:
        LOG.warning("scikit-learn not installed - falling back to Otsu")
        return float(np.expm1(threshold_otsu(np.log1p(np.clip(values, 0, None)))))

    v = np.log1p(np.clip(values, 0, None)).reshape(-1, 1)
    gm = GaussianMixture(n_components=2, random_state=0, n_init=3).fit(v)
    lo, hi = np.argsort(gm.means_.ravel())
    grid = np.linspace(v.min(), v.max(), 2000).reshape(-1, 1)
    resp = gm.predict_proba(grid)[:, hi]
    crossing = np.where(resp >= 0.5)[0]
    thr = grid[crossing[0], 0] if crossing.size else np.median(v)
    return float(np.expm1(thr))


def compute_threshold(values: np.ndarray, cfg: Config,
                      control: Optional[np.ndarray] = None,
                      marker: str = "") -> float:
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size < 20:
        return float(np.inf)

    method = cfg.threshold_method
    if method == "manual":
        if marker not in cfg.manual_thresholds:
            raise ValueError(f"No manual threshold given for '{marker}'")
        return float(cfg.manual_thresholds[marker])
    if method == "control":
        if control is None or control.size < 20:
            LOG.warning("No control cells for %s - falling back to Otsu", marker)
        else:
            c = np.log1p(np.clip(control, 0, None))
            return float(np.expm1(c.mean() + cfg.control_sd * c.std()))
    if method == "quantile":
        return float(np.quantile(v, cfg.quantile))
    if method == "gmm":
        return _gmm_threshold(v)
    return float(np.expm1(threshold_otsu(np.log1p(np.clip(v, 0, None)))))


def dapi_column(df: pd.DataFrame, stat: str) -> Tuple[str, str]:
    """Return (column_name, channel_name) for the nuclear channel."""
    names = pd.unique(df["dapi_channel"].dropna())
    if len(names) == 0:
        raise ValueError("No 'dapi_channel' column - cannot identify the nuclear channel")
    if len(names) > 1:
        LOG.warning("Multiple DAPI channel names present %s - using '%s'", names, names[0])
    return f"{stat}_{names[0]}", str(names[0])


def compute_dead_threshold(values: np.ndarray, cfg: Config) -> float:
    """Upper cutoff on DAPI intensity above which nuclei are called dead/dying."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size < 20:
        return float(np.inf)

    if cfg.dead_method == "manual":
        if cfg.dead_threshold is None:
            raise ValueError("--dead-method manual requires --dead-threshold")
        return float(cfg.dead_threshold)
    if cfg.dead_method == "quantile":
        return float(np.quantile(v, cfg.dead_quantile))
    if cfg.dead_method == "otsu":
        return float(np.expm1(threshold_otsu(np.log1p(np.clip(v, 0, None)))))

    # default: robust median + k * MAD-based SD, on a log scale
    lv = np.log1p(np.clip(v, 0, None))
    med = float(np.median(lv))
    sigma = 1.4826 * float(np.median(np.abs(lv - med)))
    if sigma <= 0:
        return float(np.inf)
    return float(np.expm1(med + cfg.dead_k * sigma))


def flag_dead_cells(df: pd.DataFrame, cfg: Config) -> Tuple[pd.DataFrame, Dict]:
    """Add a boolean 'dead' column for hyper-bright (condensed) nuclei.

    Apoptotic / pyknotic nuclei concentrate the same amount of DNA into a much
    smaller volume, so they show up as a small, very bright DAPI population that
    also tends to be non-specifically bright in every other channel. Left in,
    they inflate double- and triple-positive fractions.
    """
    col, chan = dapi_column(df, cfg.dead_stat)
    if col not in df.columns:
        raise ValueError(f"Column '{col}' not found - check --dead-stat")

    scope_key = {"screen": None, "plate": "plate",
                 "well": ["plate", "well"], "image": "image_id"}[cfg.dead_scope]
    info: Dict[str, Dict] = {}
    df = df.copy()

    def _apply(sub: pd.DataFrame, key: str) -> pd.DataFrame:
        thr = compute_dead_threshold(sub[col].to_numpy(), cfg)
        dead = sub[col] > thr
        area_cut = None
        if cfg.dead_require_small and "area" in sub.columns:
            area_cut = float(np.quantile(sub["area"].to_numpy(), cfg.dead_area_quantile))
            dead = dead & (sub["area"] <= area_cut)
        sub["dead"] = dead
        info[key] = {"dapi_column": col, "threshold": thr, "area_cut": area_cut,
                     "n_dead": int(dead.sum()), "n_cells": int(len(sub))}
        return sub

    if scope_key is None:
        df = _apply(df, "screen")
    else:
        df = pd.concat([_apply(sub, str(key))
                        for key, sub in df.groupby(scope_key, sort=False)],
                       ignore_index=True)
    info["_channel"] = chan
    info["_column"] = col
    return df, info


def gate(df: pd.DataFrame, markers: Sequence[str], cfg: Config,
         stat: str) -> Tuple[pd.DataFrame, Dict]:
    """Add per-marker positivity + population class. Returns (df, thresholds)."""
    scope_key = {"screen": None, "plate": "plate",
                 "well": ["plate", "well"], "image": "image_id"}[cfg.threshold_scope]

    thresholds: Dict[str, Dict] = {}
    df = df.copy()

    def _apply(sub: pd.DataFrame, key: str) -> pd.DataFrame:
        thresholds[key] = {}
        ctrl_mask = sub["well"].isin([w.upper() for w in cfg.control_wells]) \
            if cfg.control_wells else None
        for m in markers:
            col = f"{stat}_{m}"
            fcol = f"foci_n_{m}"
            if cfg.foci_as_marker and fcol in sub.columns:
                sub[f"pos_{m}"] = sub[fcol] >= cfg.foci_min_count
                thresholds[key][m] = f"foci_count>={cfg.foci_min_count}"
                continue
            ctrl = sub.loc[ctrl_mask, col].to_numpy() if ctrl_mask is not None else None
            thr = compute_threshold(sub[col].to_numpy(), cfg, ctrl, marker=m)
            thresholds[key][m] = thr
            sub[f"pos_{m}"] = sub[col] > thr
        return sub

    if scope_key is None:
        df = _apply(df, "screen")
    else:
        parts = []
        for key, sub in df.groupby(scope_key, sort=False):
            parts.append(_apply(sub, str(key)))
        df = pd.concat(parts, ignore_index=True)

    pos_cols = [f"pos_{m}" for m in markers]
    df["n_positive"] = df[pos_cols].sum(axis=1).astype(int)
    df["population"] = df["n_positive"].map(
        lambda n: POSITIVITY_NAMES.get(n, f"{n}-positive"))
    df["combination"] = df[pos_cols].apply(
        lambda r: "+".join([m for m, v in zip(markers, r) if v]) + "+" if r.any()
        else "all-negative", axis=1)
    return df, thresholds


# --------------------------------------------------------------------------- #
# Flow-cytometry style plotting
# --------------------------------------------------------------------------- #


def get_cmap(name: str, n: int):
    """Colormap with n discrete colours, across matplotlib versions.
    plt.get_cmap was deprecated in 3.7 and removed in 3.9."""
    n = max(int(n), 2)
    try:
        return matplotlib.colormaps[name].resampled(n)      # matplotlib >= 3.6
    except (AttributeError, KeyError):
        return plt.cm.get_cmap(name, n)                     # older matplotlib


def pick_cofactor(values: np.ndarray, cfg: Config) -> float:
    if cfg.cofactor:
        return float(cfg.cofactor)
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return 1.0
    return float(max(1.0, 0.5 * np.percentile(np.abs(v), 25)))


def asinh(x, cofactor: float):
    return np.arcsinh(np.asarray(x, dtype=float) / cofactor)


def set_asinh_axis(ax, which: str, cofactor: float, vmax: float, label: str):
    decades, d = [0], 10.0
    while d <= max(vmax, 10) * 10:
        decades.append(d)
        d *= 10
    pos = asinh(decades, cofactor)
    labels = ["0"] + [f"$10^{{{int(np.log10(t))}}}$" for t in decades[1:]]
    if which == "x":
        ax.set_xticks(pos)
        ax.set_xticklabels(labels)
        ax.set_xlabel(f"{label}  (asinh, cf={cofactor:.3g})")
    else:
        ax.set_yticks(pos)
        ax.set_yticklabels(labels)
        ax.set_ylabel(f"{label}  (asinh, cf={cofactor:.3g})")


def density_color(x: np.ndarray, y: np.ndarray, bins: int = 128) -> np.ndarray:
    H, xe, ye = np.histogram2d(x, y, bins=bins)
    ix = np.clip(np.digitize(x, xe) - 1, 0, bins - 1)
    iy = np.clip(np.digitize(y, ye) - 1, 0, bins - 1)
    return H[ix, iy]


def subsample(df: pd.DataFrame, n: int, seed: int = 0) -> pd.DataFrame:
    return df.sample(n, random_state=seed) if len(df) > n else df


def plot_histograms(df, markers, thresholds, cofactors, stat, out: Path, title: str):
    n = len(markers)
    fig, axes = plt.subplots(1, n, figsize=(4.2 * n, 3.6), squeeze=False)
    for ax, m in zip(axes[0], markers):
        v = df[f"{stat}_{m}"].to_numpy()
        cf = cofactors[m]
        thr = thresholds[m]
        ax.hist(asinh(v, cf), bins=200, color="#3b6ea5", alpha=0.85)
        ax.axvline(asinh(thr, cf), color="crimson", lw=1.4, ls="--")
        pct = 100.0 * np.mean(v > thr)
        ax.set_title(f"{m}\n{pct:.1f}% positive  (thr={thr:.4g})", fontsize=10)
        set_asinh_axis(ax, "x", cf, np.nanmax(v), f"{stat} {m}")
        ax.set_ylabel("cells")
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_biaxial(df, mx, my, thresholds, cofactors, stat, out: Path, title: str,
                 max_points: int):
    sub = subsample(df, max_points)
    x = asinh(sub[f"{stat}_{mx}"].to_numpy(), cofactors[mx])
    y = asinh(sub[f"{stat}_{my}"].to_numpy(), cofactors[my])
    tx = asinh(thresholds[mx], cofactors[mx])
    ty = asinh(thresholds[my], cofactors[my])

    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    if len(sub) > 5:
        c = density_color(x, y)
        order = np.argsort(c)
        ax.scatter(x[order], y[order], c=c[order], s=4, cmap="viridis",
                   linewidths=0, rasterized=True)
    ax.axvline(tx, color="crimson", lw=1.1, ls="--")
    ax.axhline(ty, color="crimson", lw=1.1, ls="--")

    full = df
    fx = full[f"{stat}_{mx}"] > thresholds[mx]
    fy = full[f"{stat}_{my}"] > thresholds[my]
    tot = max(len(full), 1)
    quads = {
        "UR": 100 * (fx & fy).sum() / tot,
        "UL": 100 * (~fx & fy).sum() / tot,
        "LL": 100 * (~fx & ~fy).sum() / tot,
        "LR": 100 * (fx & ~fy).sum() / tot,
    }
    for key, (hx, hy, ha, va) in {
        "UR": (0.98, 0.98, "right", "top"),
        "UL": (0.02, 0.98, "left", "top"),
        "LL": (0.02, 0.02, "left", "bottom"),
        "LR": (0.98, 0.02, "right", "bottom"),
    }.items():
        ax.text(hx, hy, f"{quads[key]:.1f}%", transform=ax.transAxes,
                ha=ha, va=va, fontsize=10, color="black",
                bbox=dict(fc="white", ec="none", alpha=0.7, pad=1.5))

    set_asinh_axis(ax, "x", cofactors[mx], df[f"{stat}_{mx}"].max(), f"{stat} {mx}")
    set_asinh_axis(ax, "y", cofactors[my], df[f"{stat}_{my}"].max(), f"{stat} {my}")
    ax.set_title(f"{title}\nn = {len(full):,} cells", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_dead_gate(df, col: str, chan: str, thr: float, area_cut: Optional[float],
                   cfg: Config, out: Path, title: str):
    """DAPI histogram + area-vs-DAPI scatter showing what the dead filter removes."""
    v = df[col].to_numpy()
    cf = pick_cofactor(v, cfg)
    dead = df["dead"].to_numpy()
    pct = 100.0 * dead.mean() if dead.size else 0.0

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    ax = axes[0]
    ax.hist(asinh(v, cf), bins=200, color="#3b6ea5", alpha=0.85)
    ax.axvline(asinh(thr, cf), color="crimson", lw=1.4, ls="--")
    ax.set_title(f"{chan} intensity - {pct:.1f}% flagged dead/dying\n"
                 f"(thr={thr:.4g}, method={cfg.dead_method})", fontsize=10)
    set_asinh_axis(ax, "x", cf, np.nanmax(v), f"{cfg.dead_stat} {chan}")
    ax.set_ylabel("nuclei")
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    sub = subsample(df, cfg.max_points)
    y = asinh(sub[col].to_numpy(), cf)
    x = sub["area"].to_numpy() if "area" in sub.columns else np.zeros(len(sub))
    d = sub["dead"].to_numpy()
    ax.scatter(x[~d], y[~d], s=4, c="#9bb7d4", linewidths=0, rasterized=True,
               label="live")
    ax.scatter(x[d], y[d], s=5, c="crimson", linewidths=0, rasterized=True,
               label="dead / dying")
    ax.axhline(asinh(thr, cf), color="crimson", lw=1.1, ls="--")
    if area_cut is not None:
        ax.axvline(area_cut, color="crimson", lw=1.1, ls=":")
    set_asinh_axis(ax, "y", cf, np.nanmax(v), f"{cfg.dead_stat} {chan}")
    ax.set_xlabel("nuclear area (px)")
    ax.legend(frameon=False, fontsize=8, markerscale=2.5, loc="upper right")
    ax.set_title("condensed nuclei sit top-left", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_foci(df, chan: str, stat: str, cfg: Config, out: Path, title: str):
    """Foci-count distribution and foci-vs-intensity agreement."""
    n = df[f"foci_n_{chan}"].to_numpy()
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.0))

    ax = axes[0]
    top = int(min(np.percentile(n, 99.5) + 1, 60)) if n.size else 1
    ax.hist(np.clip(n, 0, top), bins=np.arange(-0.5, top + 1.5), color="#6a4c93")
    ax.axvline(cfg.foci_min_count - 0.5, color="crimson", lw=1.4, ls="--")
    pct = 100.0 * np.mean(n >= cfg.foci_min_count) if n.size else 0.0
    ax.set_xlabel(f"{chan} foci per nucleus (clipped at {top})")
    ax.set_ylabel("nuclei")
    ax.set_title(f"{pct:.1f}% with >= {cfg.foci_min_count} foci\n"
                 f"median {np.median(n):.0f}, mean {np.mean(n):.1f}", fontsize=10)

    ax = axes[1]
    icol = f"{stat}_{chan}"
    if icol in df.columns:
        sub = subsample(df, cfg.max_points)
        cf = pick_cofactor(df[icol].to_numpy(), cfg)
        y = asinh(sub[icol].to_numpy(), cf)
        jitter = np.random.default_rng(0).uniform(-.3, .3, len(sub))
        x = sub[f"foci_n_{chan}"].to_numpy() + jitter
        c = density_color(x, y)
        o = np.argsort(c)
        ax.scatter(x[o], y[o], c=c[o], s=4, cmap="viridis", linewidths=0,
                   rasterized=True)
        set_asinh_axis(ax, "y", cf, df[icol].max(), f"{stat} {chan}")
        ax.set_xlim(-1, top)
        ax.set_xlabel("foci per nucleus")
        ax.set_title("foci count vs nuclear intensity", fontsize=10)

    ax = axes[2]
    col = f"foci_area_fraction_{chan}"
    if col in df.columns:
        ax.hist(100 * df[col].fillna(0).to_numpy(), bins=100, color="#2a9d8f")
        ax.set_xlabel("% of nuclear area occupied by foci")
        ax.set_ylabel("nuclei")
        ax.set_title("foci burden", fontsize=10)

    for a in axes:
        a.spines[["top", "right"]].set_visible(False)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_population_bars(df, out: Path, title: str, by: str = "well",
                         counts: bool = False):
    """Stacked bars of population classes per group.

    counts=False plots percentages (composition); counts=True plots absolute
    cell numbers, so total bar height is the well's cell count. Both are worth
    looking at: a well can have a healthy-looking composition while having lost
    most of its cells.
    """
    order = [v for k, v in sorted(POSITIVITY_NAMES.items())]
    tab = (df.groupby([by, "population"]).size()
             .unstack(fill_value=0))
    tab = tab[[c for c in order if c in tab.columns] +
              [c for c in tab.columns if c not in order]]
    frac = tab if counts else 100 * tab.div(tab.sum(axis=1), axis=0)

    fig, ax = plt.subplots(figsize=(max(6, 0.32 * len(frac)), 4.2))
    bottom = np.zeros(len(frac))
    cmap = get_cmap("cividis", len(frac.columns))
    for i, col in enumerate(frac.columns):
        ax.bar(frac.index, frac[col], bottom=bottom, label=col, color=cmap(i))
        bottom += frac[col].to_numpy()
    ax.set_ylabel("cells" if counts else "% of cells")
    ax.set_xlabel(by)
    if not counts:
        ax.set_ylim(0, 100)
    else:
        totals = tab.sum(axis=1)
        ax.axhline(totals.median(), color="crimson", lw=1.0, ls=":",
                   label=f"median {totals.median():.0f}")
    ax.set_title(title, fontsize=11)
    ax.tick_params(axis="x", rotation=90, labelsize=7)
    ax.legend(frameon=False, fontsize=8, ncol=len(frac.columns),
              loc="upper center", bbox_to_anchor=(0.5, -0.28))
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_combinations_by_well(df, out: Path, title: str, by: str = "well",
                              counts: bool = False):
    """Stacked bars per group showing WHICH markers, not just how many."""
    tab = df.groupby([by, "combination"]).size().unstack(fill_value=0)
    order = ["all-negative"] + sorted(
        [c for c in tab.columns if c != "all-negative"],
        key=lambda c: (c.count("+"), c))
    tab = tab[[c for c in order if c in tab.columns]]
    frac = tab if counts else 100 * tab.div(tab.sum(axis=1), axis=0)

    fig, ax = plt.subplots(figsize=(max(6, 0.34 * len(frac)), 4.6))
    bottom = np.zeros(len(frac))
    cmap = get_cmap("turbo", len(frac.columns))
    for i, col in enumerate(frac.columns):
        colour = "0.85" if col == "all-negative" else cmap(i)
        ax.bar(frac.index, frac[col], bottom=bottom, label=col, color=colour)
        bottom += frac[col].to_numpy()
    ax.set_ylabel("cells" if counts else "% of cells")
    ax.set_xlabel(by)
    if not counts:
        ax.set_ylim(0, 100)
    ax.set_title(title, fontsize=11)
    ax.tick_params(axis="x", rotation=90, labelsize=7)
    ax.legend(frameon=False, fontsize=7, ncol=min(len(frac.columns), 4),
              loc="upper center", bbox_to_anchor=(0.5, -0.30))
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_combinations(df, out: Path, title: str):
    counts = df["combination"].value_counts()
    frac = 100 * counts / counts.sum()
    fig, ax = plt.subplots(figsize=(max(5, 0.6 * len(frac)), 4.0))
    ax.bar(frac.index, frac.to_numpy(), color="#4c7a34")
    for i, v in enumerate(frac.to_numpy()):
        ax.text(i, v, f"{v:.1f}", ha="center", va="bottom", fontsize=8)
    ax.set_ylabel("% of cells")
    ax.set_title(title, fontsize=11)
    ax.tick_params(axis="x", rotation=45, labelsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_platemap(df, value_col: str, out: Path, title: str,
                  cbar_label: str = "%", cmap: str = "magma",
                  annotate: bool = False):
    """Heatmap of a per-well value laid out on the physical plate."""
    piv = df.pivot_table(index="well_row", columns="well_col",
                         values=value_col, aggfunc="mean")
    fig, ax = plt.subplots(figsize=(0.42 * max(piv.shape[1], 8) + 2,
                                    0.42 * max(piv.shape[0], 6) + 1.6))
    im = ax.imshow(piv.to_numpy(), cmap=cmap, aspect="equal")
    if annotate and piv.size <= 96:
        vmid = np.nanmean([np.nanmin(piv.to_numpy()), np.nanmax(piv.to_numpy())])
        for r in range(piv.shape[0]):
            for c in range(piv.shape[1]):
                v = piv.to_numpy()[r, c]
                if np.isfinite(v):
                    ax.text(c, r, f"{v:.0f}", ha="center", va="center",
                            fontsize=6,
                            color="white" if v < vmid else "black")
    ax.set_xticks(range(piv.shape[1]))
    ax.set_xticklabels([str(c + 1) for c in piv.columns], fontsize=7)
    ax.set_yticks(range(piv.shape[0]))
    ax.set_yticklabels([chr(ord("A") + r) for r in piv.index], fontsize=7)
    ax.set_title(title, fontsize=10)
    fig.colorbar(im, ax=ax, shrink=0.8, label=cbar_label)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def save_foci_crop_montage(plane: np.ndarray, nuc_labels: np.ndarray,
                           foci_labels: np.ndarray, out: Path, title: str,
                           n: int = 12, pad: int = 4, seed: int = 0,
                           line_width: int = 2):
    """Zoomed crops of individual nuclei with detected foci circled.

    At 20x a focus is 1-3 pixels wide and completely invisible in a full-field
    overlay downsampled to fit a PNG. This is the only view that lets you
    actually check foci parameters at low magnification.
    """
    ids = np.unique(nuc_labels)
    ids = ids[ids != 0]
    if ids.size == 0:
        return
    rng = np.random.default_rng(seed)
    pick = rng.choice(ids, size=int(min(n, ids.size)), replace=False)

    cols = min(4, len(pick))
    rows = int(np.ceil(len(pick) / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(2.6 * cols, 2.8 * rows),
                             squeeze=False)
    objs = ndi.find_objects(nuc_labels)
    for ax, lab in zip(axes.ravel(), pick):
        sl = objs[int(lab) - 1]
        if sl is None:
            ax.axis("off")
            continue
        y0 = max(sl[0].start - pad, 0); y1 = min(sl[0].stop + pad, plane.shape[0])
        x0 = max(sl[1].start - pad, 0); x1 = min(sl[1].stop + pad, plane.shape[1])
        crop = plane[y0:y1, x0:x1]
        fcrop = foci_labels[y0:y1, x0:x1]
        ncrop = nuc_labels[y0:y1, x0:x1] == lab

        ax.imshow(normalize_percentile(crop, 1, 99.5), cmap="gray",
                  interpolation="nearest")
        ax.contour(ncrop.astype(float), levels=[0.5], colors="#4da6ff",
                   linewidths=0.5 * line_width)
        fids = np.unique(fcrop[ncrop])
        fids = fids[fids != 0]
        for fid in fids:
            ys, xs = np.where(fcrop == fid)
            ax.plot(xs.mean(), ys.mean(), "o", mfc="none", mec="#ff3b30",
                    ms=9, mew=0.6 * line_width)
        ax.set_title(f"nucleus {lab}: {len(fids)} foci", fontsize=8)
        ax.axis("off")
    for ax in axes.ravel()[len(pick):]:
        ax.axis("off")
    fig.suptitle(f"{title}  (red = detected foci, blue = nuclear outline)",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(out, dpi=170)
    plt.close(fig)


def save_overlay(dapi: np.ndarray, labels: np.ndarray, out: Path, title: str,
                 foci_labels: Optional[np.ndarray] = None, line_width: int = 2):
    """Nuclear outlines (red) and, if given, detected foci (green) on DAPI.

    line_width is the approximate outline thickness in pixels. On a 2048x2048
    field rendered into a ~1000 px PNG a 1 px outline is close to invisible,
    so the default is deliberately thicker than one pixel.
    """
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.imshow(normalize_percentile(dapi), cmap="gray")

    edges = find_boundaries(labels, mode="outer")
    grow = max(0, int(round(line_width)) - 1)
    if grow:
        edges = ndi.binary_dilation(edges, structure=np.ones((3, 3), bool),
                                    iterations=grow)
    overlay = np.zeros((*edges.shape, 4))
    overlay[edges] = (1, 0.3, 0.1, 1)

    if foci_labels is not None and foci_labels.max() > 0:
        fmask = foci_labels > 0
        if grow:                       # foci are 1-3 px at 20x; grow to be seen
            fmask = ndi.binary_dilation(fmask, structure=np.ones((3, 3), bool),
                                        iterations=max(1, int(round(grow / 2))))
        overlay[fmask] = (0.2, 1.0, 0.4, 1)

    ax.imshow(overlay)
    n_foci = int(foci_labels.max()) if foci_labels is not None else 0
    extra = f", {n_foci} foci" if n_foci else ""
    ax.set_title(f"{title}  |  {len(np.unique(labels)) - 1} nuclei{extra}", fontsize=9)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# OMERO write-back
# --------------------------------------------------------------------------- #


def upload_masks_as_rois(conn, image, labels: np.ndarray, z: int = 0, t: int = 0):
    """Store each nucleus as an OMERO Mask ROI (one ROI per image)."""
    import omero
    from omero.rtypes import rint, rdouble, rstring

    roi = omero.model.RoiI()
    roi.setImage(image._obj)
    ids = np.unique(labels)
    ids = ids[ids != 0]
    for lab in ids:
        ys, xs = np.where(labels == lab)
        y0, y1, x0, x1 = ys.min(), ys.max() + 1, xs.min(), xs.max() + 1
        sub = (labels[y0:y1, x0:x1] == lab)
        m = omero.model.MaskI()
        m.setX(rdouble(float(x0)))
        m.setY(rdouble(float(y0)))
        m.setWidth(rdouble(float(x1 - x0)))
        m.setHeight(rdouble(float(y1 - y0)))
        m.setTheZ(rint(z))
        m.setTheT(rint(t))
        m.setTextValue(rstring(f"nucleus_{lab}"))
        m.setBytes(np.packbits(sub.astype(np.uint8)).tobytes())
        roi.addShape(m)
    conn.getUpdateService().saveAndReturnObject(roi)


def attach_file(conn, obj, path: Path, ns="omero_nuclei_flow"):
    mime = {"csv": "text/csv", "png": "image/png", "gz": "application/gzip",
            "json": "application/json"}.get(path.suffix.lstrip("."), "application/octet-stream")
    ann = conn.createFileAnnfromLocalFile(str(path), mimetype=mime, ns=ns)
    obj.linkAnnotation(ann)


# --------------------------------------------------------------------------- #
# Summaries
# --------------------------------------------------------------------------- #


EXPLORER_NAME = "gating_explorer.html"


def write_explorer(df: pd.DataFrame, outdir: Path, tag: str,
                   embed: bool = True) -> Optional[Path]:
    """Copy the interactive explorer next to the results, data optionally inlined.

    Looks for gating_explorer.html beside this script. Embedding makes the
    output a single file a student can double-click - no file dialog, no server,
    no Python. Above ~80 MB of CSV we skip embedding, because the browser has to
    hold the string and the parsed columns at once.
    """
    src = Path(__file__).resolve().parent / EXPLORER_NAME
    if not src.exists():
        LOG.warning("Explorer template not found next to the script (%s) - "
                    "skipping --export-explorer", src)
        return None
    html = src.read_text(encoding="utf-8")
    dest = outdir / EXPLORER_NAME

    if embed:
        csv = df.to_csv(index=False)
        mb = len(csv) / 1e6
        if mb > 80:
            LOG.warning("Per-nucleus table is %.0f MB - writing the explorer "
                        "without embedded data; open the CSV from it instead.", mb)
            embed = False
        else:
            payload = json.dumps(csv)
            html = html.replace(
                "if(window.EMBEDDED_CSV){",
                f"window.EMBEDDED_CSV = {payload};\n"
                f"window.EMBEDDED_NAME = {json.dumps(safe(tag) + '.csv')};\n"
                "if(window.EMBEDDED_CSV){", 1)
            LOG.info("Explorer written with %.1f MB of data embedded", mb)
    dest.write_text(html, encoding="utf-8")
    if not embed:
        LOG.info("Explorer written (open it and drop the per-nucleus CSV in)")
    return dest


def write_table(df: pd.DataFrame, stem: Path, fmt: str) -> Optional[Path]:
    """Write a dataframe in the requested format. Returns the path written."""
    fmt = (fmt or "csv.gz").lower()
    if fmt == "none":
        LOG.info("Skipping %s (--per-nucleus-format none)", stem.name)
        return None
    if fmt == "xlsx" and len(df) > 1_048_575:
        LOG.warning("%d rows exceeds the Excel limit - writing csv.gz instead", len(df))
        fmt = "csv.gz"
    try:
        if fmt == "csv":
            p = stem.with_suffix(".csv"); df.to_csv(p, index=False)
        elif fmt == "parquet":
            p = stem.with_suffix(".parquet"); df.to_parquet(p, index=False)
        elif fmt == "xlsx":
            p = stem.with_suffix(".xlsx"); df.to_excel(p, index=False)
        else:
            p = Path(str(stem) + ".csv.gz")
            df.to_csv(p, index=False, compression="gzip")
    except ImportError as exc:
        LOG.warning("Cannot write %s (%s) - falling back to csv.gz", fmt, exc)
        p = Path(str(stem) + ".csv.gz")
        df.to_csv(p, index=False, compression="gzip")
    LOG.info("Wrote %s (%d rows)", p.name, len(df))
    return p


def summarise(df: pd.DataFrame, markers: Sequence[str], by: Sequence[str]) -> pd.DataFrame:
    g = df.groupby(list(by), sort=False)
    out = g.size().rename("n_cells").to_frame()
    for m in markers:
        out[f"pct_{m}_positive"] = 100 * g[f"pos_{m}"].mean()
        out[f"n_{m}_positive"] = g[f"pos_{m}"].sum().astype(int)

    # how many markers (negative / single / double / ...), observed classes only
    seen = set(df["population"].unique())
    ordered = [v for _, v in sorted(POSITIVITY_NAMES.items()) if v in seen]
    ordered += sorted(p for p in seen if p not in POSITIVITY_NAMES.values())
    for name in ordered:
        out[f"pct_{name}"] = 100 * g["population"].apply(lambda s, n=name: (s == n).mean())
        out[f"n_{name}"] = g["population"].apply(lambda s, n=name: int((s == n).sum()))

    # WHICH markers, not just how many: one column per observed combination,
    # ordered by degree then alphabetically (GFP+, RFP+, GFP+RFP+, GFP+RFP+Cy5+)
    combos = [c for c in df["combination"].unique() if c != "all-negative"]
    for combo in sorted(combos, key=lambda c: (c.count("+"), c)):
        out[f"pct_{combo}"] = 100 * g["combination"].apply(
            lambda s, c=combo: (s == c).mean())
        out[f"n_{combo}"] = g["combination"].apply(
            lambda s, c=combo: int((s == c).sum()))

    out["pct_any_positive"] = 100 * g["n_positive"].apply(lambda s: (s > 0).mean())
    for chan in foci_channel_columns(df):
        out[f"mean_foci_{chan}"] = g[f"foci_n_{chan}"].mean()
        out[f"median_foci_{chan}"] = g[f"foci_n_{chan}"].median()
        if f"foci_pos_{chan}" in df.columns:
            out[f"pct_foci_positive_{chan}"] = 100 * g[f"foci_pos_{chan}"].mean()
    return out.reset_index()


# --------------------------------------------------------------------------- #
# Main analysis
# --------------------------------------------------------------------------- #


def process_plate(conn, plate, cfg: Config, segmenter: Segmenter,
                  inherited_kv: Optional[Dict[str, str]] = None
                  ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rows: List[pd.DataFrame] = []
    foci_rows: List[pd.DataFrame] = []
    plate_name = plate.getName()
    overlays_done = 0
    overlay_dir = cfg.outdir / safe(plate_name) / "plots" / "overlays"

    # ---- OMERO key-value pairs: screen -> plate -> well -> image ----------
    base_kv: Dict[str, str] = dict(inherited_kv or {})
    if cfg.fetch_kv and "screen" in cfg.kv_levels and not base_kv:
        # plate-only runs (--plate-id) still inherit their screen's annotations
        try:
            parent = plate.getParent()
            if parent is not None:
                base_kv.update(map_annotations(parent))
                LOG.debug("  inherited key-values from screen '%s'", parent.getName())
        except Exception as exc:                      # noqa: BLE001
            LOG.debug("  no parent screen for plate %s: %s", plate_name, exc)
    if cfg.fetch_kv and "plate" in cfg.kv_levels:
        base_kv.update(map_annotations(plate))
    well_kv_cache: Dict[str, Dict[str, str]] = {}
    well_meta_rows: List[Dict[str, str]] = []
    foci_params: Optional[FociParams] = None
    crops_done = 0
    chan_cache: Dict[int, List[str]] = {}
    strip_height = cfg.server_tile_height

    for label, row, col, fld, image, well in iter_images(plate, cfg):
        try:
            if cfg.fetch_kv and label not in well_kv_cache:
                wkv = dict(base_kv)
                if "well" in cfg.kv_levels:
                    wkv.update(map_annotations(well))
                well_kv_cache[label] = wkv
                well_meta_rows.append({"plate": plate_name, "well": label,
                                       "well_row": row, "well_col": col,
                                       **kv_to_columns(wkv, cfg.kv_prefix)})

            # Channel metadata is identical across a plate; querying it once
            # per image is thousands of needless server round trips.
            size_c = image.getSizeC()
            if size_c not in chan_cache:
                chan_cache[size_c] = channel_names(image)
                LOG.info("  channels (%d): %s", size_c,
                         ", ".join(chan_cache[size_c]))
            names = chan_cache[size_c]
            dapi_idx = resolve_dapi(names, cfg.dapi_channel)

            try:
                stack = read_stack(image, cfg, tile_height=strip_height)
            except Exception as exc:                  # noqa: BLE001
                if not _is_server_oom(exc):
                    raise
                strip_height = strip_height or 512
                LOG.warning("OMERO server ran out of Java heap reading image %s. "
                            "Retrying in %d-row strips, and using strips for the "
                            "rest of this plate.", image.getId(), strip_height)
                stack = read_stack(image, cfg, tile_height=strip_height)
            labels = filter_labels(segmenter(stack[..., dapi_idx]), cfg)
            n = int(labels.max())
            LOG.info("  %s field %d (image %d): %d nuclei", label, fld, image.getId(), n)
            if n == 0:
                continue

            meas = measure(labels, stack, names, cfg)
            if meas.empty:
                continue

            foci_overlay = None
            if cfg.foci_channels:
                if foci_params is None:
                    pixel_um = None
                    try:
                        pixel_um = image.getPixelSizeX()
                    except Exception:                 # noqa: BLE001
                        pass
                    foci_params = derive_foci_params(cfg, pixel_um)
                    if pixel_um:
                        LOG.info("  pixel size %.4f um/px -> foci params: %s",
                                 pixel_um, foci_params.describe())
                    else:
                        LOG.warning("  no pixel size in OMERO metadata - using "
                                    "raw pixel parameters (%s). Pass "
                                    "--pixel-size-um for auto-scaling.",
                                    foci_params.describe())
                for fidx in resolve_channels(names, cfg.foci_channels):
                    fchan = names[fidx]
                    plane = stack[..., fidx]
                    if cfg.background == "median" and (labels == 0).any():
                        plane = plane - float(np.median(plane[labels == 0]))
                    flab, props = detect_foci(plane, labels, cfg, foci_params)
                    if foci_overlay is None:
                        foci_overlay = flab
                        if cfg.save_foci_crops and crops_done < cfg.save_foci_crops:
                            overlay_dir.mkdir(parents=True, exist_ok=True)
                            save_foci_crop_montage(
                                plane, labels, flab,
                                overlay_dir / f"foci_crops_{safe(label)}_f{fld}.png",
                                f"{plate_name} {label} field {fld} - {fchan}",
                                line_width=cfg.overlay_line_width)
                            crops_done += 1
                    meas = summarise_foci(props, meas, fchan)
                    if not props.empty:
                        props = props.copy()
                        props.insert(0, "plate", plate_name)
                        props.insert(1, "well", label)
                        props.insert(2, "field", fld)
                        props.insert(3, "image_id", image.getId())
                        props.insert(4, "channel", fchan)
                        foci_rows.append(props)
                    LOG.debug("    %s: %d foci in %d nuclei", fchan, len(props), n)

            meas.insert(0, "plate", plate_name)
            meas.insert(1, "well", label)
            meas.insert(2, "well_row", row)
            meas.insert(3, "well_col", col)
            meas.insert(4, "field", fld)
            meas.insert(5, "image_id", image.getId())
            meas.insert(6, "image_name", image.getName())
            meas.insert(7, "dapi_channel", names[dapi_idx])

            if cfg.fetch_kv:
                kv = dict(well_kv_cache.get(label, {}))
                if "image" in cfg.kv_levels:
                    kv.update(map_annotations(image))
                at = 8
                for c, v in kv_to_columns(kv, cfg.kv_prefix).items():
                    if c in meas.columns:
                        LOG.warning("Key-value column %s collides with a "
                                    "measurement column - skipping", c)
                        continue
                    meas.insert(at, c, v)
                    at += 1

            rows.append(meas)

            if cfg.save_overlays and overlays_done < cfg.save_overlays:
                overlay_dir.mkdir(parents=True, exist_ok=True)
                save_overlay(stack[..., dapi_idx], labels,
                             overlay_dir / f"{safe(label)}_f{fld}_{image.getId()}.png",
                             f"{plate_name} {label} field {fld}",
                             foci_labels=foci_overlay,
                             line_width=cfg.overlay_line_width)
                overlays_done += 1

            if cfg.upload_rois:
                upload_masks_as_rois(conn, image, labels)

        except Exception as exc:                      # noqa: BLE001
            if _is_server_oom(exc):
                LOG.error("  image %s (%s): the OMERO SERVER ran out of Java "
                          "heap. This is a server-side limit, not your machine. "
                          "See the 'Java heap space' section of the README - the "
                          "quickest client-side lever is --server-tile-height 256.",
                          image.getId(), label)
            else:
                LOG.exception("  failed on image %s (%s): %s",
                              image.getId(), label, exc)
        finally:
            # release any server-side rendering engine held for this image
            try:
                image._closeRE()
            except Exception:                         # noqa: BLE001
                pass

    nuclei = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    foci = pd.concat(foci_rows, ignore_index=True) if foci_rows else pd.DataFrame()
    well_meta = pd.DataFrame(well_meta_rows) if well_meta_rows else pd.DataFrame()
    if not nuclei.empty:
        nuclei = coerce_numeric_kv(nuclei, cfg.kv_prefix)
    if not well_meta.empty:
        well_meta = coerce_numeric_kv(well_meta, cfg.kv_prefix)
        n_keys = len([c for c in well_meta.columns if c.startswith(cfg.kv_prefix)])
        LOG.info("  key-value pairs: %d keys across %d wells", n_keys, len(well_meta))
    return nuclei, foci, well_meta


def safe(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(name)).strip("_") or "unnamed"


def detect_markers(df: pd.DataFrame, stat: str, cfg: Config) -> List[str]:
    all_ch = [c[len(stat) + 1:] for c in df.columns if c.startswith(f"{stat}_")]
    dapi = set(df["dapi_channel"].unique()) if "dapi_channel" in df else set()
    markers = [c for c in all_ch if c not in dapi]
    if cfg.marker_channels:
        wanted = {m.lower() for m in cfg.marker_channels}
        markers = [m for m in markers if m.lower() in wanted]
    return markers


def marker_subsets(markers: Sequence[str], cfg: Config) -> List[Tuple[str, ...]]:
    """Which marker combinations to analyse as self-contained panels.

    With three markers on a plate, wells stained with only two of them are
    misdescribed by the plate-wide three-marker classification: a cell that is
    GFP+RFP+ and simply wasn't stained for Cy5 is reported as 'double' out of a
    possible three, when for that panel it is the maximum. Re-gating each
    subset independently makes each one a closed universe whose populations sum
    to 100%.
    """
    if cfg.subset_panels:
        out = []
        for spec in cfg.subset_panels:
            wanted = {s.strip().lower() for s in spec.split("+") if s.strip()}
            chosen = tuple(m for m in markers if m.lower() in wanted)
            missing = wanted - {m.lower() for m in chosen}
            if missing:
                LOG.warning("Panel '%s': no marker named %s (available: %s)",
                            spec, ", ".join(sorted(missing)), ", ".join(markers))
            if len(chosen) < 1:
                LOG.warning("Panel '%s' matched no markers in %s - skipping",
                            spec, list(markers))
                continue
            out.append(chosen)
        return out

    sizes = sorted({int(s) for s in cfg.subset_sizes if 1 <= int(s) <= len(markers)})
    out = []
    for n in sizes:
        if n >= len(markers):
            continue                      # that's the full analysis already
        out.extend(combinations(markers, n))
    return out


def analyse_subset(df: pd.DataFrame, subset: Sequence[str], cfg: Config,
                   outdir: Path, tag: str, stat: str) -> pd.DataFrame:
    """Re-gate and summarise using only `subset` of the markers.

    Positivity per marker is unchanged - the thresholds are the same - but the
    population classification (negative / single / double / ...) is recomputed
    over this subset alone, so its classes sum to 100% of the cells.
    """
    sub = df.copy()
    pos_cols = [f"pos_{m}" for m in subset]
    sub["n_positive"] = sub[pos_cols].sum(axis=1).astype(int)
    sub["population"] = sub["n_positive"].map(
        lambda n: POSITIVITY_NAMES.get(n, f"{n}-positive"))
    sub["combination"] = sub[pos_cols].apply(
        lambda r: "+".join([m for m, v in zip(subset, r) if v]) + "+" if r.any()
        else "all-negative", axis=1)

    name = "_".join(safe(m) for m in subset)
    sdir = outdir / "panels" / name
    plots = sdir / "plots"
    plots.mkdir(parents=True, exist_ok=True)
    label = " + ".join(subset)
    ptag = f"{tag} | panel: {label}"

    per_well = summarise(sub, list(subset), ["plate", "well", "well_row", "well_col"])
    kv_cols = [c for c in sub.columns if c.startswith(cfg.kv_prefix)]
    if kv_cols:
        per_well = per_well.merge(
            sub.groupby(["plate", "well"], as_index=False)[kv_cols].first(),
            on=["plate", "well"], how="left")
    per_well.to_csv(sdir / "per_well_summary.csv", index=False)
    summarise(sub, list(subset), ["plate", "well", "field", "image_id"]).to_csv(
        sdir / "per_image_summary.csv", index=False)

    plot_population_bars(sub, plots / "populations_by_well.png",
                         f"{ptag} (%)")
    plot_population_bars(sub, plots / "populations_by_well_counts.png",
                         f"{ptag} (counts)", counts=True)
    plot_combinations(sub, plots / "marker_combinations.png",
                      f"{ptag} - combinations")
    plot_combinations_by_well(sub, plots / "combinations_by_well.png",
                              f"{ptag} - which markers (%)")
    plot_combinations_by_well(sub, plots / "combinations_by_well_counts.png",
                              f"{ptag} - which markers (counts)", counts=True)

    plot_platemap(per_well, "n_cells", plots / "platemap_n_cells.png",
                  f"{ptag} - total cells", cbar_label="cells",
                  cmap="viridis", annotate=True)
    for pop in [v for _, v in sorted(POSITIVITY_NAMES.items())]:
        pcol, ncol = f"pct_{pop}", f"n_{pop}"
        plabel = pop if pop == "negative" else f"{pop} positive"
        if pcol in per_well.columns and per_well[pcol].sum() > 0:
            plot_platemap(per_well, pcol, plots / f"platemap_pct_{pop}.png",
                          f"{ptag} - % {plabel}")
        if ncol in per_well.columns and per_well[ncol].sum() > 0:
            plot_platemap(per_well, ncol, plots / f"platemap_n_{pop}.png",
                          f"{ptag} - {plabel} cells", cbar_label="cells",
                          cmap="viridis", annotate=True)

    if len(subset) == 2:
        cofactors = {m: pick_cofactor(sub[f"{stat}_{m}"].to_numpy(), cfg)
                     for m in subset}
        thr = {}
        for m in subset:
            p, n = sub.loc[sub[f"pos_{m}"], f"{stat}_{m}"], \
                   sub.loc[~sub[f"pos_{m}"], f"{stat}_{m}"]
            thr[m] = float(0.5 * (p.quantile(0.05) + n.quantile(0.95))) \
                if len(p) and len(n) else float(p.min() if len(p) else np.inf)
        plot_biaxial(sub, subset[0], subset[1], thr, cofactors, stat,
                     plots / f"biaxial_{safe(subset[0])}_vs_{safe(subset[1])}.png",
                     ptag, cfg.max_points)

    counts = sub["population"].value_counts()
    LOG.info("[%s] panel %-28s %s", tag, label,
             " | ".join(f"{k}: {100*v/len(sub):.1f}%" for k, v in counts.items()))
    return per_well


def analyse(df: pd.DataFrame, cfg: Config, outdir: Path, tag: str,
            conn=None, omero_obj=None) -> pd.DataFrame:
    """Gate, summarise and plot one plate (or a pooled screen)."""
    stat = cfg.intensity_stat
    markers = detect_markers(df, stat, cfg)
    if not markers:
        LOG.warning("[%s] no marker channels besides DAPI - nothing to gate", tag)
        return df
    LOG.info("[%s] markers: %s", tag, ", ".join(markers))

    plots = outdir / "plots"
    plots.mkdir(parents=True, exist_ok=True)

    # ---- optional dead / dying (bright DAPI) removal, BEFORE marker gating --
    dead_info: Optional[Dict] = None
    dead_per_well: Optional[pd.DataFrame] = None
    if cfg.filter_dead:
        n_total = len(df)
        df, dead_info = flag_dead_cells(df, cfg)
        col, chan = dead_info["_column"], dead_info["_channel"]
        scoped = [v for k, v in dead_info.items() if not k.startswith("_")]
        pooled_dead_thr = float(np.median([s["threshold"] for s in scoped]))
        area_cuts = [s["area_cut"] for s in scoped if s["area_cut"] is not None]
        pooled_area_cut = float(np.median(area_cuts)) if area_cuts else None
        plot_dead_gate(df, col, chan, pooled_dead_thr, pooled_area_cut, cfg,
                       plots / "dead_cell_gate.png",
                       f"{tag} - dead/dying (bright {chan}) filter")

        dead_per_well = (df.groupby(["plate", "well"])["dead"].mean() * 100) \
            .rename("pct_dead_dapi").reset_index()
        n_dead = int(df["dead"].sum())
        LOG.info("[%s] bright-%s dead/dying: %d / %d cells (%.1f%%) - action: %s",
                 tag, chan, n_dead, n_total, 100 * n_dead / max(n_total, 1),
                 cfg.dead_action)

        if cfg.dead_action == "remove" and n_dead:
            df[df["dead"]].to_csv(outdir / "excluded_dead_cells.csv.gz",
                                  index=False, compression="gzip")
            df = df.loc[~df["dead"]].reset_index(drop=True)
            if df.empty:
                LOG.warning("[%s] every cell was flagged as dead - check "
                            "--dead-k / --dead-method", tag)
                return df

    df, thresholds = gate(df, markers, cfg, stat)

    # thresholds are per scope-key; use the pooled median for plate-level plots.
    # Markers gated on foci counts have no numeric intensity threshold, so we
    # project their gate onto the intensity axis purely for the plot lines.
    def _pooled(m: str) -> float:
        vals = [v[m] for v in thresholds.values()
                if isinstance(v.get(m), (int, float)) and not isinstance(v.get(m), bool)]
        if vals:
            return float(np.median(vals))
        col, pcol = f"{stat}_{m}", f"pos_{m}"
        pos, neg = df.loc[df[pcol], col], df.loc[~df[pcol], col]
        if len(pos) and len(neg):
            return float(0.5 * (pos.quantile(0.05) + neg.quantile(0.95)))
        return float(pos.min()) if len(pos) else float("inf")

    pooled_thr = {m: _pooled(m) for m in markers}
    cofactors = {m: pick_cofactor(df[f"{stat}_{m}"].to_numpy(), cfg) for m in markers}

    plot_histograms(df, markers, pooled_thr, cofactors, stat,
                    plots / "histograms_gated.png", f"{tag} - marker distributions")
    for mx, my in combinations(markers, 2):
        plot_biaxial(df, mx, my, pooled_thr, cofactors, stat,
                     plots / f"biaxial_{safe(mx)}_vs_{safe(my)}.png",
                     f"{tag} - {mx} vs {my}", cfg.max_points)
    plot_population_bars(df, plots / "populations_by_well.png",
                         f"{tag} - population breakdown per well (%)")
    plot_population_bars(df, plots / "populations_by_well_counts.png",
                         f"{tag} - cells per population per well (counts)",
                         counts=True)
    plot_combinations(df, plots / "marker_combinations.png",
                      f"{tag} - marker combinations")
    plot_combinations_by_well(df, plots / "combinations_by_well.png",
                              f"{tag} - which markers, per well (%)")
    plot_combinations_by_well(df, plots / "combinations_by_well_counts.png",
                              f"{tag} - which markers, per well (counts)",
                              counts=True)

    foci_chans = foci_channel_columns(df)
    for chan in foci_chans:
        df[f"foci_pos_{chan}"] = df[f"foci_n_{chan}"] >= cfg.foci_min_count
        plot_foci(df, chan, stat, cfg, plots / f"foci_{safe(chan)}.png",
                  f"{tag} - {chan} foci")
        LOG.info("[%s] %s foci: mean %.2f/nucleus, %.1f%% with >= %d",
                 tag, chan, df[f"foci_n_{chan}"].mean(),
                 100 * df[f"foci_pos_{chan}"].mean(), cfg.foci_min_count)

    per_well = summarise(df, markers, ["plate", "well", "well_row", "well_col"])
    per_image = summarise(df, markers, ["plate", "well", "field", "image_id"])

    if dead_per_well is not None:
        per_well = per_well.merge(dead_per_well, on=["plate", "well"], how="left")
        plot_platemap(per_well, "pct_dead_dapi",
                      plots / "platemap_pct_dead_dapi.png",
                      f"{tag} - % dead/dying (bright DAPI)")

    # ---- carry OMERO key-value pairs into the summary tables ---------------
    kv_cols = [c for c in df.columns if c.startswith(cfg.kv_prefix)]
    if kv_cols:
        per_well = per_well.merge(
            df.groupby(["plate", "well"], as_index=False)[kv_cols].first(),
            on=["plate", "well"], how="left")
        per_image = per_image.merge(
            df.groupby(["image_id"], as_index=False)[kv_cols].first(),
            on="image_id", how="left")
        LOG.info("[%s] carried %d key-value column(s) into the summaries: %s",
                 tag, len(kv_cols), ", ".join(kv_cols[:8]) +
                 (" ..." if len(kv_cols) > 8 else ""))

    # ---- optional grouping by an experimental condition --------------------
    group_col = None
    if cfg.group_by:
        wanted = [cfg.group_by, kv_column_name(cfg.group_by, cfg.kv_prefix)]
        group_col = next((c for c in wanted if c and c in df.columns), None)
        if group_col is None:
            LOG.warning("[%s] --group-by '%s' not found. Available: %s",
                        tag, cfg.group_by, ", ".join(kv_cols) or "(none)")
    if group_col:
        per_group = summarise(df, markers, ["plate", group_col])
        per_group.to_csv(outdir / f"per_{safe(group_col)}_summary.csv", index=False)
        plot_population_bars(df, plots / f"populations_by_{safe(group_col)}.png",
                             f"{tag} - populations by {group_col} (%)",
                             by=group_col)
        plot_population_bars(df,
                             plots / f"populations_by_{safe(group_col)}_counts.png",
                             f"{tag} - cells per population by {group_col}",
                             by=group_col, counts=True)
        plot_combinations_by_well(df, plots / f"combinations_by_{safe(group_col)}.png",
                                  f"{tag} - which markers, by {group_col} (%)",
                                  by=group_col)
        LOG.info("[%s] grouped summary by %s (%d groups)",
                 tag, group_col, df[group_col].nunique())

    # ---- plate heatmaps: percentages AND absolute cell numbers -------------
    plot_platemap(per_well, "n_cells", plots / "platemap_n_cells.png",
                  f"{tag} - total cells per well", cbar_label="cells",
                  cmap="viridis", annotate=True)

    for m in markers:
        plot_platemap(per_well, f"pct_{m}_positive",
                      plots / f"platemap_pct_{safe(m)}_positive.png",
                      f"{tag} - % {m} positive")
        if f"n_{m}_positive" in per_well.columns:
            plot_platemap(per_well, f"n_{m}_positive",
                          plots / f"platemap_n_{safe(m)}_positive.png",
                          f"{tag} - {m}-positive cells per well",
                          cbar_label="cells", cmap="viridis", annotate=True)

    for pop in [v for _, v in sorted(POSITIVITY_NAMES.items())]:
        col = f"pct_{pop}"
        label = pop if pop == "negative" else f"{pop} positive"
        if col in per_well.columns and per_well[col].sum() > 0:
            plot_platemap(per_well, col, plots / f"platemap_pct_{pop}.png",
                          f"{tag} - % {label}")
        ncol = f"n_{pop}"
        if ncol in per_well.columns and per_well[ncol].sum() > 0:
            plot_platemap(per_well, ncol, plots / f"platemap_n_{pop}.png",
                          f"{tag} - {pop} cells per well",
                          cbar_label="cells", cmap="viridis", annotate=True)

    for combo in sorted({c for c in df["combination"].unique() if c.count("+") > 1}):
        col = f"pct_{combo}"
        if col in per_well.columns and per_well[col].sum() > 0:
            plot_platemap(per_well, col, plots / f"platemap_pct_{safe(combo)}.png",
                          f"{tag} - % {combo}")
    for chan in foci_chans:
        plot_platemap(per_well, f"mean_foci_{chan}",
                      plots / f"platemap_mean_foci_{safe(chan)}.png",
                      f"{tag} - mean {chan} foci per nucleus")

    if cfg.per_well_plots:
        wdir = plots / "per_well"
        wdir.mkdir(exist_ok=True)
        for well, sub in df.groupby("well"):
            if len(sub) < 20:
                continue
            plot_histograms(sub, markers, pooled_thr, cofactors, stat,
                            wdir / f"{safe(well)}_histograms.png", f"{tag} {well}")
            for mx, my in combinations(markers, 2):
                plot_biaxial(sub, mx, my, pooled_thr, cofactors, stat,
                             wdir / f"{safe(well)}_biaxial_{safe(mx)}_vs_{safe(my)}.png",
                             f"{tag} {well} - {mx} vs {my}", cfg.max_points)

    # ---- self-contained marker panels -------------------------------------
    if cfg.subset_analyses and len(markers) > 1:
        panels = marker_subsets(markers, cfg)
        if panels:
            LOG.info("[%s] %d panel analyses (each sums to 100%% on its own "
                     "markers): %s", tag, len(panels),
                     "; ".join("+".join(p) for p in panels))
            index = []
            for subset in list(panels) + [tuple(markers)]:
                is_full = len(subset) == len(markers)
                if is_full:
                    # The full set IS the top-level analysis - reuse its
                    # numbers rather than regenerating identical plots, but
                    # include the row so the index covers every level.
                    pw = per_well
                    location = "."
                else:
                    pw = analyse_subset(df, subset, cfg, outdir, tag, stat)
                    location = f"panels/{'_'.join(safe(m) for m in subset)}"
                row = {"panel": "+".join(subset), "n_markers": len(subset),
                       "n_cells": int(pw["n_cells"].sum()),
                       "results_in": location}
                # every class up to this panel's size exists, even if empty,
                # so the columns sum to 100% and don't come out as NaN
                for k in range(len(subset) + 1):
                    pop = POSITIVITY_NAMES.get(k, f"{k}-positive")
                    n = int(pw[f"n_{pop}"].sum()) if f"n_{pop}" in pw.columns else 0
                    row[f"n_{pop}"] = n
                    row[f"pct_{pop}"] = round(100 * n / max(row["n_cells"], 1), 3)
                index.append(row)
            idx = pd.DataFrame(index)
            (outdir / "panels").mkdir(parents=True, exist_ok=True)
            idx.to_csv(outdir / "panels" / "panel_index.csv", index=False)

    cells_path = write_table(df, outdir / "per_nucleus_measurements",
                             cfg.per_nucleus_format)
    per_well.to_csv(outdir / "per_well_summary.csv", index=False)
    per_image.to_csv(outdir / "per_image_summary.csv", index=False)
    with open(outdir / "thresholds.json", "w") as fh:
        json.dump({"stat": stat, "markers": markers, "method": cfg.threshold_method,
                   "scope": cfg.threshold_scope, "by_scope": thresholds,
                   "pooled": pooled_thr, "cofactors": cofactors,
                   "dead_filter": dead_info}, fh, indent=2)

    LOG.info("[%s] %d cells | %s", tag, len(df),
             " | ".join(f"{p}: {100*(df['population'] == p).mean():.1f}%"
                        for p in df["population"].value_counts().index))

    if cfg.export_explorer:
        write_explorer(df, outdir, tag, embed=not cfg.explorer_no_embed)

    if conn is not None and cfg.attach_results and omero_obj is not None:
        targets = [p for p in [cells_path, outdir / "per_well_summary.csv",
                               outdir / "thresholds.json"] if p is not None]
        for p in targets + sorted(plots.glob("*.png")):
            try:
                attach_file(conn, omero_obj, p)
            except Exception:                          # noqa: BLE001
                LOG.exception("Could not attach %s", p)
    return df


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #


def check_environment() -> int:
    """Print versions and GPU status, and flag known incompatible combinations."""
    import platform
    try:
        from importlib.metadata import version, PackageNotFoundError
    except ImportError:                                # py<3.8
        from importlib_metadata import version, PackageNotFoundError  # type: ignore

    print(f"python        {platform.python_version()}"
          f"   ({platform.system()} {platform.machine()})")

    dists = ["numpy", "pandas", "scipy", "scikit-image", "matplotlib",
             "scikit-learn", "pyarrow", "openpyxl", "omero-py", "zeroc-ice",
             "stardist", "csbdeep", "tensorflow", "keras", "tf-keras",
             "protobuf", "cellpose", "torch", "setuptools"]
    found: Dict[str, Optional[str]] = {}
    for d in dists:
        try:
            found[d] = version(d)
        except PackageNotFoundError:
            found[d] = None
    width = max(len(d) for d in dists)
    for d in dists:
        v = found[d]
        print(f"{d:<{width}}  {v if v else '-- not installed'}")

    def major_minor(v: Optional[str]) -> Tuple[int, int]:
        if not v:
            return (0, 0)
        parts = re.findall(r"\d+", v)
        return (int(parts[0]), int(parts[1])) if len(parts) > 1 else (int(parts[0]), 0)

    # ---- GPU status --------------------------------------------------------
    print()
    if found["tensorflow"]:
        try:
            import tensorflow as tf
            gpus = tf.config.list_physical_devices("GPU")
            print(f"TensorFlow GPUs visible: {len(gpus)}"
                  + (f"  {[g.name for g in gpus]}" if gpus else ""))
        except Exception as exc:                       # noqa: BLE001
            print(f"TensorFlow present but failed to import: {exc}")
    if found["torch"]:
        try:
            import torch
            print(f"PyTorch CUDA available: {torch.cuda.is_available()}")
        except Exception as exc:                       # noqa: BLE001
            print(f"PyTorch present but failed to import: {exc}")

    # ---- known-bad combinations -------------------------------------------
    problems: List[str] = []
    tf_v, np_v = major_minor(found["tensorflow"]), major_minor(found["numpy"])
    if found["numpy"] and np_v[0] >= 2 and (found["stardist"] or found["tensorflow"]):
        problems.append("numpy 2.x breaks StarDist/TensorFlow 2.10 - pin numpy<2 "
                        "(1.23.5 for the TF 2.10 stack)")
    if found["tensorflow"] and tf_v < (2, 11) and np_v >= (1, 24):
        problems.append("TensorFlow <2.11 needs numpy<1.24 - pin numpy==1.23.5")
    if found["tensorflow"] and tf_v < (2, 11) and major_minor(found["protobuf"]) >= (3, 20):
        problems.append("TensorFlow <2.11 needs protobuf<3.20 - pin protobuf==3.19.6")
    if (platform.system() == "Windows" and found["tensorflow"]
            and tf_v >= (2, 11)):
        problems.append("TensorFlow 2.11+ has no native-Windows GPU support - use "
                        "tensorflow==2.10.1 with CUDA 11.2 / cuDNN 8.1, run under "
                        "WSL2, or switch to --model cellpose (PyTorch)")
    # setuptools 82 (Feb 2026) removed pkg_resources, which StarDist still imports
    try:
        import warnings as _w
        with _w.catch_warnings():
            _w.simplefilter("ignore")
            import pkg_resources  # noqa: F401
        has_pkg_resources = True
    except Exception:                                  # noqa: BLE001
        has_pkg_resources = False
    # Keras 3 vs StarDist (TF >= 2.16 installs Keras 3 by default)
    if found["stardist"] and tf_v >= (2, 16):
        if not found["tf-keras"]:
            problems.append("TensorFlow >=2.16 ships Keras 3, which StarDist "
                            "does not support - run: pip install tf-keras "
                            "(the script then sets TF_USE_LEGACY_KERAS=1 for you)")
        elif os.environ.get("TF_USE_LEGACY_KERAS") != "1":
            print("  note: tf-keras is installed; TF_USE_LEGACY_KERAS will be "
                  "set automatically when StarDist loads")

    if found["stardist"] and not has_pkg_resources:
        problems.append("pkg_resources is missing (setuptools 82 removed it) but "
                        "StarDist imports it - run: pip install \"setuptools<81\"")
    elif not has_pkg_resources and major_minor(found["setuptools"]) >= (81, 0):
        problems.append("setuptools >=81 has dropped pkg_resources; if you later "
                        "install StarDist, pin setuptools<81 first")

    if not found["omero-py"]:
        problems.append("omero-py missing - the script cannot reach OMERO "
                        "(--from-csv still works)")
    if not (found["stardist"] or found["cellpose"]):
        problems.append("no segmentation backend - install stardist or cellpose")

    print()
    if problems:
        for p in problems:
            print(f"  [!] {p}")
        return 1
    print("  environment looks consistent")
    return 0


def parse_args(argv=None) -> Config:
    p = argparse.ArgumentParser(
        description="Nuclear segmentation + multi-channel positivity gating for "
                    "ImageXpress screens stored in OMERO.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    g = p.add_argument_group("OMERO connection")
    g.add_argument("--host", default="localhost")
    g.add_argument("--port", type=int, default=4064)
    g.add_argument("--user")
    g.add_argument("--password", help="prefer the OMERO_PASSWORD env var")
    g.add_argument("--insecure", action="store_true")
    g.add_argument("--group", default="-1")

    g = p.add_argument_group("Data selection")
    g.add_argument("--screen-id", type=int)
    g.add_argument("--plate-id", action="append", default=[], dest="plate_ids",
                   help="plate ID; repeatable AND comma-separated, e.g. "
                        "--plate-id 101,102 --plate-id 103")
    g.add_argument("--wells", help="comma separated, e.g. A1,A2,B3")
    g.add_argument("--max-fields", type=int)
    g.add_argument("--dapi-channel", help="channel index or name fragment")
    g.add_argument("--z-mode", choices=["mid", "mip"], default="mid")
    g.add_argument("--timepoint", type=int, default=0)
    g.add_argument("--server-tile-height", type=int, default=0,
                   help="read planes in strips of N rows instead of whole "
                        "planes; use 256 or 512 if the OMERO server throws "
                        "'Java heap space'")

    g = p.add_argument_group("Segmentation")
    g.add_argument("--model", choices=["stardist", "cellpose"], default="stardist")
    g.add_argument("--stardist-model", default="2D_versatile_fluo")
    g.add_argument("--cellpose-model", default="nuclei")
    g.add_argument("--prob-thresh", type=float)
    g.add_argument("--nms-thresh", type=float)
    g.add_argument("--diameter", type=float)
    g.add_argument("--n-tiles", type=int)
    g.add_argument("--no-gpu", action="store_true")
    g.add_argument("--gpu-memory-limit", type=int,
                   help="cap TensorFlow's GPU memory in MB, e.g. 4096")
    g.add_argument("--gc-every", type=int, default=50,
                   help="garbage-collect every N images (0 = never)")
    g.add_argument("--reload-model-every", type=int, default=0,
                   help="rebuild the segmentation model every N images; the "
                        "reliable cure for a slow GPU memory leak, e.g. 200")
    g.add_argument("--min-area", type=int, default=30)
    g.add_argument("--max-area", type=int, default=100000)
    g.add_argument("--min-solidity", type=float, default=0.0)
    g.add_argument("--keep-border", action="store_true")

    g = p.add_argument_group("Dead / dying cell filter (bright DAPI, optional)")
    g.add_argument("--filter-dead", action="store_true",
                   help="exclude hyper-bright (condensed/pyknotic) nuclei")
    g.add_argument("--dead-method", choices=["mad", "quantile", "otsu", "manual"],
                   default="mad")
    g.add_argument("--dead-stat", choices=["mean", "median", "integrated"],
                   default="mean", help="DAPI statistic used for the cutoff")
    g.add_argument("--dead-k", type=float, default=3.0,
                   help="robust SDs above the median (--dead-method mad)")
    g.add_argument("--dead-quantile", type=float, default=0.99)
    g.add_argument("--dead-threshold", type=float,
                   help="fixed cutoff for --dead-method manual")
    g.add_argument("--dead-scope", choices=["screen", "plate", "well", "image"],
                   default="plate")
    g.add_argument("--dead-require-small", action="store_true",
                   help="only call dead if the nucleus is ALSO small (condensed)")
    g.add_argument("--dead-area-quantile", type=float, default=0.25)
    g.add_argument("--dead-action", choices=["remove", "flag"], default="remove",
                   help="'flag' keeps the cells and only adds a 'dead' column")

    g = p.add_argument_group("Measurement / gating")
    g.add_argument("--background", choices=["none", "median"], default="median")
    g.add_argument("--intensity-stat", choices=["mean", "median", "integrated"],
                   default="mean")
    g.add_argument("--threshold-method",
                   choices=["otsu", "gmm", "quantile", "control", "manual"],
                   default="otsu")
    g.add_argument("--threshold-scope", choices=["screen", "plate", "well", "image"],
                   default="plate")
    g.add_argument("--quantile", type=float, default=0.99)
    g.add_argument("--control-wells", help="negative-control wells, e.g. A1,A2")
    g.add_argument("--control-sd", type=float, default=3.0)
    g.add_argument("--manual-thresholds", help="e.g. 'GFP=1200,RFP=800'")
    g.add_argument("--marker-channels", help="restrict gating to these channels")
    g.add_argument("--subset-analyses", action="store_true",
                   help="additionally analyse marker subsets as self-contained "
                        "panels, each summing to 100%% (useful when different "
                        "wells on a plate use different staining panels)")
    g.add_argument("--subset-sizes", default="1,2",
                   help="panel sizes to generate, e.g. '1,2' for every single "
                        "marker and every pair")
    g.add_argument("--panel", action="append", default=[], dest="subset_panels",
                   help="an explicit panel instead of all combinations, e.g. "
                        "--panel GFP+RFP --panel Cy5; repeatable")

    g = p.add_argument_group("Foci detection (gamma-H2AX / 53BP1, optional)")
    g.add_argument("--foci-channel", dest="foci_channels", action="append",
                   default=[], help="channel name fragment or index; repeatable "
                                    "(e.g. --foci-channel gH2AX)")
    g.add_argument("--foci-method", choices=["log", "tophat", "hmax"], default="log")
    g.add_argument("--foci-diameter-um", type=float, default=0.8,
                   help="expected physical focus diameter; with the pixel size "
                        "this sets the scales, top-hat radius and size filters")
    g.add_argument("--pixel-size-um", type=float,
                   help="override the pixel size from OMERO metadata")
    g.add_argument("--no-foci-auto-scale", action="store_true",
                   help="use the raw pixel parameters below instead of "
                        "deriving them from the pixel size")
    g.add_argument("--foci-upsample", type=int, default=0,
                   help="0 = auto (2-3x when foci are only 1-3 px wide)")
    g.add_argument("--no-foci-bg-clip", action="store_true",
                   help="don't sigma-clip the per-nucleus background estimate")
    g.add_argument("--save-foci-crops", type=int, default=0,
                   help="save N zoomed montages of nuclei with foci circled - "
                        "essential for checking foci at 20x")
    g.add_argument("--foci-tophat-radius", type=float, default=5.0,
                   help="background-flattening radius in px (> largest focus)")
    g.add_argument("--foci-min-sigma", type=float, default=1.0)
    g.add_argument("--foci-max-sigma", type=float, default=4.0)
    g.add_argument("--foci-num-sigma", type=int, default=5)
    g.add_argument("--foci-log-threshold", type=float, default=0.02,
                   help="LoG detector sensitivity - lower finds more foci")
    g.add_argument("--foci-k", type=float, default=3.0,
                   help="robust SDs above each nucleus's own background")
    g.add_argument("--foci-abs-threshold", type=float,
                   help="fixed top-hat cutoff instead of the per-nucleus one")
    g.add_argument("--foci-h", type=float, default=0.0, help="h-maxima depth")
    g.add_argument("--foci-min-area", type=int, default=3)
    g.add_argument("--foci-max-area", type=int, default=400)
    g.add_argument("--foci-min-count", type=int, default=3,
                   help="nuclei with >= N foci are called foci-positive")
    g.add_argument("--foci-as-marker", action="store_true",
                   help="use the foci count, not intensity, to call that channel "
                        "+/- in the single/double/triple classification")

    g = p.add_argument_group("OMERO key-value pairs (map annotations)")
    g.add_argument("--no-key-values", action="store_true",
                   help="don't fetch OMERO key-value pairs")
    g.add_argument("--kv-levels", default="screen,plate,well",
                   help="which levels to harvest, most specific wins; "
                        "add 'image' for per-field annotations")
    g.add_argument("--kv-prefix", default="kv_",
                   help="prefix for the added columns (keeps them selectable "
                        "and collision-free)")
    g.add_argument("--group-by",
                   help="key-value key to group summaries and plots by, "
                        "e.g. --group-by Compound")

    g = p.add_argument_group("Plots / output")
    g.add_argument("--per-nucleus-format",
                   choices=["csv", "csv.gz", "parquet", "xlsx", "none"],
                   default="csv.gz",
                   help="format for the per-nucleus table ('none' to skip)")
    g.add_argument("--save-foci-table", action="store_true",
                   help="also write one row per detected focus")
    g.add_argument("--export-explorer", action="store_true",
                   help="write gating_explorer.html next to the results, with "
                        "the per-nucleus data embedded, for interactive re-gating")
    g.add_argument("--explorer-no-embed", action="store_true",
                   help="write the explorer without inlining the data "
                        "(smaller file; drop the CSV in by hand)")
    g.add_argument("--cofactor", type=float, help="asinh cofactor (auto if unset)")
    g.add_argument("--max-points", type=int, default=50000)
    g.add_argument("--per-well-plots", action="store_true")
    g.add_argument("--save-overlays", type=int, default=0,
                   help="save N segmentation QC overlays per plate")
    g.add_argument("--overlay-line-width", type=int, default=2,
                   help="outline thickness in px for QC overlays and foci crops")
    g.add_argument("--outdir", type=Path, default=Path("./results"))
    g.add_argument("--pool-screen", action="store_true",
                   help="also produce a screen-level pooled analysis")
    g.add_argument("--upload-rois", action="store_true")
    g.add_argument("--attach-results", action="store_true")
    g.add_argument("--from-csv", type=Path,
                   help="skip OMERO; re-gate and re-plot an existing measurement CSV")
    g.add_argument("-v", "--verbose", action="store_true")
    g.add_argument("--check-env", action="store_true",
                   help="print package versions and GPU status, then exit")

    a = p.parse_args(argv)

    manual = {}
    if a.manual_thresholds:
        for part in a.manual_thresholds.split(","):
            k, v = part.split("=")
            manual[k.strip()] = float(v)

    plate_ids: List[int] = []
    for part in a.plate_ids:
        for piece in str(part).split(","):
            piece = piece.strip()
            if not piece:
                continue
            try:
                plate_ids.append(int(piece))
            except ValueError:
                p.error(f"--plate-id expects integers, got '{piece}'. Plate IDs "
                        f"are the numbers shown in OMERO.web, e.g. --plate-id 101,102")

    kv_prefix = a.kv_prefix or "kv_"
    if not a.kv_prefix:
        LOG.warning("Empty --kv-prefix is not supported (columns must stay "
                    "identifiable) - using 'kv_'")

    cfg = Config(
        host=a.host, port=a.port, user=a.user, password=a.password,
        secure=not a.insecure, group=a.group,
        screen_id=a.screen_id, plate_ids=plate_ids,
        wells=a.wells.split(",") if a.wells else None,
        max_fields=a.max_fields, dapi_channel=a.dapi_channel,
        z_mode=a.z_mode, timepoint=a.timepoint,
        server_tile_height=max(0, a.server_tile_height),
        model=a.model, stardist_model=a.stardist_model,
        cellpose_model=a.cellpose_model, prob_thresh=a.prob_thresh,
        nms_thresh=a.nms_thresh, diameter=a.diameter, n_tiles=a.n_tiles,
        gpu=not a.no_gpu, gpu_memory_limit=a.gpu_memory_limit,
        gc_every=max(0, a.gc_every), reload_model_every=max(0, a.reload_model_every), min_area=a.min_area, max_area=a.max_area,
        exclude_border=not a.keep_border, min_solidity=a.min_solidity,
        filter_dead=a.filter_dead, dead_method=a.dead_method,
        dead_stat=a.dead_stat, dead_k=a.dead_k, dead_quantile=a.dead_quantile,
        dead_threshold=a.dead_threshold, dead_scope=a.dead_scope,
        dead_require_small=a.dead_require_small,
        dead_area_quantile=a.dead_area_quantile, dead_action=a.dead_action,
        background=a.background, intensity_stat=a.intensity_stat,
        foci_channels=[c for part in a.foci_channels for c in part.split(",")],
        foci_method=a.foci_method, foci_tophat_radius=a.foci_tophat_radius,
        foci_min_sigma=a.foci_min_sigma, foci_max_sigma=a.foci_max_sigma,
        foci_num_sigma=a.foci_num_sigma, foci_log_threshold=a.foci_log_threshold,
        foci_k=a.foci_k, foci_abs_threshold=a.foci_abs_threshold, foci_h=a.foci_h,
        foci_min_area=a.foci_min_area, foci_max_area=a.foci_max_area,
        foci_min_count=a.foci_min_count, foci_as_marker=a.foci_as_marker,
        pixel_size_um=a.pixel_size_um, foci_diameter_um=a.foci_diameter_um,
        foci_auto_scale=not a.no_foci_auto_scale, foci_upsample=a.foci_upsample,
        foci_bg_clip=not a.no_foci_bg_clip, save_foci_crops=a.save_foci_crops,
        per_nucleus_format=a.per_nucleus_format,
        save_foci_table=a.save_foci_table,
        export_explorer=a.export_explorer, explorer_no_embed=a.explorer_no_embed,
        fetch_kv=not a.no_key_values,
        kv_levels=[s.strip().lower() for s in a.kv_levels.split(",") if s.strip()],
        kv_prefix=kv_prefix, group_by=a.group_by,
        threshold_method=a.threshold_method, threshold_scope=a.threshold_scope,
        quantile=a.quantile,
        control_wells=a.control_wells.split(",") if a.control_wells else [],
        control_sd=a.control_sd, manual_thresholds=manual,
        marker_channels=a.marker_channels.split(",") if a.marker_channels else None,
        subset_analyses=a.subset_analyses or bool(a.subset_panels),
        subset_sizes=[int(s) for s in str(a.subset_sizes).split(",") if s.strip()],
        subset_panels=a.subset_panels,
        cofactor=a.cofactor, max_points=a.max_points,
        per_well_plots=a.per_well_plots, save_overlays=a.save_overlays,
        overlay_line_width=max(1, a.overlay_line_width),
        outdir=a.outdir, upload_rois=a.upload_rois,
        attach_results=a.attach_results, from_csv=a.from_csv,
    )
    cfg._verbose = a.verbose          # type: ignore[attr-defined]
    cfg._pool_screen = a.pool_screen  # type: ignore[attr-defined]
    return cfg


def main(argv=None) -> int:
    if argv is None:
        argv = sys.argv[1:]
    if "--check-env" in argv:
        return check_environment()

    cfg = parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if getattr(cfg, "_verbose", False) else logging.INFO,
        format="%(asctime)s %(levelname)-7s %(message)s", datefmt="%H:%M:%S")
    cfg.outdir.mkdir(parents=True, exist_ok=True)
    with open(cfg.outdir / "run_config.json", "w") as fh:
        json.dump({k: (str(v) if isinstance(v, Path) else v)
                   for k, v in asdict(cfg).items() if k != "password"}, fh, indent=2)

    # ---- offline re-gating mode -------------------------------------------
    if cfg.from_csv:
        LOG.info("Re-gating %s", cfg.from_csv)
        df = pd.read_csv(cfg.from_csv)
        analyse(df, cfg, cfg.outdir, tag=cfg.from_csv.parent.name or "regated")
        LOG.info("Done -> %s", cfg.outdir)
        return 0

    if cfg.screen_id is None and not cfg.plate_ids:
        LOG.error("Provide --screen-id or --plate-id (or --from-csv)")
        return 2

    conn = connect(cfg)
    segmenter = Segmenter(cfg)
    all_plates: List[pd.DataFrame] = []
    try:
        screen_kv: Dict[str, str] = {}
        if cfg.fetch_kv and "screen" in cfg.kv_levels and cfg.screen_id is not None:
            screen = conn.getObject("Screen", cfg.screen_id)
            if screen is not None:
                screen_kv = map_annotations(screen)
                if screen_kv:
                    LOG.info("Screen-level key-value pairs: %s",
                             ", ".join(sorted(screen_kv)))

        n_targets = len(cfg.plate_ids) + (1 if cfg.screen_id is not None else 0)
        if (cfg.threshold_scope == "screen" and n_targets > 1
                and not getattr(cfg, "_pool_screen", False)):
            LOG.warning("--threshold-scope screen with several plates: each "
                        "plate's own folder is still gated on that plate's data "
                        "alone. Add --pool-screen for a genuinely common "
                        "threshold across plates (written to _pooled/).")

        for plate in iter_plates(conn, cfg):
            LOG.info("Plate %s (%s)", plate.getId(), plate.getName())
            pdir = cfg.outdir / safe(plate.getName())
            pdir.mkdir(parents=True, exist_ok=True)
            raw, foci, well_meta = process_plate(conn, plate, cfg, segmenter,
                                                 inherited_kv=screen_kv)
            if not well_meta.empty:
                well_meta.to_csv(pdir / "well_metadata.csv", index=False)
            if raw.empty:
                LOG.warning("No cells measured on plate %s", plate.getName())
                continue
            if cfg.save_foci_table:
                if foci.empty:
                    LOG.warning("No foci detected on plate %s", plate.getName())
                else:
                    fmt = cfg.per_nucleus_format
                    write_table(foci, pdir / "per_focus_measurements",
                                "csv.gz" if fmt == "none" else fmt)
            gated = analyse(raw, cfg, pdir, tag=plate.getName(),
                            conn=conn, omero_obj=plate)
            all_plates.append(gated)

        if getattr(cfg, "_pool_screen", False) and len(all_plates) > 1:
            pooled = pd.concat(all_plates, ignore_index=True)
            sdir = cfg.outdir / "_pooled"
            sdir.mkdir(parents=True, exist_ok=True)
            pooled = pooled.drop(columns=[c for c in pooled.columns
                                          if c.startswith("pos_")] +
                                 ["n_positive", "population", "combination", "dead"],
                                 errors="ignore")
            label = ("screen (pooled)" if cfg.screen_id is not None
                     else f"{len(all_plates)} plates (pooled)")
            analyse(pooled, cfg, sdir, tag=label)
        elif getattr(cfg, "_pool_screen", False):
            LOG.info("--pool-screen needs more than one plate - skipping")
    finally:
        conn.close()

    LOG.info("Done -> %s", cfg.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
