# NucleoGate

**Flow-cytometry-style analysis of high-content microscopy: nuclear
segmentation, multi-marker gating, foci counting and dose-response, straight
from OMERO.**

---

## Start here: the workbench

**`nucleoflow_workbench.html` — one file, nothing to install.** Double-click it,
drop in a results file, and you have the whole analysis in two tabs:

| Tab | What you do there |
|---|---|
| **1 · Gating** | filter out debris and dying cells, drag the marker gates, see the populations per well |
| **2 · Conditions & IC50** | group replicate wells, compare conditions, run two-way ANOVA, fit IC50 curves — per compound *and* per cell population |

Both tabs work on the **same cells**, so a gate you move on tab 1 changes the
IC50s on tab 2. That link is the point: splitting a dose-response by cell
population needs a positivity call for every cell, and tab 1 is where those come
from. Nothing leaves your computer — the browser reads the file locally.

If you are handing one thing to a student, hand them this. See
[`WORKBENCH_GUIDE.md`](WORKBENCH_GUIDE.md).

## What produces the data

**`omero_nuclei_flow.py`** takes ImageXpress high-content images stored in OMERO
as **Screens**, finds every nucleus with StarDist or Cellpose, measures its
brightness in every channel, decides which cells are positive for each marker,
counts γ-H2AX-style foci if you want it to, and writes the tables and plots the
workbench then reads.

In one sentence: **it turns a plate of images into the kind of table and plots
you'd get from a flow cytometer, except every event is a nucleus you can go back
and look at.**

```bash
conda env create -f environment.yml && conda activate omeroflow
python omero_nuclei_flow.py --check-env
python omero_nuclei_flow.py --host omero.example.org --user you \
    --screen-id 1234 --export-explorer --outdir ./results
```

`--export-explorer` writes a copy of the browser tool next to the results with
your data already inside it, so there is nothing to load by hand.

Full instructions, written for someone who has never opened a terminal, start at
[section 2](#2-opening-a-terminal).

## The standalone tools

The two halves of the workbench also exist on their own, for when you only need
one:

- **`gating_explorer.html`** — gating and populations only.
- **`condition_explorer.html`** — conditions, ANOVA and IC50 only. Use this if
  your file already carries the gates you want, or if you are working from a
  `per_well_summary.csv`.

---

**Files in this repository**

| File | What it's for |
|---|---|
| **`nucleoflow_workbench.html`** | **the main tool** — gating and analysis in one file, two tabs, one dataset |
| `WORKBENCH_GUIDE.md` | guide to the combined workflow |
| `omero_nuclei_flow.py` | the pipeline that turns OMERO images into the tables the workbench reads |
| `environment.yml` | pinned conda environment (StarDist / TensorFlow 2.10 / Python 3.10) |
| `environment-cellpose.yml` | alternative environment (Cellpose / PyTorch, no version ceiling) |
| `requirements.txt` | the same pins for pip, if you can't use conda |
| `gating_explorer.html` | the gating half on its own |
| `EXPLORER_GUIDE.md` | student-facing guide to the gating explorer |
| `condition_explorer.html` | the analysis half on its own |
| `CONDITION_EXPLORER_GUIDE.md` | guide to the condition explorer |
| `build_workbench.py` | rebuilds the workbench from the two standalone files |
| `WINDOWS_LATEST_TENSORFLOW.md` | short guide to running current TensorFlow on Windows via WSL2 |
| `README.md` | this file |

**Repository topics** — `omero` · `high-content-screening` · `microscopy` ·
`image-analysis` · `image-cytometry` · `cell-segmentation` · `stardist` ·
`cellpose` · `flow-cytometry` · `gating` · `high-throughput-screening` ·
`bioimage-analysis` · `python` · `imagexpress`

## Table of contents

**Part 1 — Getting it running (assumes zero programming experience)**
1. [What the script actually does](#1-what-the-script-actually-does)
2. [Opening a terminal](#2-opening-a-terminal)
3. [Installing everything](#3-installing-everything)
4. [Finding your Screen or Plate ID in OMERO](#4-finding-your-screen-or-plate-id-in-omero)
5. [Your first run (a small pilot)](#5-your-first-run-a-small-pilot)
6. [The three things to check before you trust anything](#6-the-three-things-to-check-before-you-trust-anything)
7. [The full run](#7-the-full-run)
8. [Reading the log as it runs](#8-reading-the-log-as-it-runs)
9. [When something goes wrong](#9-when-something-goes-wrong)

**Part 2 — What the numbers mean (for biologists)**
10. [Every output file, explained](#10-every-output-file-explained)
11. [The per-nucleus table, column by column](#11-the-per-nucleus-table-column-by-column)
    · [OMERO key-value pairs](#omero-key-value-pairs-your-experimental-metadata)
12. [How positivity is decided, and why it's the weak point](#12-how-positivity-is-decided-and-why-its-the-weak-point)
13. [Which intensity statistic to use](#13-which-intensity-statistic-to-use)
14. [The dead-cell filter](#14-the-dead-cell-filter)
15. [Foci counting](#15-foci-counting)
16. [Statistics: your n is not what you think](#16-statistics-your-n-is-not-what-you-think)
17. [Controls worth running](#17-controls-worth-running)

**Part 3 — How it works (for computer scientists and method-minded readers)**
18. [Architecture and dataflow](#18-architecture-and-dataflow)
19. [Segmentation](#19-segmentation)
20. [Measurement](#20-measurement)
21. [Thresholding algorithms](#21-thresholding-algorithms)
22. [Dead-cell detection](#22-dead-cell-detection)
23. [Foci detection](#23-foci-detection)
24. [Why asinh axes](#24-why-asinh-axes)
25. [Implementation notes and performance](#25-implementation-notes-and-performance)
26. [Extending the script](#26-extending-the-script)

**Part 4 — Reference**
27. [Complete flag reference](#27-complete-flag-reference)
28. [Recipes](#28-recipes)
29. [Known limitations](#29-known-limitations)
30. [Glossary (both directions)](#30-glossary-both-directions)

---
---

# Part 1 — Getting it running

*This part assumes you have never written a line of code. Everything you need to
type is in a grey box. Type it (or copy-paste it) and press Enter.*

## 1. What the script actually does

Imagine you've stained a 384-well plate with DAPI plus three antibodies, and the
ImageXpress has taken 4 images per well. That's 1,536 images and maybe half a
million cells. You want to know: what fraction of cells are positive for marker
A, for B, for both, for all three, and does that change with your treatment?

The script does this, in order:

1. **Finds the images.** It logs into OMERO and walks down through your Screen →
   Plates → Wells → fields → images. You never export anything by hand.
2. **Finds the nuclei.** It runs a deep-learning segmentation model on the DAPI
   channel, which draws an outline around every nucleus. Think of it as an
   automatic, very patient person tracing nuclei in ImageJ.
3. **Measures.** For each nucleus it records the average brightness inside that
   outline, *in every channel*, plus size and shape.
4. **Optionally throws out dead cells** whose DAPI is abnormally bright (dying
   cells have condensed, super-bright chromatin, and they're non-specifically
   bright in every other channel too — they'd fake being positive for everything).
5. **Optionally counts foci** — discrete bright dots inside the nucleus, like
   γ-H2AX or 53BP1 DNA-damage foci.
6. **Decides positive vs negative** for each marker by finding a cutoff in the
   distribution of brightnesses.
7. **Counts up** how many markers each cell is positive for: 0 = negative,
   1 = single positive, 2 = double positive, 3 = triple positive, and so on. The
   number of categories adapts automatically to how many channels your screen has.
8. **Writes out** a big table (one row per nucleus), summary tables (one row per
   well), and a folder of plots.

Nothing is modified in OMERO unless you explicitly ask for it.

## 2. Opening a terminal

The terminal is a window where you type commands instead of clicking. It looks
intimidating and is not.

- **macOS** — press `Cmd + Space`, type `Terminal`, press Enter.
- **Windows** — you'll install something called Miniforge in the next step, which
  adds a **"Miniforge Prompt"** to your Start menu. Use *that*, not the normal
  Command Prompt.
- **Linux** — `Ctrl + Alt + T`, or search for "Terminal".

Two things worth knowing:

- **You can paste.** `Cmd+V` on Mac, right-click on Windows, `Ctrl+Shift+V` on Linux.
- **A command spread over several lines** ends each line with a backslash `\`.
  That's just "the command continues on the next line". Copy the whole block at
  once. **On Windows this doesn't work** — I've given single-line versions where
  it matters.

## 3. Installing everything

The script needs Python plus a handful of scientific libraries. **Use conda, not
pip.** Two of the dependencies are genuinely awkward otherwise: `omero-py` needs
a C++ library called Ice, and the segmentation model needs a specific, fairly
old TensorFlow. Conda handles both; pip will fight you.

Everything is pinned in `environment.yml`, so you don't have to choose versions.

### Step 3a — install Miniforge

Miniforge is a minimal conda. Download the installer for your system from
<https://github.com/conda-forge/miniforge#download> and run it. Accept the
defaults. On Windows, tick the box that offers to add it to your Start menu.

Close and reopen your terminal afterwards.

### Step 3b — create the environment

Put `omero_nuclei_flow.py` and `environment.yml` in the same folder, `cd` into
it, and run:

```bash
conda env create -f environment.yml
conda activate omeroflow
```

`cd` means "change directory" — it's how you tell the terminal which folder to
work in. This takes several minutes.

**You have to run `conda activate omeroflow` every time you open a new
terminal.** If a command later fails with "command not found" or "no module
named", 90% of the time it's because you forgot this line.

### Step 3c — check it worked

```bash
python omero_nuclei_flow.py --check-env
```

This prints every package version, whether your GPU is visible, and — the useful
part — warns about specific known-bad version combinations:

```
python        3.10.14   (Windows AMD64)
numpy         1.23.5
tensorflow    2.10.1
protobuf      3.19.6
stardist      0.8.5
omero-py      5.19.1
...
TensorFlow GPUs visible: 1  ['/physical_device:GPU:0']

  environment looks consistent
```

If it reports problems, fix those before going further — a broken pin here shows
up later as a confusing crash or as segmentation silently running on the CPU at
a tenth the speed.

### Why the versions are pinned the way they are

This is worth understanding, because it explains why you cannot simply "install
the latest".

TensorFlow 2.10 was the last release that supported GPU on native Windows; from 2.11 onward you need WSL2 or a CPU-only install. StarDist runs
on TensorFlow, so on Windows with an NVIDIA GPU you are on TensorFlow 2.10 —
and TF 2.10 supports **Python 3.7–3.10 only**, which is where the Python 3.10
pin comes from. From there the rest follows:

| Package | Pin | Why |
|---|---|---|
| python | 3.10 | highest version TensorFlow 2.10 supports |
| tensorflow | 2.10.1 | last native-Windows GPU release |
| numpy | 1.23.5 | TF <2.11 breaks on numpy ≥1.24; numpy 2.x breaks StarDist |
| protobuf | 3.19.6 | TF <2.11 is binary-incompatible with protobuf ≥3.20 |
| cudatoolkit / cudnn | 11.2 / 8.1 | the versions TF 2.10 is built against — newer CUDA will not work |
| scikit-image | 0.21 | ≥0.19 needed for the region properties used here |
| zeroc-ice | 3.6.5 | what omero-py talks to the server with |
| setuptools | <81 | setuptools 82 (Feb 2026) **removed** `pkg_resources`, which StarDist still imports |

Break any one of these and you typically get a silent fallback to CPU, or an
`ImportError` deep inside TensorFlow that says nothing about versions.

The setuptools pin is the newest of these and catches people out: `pkg_resources`
was part of setuptools for twenty years and was removed in February 2026, so a
freshly created environment picks up a setuptools that StarDist cannot import.
The symptom is `ModuleNotFoundError: No module named 'pkg_resources'` raised from
`stardist/bioimageio_utils.py` — confusing, because nothing about your code or
StarDist changed. If you hit it in an environment created before this pin
existed, `pip install "setuptools<81"` fixes it in place.

Note that **the analysis script itself is not the constraint** — it runs fine on
numpy 2.x and current pandas/scikit-image. The ceiling comes entirely from
StarDist's TensorFlow dependency.

### Windows with an NVIDIA GPU

`environment.yml` already includes `cudatoolkit=11.2` and `cudnn=8.1.0`, so
after `conda env create` it should just work. Confirm with `--check-env` that
`TensorFlow GPUs visible: 1`. If it says 0:

- Check your NVIDIA driver is current (the driver is separate from CUDA).
- Make sure nothing has upgraded TensorFlow past 2.10 — `pip list | findstr tensorflow`.
- Confirm you have the Microsoft Visual C++ Redistributable installed.

### Escaping the Python 3.10 ceiling

Being pinned to an end-of-life TensorFlow to keep GPU support on Windows is not
a great long-term position. There are two ways out, and both are supported:

**Option A — use Cellpose instead of StarDist.** Cellpose runs on PyTorch, which
still ships current native-Windows GPU builds. `environment-cellpose.yml` gives
you a modern stack (Python 3.11, no TensorFlow at all):

```bash
conda env create -f environment-cellpose.yml
conda activate omeroflow-cp
python omero_nuclei_flow.py --model cellpose ...
```

Worth comparing on a few fields regardless — see
[section 19](#19-segmentation) for how the two differ on dense nuclei.

**Option B — run under WSL2.** Windows Subsystem for Linux gives you a real
Linux environment with current TensorFlow GPU support, so you can use recent
Python and StarDist together. More setup, but no version ceiling. Step-by-step
instructions, including the Keras 3 and VPN traps, are in
**`WINDOWS_LATEST_TENSORFLOW.md`**.

### macOS

TensorFlow's Windows constraint doesn't apply, but StarDist on Apple Silicon is
fiddly. Either use `environment-cellpose.yml` (recommended — PyTorch has good
Apple Silicon support), or follow StarDist's own Apple Silicon instructions at
<https://github.com/stardist/stardist>, which use `tensorflow-macos` and
`tensorflow-metal`.

### If conda can't solve `zeroc-ice` on your platform

This is the one dependency that occasionally has no build for a given
Python/OS combination. Glencoe Software publishes prebuilt Ice wheels, and the
OMERO docs list the current URLs per platform and Python version:
<https://omero.readthedocs.io/en/stable/developers/Python.html>. Install the
wheel first, then `pip install omero-py`.

## 4. Finding your Screen or Plate ID in OMERO

Open OMERO.web in your browser and click on your Screen in the left-hand tree.
Either:

- look at the right-hand panel — it shows something like `Screen ID: 1234`, or
- look at the URL, which will contain `?show=screen-1234`.

That number is what goes after `--screen-id`.

**You do not have to run a whole screen.** Clicking a Plate instead of a Screen
gives you a Plate ID the same way, and everything in this document works
identically with `--plate-id`:

```bash
python omero_nuclei_flow.py --host ... --user ... --plate-id 5678 --outdir ./one_plate
```

Repeat the flag for several plates: `--plate-id 5678 --plate-id 5679`. Mixing is
allowed too — `--screen-id 1234 --plate-id 9999` processes the screen's plates
plus that extra one. A plate ID that doesn't exist is reported and skipped rather
than aborting the run.

Plates run individually still inherit the key-value annotations of the screen
they belong to, so you lose no metadata by working plate by plate. Results are
written to one folder per plate either way, so a screen run and four separate
plate runs produce the same files.

**One caveat when comparing plates.** Each plate's own results folder is always
gated on that plate's data alone — that's what `--threshold-scope plate` (the
default) means, and it's what you want when plates were stained and imaged
separately. If you need a *single* threshold applied across several plates so
the percentages are directly comparable, add `--pool-screen`: the `_pooled/`
folder re-gates all the plates together. Asking for `--threshold-scope screen`
without `--pool-screen` doesn't do this, and the script now warns you about it.

Which to use:

| Situation | Use |
|---|---|
| Whole experiment, one go | `--screen-id` |
| A named set of plates | `--plate-id 101,102,103` |
| Plates arriving over several days | `--plate-id`, one run each as they land |
| Re-running one plate that failed QC | `--plate-id` |
| Plate not attached to any screen | `--plate-id` (the only option) |
| You want to use several machines or GPUs | `--plate-id`, one process per plate |

You also need your OMERO **server address** (e.g. `omero.myinstitute.org`) and
your **username**. Ask whoever runs your OMERO if you don't know the address.

## 5. Your first run (a small pilot)

**Do not start by running a whole screen.** Run two wells first, look at the
pictures, then scale up. A pilot takes two minutes; a bad full run wastes a day.

First, put your password into the terminal's memory so it doesn't sit in your
command history:

```bash
export OMERO_PASSWORD='your-password-here'
```

Windows (Miniforge Prompt): `set OMERO_PASSWORD=your-password-here`

(If you skip this, the script will just ask you for the password when it starts,
and the typing will be invisible. That's normal.)

Now the pilot:

```bash
python omero_nuclei_flow.py \
  --host omero.myinstitute.org \
  --user your-username \
  --screen-id 1234 \
  --wells A1,A2 \
  --max-fields 2 \
  --save-overlays 4 \
  --outdir ./pilot
```

Windows single-line:

```
python omero_nuclei_flow.py --host omero.myinstitute.org --user your-username --screen-id 1234 --wells A1,A2 --max-fields 2 --save-overlays 4 --outdir ./pilot
```

What each line does:

| Line | Meaning |
|---|---|
| `--host` | your OMERO server address |
| `--user` | your OMERO username |
| `--screen-id 1234` | which screen to analyse |
| `--wells A1,A2` | only these two wells |
| `--max-fields 2` | only the first 2 images per well |
| `--save-overlays 4` | save 4 pictures showing what it segmented |
| `--outdir ./pilot` | put results in a folder called `pilot` |

## 6. The three things to check before you trust anything

Open the `pilot` folder. Inside is a folder named after your plate. Look at these
three things, in this order. **This is the most important section in this
document.**

### (a) `plots/overlays/*.png` — did it find the nuclei?

Red outlines drawn on the DAPI image. You're asking: is every nucleus outlined
exactly once? Common failures and their fixes:

| What you see | What it means | Fix |
|---|---|---|
| Two touching nuclei share one outline | under-segmentation | `--nms-thresh 0.5` (allows more overlap), or try `--model cellpose` |
| One nucleus split into pieces | over-segmentation | `--prob-thresh 0.6` (be stricter) |
| Faint nuclei missed entirely | detector too strict | `--prob-thresh 0.3` (lower = more objects) |
| Debris outlined as cells | junk passing the size filter | raise `--min-area` |
| Nothing at all found | wrong channel used as DAPI | `--dapi-channel DAPI` or `--dapi-channel 0` |

If foci counting is on, detected foci appear in green inside the outlines.

Outlines are 2 px by default. A 2048x2048 field rendered into a ~1000 px PNG
halves everything, so a 1 px outline all but disappears — raise it with
`--overlay-line-width 4` if you're checking segmentation on a large screen, or
drop to 1 if thick outlines hide the nuclear edge you're trying to judge. The
same setting scales the contours and foci circles in the zoomed crop montages.

### (b) `plots/histograms_gated.png` — is the cutoff in the right place?

One panel per marker, showing the distribution of nuclear brightness with a red
dashed line where the positive/negative cutoff landed. You want the red line
sitting in the **valley between two humps**. If your data has one hump and the
line is stuck in the middle of it, the automatic threshold is meaningless —
jump to [section 12](#12-how-positivity-is-decided-and-why-its-the-weak-point).

### (c) `plots/dead_cell_gate.png` — only if you used `--filter-dead`

Shows which nuclei got thrown out as dead/dying. See
[section 14](#14-the-dead-cell-filter).

## 7. The full run

Once the pilot looks right, run everything. A realistic full command:

```bash
python omero_nuclei_flow.py \
  --host omero.myinstitute.org \
  --user your-username \
  --screen-id 1234 \
  --dapi-channel DAPI \
  --threshold-method control --control-wells A1,A2 \
  --filter-dead --dead-require-small \
  --foci-channel gH2AX --save-foci-table \
  --per-nucleus-format parquet \
  --group-by Compound \
  --save-overlays 5 \
  --pool-screen \
  --outdir ./results
```

This says: use the DAPI channel by name; set positive/negative cutoffs from your
negative-control wells rather than guessing from the data; drop dead cells that
are both bright and condensed; count γ-H2AX foci and save every individual focus;
write the big table in the fast `parquet` format; also summarise by the
`Compound` key-value annotation; save some QC overlays; and produce a screen-wide
pooled analysis on top of the per-plate ones.

Expect roughly **1–3 seconds per image on a GPU** and **10–40 seconds per image
on a CPU**. A 384-well plate at 4 fields per well is 1,536 images: about an hour
on a GPU, most of a day on a laptop. Start it and go and do something else.

## 8. Reading the log as it runs

The script prints its progress. A normal run looks like:

```
10:04:12 INFO    Connected to OMERO as jdoe
10:04:13 INFO    Screen 1234: Cisplatin dose response
10:04:13 INFO    Plate 5678 (Plate_01)
10:04:19 INFO      A1 field 0 (image 90123): 412 nuclei
10:04:24 INFO      A1 field 1 (image 90124): 388 nuclei
...
10:52:01 INFO    [Plate_01] markers: gH2AX, RFP, Cy5
10:52:03 INFO    [Plate_01] bright-DAPI dead/dying: 4021 / 61344 cells (6.6%) - action: remove
10:52:09 INFO    [Plate_01] gH2AX foci: mean 7.34/nucleus, 68.2% with >= 3
10:52:11 INFO    Wrote per_nucleus_measurements.parquet (57323 rows)
10:52:12 INFO    [Plate_01] 57323 cells | single: 41.2% | negative: 33.8% | double: 19.7% | triple: 5.3%
```

Things worth reacting to:

- **Nuclei per field wildly different between wells** — could be biology (cell
  death) or could be focus/staining failure. Check the overlays for the low wells.
- **A very high dead fraction** (>20%) — either genuine toxicity, or your
  `--dead-k` is too aggressive. Look at `dead_cell_gate.png`.
- **`failed on image ... `** — one image errored; the run continues. If it's
  only a few images out of thousands, it's usually a corrupt or blank field.

## 9. When something goes wrong

| Message | What it means | What to do |
|---|---|---|
| `command not found: python` | environment not active | `conda activate omeroflow` |
| `No module named 'stardist'` | installed into the wrong environment | `conda activate omeroflow`, then re-check with `--check-env` |
| Segmentation is slow and `--check-env` says `GPUs visible: 0` | TensorFlow can't see the GPU | on Windows, TF must be exactly 2.10.x with CUDA 11.2 / cuDNN 8.1 — re-create the environment from `environment.yml` |
| `ImportError` from deep inside TensorFlow | a pin got broken, usually numpy or protobuf | `--check-env` names the offender |
| `ModuleNotFoundError: No module named 'pkg_resources'` (from inside `stardist`) | setuptools 82 removed `pkg_resources`; StarDist 0.8.5 still imports it | `pip install "setuptools<81"` — nothing else needs reinstalling. Pinned in `environment.yml`. |
| `Java heap space` / `omero.InternalException` while reading planes | the **server** ran out of memory, not your machine | `--server-tile-height 256`, and ask your admin to raise `omero.jvmcfg.percent.blitz`. See below. |
| `Could not connect to OMERO at ...` | wrong address, or you're off the network | check the host, connect to your institution's VPN, try `--port 4064` |
| `Screen 1234 not found (check group/permissions)` | wrong ID, or the data belongs to another group | double-check the ID in OMERO.web; ask the owner to share it |
| `Channel matching 'DAPI' not found in [...]` | your channels are named something else | the error prints the real names — use one of those, or an index like `--dapi-channel 0` |
| `Foci channel 'gH2AX' not found in [...]` | same, for the foci channel | use the printed name |
| `CUDA out of memory` on the very first image | the field is too big for the card | `--n-tiles 2` (or 4), or `--no-gpu` |
| `Dst tensor is not initialized` / OOM **after a while**, even with tiling | this is GPU memory *accumulating* across predict calls, not one image being too big — more tiling won't help | `--reload-model-every 200`. See below. |
| Everything is called positive, or nothing is | the automatic threshold failed | see [section 12](#12-how-positivity-is-decided-and-why-its-the-weak-point) — use `--threshold-method control` or `manual` |
| `every cell was flagged as dead` | `--dead-k` far too low | raise it, or drop `--filter-dead` and inspect `dead_cell_gate.png` first |
| It's unbearably slow | running on CPU | use a GPU machine, or reduce with `--max-fields 2` |
| `Cannot write parquet ... falling back to csv.gz` | `pyarrow` missing | harmless; or `conda install -c conda-forge pyarrow` |

### `Java heap space` — the OMERO *server* running out of memory

A traceback ending in

```
omero.InternalException ... (java.lang.OutOfMemoryError): Java heap space
    at loci.formats.in.CellWorxReader.openBytes(...)
    at ome.io.bioformats.BfPixelsWrapper.getWholePlane(...)
```

is **the server**, not your machine. Nothing you change locally about GPUs or
Python fixes it directly. Two things are going on:

1. `BfPixelsWrapper` in the stack means the pixel data is being read **through
   Bio-Formats from your original ImageXpress files** on every request, rather
   than from OMERO's own pixel store. That happens with in-place imports.
   Each read re-opens the file and allocates a whole plane in the Java heap.
2. OMERO's Blitz JVM defaults to **15% of system memory**, which is often only
   a few hundred MB and is easily exhausted by 2048x2048 planes.

**What the script now does automatically**

- Channel metadata is read with `getChannels(noRE=True)`. Plain `getChannels()`
  spins up the server-side *rendering engine*, which loads pixel data purely to
  compute display defaults — a whole-plane allocation for information already
  present in the channel metadata. This is a common cause of exactly this error.
- Channel names are cached **per plate** instead of queried per image, removing
  thousands of round trips on a 384-well plate.
- The rendering engine is released after every image.
- On a server heap error it retries that image in horizontal strips, and keeps
  using strips for the rest of the plate.

**Your quickest lever**

```bash
--server-tile-height 256
```

This reads each plane in 256-row strips via `getTile()` rather than one
`getPlane()` call, so the server allocates a small buffer per request instead of
a whole plane. The reassembled image is byte-identical; it is purely a transfer
strategy. Slightly more round trips, much lower server memory.

**Ask your OMERO admin for more heap.** The real fix is server-side:

```bash
omero config set omero.jvmcfg.percent.blitz 40   # default is 15
# or explicitly:
omero config set omero.jvmcfg.strategy.blitz manual
omero config set omero.jvmcfg.heap_size.blitz 8000   # MB
omero admin restart
```

`omero admin jvmcfg` shows the settings that will apply.

**The structural fix.** If the plates were imported in place, the server re-reads
your original files through Bio-Formats forever. Importing them normally, so
OMERO generates its own pixel data, removes that whole code path and makes reads
dramatically cheaper. Worth raising with whoever manages the imports if this
recurs across projects.

**Don't run plates in parallel against a struggling server.** Elsewhere this
guide suggests parallel per-plate processes for speed — that advice assumes the
server has headroom. If you are hitting Java heap errors, parallel runs multiply
the pressure. Run them sequentially until the heap is raised.

### Running out of GPU memory partway through a screen

If segmentation works for hundreds of images and then fails with
`InternalError: Failed copying input tensor ... Dst tensor is not initialized`,
that is TensorFlow's obscure way of saying **the GPU ran out of memory**. The
telling detail is *when*: a field that is genuinely too large fails on the first
image. Failing later means memory is accumulating across `predict()` calls, and
increasing `--n-tiles` will not help — tiling caps the peak memory of one
forward pass, not the slow creep across thousands of them.

The script defends against this in three ways, automatically:

1. **GPU memory growth is enabled** before the model loads. TensorFlow otherwise
   reserves nearly all VRAM up front and then fragments it.
2. **Periodic garbage collection**, every `--gc-every` images (default 50).
3. **Recovery instead of a lost image.** On an out-of-memory error it frees
   memory and retries with double the tiling; if that still fails it segments
   that one field on the CPU and carries on. You lose seconds, not the run.

If it still happens, the definitive fix is to rebuild the model periodically:

```bash
--reload-model-every 200
```

This clears the Keras session and reloads StarDist every 200 images, which
resets any accumulated state. It costs a few seconds each time — negligible
against a screen. `--gpu-memory-limit 4096` (MB) is also worth trying: capping
TensorFlow rather than letting it take everything often behaves better on a card
that is also driving your display.

Two structural alternatives, if you'd rather not fight it:

- **Run plate by plate as separate processes** (`--plate-id`). Memory is
  reclaimed by the OS when each process exits, so leaks cannot accumulate across
  a whole screen.
- **Use `--model cellpose`.** PyTorch's allocator handles long-running inference
  loops better than TensorFlow 2.10's, and you are not on an end-of-life stack.

**A note on repeating an analysis.** Re-segmenting a whole screen just to change
a cutoff is a waste of a day. Once the per-nucleus table exists you can re-gate
and re-plot it in seconds, entirely offline:

```bash
python omero_nuclei_flow.py \
  --from-csv results/Plate_01/per_nucleus_measurements.csv.gz \
  --threshold-method gmm \
  --filter-dead \
  --outdir ./results_regated
```

This works for anything computed *from* the table — thresholds, the dead filter,
foci-vs-intensity marker choice, plots. It cannot redo anything that needs the
pixels, i.e. segmentation and foci *detection* (already-counted foci are reused).

---
---

# Part 2 — What the numbers mean

## 10. Every output file, explained

```
results/
├── run_config.json                     every setting used, for your methods section
├── Plate_01/
│   ├── per_nucleus_measurements.csv.gz one row per nucleus  ← the main output
│   ├── per_focus_measurements.csv.gz   one row per focus    (--save-foci-table)
│   ├── excluded_dead_cells.csv.gz      the nuclei that were removed (--filter-dead)
│   ├── per_well_summary.csv            one row per well     ← use this for stats
│   ├── panels/                         self-contained panel analyses (--subset-analyses)
│   │   ├── panel_index.csv             all panels in one table
│   │   └── GFP_RFP/                    per-panel summaries + plots, summing to 100%
│   ├── per_<key>_summary.csv           one row per condition (--group-by)
│   ├── well_metadata.csv               the OMERO key-value pairs, per well
│   ├── per_image_summary.csv           one row per field    (QC)
│   ├── thresholds.json                 the cutoffs actually applied
│   └── plots/
│       ├── overlays/                   segmentation QC images
│       │   └── foci_crops_*.png        zoomed nuclei, foci circled (--save-foci-crops)
│       ├── dead_cell_gate.png          what the dead filter removed
│       ├── histograms_gated.png        distributions + cutoffs
│       ├── biaxial_A_vs_B.png          flow-style quadrant plots
│       ├── foci_gH2AX.png              foci counts and burden
│       ├── populations_by_<key>.png    the same, grouped by condition (--group-by)
│       ├── populations_by_well.png     stacked bars, % composition per well
│       ├── populations_by_well_counts.png   stacked bars, absolute cell numbers
│       ├── combinations_by_well.png    stacked bars: WHICH markers, %
│       ├── combinations_by_well_counts.png  the same, absolute cell numbers
│       ├── marker_combinations.png     % of each combination, plate-wide
│       └── platemap_*.png              plate-shaped heatmaps
└── _pooled/                            the same, pooled across plates (--pool-screen)
```

**Which file do I actually open?** `per_well_summary.csv` — in Excel, Prism or R.
That's your data for figures and statistics. The per-nucleus table is for when
you want to look at distributions or re-gate.

## 11. The per-nucleus table, column by column

One row per nucleus. Groups of columns:

**Where the cell came from**

| Column | Meaning |
|---|---|
| `plate`, `well` | plate name, well label (`A1`) |
| `well_row`, `well_col` | 0-based row/column, used to draw plate maps |
| `field` | which image within the well |
| `image_id`, `image_name` | the OMERO image — paste the ID into OMERO.web to see it |
| `label` | the nucleus's ID **within that image** (not unique across the plate) |
| `dapi_channel` | which channel was used for segmentation |

**Shape**

| Column | Meaning |
|---|---|
| `area` | nuclear area in pixels |
| `y`, `x` | centre of the nucleus in the image |
| `eccentricity` | 0 = circle, →1 = elongated |
| `solidity` | area ÷ convex-hull area; <0.9 suggests a ragged or merged object |
| `perimeter` | outline length |

**Intensity — one set per channel**

| Column | Meaning |
|---|---|
| `mean_<channel>` | average pixel value inside the nucleus (after background subtraction) |
| `median_<channel>` | the middle pixel value; less swayed by a couple of bright specks |
| `integrated_<channel>` | sum of all pixels = mean × area, i.e. *total* amount |
| `bg_<channel>` | the background value that was subtracted for that image |

**Foci — one set per foci channel** (see [section 15](#15-foci-counting))

**Calls**

| Column | Meaning |
|---|---|
| `dead` | flagged as dead/dying (only present with `--filter-dead`) |
| `pos_<marker>` | True/False positive for that marker |
| `n_positive` | how many markers this cell is positive for |
| `population` | how *many* markers: `negative` / `single` / `double` / `triple` / … |
| `combination` | *which* markers: `GFP+RFP+`, `GFP+Cy5+`, `all-negative`, … |
| `foci_pos_<channel>` | has at least `--foci-min-count` foci |

### Interactive re-gating in the browser

`gating_explorer.html` is a self-contained web page for post-hoc analysis — made
for handing to students. Open it, drop in a `per_nucleus_measurements.csv`, and
you can click wells in or out of the analysis, switch markers on and off, and
drag the gates, with every plot, percentage and table recomputing live. No
install, no server, and the data stays on the machine.

```bash
--export-explorer          # writes gating_explorer.html next to the results,
                           # with the per-nucleus data already inside it
--explorer-no-embed        # smaller file; drop the CSV in by hand
```

With `--export-explorer` the student gets one file to double-click. Embedding is
skipped automatically above ~80 MB of CSV.

It also has an upstream **cell filter**: window gates on nuclear intensity, area,
solidity, eccentricity and perimeter, applied before any marker gating — the
FSC/SSC pre-gate idea. Alongside it is a morphology-versus-intensity density plot
with the kept cells in colour and the excluded ones in grey, which is the quickest
way to see whether debris and clumps are inflating your positives. On a test plate
seeded with 8% debris and 6% clumps, filtering on area and solidity took the
triple-positive fraction from 2.0% to 0.0%.

Three things worth knowing:

- **It opens showing your numbers.** If the file has `pos_<marker>` columns, the
  explorer recovers the thresholds that produced them, so the first view
  reproduces the pipeline's report exactly. Verified against this pipeline: the
  per-well percentages match to the digit before anything is touched.
- **Switching the intensity statistic keeps your work.** Shape windows are
  statistic-independent and carry over untouched; the nuclear-stain window and
  any hand-set gates are remapped by percentile, so a gate at 38.9% positive on
  `mean` lands at 38.6% on `integrated` despite a ~150x change of scale. A notice
  says what was kept and what was rescaled.
- **Switching a marker off re-derives the populations** over the remaining ones,
  always summing to 100% — the same logic as `--subset-analyses`, done live.
  That's the quickest way to handle a plate with mixed staining panels.
- **The flow-cytometry plots are live, and sit beside the gates.** Every marker
  pair gets a density biaxial with the gate crosshair and quadrant percentages,
  redrawn as the gate is dragged — drag on the left, watch the quadrants move on
  the right. The density cloud is cached and only the overlay repaints, and
  quadrant counts come from a per-well histogram over marker bitmasks rather
  than a fresh pass over the cells — so a drag frame costs ~17 ms on a
  quarter-million nuclei instead of ~79 ms.
- **Export summary** writes a per-well CSV containing only the included wells —
  excluded ones are left out rather than flagged, so a column of zeros can't be
  averaged over by mistake. Trailing comment lines record the gates used, the
  markers included, the wells dropped and the cell count, so the analysis
  documents itself. Read it with `pd.read_csv(path, comment="#")`.

Performance: a 96-well plate at 250,000 nuclei parses in about a second and
holds ~60 fps while a gate is dragged, with the biaxial plots, histograms,
population bars and quadrant statistics all updating. Roughly a million nuclei
is the practical ceiling.

`EXPLORER_GUIDE.md` is written for students, and includes a section on the ways
interactive gating invites you to fool yourself.

### One file, two tabs: the workbench

`nucleoflow_workbench.html` contains both explorers as tabs over a **single
shared dataset**. It exists because splitting IC50s by cell population needs a
positivity call per cell, and on its own the condition explorer can only use
whatever `pos_*` columns the pipeline stored — with no way to change them.

In the workbench, tab 1 creates those columns from live gates. Drag the GFP gate
from 282 to 2000 and the marker goes from 26.8% to 13.2% positive, and the
per-population IC50s on tab 2 change with it. Press **Apply gates →** to push
the current calls through.

A `per_well_summary.csv` still loads; the gating tab simply disables itself,
since a summary cannot be re-gated. Both standalone files remain unchanged for
when you only need one half.

See `WORKBENCH_GUIDE.md`.

### Comparing replicate wells: the Condition Explorer

`condition_explorer.html` is a second browser tool for the step after gating:
grouping replicate wells into conditions and testing between them. Load a
`per_well_summary.csv`, the Gating Explorer's export, or a per-nucleus file;
group wells automatically from an OMERO key–value column or by painting them on
the plate; then read the result.

It is built around one rule: **n is the number of wells, never the number of
cells.** Pooling cells across a condition and running a t-test is
pseudo-replication and produces meaningless p-values.

- **SuperPlot** — every well drawn as a dot, with mean and 95% CI. Bars hide the
  spread that matters at n = 3–6.
- **Estimation plot** — each condition minus the reference with a bootstrap 95%
  CI, on the measured scale, so you see how big the difference is rather than
  only that it isn't zero.
- **Tests** — Welch's t / Welch's ANOVA by default (no equal-variance
  assumption), with Student's and rank-based alternatives, Hedges' *g* with CI,
  and Holm or Benjamini–Hochberg correction across all pairs.
- **Dose–response / IC50** for multi-compound plates, with the compound and
  concentration roles labelled explicitly and a concentration-unit selector
  (guessed from the column name, applied to the axis, tables and export): a
  four-parameter logistic fit per compound with 95% CI, Hill slope, R², log-axis
  curve plot, optional
  normalisation to vehicle or % inhibition, constrained fits, and a selectivity
  index against a second readout. IC50s that are extrapolated beyond the tested
  doses, poorly fitted, or resting on an implausible Hill slope are flagged
  rather than quietly reported.
- **Assay quality** — Z′-factor, SSMD, signal window, signal/background and CV
  per control, from a designated positive and negative control condition.
- **Curve comparison** — an extra sum-of-squares F test asking whether two
  compounds share an IC50, which is the correct test; overlapping confidence
  intervals are not one.
- **IC50 per cell population within a well** — split cells on any positivity
  call, `dead`, or a population class (needs a per-nucleus file), fit a curve to
  each, and test whether the populations respond at different doses. A well
  average can hide a tenfold difference between the populations inside it.
- **Two-way ANOVA** for crossed designs (genotype × dose, siRNA × treatment),
  set up either from metadata columns or entirely by hand — add, rename and
  delete levels on each factor, then paint wells:
  both main effects plus the interaction, with an interaction plot, a cell-means
  grid, and Type II or Type III sums of squares computed by nested model
  comparison so unbalanced plates are handled correctly, plus post-hoc
  comparisons using the pooled error term — simple effects when the interaction
  is significant, marginal means when it isn't. Empty cells are
  detected and the interaction dropped rather than silently mis-estimated.
- **Warnings** when n < 3, when n = 3 throughout, when wells are unassigned, and
  when many comparisons run uncorrected.

Every test is validated against `scipy.stats` — Welch, Student, Mann–Whitney,
ANOVA, Kruskal–Wallis, Hedges' *g*, CIs, Holm and BH — agreeing to within
3.6 × 10⁻¹⁴. The two-way ANOVA is validated separately against both the textbook
formulas and an independent numpy least-squares reference, across balanced and
unbalanced designs, to within 2 × 10⁻¹⁴. The curve fitter is validated against
`scipy.optimize.curve_fit` on six dose-response shapes — IC50, Hill, R² and the
confidence bounds all agreeing to within 1.6 × 10⁻⁸. The curve-comparison F test
and the Z′/SSMD metrics are validated the same way against scipy and numpy
references.

See `CONDITION_EXPLORER_GUIDE.md`.

### Analysing marker subsets as self-contained panels

With three markers, the plate-wide classification treats every cell as one of
negative / single / double / triple **out of three**. That is wrong for a well
stained with only two of those markers: a GFP+RFP+ cell that was never stained
for Cy5 gets called "double" out of a possible three, when for its actual panel
it is the maximum. This matters whenever different wells on one plate use
different staining panels.

`--subset-analyses` re-runs the classification for each marker subset
separately, so each is a **closed universe**:

```bash
--subset-analyses                      # every single marker and every pair
--subset-analyses --subset-sizes 2     # pairs only
--panel GFP+RFP --panel Cy5            # only the panels you actually used
```

For each panel you get a full set of outputs under `panels/<panel>/`:

- 1-marker panel: `pct_negative + pct_single = 100%`
- 2-marker panel: `pct_negative + pct_single + pct_double = 100%`

**The full three-marker analysis is already there** — it's the top-level output
in the plate folder (`per_well_summary.csv`, `plots/`), which is exactly the
negative/single/double/triple classification over all three markers, summing to
100%. It isn't duplicated into `panels/` because that would mean writing every
plot twice. `--subset-sizes` therefore only generates panels *smaller* than your
full marker set; asking for `--subset-sizes 1,2,3` with three markers gives the
same result as `1,2`.

Each panel folder contains its own `per_well_summary.csv`,
`per_image_summary.csv`, stacked bars (% and counts), plate heatmaps (% and
counts) and — for a pair — the biaxial quadrant plot for those two markers.
`panels/panel_index.csv` collects **every level, including the full set**, in
one table, with a `results_in` column pointing at where each one's detailed
output lives:

```
      panel  n_markers  results_in      pct_negative  pct_single  pct_double  pct_triple
        GFP          1  panels/GFP            54.513      45.487
    GFP+RFP          2  panels/GFP_RFP        35.737      48.462      15.800
    GFP+Cy5          2  panels/GFP_Cy5        40.800      47.800      11.400
GFP+RFP+Cy5          3  .                     26.650      45.413      24.050       3.888
```

The `.` means the plate folder itself. Every row sums to 100% across the classes
that apply to it. Blank cells are correct, not missing data: a one-marker panel
has no double-positive class.

**What is and isn't recomputed.** The per-marker thresholds are *unchanged* —
`pos_GFP` means the same thing in every panel, so a cell's positivity never
depends on which panel you look at. Only the population classification and its
percentages are recalculated over the chosen subset. That is what makes the
panels comparable to each other and to the full analysis.

**Practical note for mixed-panel plates.** The subsets are computed over all
wells, because the script has no way to know which wells got which stain. To
restrict a panel to the wells that actually used it, either run those wells
separately (`--wells A1,A2,...`), or filter `panels/<panel>/per_well_summary.csv`
afterwards — if the panel is recorded as an OMERO key-value pair, it will be
sitting right there as a `kv_` column.

### Percentages and counts are both plotted

Every population breakdown comes in two forms, because they answer different
questions and the difference is often the point:

- **Percentage plots** (`populations_by_well.png`,
  `platemap_pct_<population>.png`) show *composition* — what fraction of the
  surviving cells are double positive.
- **Count plots** (`populations_by_well_counts.png`, `platemap_n_cells.png`,
  `platemap_n_<population>.png`) show *absolute cell numbers*, so the height of
  a stacked bar is that well's total cell count.

A well can look excellent by percentage and terrible by count. If a treatment
kills 80% of the cells and the survivors are enriched for your marker, the
percentage rises impressively while the count collapses — that is usually
selection, not induction. Reading the two side by side is the cheapest way to
catch it.

`platemap_n_cells.png` is also the fastest plate-level QC you have: wells that
lost cells, an edge effect, a dispensing failure or a row of out-of-focus fields
all show up immediately as a spatial pattern in cell number. Count heatmaps are
annotated with the value in each well (up to 96 wells) and use a different
colour map from the percentage maps so the two are not confused at a glance.

The per-well and per-image summary CSVs carry both: `pct_double` alongside
`n_double`, `pct_GFP_positive` alongside `n_GFP_positive`, and `n_cells`.

### OMERO key-value pairs: your experimental metadata

Anything annotated in OMERO as key-value pairs — compound, concentration, cell
line, siRNA, timepoint, whatever your lab records — is harvested automatically
and **added as columns to every exported table**. This is on by default; disable
it with `--no-key-values`.

You therefore don't need a separate plate-map spreadsheet to make sense of the
results. `per_well_summary.csv` arrives looking like this:

```
well  kv_Compound  kv_Dose_uM  kv_Cell_line  n_cells  pct_GFP_positive  pct_double
A1    DMSO                  0  U2OS             2000             15.10        4.05
B2    Cisplatin            20  U2OS             2000             76.65       22.15
```

**How it works.** Annotations are collected from the screen, the plate and the
well, in that order, with more specific levels overriding more general ones — so
a cell line recorded once on the plate applies to every well, while a compound
recorded per well stays per well. Add `image` to `--kv-levels` if you annotate
individual fields.

**Column naming.** Keys are prefixed (`--kv-prefix`, default `kv_`) and cleaned
into safe names: `Dose (uM)` becomes `kv_Dose_uM`. The prefix means metadata
never collides with a measurement column and you can select it in one go —
`df.filter(like="kv_")` in pandas, `starts_with("kv_")` in dplyr. Values that are
entirely numeric are converted to real numbers, so concentrations sort and plot
correctly instead of ordering as "10" < "5". Repeated keys (OMERO permits them)
are joined with `; `.

**Grouping by condition.** `--group-by Compound` additionally writes
`per_kv_Compound_summary.csv` — the same statistics computed per condition rather
than per well — plus `populations_by_kv_Compound.png` and
`combinations_by_kv_Compound.png`. Use the key name as it appears in OMERO; the
prefix is added for you.

A caveat worth stating: pooling cells across the wells of a condition, which is
what `--group-by` does, is convenient for plotting but it is *not* how you should
compute statistics. See [section 16](#16-statistics-your-n-is-not-what-you-think)
— for that, take the per-well rows and average them within condition, so each
well contributes once.

`well_metadata.csv` is also written per plate: a plain plate-map of well →
key-values, including wells where no cells were measured, which is often the
quickest way to spot that an annotation is missing.

### "Double positive" — but which two?

With three markers there are three different double-positive populations, and
lumping them together usually hides the interesting part. Both readouts are
produced everywhere:

- **`population`** answers *how many* — `single`, `double`, `triple`.
- **`combination`** answers *which* — `GFP+RFP+`, `GFP+Cy5+`, `RFP+Cy5+`,
  `GFP+RFP+Cy5+`, `all-negative`.

`per_well_summary.csv` carries both, so you can run statistics on either. For a
three-marker screen you get `pct_single`, `pct_double`, `pct_triple` **and**
`pct_GFP+`, `pct_RFP+`, `pct_Cy5+`, `pct_GFP+RFP+`, `pct_GFP+Cy5+`,
`pct_RFP+Cy5+`, `pct_GFP+RFP+Cy5+`, one column per combination actually observed.
The individual double columns sum to `pct_double`, which is a useful arithmetic
check that gating and counting agree.

Three plots show the same split: `combinations_by_well.png` (stacked bars per
well, by combination), `marker_combinations.png` (plate-wide percentages), and
`platemap_pct_<combination>.png` for each multi-marker combination. The biaxial
quadrant plots give you the pairwise view — but note that the upper-right
quadrant of a `GFP vs RFP` plot contains cells that are GFP+RFP+ *whether or not*
they are also Cy5+, which is exactly the ambiguity the `combination` column
removes.

## 12. How positivity is decided, and why it's the weak point

Everything downstream — your single/double/triple percentages, your quadrant
plots, your conclusions — rests on one number per marker: the cutoff. Segmentation
errors of a few percent barely move your result. **A misplaced threshold can
change it by a factor of two.** Treat this as the part to get right.

The script offers five ways to set it:

| `--threshold-method` | How it decides | When to use it |
|---|---|---|
| `otsu` *(default)* | finds the split that best separates the data into two groups | quick look; only valid if the histogram really has two humps |
| `gmm` | fits two overlapping bell curves and cuts where one becomes more likely than the other | two humps that overlap heavily |
| `control` | mean + *k* SD of your negative-control wells | **best option if you have proper controls** |
| `quantile` | "the top 1% are positive" | when you know the expected frequency |
| `manual` | you supply the number | when you've picked it by eye and want it fixed |

**The failure mode to watch for.** Otsu and GMM both *assume* there are two
populations. If a marker is present in 99% of cells, or in none, they will still
dutifully split the single hump down the middle and report ~50% positive. This is
the single most common way to get a confidently wrong answer out of a pipeline
like this. That's exactly what `histograms_gated.png` is for.

If you have unstained or secondary-only wells, use them:

```bash
--threshold-method control --control-wells A1,A2 --control-sd 3
```

This sets the cutoff at 3 standard deviations above the control mean (computed on
a log scale, which is the right scale for fluorescence), which is the same logic
as drawing a gate on an unstained sample in flow cytometry.

**`--threshold-scope` — how much data goes into one cutoff.**

| Scope | Effect |
|---|---|
| `screen` | one cutoff for everything; most comparable, most sensitive to plate effects |
| `plate` *(default)* | one per plate; a sensible compromise |
| `well` | one per well — **usually wrong**, because it forces every well to have similar positive fractions and erases your treatment effect |
| `image` | one per field; only for severe field-to-field illumination drift |

The tradeoff is bias versus variance: narrower scope adapts to local artefacts
but absorbs real biological differences into the cutoff. Default to `plate`, and
only narrow it if you can see an artefact that justifies it.

## 13. Which intensity statistic to use

`--intensity-stat` changes which number gets gated:

- **`mean`** *(default)* — concentration-like. Independent of nuclear size. The
  right default for most stains.
- **`median`** — the same but immune to a few saturated pixels or a bright speck
  of debris inside the nucleus. Worth trying if your images are noisy.
- **`integrated`** — total signal (mean × area). Use this when *amount* matters
  rather than concentration. **Important:** integrated DAPI is proportional to
  DNA content, so it doubles between G1 and G2. If you gate DAPI on integrated
  intensity you're gating on cell cycle, which is either a feature (you can pull
  out a 2N/4N distribution from the DAPI histogram) or an accident.

## 14. The dead-cell filter

Off unless you pass `--filter-dead`.

**The biology.** Apoptotic and pyknotic nuclei condense the same amount of DNA
into a much smaller volume. Per-pixel DAPI intensity therefore goes *up* sharply
while area goes *down*. These cells are also often leaky and non-specifically
bright in every other channel, so they masquerade as double and triple positives.
In a toxicity screen this can be a large fraction of your "interesting" cells.

**What the filter does.** It computes a cutoff on DAPI intensity and flags
everything above it. `--dead-require-small` additionally requires the nucleus to
be in the smallest quarter by area, which makes it much more specific — bright
*and* condensed, rather than merely bright.

Importantly, this runs **before** the marker thresholds are computed, so dead
cells don't drag your Otsu cutoffs upward either.

| Flag | Meaning |
|---|---|
| `--dead-method mad` *(default)* | median + `--dead-k` robust SDs of log DAPI |
| `--dead-method quantile` | top `--dead-quantile` fraction (default 0.99) |
| `--dead-method otsu` | two-group split; only if the dead population is large |
| `--dead-method manual --dead-threshold 8000` | fixed value |
| `--dead-stat` | which DAPI statistic; keep `mean` (see the G2 warning above) |
| `--dead-action flag` | keep the cells, just add a `dead` column |

**Check `dead_cell_gate.png`.** The right-hand panel plots area against DAPI
intensity; condensed dying nuclei form a distinct cloud in the top-left. If your
red points aren't a separate cloud, the filter is cutting into your main
population and you should raise `--dead-k`.

**Mitotic cells look similar.** Condensed mitotic chromatin is also small and
bright, so this filter removes mitotic cells too. Usually that's acceptable. If
mitotic index matters to you, use `--dead-action flag` and separate the two
yourself — phospho-H3 staining is the honest way to tell them apart.

Also: `pct_dead_dapi` in the per-well summary and `platemap_pct_dead_dapi.png`
are real data about your treatment's toxicity, not just a QC number. Look at them.

## 15. Foci counting

Off unless you pass `--foci-channel`.

**Why mean intensity isn't enough.** For a punctate marker like γ-H2AX, a cell
with 20 tight foci and a cell with diffuse pan-nuclear staining can have exactly
the same mean nuclear intensity. Counting discrete objects gives you a different,
usually more meaningful readout — and comparing the two is a useful sanity check
(the middle panel of `foci_<channel>.png` does exactly this).

```bash
--foci-channel gH2AX --foci-min-count 3 --save-foci-table
```

The flag is repeatable, so `--foci-channel gH2AX --foci-channel 53BP1` scores both.

**Columns added per nucleus**, for each foci channel:

| Column | Meaning |
|---|---|
| `foci_n_<ch>` | number of foci |
| `foci_area_total_<ch>` | combined area of all foci |
| `foci_area_fraction_<ch>` | fraction of the nucleus covered by foci |
| `foci_area_mean_<ch>` | average focus size |
| `foci_mean_intensity_<ch>` | average brightness within foci |
| `foci_max_intensity_<ch>` | brightest focus pixel |
| `foci_integrated_intensity_<ch>` | total signal in foci |
| `foci_pos_<ch>` | at least `--foci-min-count` foci |

**Counts saturate.** Validated on synthetic images with a known number of foci
per nucleus, detected counts correlate with truth at **r = 0.97**, with a mild
undercount at high density where foci physically overlap and merge. Above roughly
30–40 foci per nucleus you should stop reporting counts and switch to
`foci_area_fraction` or `foci_integrated_intensity`, which keep increasing after
counting saturates. This is a property of foci counting in general, not of this
implementation — the same is true in ImageJ or CellProfiler.

### Low magnification (20x) — read this before tuning anything

**The detector scales itself to your pixel size, using the physical focus size.**
A gamma-H2AX focus is ~0.8 um across. At 60x (~0.1 um/px) that's 8 pixels; at 20x
binned (~0.65 um/px) it's barely **1 pixel**. Parameters that work at 60x fail
badly at 20x, and not in an obvious way — you get plausible-looking counts that
are systematically 25-30% low.

By default the script reads the pixel size from OMERO metadata and derives the
LoG scales, top-hat radius, size filters and an upsampling factor from it. It
logs what it chose:

```
pixel size 0.6500 um/px -> foci params: sigma 0.78-2.35 px, tophat r=3,
                           area 1-10 px, upsample x3 [auto from 0.650 um/px]
```

Measured on synthetic fields with a known number of foci per nucleus:

| Condition | Focus size | Recall before | Recall after |
|---|---|---|---|
| 20x binned, 0.65 um/px, dim foci | ~1.2 px | 0.71 | **0.95** |
| 20x binned, 0.65 um/px, bright | ~1.2 px | 0.71 | **1.04** |
| 20x fine, 0.325 um/px, dim | ~2.5 px | 0.93 | 0.88 |
| 60x, 0.108 um/px, dim | ~7.4 px | 0.95 | 0.95 |

Three changes produce this: scales and top-hat radius derived from the focus
size; **2-3x upsampling** when foci are only 1-3 px wide, so LoG scale selection
operates where it's numerically stable; and a **sigma-clipped background**, so
that in heavily damaged nuclei the foci stop inflating the very background they
are compared against. Correlation with truth is 0.94-0.96 across all conditions
rather than excellent at 60x and poor at 20x.

If OMERO has no pixel size, the script warns and falls back to raw pixel
parameters. Supply it yourself with `--pixel-size-um 0.65`. If your foci aren't
~0.8 um, set `--foci-diameter-um`. To go back to hand-tuned pixel values, pass
`--no-foci-auto-scale`.

**Always check with `--save-foci-crops 5`.** This writes zoomed montages of
individual nuclei with each detected focus circled in red. At 20x this is the
*only* view that shows what the detector is doing — a 2048x2048 full-field
overlay squashed into a PNG cannot show a 1-pixel focus, so the normal overlays
are useless for this purpose. Look for foci circled twice (over-splitting),
obvious foci not circled (raise sensitivity), or circles on empty background
(noise being counted).

**Honest limits at 20x.** With ~0.65 um pixels you are below the diffraction
limit for these objects: two foci 1 um apart cannot be separated, period. No
software fixes that. If you need accurate counts at high damage levels, image at
higher magnification; if you must stay at 20x, prefer
`foci_integrated_intensity` or `foci_area_fraction`, which degrade gracefully
where counts saturate. Counts at 20x are best treated as comparative between
conditions imaged identically, not as absolute numbers of double-strand breaks.

**Tuning.** In the order you should reach for them:

1. `--foci-log-threshold` (default 0.02) — detector sensitivity, lower finds more
2. `--foci-k` (default 3) — how far above each nucleus's own background counts
3. `--foci-min-area` / `--foci-max-area` — size window in pixels
4. `--foci-tophat-radius` (default 5) — must be **larger than your biggest focus**

Always confirm with the overlays (`--save-overlays 5`), where foci are drawn in
green. Foci parameters are the easiest thing in this pipeline to fool yourself
with; five minutes looking at overlays is worth more than any amount of argument
about defaults.

**Foci as the positivity criterion.** By default the foci channel is *also*
gated on mean intensity like any other marker. `--foci-as-marker` switches that
channel's +/− call to "has ≥ `--foci-min-count` foci", which then feeds into the
single/double/triple classification. Both are defensible; state which you used.

## 16. Statistics: your n is not what you think

The per-nucleus table might have 500,000 rows. **That is not n = 500,000.**

Cells within a well are not independent — they share a treatment, a well, a
staining, an illumination pattern. Treating cells as replicates ("pseudo-
replication") produces p-values that are essentially arbitrary; with half a
million cells you can reach p < 10⁻¹⁰ for a difference of no consequence.

The unit of replication is whatever you randomised: usually the **well**, and
across independent experiments, the **plate**. So:

- Use `per_well_summary.csv`, one row per well, for statistics — it carries both
  the `pct_double`-style columns and the per-combination `pct_GFP+RFP+`-style ones,
  plus your OMERO key-value metadata, so it groups by condition directly.
- `per_<key>_summary.csv` from `--group-by` pools *cells* across a condition. It's
  for plots and quick looks, not for p-values — average the per-well rows instead.
- Use `per_image_summary.csv` for QC (spotting a bad field), not for stats —
  fields within a well aren't independent either.
- Report the cell count as what it is: how precisely each well was measured, not
  your sample size.
- Wells with very few cells are noisy; consider dropping wells below some
  `n_cells` (the column is right there) and say so in your methods.
- If your effect only appears when you pool cells and vanishes when you analyse
  well means, you don't have an effect.

The per-cell distributions are still genuinely useful — for showing bimodality,
for gating, for spotting subpopulations. Just don't run a t-test on them.

## 17. Controls worth running

- **Unstained / secondary-only wells** — the only principled way to set a
  positive/negative cutoff. Feed them to `--control-wells`.
- **Single-stain wells** — needed to check bleed-through. If channels overlap
  spectrally, your double-positive fraction is inflated and no amount of
  thresholding fixes it; compensate or unmix upstream.
- **A known-positive control** — confirms the assay worked at all on that plate.
- **A dose response**, if you can — a monotonic trend across doses is far more
  convincing than a single treated-vs-untreated comparison, and it validates your
  threshold at the same time.

---
---

# Part 3 — How it works

*This part is for readers who want the algorithms, the assumptions behind them,
and enough structure to modify the code.*

## 18. Architecture and dataflow

Single file, no package, ~1,600 lines, deliberately linear. Sections in order:
config → OMERO IO → segmentation → measurement → foci → gating → plotting →
OMERO write-back → summaries → orchestration → CLI.

```
BlitzGateway
     │  Screen ─▶ Plate ─▶ Well ─▶ WellSample ─▶ Image
     ▼
read_stack(image)                       (H, W, C) float32, one z/t
     │
     ├─▶ segmenter(stack[..., dapi])    ─▶ int32 label image
     │        └─ filter_labels()           border / area / solidity
     │
     ├─▶ measure(labels, stack)         ─▶ DataFrame, one row per nucleus
     │        └─ scipy.ndimage per-label reductions, all channels
     │
     └─▶ detect_foci(plane, labels)     ─▶ (focus label image, per-focus DataFrame)
              └─ summarise_foci()          collapsed onto nucleus rows
                       │
                       ▼
              concat over images ─▶ per-plate DataFrame
                       │
                       ▼
              flag_dead_cells()   ─▶ 'dead' column, then optional subset
                       │
                       ▼
              gate()              ─▶ pos_*, n_positive, population, combination
                       │
             ┌─────────┴──────────┐
             ▼                    ▼
        summarise()           plot_*()          ─▶ CSV / parquet / PNG
```

Two design decisions shape everything else:

1. **One long-form table is the interchange format.** Every downstream step —
   gating, summaries, all plots — consumes a tidy DataFrame with one row per
   nucleus. This is what makes `--from-csv` possible: the entire post-measurement
   half of the pipeline runs offline on the saved table, in seconds rather than
   hours. Anything that needs pixels (segmentation, foci detection) is on the
   near side of that boundary and cannot be replayed.
2. **Failures are per-image, not per-run.** The image loop catches and logs
   exceptions per image. On a 1,536-image plate, one corrupt field should not
   cost you the run.

Deterministic given the same inputs and versions, with two caveats: GMM
thresholding uses a fixed `random_state=0`, and scatter subsampling for plots
uses a fixed seed, but GPU non-determinism in the segmentation backends is
outside the script's control. `run_config.json` records every setting; the
password is explicitly excluded.

## 19. Segmentation

Two backends behind one `Segmenter.__call__(dapi) -> int32 labels` interface, so
adding a third is a small change.

**StarDist (default, `2D_versatile_fluo`).** Predicts, for every pixel, an
object probability and a set of radial distances to the object boundary along a
fixed number of directions — a **star-convex polygon**. Candidates are then
resolved by non-maximum suppression. The reason this is the default for
*dense* nuclei: NMS operates on whole candidate shapes, so two adjacent nuclei
each get their own polygon rather than being fused by a shared boundary. Nuclei
are close to star-convex, so the representation loses almost nothing. Knobs:
`--prob-thresh` (detection sensitivity) and `--nms-thresh` (permitted overlap;
raise it in very dense fields).

**Cellpose (`--model cellpose`).** Predicts spatial gradient flows and assigns
pixels to attractors by following them. More flexible on non-star-convex shapes,
and the v4 generalist (*cpsam*) is remarkable out of the box; the wrapper tries
`CellposeModel()` first and falls back to the v3 `Cellpose(model_type=...)` API.

Input normalisation is percentile-based (1st–99.8th → [0,1]), which is what both
models expect and makes the thresholds transferable across exposure settings.

`filter_labels()` then applies `clear_border` (a nucleus clipped by the field edge
has a truncated, meaningless intensity) plus an area window, via a vectorised
`np.isin` relabel rather than a Python loop over objects.

## 20. Measurement

Intensity statistics come from `scipy.ndimage` label-indexed reductions
(`ndi.mean`, `ndi.median`, `ndi.sum_labels`) rather than
`regionprops_table(intensity_image=...)`. Two reasons: the ndimage calls are a
single pass per channel over the label array, and — the practical one —
scikit-image renamed most region properties between 0.18 and 0.19
(`mean_intensity` → `intensity_mean`, and so on), so pinning to ndimage avoids a
whole class of version-skew bugs. Morphology, where the names are stable, still
comes from `regionprops_table`.

**Background.** `--background median` subtracts, per channel per image, the
median of all pixels *outside* any nucleus. Using the median rather than the mean
makes it insensitive to the bright objects that occupy part of that region, and
restricting it to non-nuclear pixels avoids the circularity of estimating
background from the signal you're about to measure. It's a single scalar per
image: a flat offset correction, not a shading correction. Vignetting or uneven
illumination needs a flat-field correction upstream — this will not fix it, and
`--threshold-scope image` only papers over it.

Complexity per image is O(H·W·C) for measurement, dominated in practice by the
segmentation forward pass. Memory is one field at a time: a 2048² image with 4
channels in float32 is 64 MB.

## 21. Thresholding algorithms

All thresholds operate on `log1p(x)` and are mapped back with `expm1`.
Fluorescence intensities are approximately log-normal within a population, and
`log1p` handles the zeros and small negatives that background subtraction
produces, where plain `log` would not.

**Otsu.** Exhaustive search for the cut minimising intra-class variance —
equivalently maximising between-class variance. Assumes exactly two classes.
Given a unimodal input it returns a number regardless, which is the failure mode
described in section 12.

**GMM.** Two-component `GaussianMixture` (`n_init=3`, `random_state=0`) on log
intensity. The threshold is where the posterior for the higher-mean component
first reaches 0.5, found by scanning a 2,000-point grid — the analytic solution
is a quadratic with two roots and awkward edge cases, and at this cost the grid
is not worth optimising. Handles overlapping populations better than Otsu because
it models the overlap explicitly. Falls back to Otsu with a warning if
scikit-learn is absent.

**Control-based.** `expm1(mean(log1p(control)) + k · sd(log1p(control)))`. The
most defensible option, and the direct analogue of gating on an unstained sample.

**MAD-based** (used by the dead filter). `median + k · 1.4826 · MAD`. The 1.4826
factor makes the MAD a consistent estimator of σ for Gaussian data. The point is
robustness: the population being detected is in the upper tail of the same
distribution being used to estimate the spread, so a plain mean/SD would be
inflated by the very cells you're trying to find. MAD has a breakdown point of
50%, so it stays put until half your cells are dying.

**Scope** (`screen`/`plate`/`well`/`image`) is a pooling choice, and it is a
bias–variance tradeoff with a trap: narrowing the scope reduces sensitivity to
batch artefacts but silently absorbs genuine biological differences into the
cutoff. Per-well thresholds force every well toward the same positive fraction —
exactly what you were trying to measure.

## 22. Dead-cell detection

```
thr = expm1( median(log1p(I_dapi)) + k · 1.4826 · MAD(log1p(I_dapi)) )
dead = I_dapi > thr   [ ∧  area ≤ quantile(area, q)  if --dead-require-small ]
```

Computed per scope group and applied **before** marker gating, so the dead
population contaminates neither the thresholds nor the counts. The conjunction
with the area criterion is what turns a generic outlier filter into something
specific to chromatin condensation: bright *and* small is a much narrower claim
than bright alone, and it is the actual morphological signature of pyknosis.

Excluded rows are written to `excluded_dead_cells.csv.gz` rather than dropped
silently — the removed population is data, and you should be able to audit it.

## 23. Foci detection

The interesting part of the pipeline. Three stages:

**1 — Background flattening.** A white top-hat, `f - (f ∘ B)` where `∘` is a
morphological opening with a disk of radius r. An opening with a disk of radius r
removes everything that cannot contain that disk, so subtracting it keeps small
bright structures and discards slowly varying background — including the diffuse
pan-nuclear signal that makes a global intensity threshold useless. The
constraint is simply r > the largest focus radius, hence `--foci-tophat-radius`.

**2 — Per-nucleus adaptive threshold.** This is the design decision that matters
most. For each nucleus i,

```
thr_i = median_i(tophat) + k · 1.4826 · MAD_i(tophat)
```

computed with label-indexed ndimage reductions and broadcast back to pixels
through a lookup table (`_label_lut`), so no Python loop over nuclei. The
estimate is **sigma-clipped** (two iterations, excluding pixels above
`median + 3 sigma`, implemented by zeroing those pixels in a copy of the label
array so the same ndimage reductions apply): at high foci burden the naive
median sits *inside* the foci distribution and the threshold runs away upward,
so without clipping the most heavily damaged nuclei — the interesting ones — are
exactly the ones that get undercounted. Foci
brightness varies by an order of magnitude between cells in the same field —
biologically real, plus exposure and staining variation. A single global cutoff
either misses every focus in dim cells or floods bright ones. Making the
threshold relative to each nucleus's *own* background makes the count a property
of that nucleus rather than of its overall brightness. The cost: the count is now
a relative measure, and a cell with genuinely uniform high signal will report
zero foci. That is usually the desired behaviour and is occasionally not; hence
`--foci-abs-threshold` for a fixed cutoff across a screen.

**3 — Splitting touching foci.** The thresholded mask is connected components,
which merges adjacent foci. Seeds come from one of:

- `log` *(default)* — Laplacian-of-Gaussian scale-space. The LoG response is
  maximised where blob size matches √2·σ, so scanning σ over
  `[--foci-min-sigma, --foci-max-sigma]` gives scale-aware detection with
  sub-object resolution. Blob centres outside any nucleus are discarded.
- `hmax` — h-maxima: regional maxima that stand at least h above their
  surroundings, suppressing noise peaks. h is auto-set to half the 99th–50th
  percentile spread if not given.
- `tophat` — no seeding, plain connected components. Fastest, merges neighbours.

Seeds then drive a watershed on the negated top-hat image, constrained to the
mask, which gives real focus boundaries (needed for area and integrated
intensity) rather than just centres.

**Assignment to nuclei.** A focus can straddle two touching nuclei. Rather than
using the centroid pixel (fragile for concave shapes), each focus is assigned to
the nucleus holding the majority of its pixels, computed by encoding
`(focus_id, nucleus_id)` pairs as a single integer key `focus·span + nucleus`,
running `np.unique(..., return_counts=True)`, and taking the argmax per focus via
`lexsort`. Fully vectorised — the alternative loop over thousands of foci per
image is the difference between seconds and minutes across a screen.

**Validation.** On synthetic fields with a known number of Gaussian foci per
nucleus: r = 0.97 between true and detected counts, |error| ≤ 1 for 62% of
nuclei, with systematic undercounting at high density where foci genuinely
overlap. Any counting method has this ceiling; the honest response is to switch
readout (area fraction, integrated intensity) rather than to tune it away.

## 24. Why asinh axes

Flow cytometry stopped using pure log axes because log is undefined at zero and
explodes near it, while compensated or background-subtracted data legitimately
contains zeros and small negatives — producing the notorious pile-up at the axis.
The standard fix is a biexponential/logicle transform: linear near zero, log-like
far from it. `asinh(x/c)` has exactly that behaviour, with one parameter and no
fitting:

```
asinh(x/c) ≈ x/c        for |x| ≪ c      (linear, symmetric about 0)
asinh(x/c) ≈ ln(2x/c)   for x ≫ c        (log-like)
```

The cofactor c sets where the transition happens; it defaults to a robust
per-channel estimate (`max(1, 0.5·P25(|x|))`) and is overridable with
`--cofactor`. Axis ticks are drawn at real decade values mapped through the
transform, so you read them as ordinary intensities.

Density colouring in the biaxial plots is a 2D histogram lookup — `histogram2d`,
then `digitize` each point back into its bin — which is O(n) and works at
hundreds of thousands of points, unlike a KDE. Points are drawn in ascending
density order so the dense core isn't buried under the sparse halo, and
`--max-points` subsamples for the scatter while the quadrant percentages are
always computed on the *full* data.

## 25. Implementation notes and performance

- **Single process, no parallelism across images.** Deliberate: the GPU is the
  bottleneck and would be contended, and OMERO connections are stateful. To use
  more cores, run separate processes per plate (`--plate-id`) and merge the
  per-well summaries afterwards.
- **Throughput** is dominated by the segmentation forward pass: roughly 1–3 s per
  2048² field on a modern GPU, 10–40 s on CPU. Foci detection adds ~0.5–2 s per
  channel. Plotting is negligible until `--per-well-plots` on a 384-well plate.
- **Memory** is bounded by one field: `H·W·C·4` bytes, plus the accumulating
  DataFrame (~200 bytes/nucleus/channel). A screen of 500k nuclei × 4 channels
  is a few hundred MB — fine, but `--per-nucleus-format parquet` writes it 5–10×
  smaller and far faster than CSV, with dtypes preserved.
- **Tiling** via `--n-tiles` splits the field for the segmentation forward pass
  when GPU memory is short.
- **The label column is per-image.** Join the focus table to the nucleus table on
  `(image_id, nucleus_label)` ↔ `(image_id, label)`. There is deliberately no
  global cell ID; `image_id` is a stable OMERO identifier and composite keys make
  the provenance explicit.
- **ROI upload** (`--upload-rois`) packs each mask with `np.packbits` into an
  OMERO `MaskI` shape, one ROI per image. It's slow on large screens and off by
  default, but being able to see the segmentation in OMERO.iviewer next to the
  original data is worth it on a pilot.

## 26. Extending the script

- **A new measurement per nucleus** — add it in `measure()`. Anything expressible
  as a label-indexed reduction is a one-liner; everything downstream picks up new
  columns automatically as long as you follow the `<stat>_<channel>` convention.
- **A cytoplasmic ring** — dilate the label image, subtract the nuclei, and run
  `measure()` again on the resulting annulus labels. The measurement code is
  agnostic about what the labels represent.
- **A new segmentation backend** — implement `__call__(dapi) -> int32 labels` in
  `Segmenter` and add a choice to `--model`.
- **A new thresholding rule** — add a branch to `compute_threshold()`; scope
  handling, JSON logging and plotting are already generic.
- **A new plot** — write `plot_x(df, ..., out, title)` and call it from
  `analyse()`. Every plot takes the tidy DataFrame, so nothing else needs to know.
- **3D / z-stacks** — currently a mid-plane or MIP (`--z-mode`). True 3D would
  mean swapping in StarDist-3D and moving the measurement to 3D label arrays;
  ndimage reductions already work on N-D input, so `measure()` needs less change
  than you'd expect.

---
---

# Part 4 — Reference

## 27. Complete flag reference

**Connection**

| Flag | Default | Meaning |
|---|---|---|
| `--host` | `localhost` | OMERO server address |
| `--port` | `4064` | OMERO port |
| `--user` | — | username |
| `--password` | — | prefer the `OMERO_PASSWORD` env var, or let it prompt |
| `--insecure` | off | disable the secure connection |
| `--group` | `-1` | OMERO group; `-1` = all groups |

**Data selection**

| Flag | Default | Meaning |
|---|---|---|
| `--screen-id` | — | screen to process |
| `--plate-id` | — | a single plate; repeatable, and combinable with `--screen-id` |
| `--wells` | all | e.g. `A1,A2,B3` |
| `--max-fields` | all | first N fields per well |
| `--dapi-channel` | auto | name fragment or index; auto-detects DAPI/Hoechst/405/w1 |
| `--z-mode` | `mid` | `mid` plane or `mip` max-projection |
| `--timepoint` | `0` | which timepoint |
| `--server-tile-height` | `0` | read planes in N-row strips; use `256` if the server throws Java heap errors |

**Segmentation**

| Flag | Default | Meaning |
|---|---|---|
| `--model` | `stardist` | `stardist` or `cellpose` |
| `--stardist-model` | `2D_versatile_fluo` | pretrained model name |
| `--cellpose-model` | `nuclei` | used by the Cellpose 3 API |
| `--prob-thresh` | model default | lower = more objects |
| `--nms-thresh` | model default | higher = allows more overlap |
| `--diameter` | auto | Cellpose expected diameter |
| `--n-tiles` | off | tile the forward pass if GPU memory is short |
| `--no-gpu` | off | force CPU |
| `--gpu-memory-limit` | — | cap TensorFlow's GPU memory in MB, e.g. `4096` |
| `--gc-every` | `50` | garbage-collect every N images (0 = never) |
| `--reload-model-every` | `0` | rebuild the model every N images — cures a slow GPU memory leak |
| `--min-area` / `--max-area` | `30` / `100000` | size window in pixels |
| `--min-solidity` | `0` | drop ragged/merged objects |
| `--keep-border` | off | keep nuclei touching the image edge |

**Dead / dying filter**

| Flag | Default | Meaning |
|---|---|---|
| `--filter-dead` | off | enable it |
| `--dead-method` | `mad` | `mad`, `quantile`, `otsu`, `manual` |
| `--dead-stat` | `mean` | DAPI statistic used |
| `--dead-k` | `3.0` | robust SDs above the median |
| `--dead-quantile` | `0.99` | for `quantile` |
| `--dead-threshold` | — | for `manual` |
| `--dead-scope` | `plate` | pooling level |
| `--dead-require-small` | off | also require small area — recommended |
| `--dead-area-quantile` | `0.25` | what counts as small |
| `--dead-action` | `remove` | `remove` or `flag` |

**Measurement and gating**

| Flag | Default | Meaning |
|---|---|---|
| `--background` | `median` | flat background subtraction, or `none` |
| `--intensity-stat` | `mean` | `mean`, `median`, `integrated` |
| `--threshold-method` | `otsu` | `otsu`, `gmm`, `quantile`, `control`, `manual` |
| `--threshold-scope` | `plate` | `screen`, `plate`, `well`, `image` |
| `--quantile` | `0.99` | for `quantile` |
| `--control-wells` | — | e.g. `A1,A2` |
| `--control-sd` | `3.0` | SDs above the control mean |
| `--manual-thresholds` | — | `'GFP=1200,RFP=800'` |
| `--marker-channels` | all non-DAPI | restrict which channels are gated |
| `--subset-analyses` | off | also analyse marker subsets as self-contained panels |
| `--subset-sizes` | `1,2` | which panel sizes to generate |
| `--panel` | — | an explicit panel, e.g. `--panel GFP+RFP`; repeatable |

**Foci**

| Flag | Default | Meaning |
|---|---|---|
| `--foci-channel` | — | name or index; repeatable |
| `--foci-method` | `log` | `log`, `tophat`, `hmax` |
| `--foci-diameter-um` | `0.8` | expected physical focus diameter — drives all derived parameters |
| `--pixel-size-um` | from OMERO | override the metadata pixel size |
| `--no-foci-auto-scale` | off | use raw pixel parameters instead of deriving them |
| `--foci-upsample` | `0` (auto) | 2-3x upsampling for 1-3 px foci |
| `--no-foci-bg-clip` | off | disable sigma-clipping of the in-nucleus background |
| `--save-foci-crops` | `0` | N zoomed montages with foci circled — **use this at 20x** |
| `--foci-tophat-radius` | `5.0` | must exceed the largest focus |
| `--foci-min-sigma` / `--foci-max-sigma` | `1.0` / `4.0` | LoG scale range |
| `--foci-num-sigma` | `5` | scales sampled |
| `--foci-log-threshold` | `0.02` | detector sensitivity; lower = more foci |
| `--foci-k` | `3.0` | SDs above each nucleus's own background |
| `--foci-abs-threshold` | — | fixed cutoff instead of per-nucleus |
| `--foci-h` | auto | h-maxima depth |
| `--foci-min-area` / `--foci-max-area` | `3` / `400` | focus size window |
| `--foci-min-count` | `3` | foci-positive cutoff |
| `--foci-as-marker` | off | count decides +/−, not intensity |

**OMERO key-value pairs (map annotations)**

| Flag | Default | Meaning |
|---|---|---|
| `--no-key-values` | off | don't fetch key-value pairs |
| `--kv-levels` | `screen,plate,well` | which levels to harvest; add `image` for per-field annotations |
| `--kv-prefix` | `kv_` | prefix for the added columns |
| `--group-by` | — | key to group summaries and plots by, e.g. `--group-by Compound` |

**Output**

| Flag | Default | Meaning |
|---|---|---|
| `--per-nucleus-format` | `csv.gz` | `csv`, `csv.gz`, `parquet`, `xlsx`, `none` |
| `--save-foci-table` | off | one row per focus |
| `--export-explorer` | off | write the interactive HTML explorer next to the results |
| `--explorer-no-embed` | off | write the explorer without inlining the data |
| `--cofactor` | auto | asinh cofactor |
| `--max-points` | `50000` | points drawn in scatter plots |
| `--per-well-plots` | off | full plot set per well |
| `--save-overlays` | `0` | N segmentation QC images per plate |
| `--overlay-line-width` | `2` | outline thickness in px for overlays and foci crops |
| `--outdir` | `./results` | where results go |
| `--pool-screen` | off | extra analysis pooling all processed plates |
| `--upload-rois` | off | push masks to OMERO as ROIs |
| `--attach-results` | off | attach CSVs/plots to the plate in OMERO |
| `--from-csv` | — | offline re-gate and re-plot |
| `-v` | off | verbose logging |
| `--check-env` | — | print package versions and GPU status, then exit |

## 28. Recipes

**Pilot before committing**
```bash
--wells A1,A2 --max-fields 2 --save-overlays 4 --outdir ./pilot
```

**Proper controls-based gating**
```bash
--threshold-method control --control-wells A1,A2,A3 --control-sd 3
```

**DNA-damage assay with dead cells removed**
```bash
--filter-dead --dead-require-small \
--foci-channel gH2AX --foci-min-count 5 --foci-as-marker --save-foci-table
```

**Foci at 20x (small foci) — always pilot with the crops first**
```bash
--foci-channel gH2AX --foci-diameter-um 0.8 --save-foci-crops 5 \
--wells A1,B1 --max-fields 2
```

**Dense nuclei that are being merged**
```bash
--nms-thresh 0.5 --prob-thresh 0.4 --save-overlays 10
```

**Mixed staining panels on one plate**
```bash
--subset-analyses                      # every single marker and every pair
--panel GFP+RFP --panel GFP+Cy5        # or just the panels you used
```

**Hand a plate to a student for interactive re-gating**
```bash
--export-explorer --per-nucleus-format csv
```

**Summarise by treatment instead of by well**
```bash
--group-by Compound          # uses the OMERO key-value pair named "Compound"
```

**Long screen that exhausts GPU memory partway through**
```bash
--reload-model-every 200 --gpu-memory-limit 4096
```

**Big screen, fast output**
```bash
--per-nucleus-format parquet --max-points 20000
```

**Re-gate yesterday's run without the server**
```bash
--from-csv results/Plate_01/per_nucleus_measurements.parquet \
--threshold-method gmm --outdir ./results_v2
```

**One plate only**
```bash
--plate-id 5678 --outdir ./one_plate
```

**Several named plates in one run** (repeatable *and* comma-separated)
```bash
--plate-id 101,102,103 --pool-screen
--plate-id 101,102 --plate-id 103        # equivalent
--screen-id 1234 --plate-id 999          # a screen plus one extra plate
```

**Run four plates in parallel** (four terminals, or `&` on one line)
```bash
python omero_nuclei_flow.py --plate-id 101 --outdir ./p101 ...
python omero_nuclei_flow.py --plate-id 102 --outdir ./p102 ...
```
Each process needs its own `--outdir`. **Only do this if the OMERO server has
headroom** — parallel clients multiply the server's memory pressure, and are a
good way to trigger `Java heap space` errors. Merge the `per_well_summary.csv`
files afterwards — they have identical columns, so in pandas:
`pd.concat([pd.read_csv(f) for f in glob("p*/**/per_well_summary.csv")])`.

## 29. Known limitations

- **Nuclear intensity only.** Cytoplasmic markers will be underestimated. See
  [section 26](#26-extending-the-script) for the ring-measurement approach.
- **No spectral compensation.** Overlapping fluorophores inflate double- and
  triple-positive fractions, and no threshold choice fixes that. Unmix upstream.
- **Flat background subtraction, not shading correction.** Vignetting needs
  flat-field correction before this script sees the data.
- **2D only.** Mid-plane or MIP; a MIP of a thick stack overestimates mean
  intensity and merges foci that are separated in z.
- **Foci counts saturate** above ~30–40 per nucleus.
- **The dead filter also removes mitotic cells.**
- **Thresholds are the dominant source of error**, not segmentation.
- **`--from-csv` cannot redo segmentation or foci detection**, only everything
  computed from the table.
- **Key-value pairs are read, never written.** If a well isn't annotated in
  OMERO, its metadata columns are empty — check `well_metadata.csv`. Annotations
  on PlateAcquisitions or WellSamples aren't harvested; screen, plate, well and
  image are.
- **Not a substitute for looking at your images.**

## 30. Glossary (both directions)

**For computer scientists — the biology**

| Term | Meaning |
|---|---|
| Screen / Plate / Well / Field | Nested containers: a screen holds plates (usually 96 or 384 wells), each well holds several imaged positions ("fields") |
| ImageXpress | A high-content microscope that images plates automatically |
| Channel | One fluorescence colour = one grayscale image of one stain |
| DAPI / Hoechst | Dyes that bind DNA, so they light up nuclei — used to find cells |
| Marker | A stain reporting a biological property; cells are called + or − for it |
| Foci | Discrete bright dots inside a nucleus where a protein has accumulated |
| γ-H2AX | A histone modification marking DNA double-strand breaks; one focus ≈ one break |
| Pyknotic / apoptotic | Dying cells, whose chromatin condenses into a small, very bright nucleus |
| Flow cytometry | Measures one cell at a time in a fluid stream; the source of the histogram/quadrant plotting conventions used here |
| Gating | Choosing a cutoff to split cells into populations |
| Bleed-through | Signal from one fluorophore leaking into another channel |
| Key-value pair | OMERO's per-object metadata (a "map annotation"), e.g. `Compound = Cisplatin`; how plate maps are usually recorded |

**For biologists — the computing**

| Term | Meaning |
|---|---|
| Segmentation | Deciding which pixels belong to which object |
| Label image | An image where every pixel holds the ID number of the object it belongs to; 0 = background |
| StarDist | A neural network that outlines nuclei as star-shaped polygons; good at separating touching ones |
| Cellpose | An alternative network that outlines cells by following predicted gradient flows |
| NMS (non-maximum suppression) | Keeps the best of several overlapping candidate outlines |
| Otsu's method | Picks the cutoff that best splits a histogram into two groups |
| GMM | Fits two overlapping bell curves and cuts where one becomes more likely |
| MAD | Median absolute deviation — a spread measure that ignores outliers, unlike SD |
| Top-hat filter | Keeps small bright things, removes smooth background |
| LoG (Laplacian of Gaussian) | A blob detector that responds most strongly at a chosen spot size |
| Watershed | Splits touching objects by "flooding" from seed points |
| asinh | A log-like axis that behaves properly at and below zero |
| DataFrame | A table with named columns — what a CSV becomes once loaded |
| parquet | A compressed table format; smaller and much faster than CSV |
| CLI flag | The `--something` options you type after the script name |

