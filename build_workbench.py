#!/usr/bin/env python3
"""
Build nucleoflow_workbench.html — the gating explorer and the condition
explorer in one file, as two tabs over a single shared dataset.

The two apps are merged without editing either one's logic:

  * ids are prefixed per pane (G- / C-) in the HTML and CSS, so nothing
    collides and <label for=...> still points at the right input;
  * each app's JS runs inside an IIFE where `document` is shadowed by a
    proxy that adds the prefix to getElementById and scopes querySelector
    to that pane's subtree. The app code is untouched;
  * CSS rules are prefixed with the pane id, so the two palettes and the
    two sets of class names cannot interfere.

The point of the merge is the shared data: the gating tab writes its
positivity calls back into the columns that the condition tab reads, so a
gate you drag in tab 1 changes the IC50s in tab 2.
"""
import re, pathlib, sys

SRC_G = pathlib.Path("/mnt/user-data/outputs/gating_explorer.html")
SRC_C = pathlib.Path("/mnt/user-data/outputs/condition_explorer.html")
OUT   = pathlib.Path("/home/claude/nucleoflow_workbench.html")


# ---------------------------------------------------------------- parsing
def split_file(html):
    css  = re.search(r"<style>(.*?)</style>", html, re.S).group(1)
    body = re.search(r"<body>(.*?)</body>", html, re.S).group(1)
    body = re.sub(r"<script.*?</script>", "", body, flags=re.S)
    js   = html[html.rindex("<script>") + 8: html.rindex("</script>")]
    return css, body, js


# ---------------------------------------------------------------- CSS
def root_vars(css):
    m = re.search(r":root\{(.*?)\}", css, re.S)
    return m.group(1) if m else ""


def strip_root(css):
    return re.sub(r":root\{.*?\}", "", css, count=1, flags=re.S)


def split_rules(css):
    """Yield (selector, body) pairs, keeping @-blocks whole."""
    out, i, n = [], 0, len(css)
    while i < n:
        brace = css.find("{", i)
        if brace < 0:
            break
        sel = css[i:brace].strip()
        depth, j = 1, brace + 1
        while j < n and depth:
            if css[j] == "{":
                depth += 1
            elif css[j] == "}":
                depth -= 1
            j += 1
        out.append((sel, css[brace + 1:j - 1]))
        i = j
    return out


def scope_css(css, pane, pfx):
    """Prefix every selector with the pane id and every #id with the pane prefix."""
    pieces = []
    for sel, body in split_rules(css):
        if sel.startswith("@"):
            if sel.startswith(("@media", "@supports")):
                pieces.append(f"{sel}{{{scope_css(body, pane, pfx)}}}")
            else:                       # @keyframes, @font-face - leave alone
                pieces.append(f"{sel}{{{body}}}")
            continue
        parts = []
        for one in sel.split(","):
            one = one.strip()
            if not one:
                continue
            one = re.sub(r"#([A-Za-z][\w-]*)", lambda m: f"#{pfx}{m.group(1)}", one)
            # a rule on the pane root itself must not be nested inside itself
            parts.append(one if one.startswith(f"#{pane}") else f"#{pane} {one}")
        body = re.sub(r"#([A-Za-z][\w-]*)", lambda m: f"#{pfx}{m.group(1)}", body)
        pieces.append(f"{', '.join(parts)}{{{body}}}")
    return "\n".join(pieces)


# ---------------------------------------------------------------- HTML
def prefix_ids(body, pfx):
    body = re.sub(r'\bid="([\w-]+)"', lambda m: f'id="{pfx}{m.group(1)}"', body)
    body = re.sub(r'\bfor="([\w-]+)"', lambda m: f'for="{pfx}{m.group(1)}"', body)
    return body


# ---------------------------------------------------------------- assemble
def main():
    gcss, gbody, gjs = split_file(SRC_G.read_text())
    ccss, cbody, cjs = split_file(SRC_C.read_text())

    merged_vars = root_vars(gcss) + "\n" + root_vars(ccss)
    gcss_s = scope_css(strip_root(gcss), "pane-gate", "G-")
    ccss_s = scope_css(strip_root(ccss), "pane-cond", "C-")

    gbody_p = prefix_ids(gbody, "G-")
    cbody_p = prefix_ids(cbody, "C-")

    shim = """
  /* --- pane isolation -------------------------------------------------
     `document` is shadowed for the app inside this closure: ids get the
     pane prefix, and element queries are scoped to the pane's subtree.
     The app's own code needs no changes. */
  const __ROOT = window.document.getElementById("%(pane)s");
  const __PFX  = "%(pfx)s";
  const document = new Proxy(window.document, {
    get(t, k){
      if(k === "getElementById") return id => __ROOT.querySelector("#" + CSS.escape(__PFX + id));
      if(k === "querySelector")  return s  => __ROOT.querySelector(s);
      if(k === "querySelectorAll") return s => __ROOT.querySelectorAll(s);
      const v = t[k];
      return typeof v === "function" ? v.bind(t) : v;
    }
  });
"""

    tail_g = """
  // expose what the workbench shell needs
  window.NF_GATE = {
    S, refresh,
    initFrom(parsed, name){
      S.name = name; S.header = parsed.header; S.cols = parsed.cols; S.rows = parsed.rows;
      detect(); detectQC(); initThresholds();
      document.getElementById("sFile").textContent = name;
      document.getElementById("drop").classList.add("hide");
      document.getElementById("app").classList.add("ready");
      buildSidebar(); buildGates(); buildQC();
      recompute(); buildBiaxial(); refresh();
    },
    parseCSV,
    markers: () => S.markersOn || [],
  };
"""
    tail_c = """
  window.NF_COND = {
    S, refresh, parseCSV,
    initFrom(parsed, name){
      S.name = name; S.header = parsed.header; S.cols = parsed.cols; S.rows = parsed.rows;
      S.groups = []; S.assign = new Map(); S.active = 0; S.ref = 0;
      detect(); buildUI();
      const auto = document.getElementById("autoCol");
      if(auto && auto.value) autoGroup(auto.value);
      document.getElementById("sFile").textContent = name +
        (S.kind === "cell" ? "  (per-nucleus)" : "  (per-well)");
      document.getElementById("drop").classList.add("hide");
      document.getElementById("app").classList.add("ready");
      refresh();
    },
    rebuild(){                       // gates changed upstream: re-read the columns
      if(!S.rows) return;
      const keepSplit = S.split;
      buildWellRows();               // the exclusion set may have changed
      buildMetrics();
      buildMetricOptions();          // the LIST changes, not just the selection
      buildSiOptions(); buildSplitOptions();
      if(splitCandidates().some(c => c[0] === keepSplit)) S.split = keepSplit;
      const sp = document.getElementById("split"); if(sp) sp.value = S.split;
      refresh();
    },
  };
"""

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>nucleoflow workbench — gate, then compare</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Space+Grotesk:wght@400;500;700&family=JetBrains+Mono:wght@400;500;700&display=swap" rel="stylesheet">
<style>
:root{{{merged_vars}}}
*{{box-sizing:border-box;margin:0;padding:0}}
html,body{{height:100%}}
body{{background:var(--void); color:var(--ink); font-family:var(--sans);
  font-size:13px; line-height:1.45; overflow:hidden; -webkit-font-smoothing:antialiased}}
button,input,select{{font:inherit;color:inherit}}
:focus-visible{{outline:2px solid var(--c488); outline-offset:2px}}

/* ---- workbench shell ---- */
#wb{{display:grid; grid-template-rows:auto 1fr; height:100vh}}
#tabs{{display:flex; align-items:center; gap:2px; padding:0 14px;
  background:var(--void); border-bottom:1px solid var(--line)}}
#tabs .brand{{font-weight:700; font-size:14px; letter-spacing:-.02em;
  margin-right:16px; white-space:nowrap; padding:9px 0}}
#tabs .brand span{{color:var(--c488)}}
.tab{{background:none; border:0; border-bottom:2px solid transparent; cursor:pointer;
  padding:11px 15px 9px; font-size:12.5px; color:var(--dim); white-space:nowrap;
  display:flex; align-items:center; gap:7px}}
.tab:hover{{color:var(--ink)}}
.tab.on{{color:var(--ink); border-bottom-color:var(--c488)}}
.tab .num{{font-family:var(--mono); font-size:9.5px; color:var(--faint);
  border:1px solid var(--line); border-radius:3px; padding:0 4px}}
.tab.on .num{{color:var(--c488); border-color:var(--c488)}}
.tab:disabled{{opacity:.4; cursor:default}}
#tabs .sp{{flex:1}}
#wbFile{{font-family:var(--mono); font-size:11px; color:var(--dim); white-space:nowrap}}
#wbFile b{{color:var(--ink); font-weight:500}}
.wbtn{{background:var(--panel-2); border:1px solid var(--line); border-radius:3px;
  padding:4px 10px; font-size:11.5px; cursor:pointer; margin-left:8px}}
.wbtn:hover{{border-color:var(--faint)}}
#panes{{position:relative; min-height:0}}
.pane{{position:absolute; inset:0; overflow:hidden; display:none}}
.pane.on{{display:block}}

/* ---- shared loader ---- */
#wbDrop{{position:fixed; inset:0; display:flex; align-items:center; justify-content:center;
  flex-direction:column; gap:22px; background:var(--void); z-index:80; padding:24px}}
#wbDrop.hide{{display:none}}
.wbdz{{border:1px dashed var(--line); border-radius:6px; padding:44px 54px;
  text-align:center; max-width:660px; transition:border-color .15s, background .15s}}
.wbdz.hot{{border-color:var(--c488); background:rgba(63,209,107,.05)}}
.wbdz h1{{font-size:24px; font-weight:700; letter-spacing:-.025em; margin-bottom:8px}}
.wbdz h1 em{{font-style:normal; color:var(--c488)}}
.wbdz p{{color:var(--dim); font-size:12.5px; margin-bottom:20px; line-height:1.65}}
.wbdz code{{font-family:var(--mono); font-size:11.5px; color:var(--c561);
  background:var(--panel); padding:1px 5px; border-radius:2px}}
.steps{{display:flex; gap:0; justify-content:center; margin-bottom:26px;
  font-family:var(--mono); font-size:10.5px; color:var(--faint)}}
.steps i{{font-style:normal; padding:0 12px; position:relative}}
.steps i+i::before{{content:"→"; position:absolute; left:-4px; color:var(--line)}}
#wbErr{{color:var(--c640); font-family:var(--mono); font-size:11.5px;
  max-width:660px; text-align:center; min-height:17px}}
#wbBusy{{position:fixed; inset:0; background:rgba(11,16,23,.85); display:none;
  align-items:center; justify-content:center; z-index:90;
  font-family:var(--mono); font-size:12px; color:var(--c488)}}
#wbBusy.on{{display:flex}}
#wbNote{{position:fixed; left:50%; bottom:20px; transform:translate(-50%,14px);
  background:var(--panel-2); border:1px solid var(--c561); border-left-width:3px;
  border-radius:3px; padding:10px 14px; max-width:620px; z-index:85; opacity:0;
  pointer-events:none; transition:opacity .18s, transform .18s; font-size:12px; line-height:1.5}}
#wbNote.show{{opacity:1; transform:translate(-50%,0)}}
#wbNote b{{color:var(--c561)}}
::-webkit-scrollbar{{width:9px;height:9px}}
::-webkit-scrollbar-track{{background:var(--void)}}
::-webkit-scrollbar-thumb{{background:var(--line);border-radius:4px}}
@media (prefers-reduced-motion:reduce){{*{{transition:none!important}}}}

/* ================= gating pane ================= */
{gcss_s}

/* ================= condition pane ================= */
{ccss_s}
</style>
</head>
<body>

<div id="wbDrop">
  <div class="steps"><i>load once</i><i>gate the cells</i><i>compare conditions</i></div>
  <div class="wbdz" id="wbDz">
    <h1>nucleoflow <em>workbench</em></h1>
    <p>Drop a <code>per_nucleus_measurements.csv</code> here.<br>
       Set the gates on the first tab, then compare conditions and fit IC50s on
       the second — <br>both tabs work on the same cells, so a gate you move
       changes the analysis downstream.<br>
       A <code>per_well_summary.csv</code> also works, but only the second tab.</p>
    <button class="wbtn" id="wbPick" style="margin:0">Choose a file</button>
    <input type="file" id="wbFileInput" accept=".csv,.gz,.tsv,.txt" hidden>
  </div>
  <div id="wbErr"></div>
</div>
<div id="wbBusy">reading…</div>
<div id="wbNote" role="status" aria-live="polite"></div>

<div id="wb">
  <div id="tabs">
    <div class="brand">nucleoflow <span>workbench</span></div>
    <button class="tab on" id="tabGate"><span class="num">1</span> Gating</button>
    <button class="tab" id="tabCond"><span class="num">2</span> Conditions &amp; IC50</button>
    <div class="sp"></div>
    <div id="wbFile">—</div>
    <button class="wbtn" id="wbSync" title="push the current gates into the analysis tab">Apply gates →</button>
    <button class="wbtn" id="wbNew">New file</button>
  </div>
  <div id="panes">
    <div class="pane on" id="pane-gate">{gbody_p}</div>
    <div class="pane" id="pane-cond">{cbody_p}</div>
  </div>
</div>

<script src="https://cdn.jsdelivr.net/npm/pako@2.1.0/dist/pako.inflate.min.js"></script>
<script>
(function(){{{shim % {"pane": "pane-gate", "pfx": "G-"}}
{gjs}
{tail_g}
}})();
</script>
<script>
(function(){{{shim % {"pane": "pane-cond", "pfx": "C-"}}
{cjs}
{tail_c}
}})();
</script>
<script>
/* ============================================================
   Workbench shell: one file, two tabs, one set of cells.
   ============================================================ */
"use strict";
const WB = {{ parsed:null, name:"", kind:"", tab:"gate" }};
const $ = id => document.getElementById(id);
const busy = on => $("wbBusy").classList.toggle("on", on);
let noteTimer = null;
function note(html, ms = 6000){{
  const el = $("wbNote"); el.innerHTML = html; el.classList.add("show");
  clearTimeout(noteTimer); noteTimer = setTimeout(()=> el.classList.remove("show"), ms);
}}

function showTab(which){{
  WB.tab = which;
  $("tabGate").classList.toggle("on", which === "gate");
  $("tabCond").classList.toggle("on", which === "cond");
  $("pane-gate").classList.toggle("on", which === "gate");
  $("pane-cond").classList.toggle("on", which === "cond");
  // each pane sizes its canvases from clientWidth, which is 0 while hidden
  const app = which === "gate" ? window.NF_GATE : window.NF_COND;
  if(app && app.S && app.S.rows) requestAnimationFrame(()=> app.refresh());
}}

/** Write the gating tab's live calls back into the shared columns, so the
 *  analysis tab sees the gates you actually set rather than whatever the
 *  pipeline stored. This is the whole reason the two live in one file. */
function applyGates(quiet){{
  const G = window.NF_GATE, C = window.NF_COND;
  if(!G || !G.S.rows || !C || !C.S.rows) return;
  const S = G.S, act = S.active;
  if(!act || !S.markersOn) return;
  const cols = S.cols, header = S.header;
  const add = (name, arr) => {{
    cols[name] = arr;
    if(!header.includes(name)) header.push(name);
  }};
  const stat = S.stat;
  const on = S.markersOn;
  const dropped = [];
  // per-marker positivity
  const pos = on.map(m => {{
    const a = new Float64Array(S.rows);
    const v = cols[`${{stat}}_${{m}}`], thr = S.thr[m];
    for(let j=0;j<act.length;j++){{ const i = act[j]; a[i] = v[i] > thr ? 1 : 0; }}
    return a;
  }});
  on.forEach((m,k)=> add("pos_"+m, pos[k]));
  // A marker switched off in the gating tab must not linger here as a stale
  // column: it would still be offered as a way to split cells, using calls
  // that no longer reflect any gate on screen.
  for(const m of S.markers){{
    if(on.includes(m)) continue;
    const key = "pos_"+m;
    if(key in cols){{
      delete cols[key];
      const at = header.indexOf(key);
      if(at >= 0) header.splice(at, 1);
      dropped.push(m);
    }}
  }}
  // population / combination, over the ACTIVE markers only
  const POP = ["negative","single","double","triple","quadruple","quintuple"];
  const popCol = new Array(S.rows).fill("");
  const comboCol = new Array(S.rows).fill("");
  for(let j=0;j<act.length;j++){{
    const i = act[j];
    let n = 0, names = [];
    for(let k=0;k<on.length;k++) if(pos[k][i]){{ n++; names.push(on[k]); }}
    popCol[i] = POP[n] || (n+"-positive");
    comboCol[i] = names.length ? names.join("+")+"+" : "all-negative";
  }}
  add("population", popCol); add("combination", comboCol);
  // foci positivity, if the gating tab computed any
  for(const c of header.filter(h => h.startsWith("foci_n_"))){{
    const ch = c.slice(7), a = new Float64Array(S.rows), v = cols[c];
    for(let j=0;j<act.length;j++){{ const i = act[j]; a[i] = v[i] >= S.foci_min_count_hint ? 1 : 0; }}
  }}
  // The cell filter carries across. A cell the gating tab threw out as debris
  // or a dying nucleus must not sit in the analysis tab's denominator: it would
  // silently drag every population percentage below 100.
  if(S.qcKeep){{
    const d = new Float64Array(S.rows);
    for(let i=0;i<S.rows;i++) d[i] = S.qcKeep[i] ? 0 : 1;
    add("excluded_by_filter", d);
    C.S.excludeCol = "excluded_by_filter";
  }} else {{
    C.S.excludeCol = "";
    if("excluded_by_filter" in cols){{
      delete cols["excluded_by_filter"];
      const at = header.indexOf("excluded_by_filter");
      if(at >= 0) header.splice(at, 1);
    }}
  }}
  // cells outside the gating tab's active set are unclassified; drop them too
  // rather than leaving them to dilute the percentages
  {{
    const un = new Float64Array(S.rows);
    let anyUn = false;
    const inAct = new Uint8Array(S.rows);
    for(let j=0;j<act.length;j++) inAct[act[j]] = 1;
    for(let i=0;i<S.rows;i++) if(!inAct[i]){{ un[i] = 1; anyUn = true; }}
    if(anyUn){{
      add("excluded_by_filter", un);
      C.S.excludeCol = "excluded_by_filter";
    }}
  }}
  // the analysis tab hides readouts and splits for markers that are off
  C.S.hiddenMarkers = S.markers.filter(m => !on.includes(m));
  C.rebuild();
  if(!quiet) note(`<b>Gates applied.</b> ${{on.length}} marker${{on.length===1?"":"s"}} `+
    `(${{on.join(", ")}}) and the population calls were written into the analysis tab.` +
    (dropped.length ? ` <b>${{dropped.join(", ")}}</b> ${{dropped.length===1?"is":"are"}} `+
      `switched off in the gating tab, so ${{dropped.length===1?"it was":"they were"}} `+
      `removed from the analysis too.` : "") +
    ` Wells you excluded in the gating tab are still present here — exclude them again `+
    `on this tab if you meant to drop them.`, 9000);
}}

/* ---------- loading, once, for both panes ---------- */
async function loadFile(file){{
  $("wbErr").textContent = ""; busy(true);
  await new Promise(r=>setTimeout(r,30));
  try{{
    const buf = await file.arrayBuffer();
    const head = new Uint8Array(buf,0,2);
    const gz = file.name.toLowerCase().endsWith(".gz") || (head[0]===0x1f && head[1]===0x8b);
    let text;
    if(gz){{
      if(typeof pako === "undefined") throw new Error(
        "This file is gzipped and the decompressor didn't load (no internet?). "+
        "Use a plain .csv, or unzip it first.");
      text = pako.inflate(new Uint8Array(buf), {{to:"string"}});
    }} else text = new TextDecoder().decode(buf);

    // both panes share these columns, so parse once with the wider type
    const parsed = window.NF_COND.parseCSV(text);
    WB.parsed = parsed; WB.name = file.name;
    const perCell = parsed.header.some(c => /^(mean|median|integrated)_/.test(c))
                 && !parsed.header.includes("n_cells");
    WB.kind = perCell ? "cell" : "well";

    let gateOK = false;
    if(perCell){{
      try {{ window.NF_GATE.initFrom(parsed, file.name); gateOK = true; }}
      catch(e){{ console.warn("gating tab:", e); }}
    }}
    window.NF_COND.initFrom(parsed, file.name);

    $("tabGate").disabled = !gateOK;
    $("wbSync").style.display = gateOK ? "" : "none";
    $("wbFile").innerHTML = `<b>${{file.name}}</b> · ${{
      parsed.rows.toLocaleString()}} rows · ${{perCell ? "per-nucleus" : "per-well"}}`;
    $("wbDrop").classList.add("hide");
    showTab(gateOK ? "gate" : "cond");
    if(gateOK){{ applyGates(true); markerSig = (window.NF_GATE.S.markersOn||[]).join("|"); }}
    else note("<b>Loaded a per-well summary.</b> The gating tab needs per-nucleus "+
      "data, so it is disabled — the analysis tab has everything it needs.", 8000);
  }}catch(e){{
    $("wbErr").textContent = e.message || String(e);
    console.error(e);
  }}finally{{ busy(false); }}
}}

$("wbPick").onclick = () => $("wbFileInput").click();
$("wbFileInput").onchange = () => $("wbFileInput").files[0] && loadFile($("wbFileInput").files[0]);
["dragenter","dragover"].forEach(ev => document.addEventListener(ev, e=>{{
  e.preventDefault(); $("wbDz").classList.add("hot"); }}));
document.addEventListener("dragleave", e=>{{ if(e.relatedTarget===null) $("wbDz").classList.remove("hot"); }});
document.addEventListener("drop", e=>{{ e.preventDefault(); $("wbDz").classList.remove("hot");
  const f = e.dataTransfer.files[0]; if(f) loadFile(f); }});
// Switching a marker on or off changes what the analysis *is*, so it flows
// through immediately. Thresholds still need "Apply gates" — those you nudge
// repeatedly, and the analysis tab should not churn underneath you.
let markerSig = "";
function watchMarkers(){{
  const G = window.NF_GATE;
  if(!G || !G.S.rows) return;
  const sig = (G.S.markersOn || []).join("|");
  if(markerSig === "") {{ markerSig = sig; return; }}
  if(sig === markerSig) return;
  markerSig = sig;
  applyGates(true);
  const off = G.S.markers.filter(m => !(G.S.markersOn||[]).includes(m));
  note(`<b>Markers changed.</b> The analysis tab now uses ` +
    `${{(G.S.markersOn||[]).join(", ") || "no markers"}}` +
    (off.length ? `; ${{off.join(", ")}} ${{off.length===1?"is":"are"}} excluded from ` +
      `readouts, splits and the population classes.` : "."), 8000);
}}
$("pane-gate").addEventListener("click", e => {{
  if(e.target.closest("[data-m]") || e.target.closest("#G-allMark"))
    setTimeout(watchMarkers, 0);
}});

$("tabGate").onclick = () => showTab("gate");
$("tabCond").onclick = () => {{ showTab("cond"); }};
$("wbSync").onclick = () => applyGates(false);
$("wbNew").onclick = () => {{ $("wbDrop").classList.remove("hide");
  $("wbFileInput").value = ""; }};
</script>
</body>
</html>
"""
    OUT.write_text(html, encoding="utf-8")
    print(f"built {OUT} — {len(html):,} chars")


if __name__ == "__main__":
    main()
