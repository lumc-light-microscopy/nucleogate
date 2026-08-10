# nucleoflow workbench — gate, then compare, in one file

`nucleoflow_workbench.html` puts the two explorers in one page as two tabs over
**one shared set of cells**. You load a file once, set the gates on tab 1, and
tab 2 analyses the cells *as you gated them* rather than as the pipeline stored
them.

It runs in the browser, needs no install, and nothing leaves your computer.

---

## Why combine them

Splitting cells into populations for a dose-response — GFP-positive versus
GFP-negative, say — needs a positivity call for every cell. On its own, the
condition explorer can only use whatever `pos_*` columns the pipeline happened
to write, and offers no way to change them. If those gates are wrong, or the
file has none at all, you are stuck.

In the workbench, tab 1 *creates* those columns from live gates:

```
tab 1: drag the GFP gate          tab 2: IC50 per population
       threshold 282  → 26.8% +          GFP-negative  3.04 µM
                                         GFP-positive  0.294 µM   (10.3x apart)

       threshold 2000 → 13.2% +          GFP-negative  2.38 µM
                                         GFP-positive  0.302 µM
```

The gate is now an analysis choice you can see the consequences of, rather than
a fixed property of the file.

---

## The workflow

**1 — Load.** Drop a `per_nucleus_measurements.csv`. Both tabs read the same
parsed columns. A `per_well_summary.csv` also loads, but the gating tab is
disabled: a summary has already averaged the cells together and cannot be
re-gated.

**2 — Gating tab.** Everything the standalone gating explorer does: the cell
filter (DAPI and morphology windows), marker gates you drag, biaxial plots,
per-well populations. Get the gates where you want them.

**3 — Apply gates →.** Writes the current calls into the shared columns:
`pos_<marker>` for each active marker, plus `population` and `combination`
computed over the markers you left switched on. This happens automatically when
the file first loads, and any time you press the button.

**4 — Conditions & IC50 tab.** Group replicate wells, run two-way ANOVA, or fit
dose-response curves — now with **Split cells by** offering the `pos_*` columns
the gating tab just made.

Go back to tab 1, move a gate, press **Apply gates** again, and tab 2 updates.

---

## What carries across, and what does not

**Carries across:** marker positivity, population and combination classes, and
which markers were switched on.

**Switching a marker off in tab 1 flows through immediately — no button press.**
That is a structural choice, not a fine-tuning knob, so it does not wait for
*Apply gates*. The marker disappears from tab 2 completely: its `pos_` column is
deleted, its intensity readouts leave the readout list, it leaves the split
options, and the population classes are recomputed over the markers that remain.
A notice says which markers were dropped.

Threshold changes still need *Apply gates*, deliberately — those you nudge
repeatedly, and the analysis tab should not churn underneath you while you are
reading it.

With three markers active you get `negative / single / double / triple` as split
levels, and `% triple` as a readout in its own right. Switch one marker off and
those become `negative / single / double`.

That matters more than it sounds. A channel with no real positive population
gets split down the middle by an automatic threshold, inventing a
"double-positive" class out of noise. Switching it off in tab 1 and re-applying
collapses the analysis back to the populations that exist:

```
Cy5 on  (no real signal):  negative / single / double / triple
Cy5 off:                   negative / single / double
```

**Does not carry across: well exclusions.** Wells you exclude in the gating tab
stay present in the analysis tab. The two tabs mean different things by
excluding a well — one is "don't show me these cells", the other is "this well
is not part of the design" — and silently coupling them would be worse than
making you say it twice. The notice after applying gates reminds you.

**Carries across: the cell filter.** Cells removed by the DAPI/morphology filter
on tab 1 — debris, clumps, dying nuclei — are dropped from the analysis tab
entirely, not just from the classification. The cell counts on tab 2 are the
filtered counts, and the header shows how many cells were removed.

This has to be all-or-nothing. An earlier version excluded those cells from the
population classes but left them in the denominator, so `negative + single +
double` came to 96% instead of 100% and the missing 4% had no visible
explanation. A cell removed by a filter has to leave the numerator *and* the
denominator.

---

## Standalone versus workbench

The two standalone files still exist and are unchanged:

| File | Use it when |
|---|---|
| `gating_explorer.html` | you only want to set gates and look at populations |
| `condition_explorer.html` | your file already has the gates you want, or you're working from a `per_well_summary.csv` |
| `nucleoflow_workbench.html` | you want to gate *and* analyse, especially splitting IC50s by cell population |

The workbench contains both apps in full — nothing is cut down. If you're
handing a single file to a student, hand them this one.

---

## Honest limits

- **One file at a time.** Both tabs share the one dataset. For several plates,
  run them one at a time, or split per plate with the Plate selector.
- **Big per-nucleus files are slower here** than in the standalone gating tool,
  because both apps hold the same columns and both redraw on tab switch. Around
  a million nuclei is still workable; beyond that, gate in the standalone tool,
  export, and analyse from the export.
- **Applying gates is a deliberate action.** It runs once on load and then only
  when you press the button, so tab 2 never shifts under you while you are
  reading it. The flip side is that you must press it after changing a gate.
- **The statistics are the same code** as the standalone condition explorer,
  validated against `scipy` — see `CONDITION_EXPLORER_GUIDE.md`, which applies
  in full to tab 2.
