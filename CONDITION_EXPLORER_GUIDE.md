# Condition Explorer — comparing replicate wells

`condition_explorer.html` groups replicate wells into conditions and compares
them properly. It's a companion to `gating_explorer.html`, not a replacement:
gate there, group and test here. Runs in the browser, nothing to install,
nothing leaves your computer.

---

## The one rule this tool is built around

**The unit of replication is the well, not the cell.**

Cells inside a well share a treatment, a staining, a dispense and an
illumination pattern. They are not independent measurements. If you feed half a
million cells into a t-test you will get p < 10⁻¹⁰ for a difference of no
consequence — that is pseudo-replication, and it is the most common statistical
error in high-content imaging.

So every number here has **n = number of wells**. The header says so, the tables
say so, and the tool warns you when n is too small to support the test you asked
for. If that makes your effect look less impressive than a per-cell test would,
the per-cell test was lying to you.

(Wells within a plate aren't fully independent either. For a proper multi-plate
experiment, the plate is arguably the unit — analyse each plate here, then
compare the plate-level means.)

---

## What to load

Any of these:

| File | Where from |
|---|---|
| `per_well_summary.csv` | the pipeline — the natural input |
| `*_regated.csv` | exported from the Gating Explorer, after you set gates by hand |
| `per_nucleus_measurements.csv` | also works; readouts are aggregated per well |

Trailing `#` comment lines are ignored, so the Gating Explorer's export loads
directly with its provenance intact.

---

## Grouping wells

**Automatically.** If the wells carry metadata — an OMERO key–value pair such as
`kv_Condition` or `kv_Compound` — pick it under *Auto-group by* and the
conditions build themselves. This is the usual case and takes one click.

Conditions are sorted naturally, so a dose series reads `1uM, 3uM, 10uM, 30uM`
rather than the alphabetical `1uM, 10uM, 30uM, 3uM`. The reference is set to
whatever looks like a control (DMSO, vehicle, untreated, 0) rather than to
whichever name sorts first.

**By hand.** Press *+ add*, then click or drag across wells on the plate to paint
them into the active condition. Click an assigned well again to remove it.
Rename a condition by typing over its label. There's no limit worth worrying
about — up to twelve conditions get distinct colours.

Wells left unassigned are excluded from every statistic, and the tool tells you
how many there are so it can't happen silently.

---

## Reading the two plots

### SuperPlot (left)

Every well is drawn as a dot. The bar is the mean, the whiskers are the 95%
confidence interval of the mean.

Dots rather than bars is deliberate. A bar chart with error bars hides the very
thing you need to see with n = 3–6: whether the conditions actually separate, or
whether one odd well is carrying the result. If two dots from one condition sit
inside the spread of another, no p-value will make that difference convincing.

### Estimation plot (right)

Each condition minus the reference, with a bootstrap 95% confidence interval, on
the units you actually measured. A dot is hollow when its interval includes zero.

This is the plot to lead with. A p-value only says "probably not exactly zero";
this says **how big** the difference is and **how uncertain** you are. A
difference of 0.4 percentage points can be highly significant and completely
irrelevant, and this plot makes that obvious where a p-value hides it.

---

## The statistics

**Per condition**: n wells, total cells, mean, SD, SEM, 95% CI, median, range.

**Omnibus** (three or more conditions), shown next to *Comparisons*:

| Test | When |
|---|---|
| Welch's ANOVA | default — does not assume equal variances |
| One-way ANOVA | assumes equal variance across conditions |
| Kruskal–Wallis | rank-based, assumes nothing about shape |

**Pairwise**: every pair, with the difference and its bootstrap CI, Hedges' *g*
and its CI, the test statistic, the raw p and the adjusted p.

Welch is the default rather than Student's t because treatments routinely change
the spread as well as the mean, and Welch costs almost nothing when the variances
do happen to match. Hedges' *g* is the standardised effect size with the
small-sample correction applied, so it's honest at n = 3.

**Multiple comparisons.** Six conditions is fifteen pairwise tests; at α = 0.05
you would expect a false positive by chance alone. Holm (family-wise) is the
default; Benjamini–Hochberg (FDR) is available when you're screening and can
tolerate some false positives; *None* exists but earns you a warning.

**Accuracy.** Every test was validated against `scipy.stats` — Welch's t,
Student's t, Mann–Whitney U, one-way ANOVA, Kruskal–Wallis, the normal and
incomplete-beta functions underneath them, Hedges' *g*, confidence intervals,
Holm and Benjamini–Hochberg. Worst relative disagreement across the whole suite:
**3.6 × 10⁻¹⁴**, i.e. floating-point noise.

---

## The warnings are the point

The tool interrupts when the design won't support the analysis:

- **n < 3 in a condition** — a spread cannot be estimated from two wells. The
  interval and p-value shown are not meaningful; treat those conditions as
  descriptive.
- **n = 3 throughout** — workable but weak. A 95% interval from three wells is
  very wide, and a non-significant result says very little. Lead with the effect
  size.
- **Unassigned wells** — fine if deliberate, worth checking if not.
- **Many comparisons, no correction** — read the adjusted column instead.

---

## Export

*Export report* writes a single CSV with four sections: the settings used, per-
condition descriptives (including which wells went into each), the omnibus test,
every pairwise comparison with effect sizes and adjusted p-values, and the
per-well values. It's a complete record of the analysis — enough for a methods
section and enough for someone else to reproduce your numbers.

---

## Honest limits

- **No mixed models.** Well-level means with equal weight per well. If your wells
  have very different cell counts, a weighted or hierarchical model would use the
  data better; this tool doesn't do that, and the cell counts are in the table so
  you can see whether it matters.
- **No plate-effect correction.** One plate at a time. If you have several,
  analyse them separately and compare the per-plate answers — consistent
  direction across plates is far stronger evidence than one pooled p-value.
- **Bootstrap CIs at n = 3–4 are approximate.** With so few wells the resampling
  can only be so informative. The interval is honest about being wide; don't read
  more into its exact endpoints than that.
- **It won't stop you fishing.** Regrouping wells until something turns
  significant will produce something significant. Decide the grouping from the
  plate layout before you look at the readout.
