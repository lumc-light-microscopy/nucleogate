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

## One factor or two

The **Design** selector at the top of the sidebar chooses between:

- **One factor** — conditions compared pairwise. This is the original mode and
  the default.
- **Two factors, crossed** — a factorial design, e.g. genotype × dose, or
  siRNA × treatment, analysed with a two-way ANOVA including the interaction.

Pick two factors and every well gets a level of each. The plate map splits each
well into two halves, factor A on the left and factor B on the right, so a
factorial layout is legible at a glance.

### Setting up a two-factor design

**Automatically** — choose a metadata column for each factor. That's usually all
it takes, since a factorial layout is normally already recorded that way.

**By hand** — each factor has the same controls as the one-factor condition
list: rename the factor by typing over its name, **+ level** to add a level,
type over a level to rename it, **×** to delete one, and click a level to make it
active. Well counts sit beside each level so you can see the design filling in.

Then click or drag wells on the plate. *Clicking a well sets* lets you paint both
factors at once, or only A, or only B — useful when one factor runs in columns
and the other in rows, which is the common layout.

You do not need any metadata for this. Switching to two-factor mode on a plate
with nothing to auto-group from seeds two levels per factor so the controls are
immediately usable; **clear** empties the assignments without losing the levels.

Adding a level guesses a sensible name by continuing whatever series is already
there — `0, 1, 2` offers `3`, and `0.1, 1, 10` offers `100` — but it is only a
placeholder, so type over it.

**swap A/B** exchanges the two. The interaction plot puts factor A on the x-axis,
so if one factor is ordered (dose, time) it usually reads better there.

### Reading the two-way output

**The interaction plot** is the one to look at first. Cell means are joined
across factor A, one line per level of B, with the individual wells shown
faintly behind. **Parallel lines mean no interaction.** Lines that converge,
diverge or cross mean the effect of one factor depends on the other — which is
exactly what the interaction term tests.

**The ANOVA table** gives the two main effects, the interaction, and the
residual, each with SS, df, MS, F and p.

**If the interaction is significant, be careful with the main effects.** "Dose
has an effect" is not a useful summary when the effect exists in the KO and not
the WT. Compare cells, not margins. The tool says so underneath the table when
this happens.

### Post-hoc comparisons

Under the ANOVA table you get pairwise comparisons, and **which set is shown
first depends on the interaction** — because that is what determines which ones
mean anything.

- **Interaction significant** → *simple effects* lead: levels of one factor
  compared inside each level of the other. This is the honest reading when the
  effect of a factor depends on the other one. On the example plate the genotype
  difference is `ns` at dose 0, then p = 0.012, 0.0007 and 5×10⁻⁷ as the dose
  climbs — the interaction made visible one row at a time.
- **Interaction not significant** → *marginal means* lead: each factor averaged
  over the other.

The other set is still there, one click away in a collapsed section, labelled
with a caution.

All comparisons use the **pooled error term** from the model (MSE and df are
printed underneath), which is the standard post-hoc: every well contributes to
the noise estimate rather than each pair being tested in isolation. Correction is
applied **within each family** of comparisons, not across all of them.

Validated against `scipy`: MSE, differences, t, p, confidence intervals and the
Holm adjustment all agree to within 5×10⁻¹⁵.

### Type II or Type III sums of squares

This only matters when the design is **unbalanced** — different numbers of wells
per cell, which happens the moment one well fails QC.

- **Balanced**: Type II and Type III are identical, and both match the textbook
  formulas. The selector makes no difference.
- **Unbalanced**: they differ, and the textbook shortcut formulas are simply
  wrong. Type II is the default and is the better choice when the interaction is
  not significant; Type III is what SPSS reports by default. Say which you used.

Both are computed here by comparing nested least-squares fits, which is correct
for balanced and unbalanced designs alike. For reference, on one of the test
designs the textbook shortcut disagreed with the correct answer by 30–45%.

**Empty cells.** If no well has a particular combination of levels, the
interaction cannot be estimated at all. The tool detects this, drops the
interaction term, fits main effects only, and tells you which cells are missing.

### Accuracy

The two-way ANOVA was validated on five designs — balanced 2×2, 3×2 and 2×4, and
unbalanced 2×2 and 3×3 — against two independent references: the exact textbook
formulas (valid for the balanced cases) and a numpy least-squares implementation
(valid for all of them). Every SS, df, F and p agreed to within **2 × 10⁻¹⁴**,
and Type II and Type III were confirmed identical on balanced designs and
different on unbalanced ones, as they should be.

## Dose–response and IC50

Choose **Dose–response** in the Design selector. The two factors take fixed
roles, stated in a banner above them and on their labels:

- **Factor A is the compound** (or condition, or cell line — whatever you want
  one curve per). Its levels can be any text.
- **Factor B is the concentration.** Its levels *must* contain a number: `10`,
  `10uM`, `10 µM`, `1e-3`, `0,5` (decimal comma) and `1,000` (thousands
  separator) all work, but `high` does not.

  **Every concentration level shows what it was read as** — `= 0.3 µM` beside
  the name, or a red `not a number` if nothing could be parsed — and the banner
  lists the whole series. If a value isn't being picked up, you can see it
  immediately rather than discovering it in the fit.

  Level names must be unique. Typing a name another level already has is
  refused with a message saying so, and your text stays in the box so you can
  correct it.

Assigning wells works exactly as it does in two-factor mode: each well is drawn
split in half — compound on the left, concentration on the right — and clicking
or dragging paints the active level of each. All the level controls (add, rename,
delete) are there too, so a plate with no metadata can be set up entirely by hand.

If the two end up the wrong way round, the banner says so and tells you to press
**swap A/B** — a common slip, since the columns are often in the other order in
the metadata.

**Concentration unit.** Pick one from the dropdown (nM, µM, mM, M, µg/mL, mg/mL,
ng/mL, %, Gy, or none). It is guessed from the column name where possible, so
`kv_Dose_uM` selects µM by itself, and `kv_Dose_nM` selects nM. The unit then
appears on the plot axis, in the IC50 and confidence-interval column headers, in
the selectivity and curve-comparison tables, and in the exported CSV column names
(`ic50_uM`, `conc_min_uM`). Changing it relabels everything at once — it does
**not** convert your numbers, so use the unit your concentrations are actually in.

A four-parameter logistic (Hill) curve is fitted per compound:

```
Y = Bottom + (Top − Bottom) / (1 + 10^((log10 C − log10 IC50) · Hill))
```

You get IC50 (or EC50 — the direction is detected from the fitted asymptotes,
not assumed), a 95% confidence interval, the Hill slope, Top, Bottom, the span
and R², plus a curve plot on a log axis with every well shown.

**Vehicle wells** (concentration 0) cannot sit on a log axis, so they are drawn
in a separate band at the left and used for normalisation, but excluded from the
fit. That is standard practice, not a limitation.

**Normalisation**: raw, % of vehicle, or % inhibition (vehicle = 0%, largest
observed effect = 100%). **Constraints**: free, or hold Top at 100 and/or Bottom
at 0 — useful for normalised data where those plateaus are known by definition
and fitting them wastes degrees of freedom.

### The three ways an IC50 misleads you

The tool flags all three in the table rather than quietly reporting a number:

- **`extrapolated`** — the IC50 falls outside the doses you actually tested. It
  is a property of the fitted curve, not a measurement. Report "> highest dose".
- **`R² low`** — the curve doesn't describe the data.
- **`Hill n`** — an implausibly steep slope, which almost always means an
  inactive compound whose curve is unconstrained. In the example plate, an inert
  compound returned an IC50 of 38 µM with Hill = 18.9 and R² = 0.34: all three
  flags fire, and the correct report is "inactive", not "IC50 = 38 µM".

**One more thing worth knowing.** The IC50 here is the *relative* IC50 — the dose
halfway between the fitted Top and Bottom. If a curve never reaches a lower
plateau, that is not the same as half the vehicle response (the *absolute*
IC50), and the two can differ several-fold. Most publications report the
relative one; say which you mean.

### Comparing cell populations inside the same well

A well average hides the possibility that two cell populations in it responded
completely differently. **Split cells by** fits a separate curve to each
population and tests whether their IC50s differ.

This needs a **per-nucleus file** — a per-well summary has already averaged the
populations together, and they cannot be separated again afterwards. The control
is disabled with an explanation if you loaded a summary.

You can split on any positivity call (`pos_GFP`, `foci_pos_gH2AX`), on `dead`,
or on the class labels — **including double- and triple-positive cells**. The
options name their classes rather than counting them, so you can see what you
are getting:

```
population: negative / single / double
marker combination: all-negative / GFP+ / GFP+RFP+
```

**Population percentages always sum to 100%.** Every cell in a well belongs to
exactly one class, and any cell excluded upstream leaves the denominator as well
as the numerator. If you ever see the classes summing to less than 100, that is
a bug — please report it.

**The classes on offer can never exceed the markers in play.** "Triple" is not
possible with two markers, so it is not offered — even if a stale `population`
column in the file claims otherwise. The ceiling is derived from the number of
active `pos_` columns rather than from the strings in the column, so a file
carrying leftover class labels from an earlier analysis cannot reintroduce them.

With three markers you therefore get `negative / single / double / triple`, and
`% triple` also appears as a readout you can fit a curve to directly.

`population` splits by *how many* markers a cell is positive for, which is the
usual way to ask "do double-positive cells respond differently?". `combination`
splits by *which* markers, if the identity matters. Levels are ordered
negative → single → double → triple so the curves read as a progression. Each level of the split gets
its own curve: solid, dashed and dotted lines distinguish populations, colour
distinguishes compounds, and both legends are drawn.

**Choose the readout carefully.** *Cells per well* is usually what you want —
"how many of this population survive at each dose" — and *% of well in this
population* is the natural companion. A marker percentage measured *within* the
population you split on is circular, so the tool warns when the combination looks
odd.

The comparison table then leads with the question you actually asked: **within
each compound, one population against the other**, followed by the same
population across compounds.

**Why it matters, with numbers.** On a test dataset where GFP-positive cells were
made ten times more sensitive than GFP-negative ones:

```
                          IC50      95% CI
Cisplatin  GFP-negative   2.98      2.45 – 3.62     (truth 3.0)
Cisplatin  GFP-positive   0.270     0.236 – 0.309   (truth 0.3)
           -> 11.0x apart, p = 2.1e-9

Etoposide  GFP-negative   0.990     0.820 – 1.20    (truth 1.0)
Etoposide  GFP-positive   1.08      0.944 – 1.24    (truth 1.0)
           -> 1.09x apart, p = 0.62, ns
```

Without the split, Cisplatin reports a single IC50 of **1.83** — a number that
describes neither population and would be reported as if it described both.

### Comparing two curves

"Is compound A more potent than B?" is not answered by checking whether their
confidence intervals overlap. That comparison is conservative and has no stated
error rate — two genuinely different IC50s can have overlapping intervals.

The **Curve comparison** table under the dose-response plot does it properly.
Each pair of compounds is fitted twice: once with its own IC50 for each, once
forced to share a single IC50. An extra sum-of-squares F test then asks whether
forcing them together costs a significant amount of fit. You get the fold
difference in potency, F, the raw p and the p adjusted across all pairs by
whichever correction you chose in the sidebar.

A compound whose own curve is flagged (extrapolated, low R², extreme Hill) will
usually come back "ns" against everything, which is the honest answer — an
unconstrained curve provides no evidence that its IC50 differs from anything.

### Selectivity index

Pick a second readout under **Selectivity vs** and the same curves are fitted to
it, giving SI = IC50(secondary) ÷ IC50(primary) per compound.

The usual use in toxicology is potency of the effect you want versus the dose
that kills the cells: set the effect as your primary readout and viability as
the secondary, and an SI of 10 means you get the effect at a tenth of the toxic
dose. Both readouts must come from the same wells or the ratio means nothing —
which they do here, since it's the same plate.

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

**Omnibus** (three or more conditions, one-factor mode), shown next to
*Comparisons*:

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

## Assay quality: Z′-factor

Pick a **negative** and a **positive** control condition in the sidebar and the
*Assay quality* card reports:

| Metric | What it tells you |
|---|---|
| **Z′-factor** | `1 − 3(σpos + σneg) / |μpos − μneg|`. Above 0.5 is the usual screening acceptance threshold; 0–0.5 is marginal; below 0 means the controls don't separate |
| **SSMD** | separation in units of combined SD; \|SSMD\| ≥ 3 is strong |
| **Signal window** | the raw difference between control means |
| **Signal / background** | their ratio |
| **CV per control** | under ~10% is tight, over ~20% is noisy |

The important thing about Z′ is that it is a property of **the assay, not your
treatment**. It answers "could this plate have detected an effect at all?" —
which is worth knowing *before* you interpret one. A plate with Z′ = 0.1 and a
non-significant treatment effect has told you nothing either way.

It is computed on well means, like everything else here, and it is
readout-specific: a plate can have an excellent Z′ on viability and a poor one on
a marker, so check the readout you are actually reporting.

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

- **Dose–response fits are unweighted.** Every well counts equally. If your
  variance grows with the response, a weighted fit would be better; export the
  values and use `drc` in R or `scipy.optimize.curve_fit` with `sigma`.
- **Confidence intervals on IC50 come from the fit's covariance matrix.** They
  are symmetric on the log scale and therefore asymmetric in concentration,
  which is correct, but they assume the model is right and the noise is normal.
  At five or six doses treat them as indicative.
- **Two-way ANOVA only; no three-way, no nesting.** Two crossed factors is what
  plate layouts usually are. If you need a third factor or a random effect for
  plate, export the per-well values and fit the model in R or statsmodels.
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
