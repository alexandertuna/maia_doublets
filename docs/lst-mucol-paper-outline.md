# Paper outline v3: pixelated double-layers for tracking at a muon collider

**Revisions in this version**
- Reweighted toward detector design: §6 is now the centerpiece (~2.25 pages), §3 compressed (~1.25).
- Staggering is a **property of the baseline geometry**, and the unstaggered geometry is **retained as a comparison** — see the correction note below. Only one efficiency definition (overall) is used anywhere.
- BIB statistics caveat removed on the assumption of ≥100 BIB events; §5.1 now carries a quantitative upper limit instead of a warning.
- Added Table 3, a summary of design requirements — the citable artifact of a design-oriented paper.

**Correction to v2.** v2 treated "drop the unstaggered result" and "use one efficiency definition" as the same decision. They are independent. What forces two definitions is *algorithmic* efficiency, whose denominator depends on which hits a muon happened to leave. Overall efficiency has a fixed denominator and can be plotted for any geometry. So the unstaggered curve can stay, with one definition throughout — it costs one extra line on a figure you are already making, and it turns staggering from an assertion into a measurement.

**Working title options** (design-oriented)

1. *Tracker design requirements for geometric track finding at a muon collider*
2. *Pixelated double-layers for timing-independent tracking at a muon collider*
3. *Geometric background rejection at a muon collider with pixelated double-layers*

Option 2 balances the design message with the concrete thing you built. Option 1 promises broader design coverage than a barrel-only, single-algorithm study delivers — reach for it only if §6 grows further.

**Target length**: 8–9 pages two-column, 9 figures, 3 tables.

---

## Abstract (~150–200 words)

1. BIB makes conventional combinatorial track finding prohibitive at a muon collider; the standard mitigation is precision (30–60 ps) timing.
2. We study an alternative: a barrel tracker built entirely from *pixelated double-layers*, reconstructed with a purely geometric, LST-inspired algorithm.
3. In simulated 10 TeV collisions, double-layer hit filtering alone gives >10× BIB hit reduction — comparable to or better than a precision timing cut — and the full chain reduces BIB objects by O(10⁶) at ≈90% efficiency for muons with pT > 2 GeV.
4. Performance depends strongly on three layout choices: radial grouping of double-layers, sensor staggering in z, and single-hit resolution. We quantify each and state the resulting requirements.
5. Conclusion: geometric rejection is a viable alternative or complement to timing, and it constrains tracker layout in specific, stateable ways.

The design framing should be visible by sentence 4 of the abstract, not just in §6.

---

## 1. Introduction (~1 page)

### 1.1 Motivation
- Muon collider physics case in two or three sentences; cite recent design/physics reports. Do not re-litigate the machine.
- Tracking is the entry point for every downstream object.

### 1.2 The tracking problem at a muon collider
- BIB origin: beam decays interacting with the shielding nozzles. Contrast with LHC pile-up: out-of-time, displaced origin, broad momentum spectrum, but far lower crossing rate (30 kHz vs 40 MHz).
- Occupancy for a representative concept (MAIA): hits/cm² per layer; O(10⁷) tracker hits/event; ~10⁶ candidate tracks before cleaning.
- Why the combinatorial Kalman filter does not scale.
- Displaced tracking must survive the cleaning (500 GeV b-hadron, βγcτ ≈ 5 cm) — motivates keeping IP requirements loose.

### 1.3 This work
- The idea: the CMS Phase-2 Outer Tracker double-layer "pT module" concept, but **pixelated** and applied across the whole barrel, reconstructed with an LST-inspired hierarchy.
- Novelty claim: first study of double-layer-based, timing-independent track finding in a muon collider environment.
- **Design framing sentence** — the load-bearing addition for this version: *we treat the algorithm as a fixed probe and the detector layout as the variable*, and report what the layout must satisfy.
- Roadmap paragraph.

State scope limits here, not in the conclusions: barrel only, single particles plus BIB, no track fit and therefore no momentum resolution or jet-environment fake rate.

---

## 2. Simulation and detector model (~1.5 pages)

### 2.1 Software and samples
- Key4hep / DD4hep / Geant4; MuonColliderSoft; pinned container.
- Signal: single muons, flat in pT and η/φ, from the IP; state N.
- Background: BIB, EU24 lattice, √s = 10 TeV; **state N ≥ 100** and the assumed time window.
- Digitization: Gaussian smearing at fixed resolution, default 10 µm. State plainly that charge sharing, clustering and sensor inefficiency are not modeled.

### 2.2 Baseline tracker geometry
This is now the paper's proposed design, so present it as such rather than as "what we happened to simulate."
- Pixelated double-layers throughout: 16 detection layers as 8 double-layers, split into inner (IT) and outer (OT) trackers. Radii, intra-double-layer gap, z-extent.
- **Double-layers are grouped radially** and **sensors are staggered in z**. Both are baseline properties; §6 shows what each buys.
- Figure 1: r–z view with a double-layer inset showing both the radial gap and the z-stagger. The stagger must be visible here, since it is no longer given its own figure.

### 2.3 Geometry variants studied
Three axes, defined once and reused as labels throughout:
- **Grouped vs equally-spaced** double-layers.
- **Staggered vs unstaggered** in z. The staggered layout is the baseline; the unstaggered one appears only in §4.3 and §6.2, as the comparison that measures what staggering buys.
- **Single-hit resolution**: 0, 10, 20, 50 µm.
- Table 1: radii and z-extents for the grouped and equally-spaced layouts.

All results in §4–§6 use the staggered baseline unless a figure explicitly compares against the unstaggered geometry. Say this once, here, and it does not need repeating per figure.

### 2.4 Scope and assumptions
Bulleted: barrel only; no precision timing assumed anywhere (the point, not a limitation); prompt particles with pT > 2 GeV; state exactly how signal and BIB samples were combined or kept separate.

---

## 3. Reconstruction algorithm (~1.25 pages)

Compressed relative to v1. It needs to be reimplementable, not exhaustive.

### 3.1 Object hierarchy
- **Mini-doublet (MD)**: hit pair across one double-layer, consistent with a track from the collision point.
- **T2** = two MDs; **T4** = two T2s; **T8** = IT T4 + OT T4 (8 MDs, 16 hits).
- Figure 2: hierarchy on an x–y event display. Fold the Δz separation plot in as a second panel — see Figure 3 below.
- One sentence mapping MD/T2/T4/T8 onto CMS LST's segments/triplets/quintuplets, so LST readers are not confused.

### 3.2 Selection requirements
- MD level: Δz and Δr consistency with the IP.
- Figure 3: Δz for signal muons vs BIB in one OT double-layer, overlaid and normalized. The single most persuasive plot in the paper; merge the talk's two panels.
- T2/T4/T8: Δz, Δr, xy and sz consistency (χ²-like) quantities, each with an equation.
- Table 2: requirements per stage with objects expected per particle (talk slide 24, near-verbatim).
- State how thresholds were chosen and that they are deliberately unoptimized.

### 3.3 Computational considerations
- Intrinsically parallel and local: nearest-neighbor linking, no global fit, no propagation.
- Include at least one measured number, even single-threaded wall-clock per event, or soften "fast" to an explicit scaling argument. Unchanged from v1 — still the second-most-likely referee question.

---

## 4. Performance in the baseline geometry (~1 page)

### 4.1 Efficiency definition
One definition only: **fraction of generated muons (pT > 2 GeV, within barrel acceptance) with at least one reconstructed object of the given type**. Defined once, in one sentence, used everywhere. No qualifier needed in results sentences.

### 4.2 Efficiency vs pT
- Figure 4: efficiency vs pT for T4 and T8, staggered baseline.
- Quote plateau values.

### 4.3 Efficiency vs η
- Figure 5: efficiency vs η for T4 and T8, **two curves per panel — staggered baseline and unstaggered**.
- The staggered curve is flat; the unstaggered one dips at module boundaries. The figure carries its own explanation, where a lone flat curve would have needed a paragraph asserting that flatness was designed in.
- Both curves are overall efficiency. No second definition is introduced anywhere in the paper.
- Quote the plateau difference in the text; §6.2 then interprets it rather than re-deriving it.

---

## 5. Background rejection (~1 page)

### 5.1 BIB yields through the chain
- Figure 6: average BIB object yield per event at hits / MDs / T2 / T4 / T8, log scale, with the **timing-equivalent hit count as a dashed line** (folding v1's separate Fig. 8 into this one).
- Numbers: ~3.8×10⁷ hits → ~2.5×10⁶ in MDs (>10× from geometry alone) → ~1.4×10⁵ in T2s → O(10²) in T4s → zero T8s.
- With N ≥ 100 BIB events, quote the zero-T8 result as a **95% CL upper limit of <3/N T8 per event**, and give the corresponding total rejection factor. This converts the paper's headline claim from an anecdote into a measurement.

### 5.2 Event display
- Figure 7: three panels (all hits / T2 / T4) with signal muon hits overlaid. Trimmed from five for space.

### 5.3 Comparison with timing-based rejection
- Read off the dashed line in Figure 6.
- Frame as complementarity plus sensor-technology maturity, not as timing being unnecessary. A timing-community referee will read this section closely; a sentence explicitly welcoming the combination costs nothing and defuses most of the objection.

---

## 6. Detector design studies (~2.25 pages) — **centerpiece**

Open with a framing paragraph: the algorithm is held fixed; each subsection varies one layout choice and reports the cost in efficiency, background, or both.

### 6.1 Radial grouping: super-layers vs equal spacing
- Figure 8: BIB yields through the chain for both layouts, side by side.
- Result: equal spacing leaves ~4× more objects at T2 and ~40× more at T4.
- Mechanism: grouped double-layers give short lever arms at the first linking step, so IP/curvature consistency requirements are tighter relative to the combinatorics.
- **Report signal efficiency for both layouts too.** Background-only comparison invites the objection that grouping simply cuts harder. If efficiencies match, that is the clean result and should be stated in one sentence.

### 6.2 Sensor staggering in z
Short subsection, but now backed by a figure rather than an assertion.
- Refer back to Figure 5 rather than repeating it — the measurement lives in §4.3, the interpretation here.
- Mechanism, with the small schematic as an inset or single-column figure: without staggering, a particle crossing a module edge in one sensor crosses the aligned edge in the other, losing the whole MD; the loss is localized in η at module boundaries.
- Quote the integrated cost (≈10 points at T8) and note that it is concentrated rather than uniform — a detector group cares that the loss is localized and therefore recoverable by layout alone, at no cost in material or channel count.
- Note the precedent in the CMS Phase-2 Outer Tracker, applied here to a pixelated double-layer-everywhere tracker.

### 6.3 Single-hit position resolution
- Figure 9: efficiency and BIB yield at T4 vs assumed resolution (0, 10, 20, 50 µm), ideally as a two-panel or twin-axis figure so the trade-off is visible at once.
- Report signal loss and background leakage separately.
- State in the caption whether selection windows were re-tuned per resolution point or held fixed — it changes how the plot must be read.
- Output a requirement: *resolution better than ~X µm is needed to retain O(10⁶) rejection.*

### 6.4 Summary of design requirements
- **Table 3**: one row per design choice — double-layer gap, radial grouping, z-staggering, hit resolution, pixel vs strip — with the requirement or recommendation and the quantitative basis (figure/section reference).
- This table is what a detector-concept working group cites. In a design-oriented paper it is arguably the single most valuable object, and it costs a third of a page.

---

## 7. Discussion and outlook (~0.5 page)

- What is and is not established: barrel only; single particles, not jets; no fit, so no momentum resolution or in-jet fake rate.
- Next steps: endcaps; full physics events; a fitting stage; GPU implementation; joint optimization of radii and selection windows; combining with modest (not 30 ps) timing.
- Design implications as a short list, pointing back at Table 3.

## 8. Conclusions (~0.25 page)

Three or four sentences, no new information.

## Back matter
Acknowledgements and funding; data/code availability (container image, geometry description, analysis code); ~30–40 references.

---

## Notes for the first draft

**On the unstaggered result.** Keeping one unstaggered curve on Figure 5 costs almost nothing and buys the third row of Table 3. The simplification you wanted — dropping algorithmic efficiency — is fully preserved: every number in the paper is overall efficiency, with one denominator, defined once in §4.1. What made the old presentation confusing was the two definitions, not the two geometries.

**Consistency, now resolved.** Since both geometries were run for every result, §4–§6 can all use the staggered baseline, with the unstaggered geometry appearing only where it is the point of the comparison. Worth one explicit sentence in §2.3 (added above) so no reader has to infer it per figure.

**A possible bonus, if it is cheap.** You also have the unstaggered version of the grouping comparison and the resolution scan. Those do not belong in the paper — four permutations would bury the message. But if staggering and grouping interact at all (e.g. staggering matters more in the equally-spaced layout, where lever arms are longer), a single sentence noting the interaction, or its absence, would preempt the obvious follow-up question and costs no figure.

**Figure and table count.** 9 figures, 3 tables, from 10 and 2 in v1 — unchanged from v2, since the unstaggered result returns as an extra curve rather than an extra figure. Comfortable for 8–9 pages; Figures 3, 5 and the §6.2 schematic are natural single-column candidates.

**Terminology.** Fix MD / T2 / T4 / T8 once; avoid "line segment" for objects that are not segments; "efficiency" now needs no qualifier, having exactly one meaning.

**Remaining referee-bait, ranked.**
1. A quantitative timing/throughput number (§3.3).
2. Fake rate or purity — how many spurious T4/T8 per real muon.
3. Sensitivity to the BIB model itself (lattice version, time window) — one acknowledging sentence is enough.
4. Whether grouping helps efficiency or only cuts harder (§6.1, addressed if you report both).

Item 1 is now the largest open gap, since the BIB statistics one closes with the additional events.
