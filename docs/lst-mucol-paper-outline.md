# Paper outline: pixelated double-layers for tracking at a muon collider

A proposed structure for a journal paper based on the LST-at-a-muon-collider results, including the post-conference additions (grouped vs equally-spaced layers, z-staggering, and the hit-resolution scan).

**Orientation.** Detector design. The algorithm is held fixed and treated as a probe; the tracker layout is the variable. The paper's contribution is not only "geometric track finding works at a muon collider" but "here is what a muon collider tracker must look like for it to work."

**Target.** 8–9 pages, two-column, 9 figures, 3 tables. JINST is the natural home given this orientation.

---

## Contents

**Front matter**
- Title options
- Abstract

**Paper body**

1. Introduction
   - 1.1 Motivation
   - 1.2 The tracking problem at a muon collider
   - 1.3 This work
2. Simulation and detector model
   - 2.1 Software and samples
   - 2.2 Baseline tracker geometry
   - 2.3 Geometry variants studied
   - 2.4 Scope and assumptions
3. Reconstruction algorithm
   - 3.1 Object hierarchy
   - 3.2 Selection requirements
   - 3.3 Computational considerations
4. Performance in the baseline geometry
   - 4.1 Efficiency definition
   - 4.2 Efficiency versus pT
   - 4.3 Efficiency versus η
5. Background rejection
   - 5.1 BIB yields through the chain
   - 5.2 Event display
   - 5.3 Comparison with timing-based rejection
6. Detector design studies
   - 6.1 Radial grouping: super-layers versus equal spacing
   - 6.2 Sensor staggering in z
   - 6.3 Single-hit position resolution
   - 6.4 Summary of design requirements
7. Discussion and outlook
8. Conclusions

**Back matter**
- Acknowledgements, data and code availability, references

**Working notes** (not part of the paper)
- Figure and table inventory
- Notes for the first draft

---

## Title options

1. *Pixelated double-layers for timing-independent tracking at a muon collider*
2. *Geometric background rejection at a muon collider with pixelated double-layers*
3. *Tracker design requirements for geometric track finding at a muon collider*

Option 1 balances the design message against the concrete thing we built. Option 3 promises broader design coverage than a barrel-only, single-algorithm study delivers; it is the right title only if §6 grows further.

---

## Abstract (~150–200 words)

Beats, in order:

1. Beam-induced background (BIB) makes conventional combinatorial track finding prohibitive at a muon collider; the standard mitigation is precision (30–60 ps) timing.
2. We study an alternative: a barrel tracker built entirely from *pixelated double-layers*, reconstructed with a purely geometric, LST-inspired algorithm.
3. In simulated 10 TeV collisions, double-layer hit filtering alone gives >10× BIB hit reduction — comparable to or better than a precision timing cut — and the full chain reduces BIB objects by O(10⁶) at ≈90% efficiency for muons with pT > 2 GeV.
4. Performance depends strongly on three layout choices: radial grouping of double-layers, sensor staggering in z, and single-hit resolution. We quantify each and state the resulting requirements.
5. Geometric rejection is therefore a viable alternative or complement to timing, and it constrains tracker layout in specific, stateable ways.

The design framing should be visible by sentence 4 of the abstract, not held back until §6.

---

## 1. Introduction (~1 page)

### 1.1 Motivation

- Muon collider physics case in two or three sentences, citing the recent design and physics reports. The machine's case does not need re-litigating here.
- Tracking is the entry point for essentially every downstream reconstructed object; if tracking is not solvable, nothing else is.

### 1.2 The tracking problem at a muon collider

- BIB origin: muon beam decays interacting with the shielding nozzles. Contrast with LHC pile-up — out-of-time, displaced origin, broad momentum spectrum, but a far lower crossing rate (30 kHz vs 40 MHz).
- Occupancy for a representative concept (MAIA): hits/cm² per layer, O(10⁷) tracker hits per event, ~10⁶ candidate tracks before cleaning.
- Why the combinatorial Kalman filter does not scale in this environment.
- Displaced tracking must survive whatever cleaning is applied — a 500 GeV b-hadron has βγcτ ≈ 5 cm. This motivates keeping IP-consistency requirements loose, and is worth stating even though this study does apply a loose IP constraint.

### 1.3 This work

- The idea: take the CMS Phase-2 Outer Tracker double-layer "pT module" concept, make it **pixelated** rather than strip-based, apply it across the whole barrel, and reconstruct with an LST-inspired hierarchy of geometric objects.
- Novelty claim: to our knowledge, the first study of double-layer-based, timing-independent track finding in a muon collider environment.
- Design framing: we hold the algorithm fixed and vary the detector layout, reporting what the layout must satisfy.
- Roadmap paragraph.

Scope limits belong here rather than in the conclusions — barrel only, single particles plus BIB, and no track fit, hence no momentum resolution or in-jet fake rate claim. Reviewers punish deferred caveats.

---

## 2. Simulation and detector model (~1.5 pages)

### 2.1 Software and samples

- Key4hep / DD4hep / Geant4 chain; MuonColliderSoft; pinned container image.
- Signal: single muons, flat in pT over the studied range and flat in η and φ, originating at the IP. State N.
- Background: BIB from the EU24 lattice at √s = 10 TeV. State N (≥100) and the assumed time window.
- Digitization: Gaussian smearing at fixed position resolution, 10 µm by default. State plainly that charge sharing, clustering and sensor inefficiency are not modeled.

### 2.2 Baseline tracker geometry

Presented as the proposed design, not merely as what happened to be simulated.

- Pixelated double-layers throughout: 16 detection layers arranged as 8 double-layers, split into an inner tracker (IT) and an outer tracker (OT). Give radii, the intra-double-layer gap, and the z-extent.
- Double-layers are **grouped radially** and sensors are **staggered in z**. Both are baseline properties; §6 quantifies what each buys.
- **Figure 1**: r–z view of the layout with an inset showing a single double-layer, with both the radial gap and the z-stagger visible.

### 2.3 Geometry variants studied

Three axes, defined once here and reused as labels throughout:

- **Grouped vs equally-spaced** double-layers.
- **Staggered vs unstaggered** in z. Staggered is the baseline; unstaggered appears only in §4.3 and §6.2, as the comparison that measures what staggering buys.
- **Single-hit resolution**: 0, 10, 20, 50 µm.
- **Table 1**: radii and z-extents for the grouped and equally-spaced layouts.

All results in §4–§6 use the staggered baseline unless a figure explicitly compares against the unstaggered geometry. Stating this once here removes the need to repeat it per figure.

### 2.4 Scope and assumptions

Short and bulleted:

- Barrel only; endcaps deferred.
- No precision timing assumed anywhere — this is the point of the study, not a limitation of it.
- Target is prompt particles with pT > 2 GeV.
- State exactly how the signal and BIB samples were combined or kept separate. Whether efficiency and BIB rates were measured in the same events matters, and a referee will ask.

---

## 3. Reconstruction algorithm (~1.25 pages)

Descriptive and geometric. It needs to be reimplementable, not exhaustive.

### 3.1 Object hierarchy

- **Mini-doublet (MD)**: a pair of hits in the two sensors of one double-layer, consistent with a track from the collision point.
- **T2**: two MDs. **T4**: two T2s. **T8**: an IT T4 combined with an OT T4 — 8 MDs, 16 hits.
- **Figure 2**: the object hierarchy drawn on an x–y event display.
- One sentence mapping MD/T2/T4/T8 onto CMS LST's segments, triplets and quintuplets, so LST readers are not confused by the naming.

### 3.2 Selection requirements

- MD level: Δz and Δr consistency with the IP.
- **Figure 3**: Δz for signal muons and for BIB in one OT double-layer, overlaid and normalized on shared axes. This is the single most persuasive plot in the paper and should not be split across panels.
- T2, T4 and T8 levels: Δz, Δr, and the xy and sz consistency (χ²-like) quantities, each defined with an equation.
- **Table 2**: requirements applied at each stage, with the number of objects expected per particle.
- State how the thresholds were chosen (e.g. as a quantile of the signal distribution) and note that the tuning is deliberately loose and unoptimized.

### 3.3 Computational considerations

- Why the hierarchy is intrinsically parallel and local: linking needs only nearest-neighbor information, with no global fit and no track propagation.
- Include at least one measured number here — even single-threaded wall-clock per event — or soften the "fast" claim to an explicit scaling argument. This is currently the largest open gap in the paper.

---

## 4. Performance in the baseline geometry (~1 page)

### 4.1 Efficiency definition

One definition, used everywhere: the **fraction of generated muons with pT > 2 GeV and within barrel acceptance for which at least one object of the given type is reconstructed**. Defined once, in one sentence. With a single definition in play, results sentences need no qualifier.

### 4.2 Efficiency versus pT

- **Figure 4**: efficiency vs pT for T4 and T8 in the staggered baseline.
- Quote the plateau values.

### 4.3 Efficiency versus η

- **Figure 5**: efficiency vs η for T4 and T8, with **two curves per panel — staggered baseline and unstaggered**.
- The staggered curve is flat; the unstaggered one dips at module boundaries. The figure then carries its own explanation, where a single flat curve would need a paragraph asserting that the flatness was designed in.
- Both curves are the efficiency of §4.1. No second definition is introduced anywhere in the paper.
- Quote the plateau difference here; §6.2 interprets it rather than re-deriving it.

---

## 5. Background rejection (~1 page)

### 5.1 BIB yields through the chain

- **Figure 6**: average BIB object yield per event at hits, MDs, T2, T4 and T8, on a log scale, with the timing-equivalent hit count overlaid as a dashed line.
- Numbers: ~3.8×10⁷ hits → ~2.5×10⁶ hits in MDs (>10× from geometry alone) → ~1.4×10⁵ in T2s → O(10²) in T4s → zero T8s observed.
- With N ≥ 100 BIB events, quote the zero-T8 result as a 95% CL upper limit of <3/N T8 per event and give the corresponding total rejection factor. This turns the paper's headline claim from an anecdote into a measurement.

### 5.2 Event display

- **Figure 7**: three panels — all hits, T2s, T4s — showing surviving BIB hits with signal muon hits overlaid. Visually compelling and worth the space at three panels; five is too many.

### 5.3 Comparison with timing-based rejection

- Read the comparison off the dashed line in Figure 6.
- Frame it as complementarity plus sensor-technology maturity: geometric double-layer filtering achieves comparable or better BIB hit reduction than a precision timing cut, using a more established technology, and the two approaches combine. Avoid any framing in which timing is unnecessary. A referee from the timing community will read this section closely, and one sentence explicitly welcoming the combination defuses most of the objection.

---

## 6. Detector design studies (~2.25 pages)

The centerpiece. Open with a framing paragraph: the algorithm is held fixed, each subsection varies one layout choice, and each reports the cost in efficiency, in background, or in both.

### 6.1 Radial grouping: super-layers versus equal spacing

- **Figure 8**: BIB yields through the chain for both layouts, side by side.
- Result: equal spacing leaves ~4× more objects at T2 and ~40× more at T4.
- Mechanism: grouped double-layers give short lever arms at the first linking step, so IP and curvature consistency requirements are tighter relative to the combinatorial background.
- Report the signal efficiency for both layouts as well. A background-only comparison invites the objection that grouping simply cuts harder. If the efficiencies match, that is the clean result and deserves an explicit sentence.

### 6.2 Sensor staggering in z

- Refer back to Figure 5 rather than repeating it — the measurement lives in §4.3 and the interpretation here.
- Mechanism, with a small schematic as an inset or single-column figure: without staggering, a particle crossing a module edge in one sensor tends to cross the aligned edge in the other, losing the whole MD, and the loss is localized in η at module boundaries.
- Quote the integrated cost (≈10 points at T8) and note that it is concentrated rather than uniform. A detector group cares that the loss is localized and therefore recoverable by layout alone, at no cost in material or channel count.
- Note the precedent in the CMS Phase-2 Outer Tracker, applied here to a pixelated, double-layer-everywhere tracker.

### 6.3 Single-hit position resolution

- **Figure 9**: efficiency and T4 BIB yield versus assumed resolution (0, 10, 20, 50 µm), as a two-panel or twin-axis figure so the trade-off is visible at once.
- Report signal loss and background leakage separately — they trade off against each other.
- State in the caption whether the selection windows were re-tuned per resolution point or held fixed. It changes how the plot must be read.
- Close with a requirement: resolution better than ~X µm is needed to retain O(10⁶) rejection.

### 6.4 Summary of design requirements

- **Table 3**: one row per design choice — double-layer gap, radial grouping, z-staggering, hit resolution, pixels vs strips — giving the requirement or recommendation and the quantitative basis, with a section or figure reference.
- This table is what a detector-concept working group will cite. In a design-oriented paper it is arguably the most valuable single object, and it costs about a third of a page.

---

## 7. Discussion and outlook (~0.5 page)

- What this does and does not establish: barrel only, single particles rather than jets or dense environments, and no track fit, so no claim is being made about momentum resolution or in-jet fake rates.
- Next steps: endcaps; full physics events; adding a fitting stage; a GPU implementation; joint optimization of layer radii and selection windows; combining with modest — not 30 ps — timing.
- Design implications as a short list, pointing back to Table 3.

## 8. Conclusions (~0.25 page)

Three or four sentences. No new information.

## Back matter

- Acknowledgements and funding.
- Data and code availability: container image, geometry description files, analysis code.
- References (~30–40).

---

## Figure and table inventory

| # | Content | Section | Status |
|---|---|---|---|
| Fig. 1 | Tracker layout, r–z, with double-layer inset showing gap and z-stagger | 2.2 | Adapt from existing |
| Fig. 2 | Object hierarchy (MD → T2 → T4 → T8) on x–y display | 3.1 | Adapt from existing |
| Fig. 3 | Δz for signal vs BIB, overlaid and normalized | 3.2 | Merge two existing panels |
| Fig. 4 | Efficiency vs pT, T4 and T8 | 4.2 | Adapt from existing |
| Fig. 5 | Efficiency vs η, staggered and unstaggered | 4.3 | **New** |
| Fig. 6 | BIB yields through the chain, with timing-equivalent line | 5.1 | Adapt from existing |
| Fig. 7 | Event display, three panels | 5.2 | Trim from existing five |
| Fig. 8 | BIB yields, grouped vs equally-spaced | 6.1 | Existing |
| Fig. 9 | Efficiency and BIB yield vs hit resolution | 6.3 | **New** |
| Tab. 1 | Layer radii and z-extents per layout | 2.3 | New, small |
| Tab. 2 | Requirements per stage, objects per particle | 3.2 | Existing |
| Tab. 3 | Summary of design requirements | 6.4 | **New** |

Figures 3, 5 and the §6.2 schematic are natural single-column candidates.

---

## Notes for the first draft

**Terminology.** Fix the object names once and use them consistently. MD / T2 / T4 / T8 are LST-adjacent but not identical to CMS LST's segments, triplets and quintuplets, so §3.1 should map between them explicitly. Avoid the phrase "line segment" for objects that are not segments. With a single efficiency definition, "efficiency" needs no qualifier anywhere.

**Geometry consistency.** Since both staggered and unstaggered results exist for everything, §4–§6 should all use the staggered baseline, with the unstaggered geometry appearing only where the comparison is the point. The sentence in §2.3 covers this for the reader.

**A possible addition, if cheap.** Unstaggered versions of the grouping comparison and the resolution scan also exist. They do not belong in the paper — four permutations would bury the message — but if staggering and grouping interact (for instance, if staggering matters more in the equally-spaced layout, where lever arms are longer), one sentence noting the interaction or its absence would preempt the obvious follow-up question at no cost in figures.

**Ranked list of what a referee is most likely to ask for.**

1. A quantitative timing or throughput measurement behind the "fast" claim (§3.3). This is now the largest gap.
2. Fake rate or purity for signal — how many spurious T4 or T8 objects are built per real muon.
3. Sensitivity to the BIB model itself, including lattice version and time window. One acknowledging sentence is likely enough.
4. Whether grouping improves efficiency or merely cuts harder (§6.1, addressed if both efficiencies are reported).
