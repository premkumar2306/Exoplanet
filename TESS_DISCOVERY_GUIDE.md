# TESS Exoplanet Discovery Guide for Citizen Scientists

**Goal**: Find real exoplanet candidates using public NASA TESS data with paper-level rigor.

**Requirements**: Python 3.8+, laptop, internet, patience.

---

## PART 1 — STRATEGIC RATIONALE

### Why TESS Light Curves Work for Citizens

**1. Data Quality and Accessibility**
- TESS observes 200,000+ stars per sector with 2-minute or 10-minute cadence
- All data is public via MAST archive within weeks of collection
- Light curves are pre-processed (flat-fielded, calibrated) by SPOC pipeline
- No telescope required - NASA did the hard work

**2. Where Automated Pipelines Fail**

The TESS SPOC pipeline finds ~90% of strong signals, but misses:

- **Long-period planets** (P > 20 days): Fewer transits per sector, lower SNR
- **Small planets around faint stars**: Below automated detection threshold (SNR < 7.1)
- **Grazing transits**: Shallow, asymmetric signals flagged as stellar activity
- **Multi-planet systems**: Pipeline optimized for single dominant signal
- **Edge cases in stitching**: Artifacts at sector boundaries confuse algorithms

**3. Human Advantages**

- Pattern recognition across noisy data
- Contextual judgment (stellar activity vs transit)
- Willingness to inspect marginal candidates pipeline ignores
- Cross-sector consistency checks

**4. What Planets Citizens Actually Find**

Recent citizen discoveries (Planet Hunters TESS):

- **TOI-813 b**: 84-day period Neptune (missed by pipeline - long period)
- **TOI-1338 b**: Circumbinary planet (complex binary signals)
- **TOI-2180 b**: 261-day giant (0.5% depth, only 3 transits)

Common traits:
- Periods: 15-200 days (pipeline sweet spot is 0.5-15 days)
- Depths: 0.2-1.5% (below auto-detection for faint stars)
- Host stars: Tmag 10-13 (not bright priorities)

**5. Realistic Success Metrics**

- **Candidate identification**: 1-5% of carefully selected targets
- **Confirmed planet**: 10-30% of vetted candidates (need follow-up)
- **Co-authorship**: Standard if you find it first and vet properly
- **Timeline**: 3-12 months from discovery to paper submission

**Not a hobby. This is real science with real publication potential.**

---

## PART 2 — TARGET SELECTION (HIGH-ODDS STRATEGY)

### Scientific Rationale for Filters

**Goal**: Identify 100 stars with maximum transit probability and minimum false positive rate.

### Selection Criteria

| Filter | Value | Justification |
|--------|-------|---------------|
| **Spectral Type** | K5-M4 (Teff: 2700-4500 K) | Smaller stars → larger transit depth for same planet size. M dwarfs have R ≈ 0.2-0.6 R☉ |
| **TESS Magnitude** | 10.0 < Tmag < 13.5 | Bright enough for good SNR, faint enough to avoid saturation and pipeline scrutiny |
| **Sectors Observed** | ≥ 4 | More data → better period coverage, lower noise, more transits |
| **Stellar Radius** | 0.3 < R < 0.7 R☉ | Confirms K/M dwarf, avoids giants (false positive risk) |
| **Distance** | d < 100 pc | Ground follow-up feasible, proper motion well-measured |
| **Contamination Ratio** | < 0.1 | Minimize blended eclipsing binaries |
| **Known TOIs** | Exclude | Avoid duplicating pipeline discoveries |
| **Known Planets** | Exclude | Focus on unexplored systems |

### Why This Works

- **Small stars**: Earth-sized planet creates 0.5-2% depth (vs 0.01% for Sun-like star)
- **Multi-sector**: 4 sectors = 108 days of data → can detect P up to ~50 days
- **Magnitude range**: Pipeline prioritizes Tmag < 10; we exploit 10-13.5
- **Contamination**: Nearby eclipsing binary can mimic planet signal

### Expected Yield

From ~500,000 TESS targets:
- K/M dwarfs with Tmag 10-13.5: ~80,000 stars
- Multi-sector (≥4): ~12,000 stars
- After quality filters: ~800 candidates
- **Random selection of 100 from this pool gives ~1-5 likely candidates**

---

## PART 3 — DATA INGESTION PIPELINE

### Workflow

1. Download all SPOC 2-minute light curves for target
2. Stitch sectors together with normalized flux
3. Remove systematic trends (not astrophysical signals)
4. Identify and sigma-clip outliers
5. Save processed light curve for analysis

### Why Each Step Matters

**Stitching**: TESS observes in ~27-day sectors. Multi-sector targets have gaps - naive concatenation creates edge artifacts.

**Normalization**: Each sector has different mean flux due to:
- CCD temperature variations
- Scattered light
- Different pixel apertures

**Detrending**: Remove:
- Spacecraft momentum dumps (every ~2.5 days)
- Pointing drift (correlated with centroid motion)
- Long-term stellar variability (not transit-related)

**But preserve**:
- Transit signals (hours to days)
- True astrophysical variation

### Critical Choices

- Use **lightkurve.flatten(window_length)**: Removes trends longer than window
- **Window = 0.5-1.0 days**: Preserves transits (hours) but removes systematics (days+)
- **Sigma clipping = 5σ**: Aggressive outlier removal can clip transits - use conservative threshold

---

## PART 4 — TRANSIT SEARCH (BLS IMPLEMENTATION)

### Box Least Squares Method

**Principle**: Fit a periodic box-shaped signal to light curve and identify period with strongest match.

**Why BLS for Transits**

- Transits are approximately box-shaped (flat-bottomed, sharp edges)
- BLS maximizes signal detection efficiency (SDE) for this morphology
- Robust to red noise (stellar variability)

### Parameter Choices

| Parameter | Value | Reason |
|-----------|-------|--------|
| **Minimum Period** | 0.5 days | Below this: hot Jupiter territory, pipeline already found |
| **Maximum Period** | 0.4 × baseline | Nyquist limit - need 2.5+ transits to constrain period |
| **Frequency Resolution** | df = 0.1 / T² | Ensures period bins finer than transit duration |
| **Duration Range** | 0.5-8 hours | Shorter = grazing; longer = blend or stellar activity |

### Signal Detection Efficiency (SDE)

**SDE = (signal depth) / (noise level)**

- **SDE > 7**: High confidence (pipeline threshold)
- **SDE 5-7**: Marginal, needs careful vetting
- **SDE < 5**: Likely noise

**Citizens target SDE 4-6 range that pipeline discards.**

### Multi-Peak Analysis

BLS returns power spectrum. Top peak = best period, but:

- **Check harmonics**: Period and 2P, 0.5P should have different depths
- **Inspect top 5 peaks**: Real planets sometimes at 2nd-3rd peak (window function aliasing)
- **1-day aliases**: TESS observes continuously, but Earth-based systematic at 1.0 day period

---

## PART 5 — VETTING CHECKLIST (CRITICAL STEP)

### False Positive Sources

1. **Eclipsing Binary (EB)**: Background or bound companion
2. **Blended EB**: Nearby star in photometric aperture
3. **Stellar Activity**: Spots rotating in/out of view
4. **Systematics**: Spacecraft artifacts aliased to orbital period
5. **Centroid Motion**: Signal from different pixel location (blend)

### Automated Vetting Tests

#### Test 1: Odd/Even Transit Depth Consistency

**Method**: Fold light curve at 2× period, compare depths of alternating transits.

**Logic**: Real planet has identical transits. EB has two different depths (primary/secondary eclipse).

**Decision**:
- Depth difference < 20%: PASS
- Depth difference > 50%: REJECT (EB)
- 20-50%: FLAG (marginal, possible grazing)

#### Test 2: Transit Shape

**Method**: Fit transit model, measure ingress/egress symmetry and bottom flatness.

**V-shaped**: Grazing geometry or blended EB
**U-shaped**: Central transit (good)
**Asymmetric**: Systematics or spot-crossing event

**Decision**:
- χ² reduced < 2.0: PASS
- Ingress/egress duration ratio 0.7-1.3: PASS
- Otherwise: FLAG

#### Test 3: Secondary Eclipse Check

**Method**: Fold at same period, phase = 0.5, check for signal.

**Logic**:
- Planets: No secondary (non-luminous) unless hot Jupiter (IR emission, depth << transit)
- EBs: Clear secondary eclipse

**Decision**:
- Secondary depth > 30% of primary: REJECT (EB)
- Secondary depth 10-30%: FLAG (possible hot planet or blend)
- No secondary: PASS

#### Test 4: Depth Plausibility

**Method**: Compare transit depth to stellar radius.

**Planet radius from depth**:
R_planet = R_star × √(depth)

**Decision**:
- R_p < 0.3 R_Jup (3.4 R⊕): PASS (terrestrial to Neptune)
- 0.3-2.0 R_Jup: PASS (gas giant)
- R_p > 2.0 R_Jup: REJECT (unphysical - likely EB)

#### Test 5: Sector-to-Sector Consistency

**Method**: Run BLS on each sector independently, check if same period detected.

**Decision**:
- Same period (±1%) in ≥75% of sectors: PASS
- Inconsistent: REJECT (systematic artifact)

#### Test 6: Centroid Analysis

**Method**: Compare pixel-level centroid during vs out-of-transit.

**Logic**: If signal from background star, centroid shifts toward that star during transit.

**Decision**:
- Centroid shift < 1 pixel: PASS
- Shift > 2 pixels: REJECT (blend)

#### Test 7: Red Noise Check

**Method**: Compare transit depth to RMS of binned light curve.

**SNR = depth / σ_binned**

**Decision**:
- SNR > 5: PASS
- SNR 3-5: FLAG (marginal)
- SNR < 3: REJECT (noise)

### Manual Inspection Checklist

For each candidate that passes automated tests:

- [ ] Visually inspect raw light curve for data quality
- [ ] Check if transits occur near sector boundaries (edge artifact)
- [ ] Look for correlated systematics (thruster fires, safe modes)
- [ ] Verify star is not in crowded field (Gaia EDR3 neighbors)
- [ ] Check SIMBAD for known variables, binaries
- [ ] Review TIC flags for contamination warnings

### Decision Tree

```
PASS all automated tests + manual OK → KEEP (strong candidate)
PASS most tests, 1-2 FLAGS → KEEP (moderate candidate, needs follow-up)
Any REJECT → DISCARD
Multiple FLAGS → DISCARD (too uncertain)
```

---

## PART 6 — OUTPUT AND RANKING

### Quality Score Formula

Assign each candidate a score (0-100):

```
Score = 30×(SNR/10) + 20×(N_transits/5) + 15×(1 - σ_depth/depth)
        + 15×(χ²_fit < 1.5) + 10×(odd_even_consistency) + 10×(no_secondary)
```

Where:
- **SNR**: Signal-to-noise ratio (cap at 10)
- **N_transits**: Number of observed transits (cap at 5)
- **σ_depth**: Uncertainty on depth measurement
- **χ²_fit**: Quality of transit model fit (binary: 15 if χ² < 1.5, else 0)
- **Odd/even consistency**: Binary (10 if consistent, else 0)
- **No secondary**: Binary (10 if none detected, else 0)

**Interpretation**:
- Score > 70: High-priority candidate
- Score 50-70: Moderate candidate
- Score < 50: Low confidence, likely false positive

### Output Table Format

Generate `candidates_ranked.csv`:

```
TIC_ID, Period_days, Depth_ppt, Radius_Rearth, N_transits, SNR, Score, Flags
```

Example:
```
TIC_12345678, 23.4, 4.2, 2.1, 4, 6.8, 78, "PASS"
TIC_87654321, 15.1, 1.8, 1.3, 5, 4.2, 65, "FLAG:secondary_marginal"
```

### Diagnostic Plots for Each Candidate

Auto-generate `top_10_candidates/TIC_XXXXXXXX/` containing:

1. **full_lightcurve.png**: All sectors stitched, with transits marked
2. **folded_lightcurve.png**: Phase-folded at detected period
3. **bls_periodogram.png**: Power spectrum with peak marked
4. **odd_even_comparison.png**: Odd vs even transits overlaid
5. **secondary_eclipse_check.png**: Phase 0.5 region
6. **summary.txt**: All parameters and vetting results

### Summary Markdown Template

```markdown
# TIC XXXXXXXX - Candidate Summary

## Stellar Properties
- TESS Magnitude: XX.X
- Spectral Type: MX / KX
- Radius: X.XX R☉
- Distance: XX pc

## Candidate Properties
- Period: XX.XX ± X.XX days
- Transit Depth: X.XX ± X.XX ppt (parts per thousand)
- Planet Radius: X.XX ± X.XX R⊕
- Number of Transits: X
- SNR: X.X

## Vetting Status
- Odd/Even Consistency: PASS / FLAG / REJECT
- Secondary Eclipse: Not detected
- Transit Shape: U-shaped, symmetric
- Centroid Shift: < 0.5 pixels
- Quality Score: XX / 100

## Notes
[Any unusual features, follow-up recommendations]
```

---

## PART 7 — SUBMISSION AND CREDIT STRATEGY

### Where to Submit

**Primary Path: Planet Hunters TESS (Zooniverse)**

- URL: https://www.zooniverse.org/projects/nora-dot-eisner/planet-hunters-tess
- Active professional astronomer review
- Established co-authorship policy

**Process**:
1. Complete vetting independently
2. Submit via "Talk" discussion board for subject
3. Tag @scientists with summary
4. Provide link to plots and analysis code
5. Wait for professional review (1-4 weeks)

**Alternative: Direct ExoFOP-TESS Submission**

- Requires IPAC account (free)
- Upload finding charts, light curves, vetting reports
- Community-accessible database
- Less direct path to co-authorship but establishes priority

### What to Include in Submission

**Required**:
- TIC ID and coordinates
- Period and epoch of first transit (BJD format)
- Transit depth and duration
- Number of observed transits
- Folded light curve plot
- BLS periodogram

**Strongly Recommended**:
- Odd/even comparison
- Secondary eclipse check
- Description of vetting performed
- Link to code repository (GitHub)

**Optional but Helpful**:
- Stellar characterization (Teff, radius from TIC)
- Predicted planet radius
- Observability windows for ground follow-up

### Co-Authorship Guidelines (Planet Hunters TESS Policy)

**You WILL be a co-author if**:
- You identified the candidate first
- You performed initial vetting
- Candidate confirmed as real planet

**Author list position**:
- Typically middle authorship (not first or last)
- First author = professional who led confirmation
- You appear in "Planet Hunters TESS Collaboration" group

**What You Must Provide**:
- Detailed discovery narrative (how you found it)
- All analysis scripts
- Vetting documentation
- Response to professional astronomer questions

**Timeline**:
- Discovery → Initial review: 2-8 weeks
- Review → Follow-up decision: 1-3 months
- Follow-up → Paper submission: 6-18 months
- Submission → Publication: 3-6 months

**Total: 1-2 years from discovery to published paper with your name on it.**

### What NOT to Claim

**Unacceptable**:
- "I discovered a planet" (before confirmation)
- Announcing on social media before professional review
- Claiming discovery if you missed vetting red flags

**Acceptable**:
- "I identified a planet candidate in TESS data"
- "I'm working with astronomers to confirm a potential exoplanet"
- Sharing after professional review confirms it's worth follow-up

### Submission Checklist Template

```
[ ] Candidate passes all automated vetting tests
[ ] Manual inspection completed
[ ] Checked SIMBAD, TIC, ExoFOP for prior reports
[ ] Generated diagnostic plots
[ ] Written discovery summary (< 300 words)
[ ] Code uploaded to GitHub (public repository)
[ ] Submitted to Planet Hunters Talk board
[ ] Tagged relevant scientists
[ ] Prepared to answer follow-up questions
[ ] Documented all analysis decisions
```

### Discovery Summary Template (for Submission)

```
Subject: Candidate Transit Signal in TIC XXXXXXXX

I identified a candidate planetary transit signal in TESS data for TIC XXXXXXXX using Box Least Squares period analysis.

Stellar Properties:
- Tmag: XX.X
- Spectral Type: [K/M]
- Radius: X.XX R☉

Signal Properties:
- Period: XX.XX days
- Depth: X.XX ppt
- Duration: X.X hours
- Implied Radius: X.X R⊕
- Number of Transits: X

Vetting Performed:
- Odd/even transit depths consistent
- No secondary eclipse detected
- Transit shape U-shaped and symmetric
- Centroid stable during transits
- Signal present in all observed sectors

Diagnostic plots and analysis code: [GitHub link]

Requesting professional review for follow-up consideration.

[Your name]
```

---

## PART 8 — REALISTIC ODDS AND EXPECTATIONS

### Hit Rate Projections

**Starting Pool**: 100 carefully selected K/M dwarfs (Tmag 10-13.5, multi-sector)

**Expected Outcomes**:

| Stage | Count | Percentage |
|-------|-------|------------|
| Stars analyzed | 100 | 100% |
| BLS detects periodic signal (SDE > 4) | 15-25 | 15-25% |
| Passes automated vetting | 3-8 | 3-8% |
| Passes manual inspection | 1-5 | 1-5% |
| Professionals confirm worth follow-up | 1-2 | 1-2% |
| Follow-up confirms planet | 0-1 | 0-1% |

**Most Likely Outcome**: 1-3 candidates worth submitting per 100 stars analyzed.

**Confirmation Rate**: ~20-40% of well-vetted citizen candidates confirm as real planets.

**Base Rate**: ~0.2-0.5 confirmed planets per 100 stars (1 in 200-500).

### Time Investment

**Per-Star Breakdown**:
- Target selection script: 2-4 hours (one-time)
- Data download per star: 2-5 minutes
- BLS analysis per star: 5-15 minutes
- Automated vetting per star: 5 minutes
- Manual inspection per candidate: 30-60 minutes
- Deep vetting per strong candidate: 2-4 hours

**Total Time Estimates**:

| Scenario | Time Required |
|----------|---------------|
| Analyze 100 stars (bulk processing) | 20-40 hours |
| Vet 5 candidates in detail | 10-20 hours |
| Prepare 1 submission-quality report | 4-8 hours |
| **Total for first discovery** | **35-70 hours** |

**Calendar Time**:
- Working 5 hours/week: 7-14 weeks (2-3 months)
- Working 10 hours/week: 4-7 weeks (1-2 months)

### Outcome Probabilities (90 days, 10 hrs/week)

**Possible Outcomes**:

1. **No viable candidates** (30% probability)
   - You analyzed sample, found signals, all failed vetting
   - Learned methods, contributed null results (still scientifically valuable)
   - Action: Analyze 100 more stars or refine selection criteria

2. **1-2 marginal candidates** (50% probability)
   - Submitted to Planet Hunters
   - Professional review: "Interesting but low priority" or "Likely false positive"
   - No follow-up scheduled
   - Action: Continue searching, incorporate feedback

3. **1 strong candidate** (18% probability)
   - Professionals confirm worth follow-up
   - Ground-based observations scheduled within 6-12 months
   - You're added to collaboration email list
   - 40% chance this confirms as planet → co-authorship

4. **Multiple strong candidates** (2% probability)
   - Lightning strikes
   - You found a productive niche or got lucky
   - High likelihood of at least one confirmation

### Realistic Success Stories

**Case Study 1: Tom Jacobs (Visual Survey)**
- Analyzed ~10,000 TESS targets over 2 years
- Found 6 confirmed planets (0.06% hit rate)
- Co-author on 6 papers
- Most notable: TIC 172900988 (warm Jupiter)

**Case Study 2: Citizen Science Team (Planet Hunters)**
- TOI-2180 b: Found by volunteer reviewing ~1000 stars
- 261-day period, 0.5% depth
- Took 8 months from discovery to confirmation
- Team member co-authored Nature paper

**Takeaway**: High volume + persistence + rigor = results.

### What Will Definitely Happen

**After 30 days**:
- You'll understand TESS data structure intimately
- You'll have analyzed 20-40 stars
- You'll have seen 5-10 false positives (eclipsing binaries, systematics)
- You'll have refined your vetting intuition

**After 60 days**:
- You'll have analyzed 50-80 stars
- You'll likely have 1-2 candidates under review
- You'll be active in Planet Hunters community
- You'll be debugging edge cases in your code

**After 90 days**:
- You'll have completed 100-star sample
- You'll have submitted at least one candidate (if selection criteria sound)
- You'll be waiting on professional feedback
- You'll be planning next target list

### What Won't Happen

**Unrealistic Expectations**:
- Finding a planet in your first 10 stars
- Discovering an Earth analog in habitable zone (too faint to detect at this SNR)
- Getting immediate confirmation (follow-up takes months to years)
- Becoming famous (even confirmed planet is niche science news)

### Expected Value Calculation

**Investment**: 70 hours over 3 months

**Potential Outcomes**:
- Learn advanced time-series analysis: 100% probability
- Contribute to citizen science: 100% probability
- Submit viable candidate: 60-70% probability
- Trigger professional follow-up: 15-20% probability
- Co-author a paper: 5-10% probability

**Is It Worth It?**

If your goal is:
- Learn real astrophysics: **Yes**
- Contribute to science: **Yes**
- Co-author a paper: **Maybe** (but better odds than lottery)
- Get rich or famous: **No**

### Failure Modes and Mitigations

**Problem**: "I analyzed 100 stars, found nothing."

**Cause**: Selection criteria too conservative or pipeline already caught everything in this regime.

**Mitigation**: Expand to 200 stars, or adjust Tmag range to 13.5-14.5 (fainter, noisier, but less picked-over).

---

**Problem**: "I keep finding eclipsing binaries."

**Cause**: Vetting not strict enough, or target selection includes close binaries.

**Mitigation**: Add Gaia RUWE < 1.2 filter (rejects unresolved binaries), stricter odd/even test.

---

**Problem**: "My candidates get rejected by professionals."

**Cause**: Insufficient vetting documentation or marginal SNR.

**Mitigation**: Include more diagnostic plots, calculate formal false alarm probability (FAP), wait for higher-SNR candidates.

---

### The Bottom Line

**If you follow this workflow rigorously for 100 stars:**

- **60-70% chance**: You'll submit at least one candidate to professionals
- **15-25% chance**: Professionals will request follow-up
- **5-10% chance**: You'll co-author a paper within 2 years

**These are real odds. This is a real scientific contribution.**

Do not expect overnight results. Do expect to learn more about data analysis, statistical inference, and astrophysics than any textbook can teach.

The planets are in the data. You just have to find them.

---

*End of Strategic Guide. Implementation code follows in separate files.*
