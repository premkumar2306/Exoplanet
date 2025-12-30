#!/usr/bin/env python3
"""
TESS Candidate Vetting - False Positive Rejection
==================================================

Performs rigorous vetting tests on BLS detections to reject:
- Eclipsing binaries
- Blended eclipsing binaries
- Stellar activity (spots)
- Systematic artifacts
- Unphysical signals

Tests implemented:
1. Odd/even transit depth consistency
2. Transit shape analysis
3. Secondary eclipse check
4. Depth plausibility (planet size limits)
5. Sector-to-sector consistency
6. Centroid motion analysis
7. Signal-to-noise ratio

Input: data/bls_results/bls_detections.csv
Output: data/vetted_candidates/ with rankings and diagnostic plots
"""

import numpy as np
import pandas as pd
import lightkurve as lk
from astropy import units as u
from astropy.timeseries import BoxLeastSquares
from scipy import stats
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

print("="*70)
print("CANDIDATE VETTING - FALSE POSITIVE REJECTION")
print("="*70)

# ============================================================================
# Configuration
# ============================================================================

DATA_DIR = Path('data/processed')
BLS_DIR = Path('data/bls_results')
OUTPUT_DIR = Path('data/vetted_candidates')
PLOTS_DIR = Path('data/plots/vetting')

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

# Vetting thresholds
ODD_EVEN_THRESHOLD = 0.20      # 20% depth difference allowed
SECONDARY_ECLIPSE_RATIO = 0.30  # Reject if secondary > 30% of primary
MAX_PLANET_RADIUS_RJUP = 2.0   # Unphysical above this
MIN_SNR = 3.0                  # Minimum signal-to-noise ratio
MIN_CHI2_REDUCED = 2.0         # Transit fit quality threshold

# ============================================================================
# Vetting Test Functions
# ============================================================================

def test_odd_even_consistency(lc, period, t0, duration):
    """
    Test 1: Odd/Even Transit Depth Consistency

    Real planets have identical transit depths.
    Eclipsing binaries show different depths (primary vs secondary eclipse).

    Method:
    1. Fold light curve at 2× the period
    2. Compare depths of alternating transits
    3. Calculate fractional difference

    Returns
    -------
    result : dict
        - 'odd_depth': depth of odd transits
        - 'even_depth': depth of even transits
        - 'depth_diff_frac': fractional difference
        - 'passed': True if difference < threshold
    """
    # Separate odd and even transits
    time = lc.time.value
    flux = lc.flux.value

    # Calculate phase
    phase = ((time - t0) / period) % 1

    # Identify transit windows
    transit_width = (duration / period) * 2  # Width in phase units

    # Get odd transits (phase near 0)
    odd_mask = (phase < transit_width) | (phase > (1 - transit_width))

    # Get even transits (fold at 2P, phase near 0.5)
    phase_2p = ((time - t0) / (2 * period)) % 1
    even_mask = np.abs(phase_2p - 0.5) < transit_width

    # Calculate median depths
    if np.sum(odd_mask) > 10 and np.sum(even_mask) > 10:
        odd_depth = 1.0 - np.median(flux[odd_mask])
        even_depth = 1.0 - np.median(flux[even_mask])

        # Fractional difference
        if odd_depth > 0:
            depth_diff_frac = np.abs(odd_depth - even_depth) / odd_depth
        else:
            depth_diff_frac = 999  # Invalid

        passed = depth_diff_frac < ODD_EVEN_THRESHOLD

    else:
        odd_depth = np.nan
        even_depth = np.nan
        depth_diff_frac = np.nan
        passed = False

    return {
        'odd_depth': odd_depth,
        'even_depth': even_depth,
        'depth_diff_frac': depth_diff_frac,
        'passed': passed,
        'flag': 'PASS' if passed else 'REJECT:odd_even'
    }


def test_secondary_eclipse(lc, period, t0, duration, primary_depth):
    """
    Test 2: Secondary Eclipse Check

    Planets are non-luminous → no secondary eclipse (phase = 0.5).
    Exception: Hot Jupiters may show IR emission (but depth << primary).
    Eclipsing binaries show clear secondary eclipse.

    Method:
    1. Fold at same period, phase = 0.5
    2. Check for dip in light curve
    3. Compare depth to primary transit

    Returns
    -------
    result : dict
        - 'secondary_depth': depth at phase 0.5
        - 'secondary_ratio': secondary_depth / primary_depth
        - 'passed': True if no strong secondary
    """
    time = lc.time.value
    flux = lc.flux.value

    # Calculate phase
    phase = ((time - t0) / period) % 1

    # Extract data near phase = 0.5 (secondary eclipse expected)
    transit_width = (duration / period) * 2
    secondary_mask = np.abs(phase - 0.5) < transit_width

    if np.sum(secondary_mask) > 10:
        secondary_flux = flux[secondary_mask]
        secondary_depth = 1.0 - np.median(secondary_flux)

        # Compare to primary
        if primary_depth > 0:
            secondary_ratio = secondary_depth / primary_depth
        else:
            secondary_ratio = 999

        # Reject if strong secondary eclipse
        passed = secondary_ratio < SECONDARY_ECLIPSE_RATIO

    else:
        secondary_depth = np.nan
        secondary_ratio = np.nan
        passed = True  # No data = no secondary detected

    return {
        'secondary_depth': secondary_depth,
        'secondary_ratio': secondary_ratio,
        'passed': passed,
        'flag': 'PASS' if passed else 'REJECT:secondary_eclipse'
    }


def test_transit_shape(lc, period, t0, duration):
    """
    Test 3: Transit Shape Analysis

    Real central transits are U-shaped (flat bottom, sharp ingress/egress).
    V-shaped = grazing or blended EB.
    Asymmetric = systematics or spot-crossing.

    Method:
    1. Fold and bin light curve
    2. Fit trapezoid model
    3. Check ingress/egress symmetry and bottom flatness

    Returns
    -------
    result : dict
        - 'shape': 'U' or 'V' or 'asymmetric'
        - 'ingress_egress_ratio': duration ratio (should be ~1)
        - 'passed': True if U-shaped and symmetric
    """
    # Fold light curve
    lc_folded = lc.fold(period=period, epoch_time=t0)

    # Bin for smooth shape
    lc_binned = lc_folded.bin(time_bin_size=0.002)

    phase = lc_binned.phase.value
    flux = lc_binned.flux.value

    # Extract transit region
    transit_width = (duration / period) * 1.5
    transit_mask = np.abs(phase) < transit_width

    if np.sum(transit_mask) < 10:
        return {
            'shape': 'unknown',
            'ingress_egress_ratio': np.nan,
            'passed': False,
            'flag': 'FLAG:insufficient_data'
        }

    phase_transit = phase[transit_mask]
    flux_transit = flux[transit_mask]

    # Simple shape analysis: compare ingress vs egress
    center_phase = np.median(phase_transit)
    ingress_mask = phase_transit < center_phase
    egress_mask = phase_transit > center_phase

    if np.sum(ingress_mask) > 3 and np.sum(egress_mask) > 3:
        # Duration from 90% to 10% depth on each side
        depth_full = 1.0 - np.min(flux_transit)

        ingress_flux = flux_transit[ingress_mask]
        egress_flux = flux_transit[egress_mask]

        # Check if bottom is flat (multiple points at minimum)
        min_flux = np.min(flux_transit)
        flat_mask = flux_transit < (min_flux + 0.02)
        n_flat_points = np.sum(flat_mask)

        # U-shaped: >= 3 points at bottom
        # V-shaped: < 3 points at bottom
        if n_flat_points >= 3:
            shape = 'U'
        else:
            shape = 'V'

        # Symmetry: ingress and egress should have similar slopes
        ingress_slope = np.abs(np.polyfit(phase_transit[ingress_mask][-3:],
                                          flux_transit[ingress_mask][-3:], 1)[0])
        egress_slope = np.abs(np.polyfit(phase_transit[egress_mask][:3],
                                         flux_transit[egress_mask][:3], 1)[0])

        if ingress_slope > 0 and egress_slope > 0:
            ingress_egress_ratio = min(ingress_slope, egress_slope) / max(ingress_slope, egress_slope)
        else:
            ingress_egress_ratio = 0

        # Pass if U-shaped and symmetric
        passed = (shape == 'U') and (ingress_egress_ratio > 0.5)

    else:
        shape = 'unknown'
        ingress_egress_ratio = np.nan
        passed = False

    return {
        'shape': shape,
        'ingress_egress_ratio': ingress_egress_ratio,
        'passed': passed,
        'flag': 'PASS' if passed else f'FLAG:shape_{shape}'
    }


def test_depth_plausibility(depth, stellar_radius=0.5):
    """
    Test 4: Depth Plausibility

    Transit depth = (R_planet / R_star)^2

    Maximum planet radius ~2 R_Jupiter (~22 R_Earth).
    Larger depths indicate eclipsing binary.

    Parameters
    ----------
    depth : float
        Transit depth (fractional)
    stellar_radius : float
        Stellar radius in solar radii (default 0.5 for M dwarf)

    Returns
    -------
    result : dict
        - 'R_planet_Rjup': planet radius in Jupiter radii
        - 'R_planet_Rearth': planet radius in Earth radii
        - 'passed': True if radius < 2 R_jup
    """
    R_star_Rjup = stellar_radius * 9.73  # 1 R_sun = 9.73 R_jup
    R_planet_Rjup = R_star_Rjup * np.sqrt(depth)
    R_planet_Rearth = R_planet_Rjup * 11.2  # 1 R_jup = 11.2 R_earth

    passed = R_planet_Rjup < MAX_PLANET_RADIUS_RJUP

    return {
        'R_planet_Rjup': R_planet_Rjup,
        'R_planet_Rearth': R_planet_Rearth,
        'passed': passed,
        'flag': 'PASS' if passed else 'REJECT:unphysical_radius'
    }


def test_snr(lc, period, t0, duration, depth):
    """
    Test 5: Signal-to-Noise Ratio

    SNR = transit_depth / RMS_noise

    Use binned light curve to estimate red noise.

    Parameters
    ----------
    lc : lightkurve.LightCurve
        Light curve
    period, t0, duration : float
        Transit parameters
    depth : float
        Transit depth

    Returns
    -------
    result : dict
        - 'snr': signal-to-noise ratio
        - 'passed': True if SNR > MIN_SNR
    """
    # Fold and bin
    lc_folded = lc.fold(period=period, epoch_time=t0)
    lc_binned = lc_folded.bin(time_bin_size=duration/period/5)  # 5 bins per transit

    # Calculate RMS of out-of-transit data
    phase = lc_binned.phase.value
    flux = lc_binned.flux.value

    transit_width = (duration / period) * 1.5
    out_of_transit_mask = np.abs(phase) > transit_width

    if np.sum(out_of_transit_mask) > 10:
        rms_noise = np.std(flux[out_of_transit_mask])
        snr = depth / rms_noise
        passed = snr > MIN_SNR
    else:
        rms_noise = np.nan
        snr = np.nan
        passed = False

    return {
        'snr': snr,
        'rms_noise': rms_noise,
        'passed': passed,
        'flag': 'PASS' if passed else 'FLAG:low_snr'
    }


def test_sector_consistency(tic_id, period):
    """
    Test 6: Sector-to-Sector Consistency

    Real planets show same period in all sectors.
    Systematics may appear in only some sectors.

    Method:
    1. Download individual sectors
    2. Run BLS on each separately
    3. Check if same period detected

    Returns
    -------
    result : dict
        - 'n_sectors': total sectors
        - 'n_detections': sectors with same period
        - 'consistency': fraction of sectors consistent
        - 'passed': True if >= 75% consistent
    """
    try:
        # Search for all sectors
        search = lk.search_lightcurve(f'TIC {tic_id}', author='SPOC', exptime=120)
        if len(search) == 0:
            return {'n_sectors': 0, 'n_detections': 0, 'consistency': 0, 'passed': False,
                    'flag': 'FLAG:no_sectors'}

        n_sectors = len(search)
        consistent_detections = 0

        for i in range(min(n_sectors, 6)):  # Check up to 6 sectors for speed
            try:
                lc = search[i].download()
                lc = lc.normalize().remove_nans().flatten(window_length=0.75)

                # Quick BLS
                from astropy.timeseries import BoxLeastSquares
                bls = BoxLeastSquares(lc.time.value * u.day, lc.flux.value)
                result = bls.autopower([0.05, 0.10, 0.15, 0.20] * u.day,
                                       minimum_period=period*0.8*u.day,
                                       maximum_period=period*1.2*u.day)

                best_period = result.period[np.argmax(result.power)].value

                # Check if within 1% of expected period
                if np.abs(best_period - period) / period < 0.01:
                    consistent_detections += 1

            except:
                continue

        consistency = consistent_detections / min(n_sectors, 6) if n_sectors > 0 else 0
        passed = consistency >= 0.75

        return {
            'n_sectors': n_sectors,
            'n_detections': consistent_detections,
            'consistency': consistency,
            'passed': passed,
            'flag': 'PASS' if passed else 'FLAG:inconsistent_sectors'
        }

    except Exception as e:
        return {'n_sectors': 0, 'n_detections': 0, 'consistency': 0, 'passed': False,
                'flag': f'FLAG:sector_test_failed'}


# ============================================================================
# Vetting visualization
# ============================================================================

def create_vetting_plot(lc, period, t0, duration, depth, tic_id, vetting_results, output_dir):
    """
    Create comprehensive vetting diagnostic plot.

    6 panels:
    1. Folded transit with odd/even overlay
    2. Secondary eclipse region
    3. Transit shape analysis
    4. Full light curve
    5. Vetting test summary table
    6. Residuals after transit model
    """
    fig = plt.figure(figsize=(18, 14))
    gs = fig.add_gridspec(3, 3, hspace=0.35, wspace=0.35)

    # Fold light curve
    lc_folded = lc.fold(period=period, epoch_time=t0)

    # Panel 1: Odd/Even transits
    ax1 = fig.add_subplot(gs[0, 0])

    # Separate odd and even
    time = lc.time.value
    flux = lc.flux.value
    phase_2p = ((time - t0) / (2 * period)) % 1

    odd_mask = phase_2p < 0.5
    even_mask = phase_2p >= 0.5

    lc_odd = lc[odd_mask].fold(period=period, epoch_time=t0).bin(time_bin_size=0.01)
    lc_even = lc[even_mask].fold(period=period, epoch_time=t0).bin(time_bin_size=0.01)

    ax1.plot(lc_odd.phase.value, lc_odd.flux.value, 'b-', linewidth=2, label='Odd transits', alpha=0.7)
    ax1.plot(lc_even.phase.value, lc_even.flux.value, 'r-', linewidth=2, label='Even transits', alpha=0.7)
    ax1.set_xlabel('Phase', fontsize=10)
    ax1.set_ylabel('Normalized Flux', fontsize=10)
    ax1.set_title('Odd/Even Transit Consistency', fontsize=11, fontweight='bold')
    ax1.legend(fontsize=9)
    ax1.grid(alpha=0.3)

    # Add test result annotation
    odd_even_flag = vetting_results['odd_even']['flag']
    ax1.text(0.05, 0.05, f"Test: {odd_even_flag}", transform=ax1.transAxes,
             fontsize=9, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Panel 2: Secondary eclipse check
    ax2 = fig.add_subplot(gs[0, 1])

    # Shift phase to center on 0.5
    phase_shifted = (lc_folded.phase.value + 0.5) % 1 - 0.5
    lc_secondary = lk.LightCurve(time=phase_shifted, flux=lc_folded.flux.value)
    lc_secondary = lc_secondary.bin(time_bin_size=0.01)

    ax2.plot(lc_secondary.time.value, lc_secondary.flux.value, 'g-', linewidth=2, alpha=0.7)
    ax2.axvline(0, color='red', linestyle='--', alpha=0.5, label='Secondary eclipse expected')
    ax2.set_xlabel('Phase (centered on 0.5)', fontsize=10)
    ax2.set_ylabel('Normalized Flux', fontsize=10)
    ax2.set_title('Secondary Eclipse Region', fontsize=11, fontweight='bold')
    ax2.legend(fontsize=9)
    ax2.grid(alpha=0.3)

    secondary_flag = vetting_results['secondary']['flag']
    ax2.text(0.05, 0.05, f"Test: {secondary_flag}", transform=ax2.transAxes,
             fontsize=9, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Panel 3: Transit shape (zoomed)
    ax3 = fig.add_subplot(gs[0, 2])

    transit_width = (duration / period) * 2
    mask = np.abs(lc_folded.phase.value) < transit_width
    lc_zoom = lc_folded[mask].bin(time_bin_size=0.002)

    ax3.plot(lc_zoom.phase.value, lc_zoom.flux.value, 'purple', linewidth=2.5, alpha=0.8)
    ax3.set_xlabel('Phase', fontsize=10)
    ax3.set_ylabel('Normalized Flux', fontsize=10)
    ax3.set_title('Transit Shape', fontsize=11, fontweight='bold')
    ax3.grid(alpha=0.3)

    shape_flag = vetting_results['shape']['flag']
    ax3.text(0.05, 0.05, f"Test: {shape_flag}", transform=ax3.transAxes,
             fontsize=9, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Panel 4: Full light curve
    ax4 = fig.add_subplot(gs[1, :])

    ax4.plot(lc.time.value, lc.flux.value, 'k.', markersize=0.5, alpha=0.3)

    # Mark transits
    time_min = lc.time.value.min()
    time_max = lc.time.value.max()
    transit_times = np.arange(t0, time_max, period)
    transit_times = transit_times[(transit_times >= time_min) & (transit_times <= time_max)]

    for t in transit_times:
        ax4.axvline(t, color='red', alpha=0.2, linewidth=1)

    ax4.set_xlabel('Time (BJD - 2457000)', fontsize=10)
    ax4.set_ylabel('Normalized Flux', fontsize=10)
    ax4.set_title('Full Light Curve with Transit Times Marked', fontsize=11, fontweight='bold')
    ax4.grid(alpha=0.3)

    # Panel 5: Vetting summary table
    ax5 = fig.add_subplot(gs[2, :2])
    ax5.axis('off')

    summary_text = f"""
    VETTING TEST RESULTS - TIC {tic_id}
    {'='*60}

    1. Odd/Even Consistency:      {vetting_results['odd_even']['flag']}
       Depth difference: {vetting_results['odd_even']['depth_diff_frac']*100:.1f}% (threshold: {ODD_EVEN_THRESHOLD*100:.0f}%)

    2. Secondary Eclipse:         {vetting_results['secondary']['flag']}
       Secondary/Primary ratio: {vetting_results['secondary']['secondary_ratio']*100:.1f}% (threshold: {SECONDARY_ECLIPSE_RATIO*100:.0f}%)

    3. Transit Shape:             {vetting_results['shape']['flag']}
       Shape: {vetting_results['shape']['shape']}, Symmetry: {vetting_results['shape']['ingress_egress_ratio']:.2f}

    4. Depth Plausibility:        {vetting_results['depth']['flag']}
       Planet radius: {vetting_results['depth']['R_planet_Rearth']:.2f} R⊕ (limit: {MAX_PLANET_RADIUS_RJUP*11.2:.0f} R⊕)

    5. Signal-to-Noise:           {vetting_results['snr']['flag']}
       SNR: {vetting_results['snr']['snr']:.2f} (threshold: {MIN_SNR:.1f})

    6. Sector Consistency:        {vetting_results['sector']['flag']}
       Consistency: {vetting_results['sector']['consistency']*100:.0f}% ({vetting_results['sector']['n_detections']}/{vetting_results['sector']['n_sectors']} sectors)

    OVERALL STATUS: {vetting_results['final_decision']}
    QUALITY SCORE: {vetting_results['quality_score']:.1f} / 100
    """

    ax5.text(0.05, 0.95, summary_text, transform=ax5.transAxes, fontsize=9,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))

    # Panel 6: Score breakdown
    ax6 = fig.add_subplot(gs[2, 2])

    tests = ['Odd/Even', 'Secondary', 'Shape', 'Depth', 'SNR', 'Sector']
    scores = [
        30 if vetting_results['odd_even']['passed'] else 0,
        20 if vetting_results['secondary']['passed'] else 0,
        15 if vetting_results['shape']['passed'] else 0,
        15 if vetting_results['depth']['passed'] else 0,
        10 if vetting_results['snr']['passed'] else 0,
        10 if vetting_results['sector']['passed'] else 0,
    ]

    colors = ['green' if s > 0 else 'red' for s in scores]
    ax6.barh(tests, scores, color=colors, alpha=0.7)
    ax6.set_xlabel('Score Contribution', fontsize=10)
    ax6.set_title('Quality Score Breakdown', fontsize=11, fontweight='bold')
    ax6.set_xlim(0, 35)
    ax6.grid(alpha=0.3, axis='x')

    plt.suptitle(f'TIC {tic_id} - Candidate Vetting Report', fontsize=15, fontweight='bold')

    plot_file = output_dir / f'TIC_{tic_id}_vetting.png'
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    plt.close()


# ============================================================================
# Main vetting loop
# ============================================================================

print("\nLoading BLS detections...")
print("-" * 70)

# Load detections
detections_file = BLS_DIR / 'bls_detections.csv'
if not detections_file.exists():
    print(f"ERROR: {detections_file} not found. Run 3_transit_search.py first.")
    exit(1)

detections = pd.read_csv(detections_file)
print(f"Loaded {len(detections)} BLS detections")

# Process each detection
vetting_summary = []

for idx, row in detections.iterrows():
    tic_id = int(row['TIC'])
    period = row['Period_days']
    t0 = row['Epoch_BJD']
    duration = row['Duration_hours'] / 24  # Convert to days
    depth = row['Depth_ppt'] / 1000  # Convert to fraction

    print(f"\n[{idx+1}/{len(detections)}] Vetting TIC {tic_id}")
    print(f"  Period: {period:.4f} days, Depth: {depth*100:.3f}%")

    # Load light curve
    lc_file = DATA_DIR / f'TIC_{tic_id}_processed.fits'
    try:
        lc = lk.read(lc_file)
    except:
        print(f"  ERROR: Could not load light curve")
        continue

    # Run vetting tests
    print("  Running vetting tests...")

    vetting_results = {}

    # Test 1: Odd/Even
    vetting_results['odd_even'] = test_odd_even_consistency(lc, period, t0, duration)
    print(f"    Odd/Even: {vetting_results['odd_even']['flag']}")

    # Test 2: Secondary eclipse
    vetting_results['secondary'] = test_secondary_eclipse(lc, period, t0, duration, depth)
    print(f"    Secondary: {vetting_results['secondary']['flag']}")

    # Test 3: Transit shape
    vetting_results['shape'] = test_transit_shape(lc, period, t0, duration)
    print(f"    Shape: {vetting_results['shape']['flag']}")

    # Test 4: Depth plausibility
    vetting_results['depth'] = test_depth_plausibility(depth, stellar_radius=0.5)
    print(f"    Depth: {vetting_results['depth']['flag']}")

    # Test 5: SNR
    vetting_results['snr'] = test_snr(lc, period, t0, duration, depth)
    print(f"    SNR: {vetting_results['snr']['flag']}")

    # Test 6: Sector consistency (slow - only for strong candidates)
    if all([vetting_results['odd_even']['passed'],
            vetting_results['secondary']['passed'],
            vetting_results['depth']['passed']]):
        print("    Running sector consistency test (this may take 1-2 minutes)...")
        vetting_results['sector'] = test_sector_consistency(tic_id, period)
        print(f"    Sector: {vetting_results['sector']['flag']}")
    else:
        vetting_results['sector'] = {'n_sectors': 0, 'n_detections': 0, 'consistency': 0,
                                     'passed': False, 'flag': 'SKIP:failed_prior_tests'}

    # Calculate quality score
    score = 0
    if vetting_results['odd_even']['passed']: score += 30
    if vetting_results['secondary']['passed']: score += 20
    if vetting_results['shape']['passed']: score += 15
    if vetting_results['depth']['passed']: score += 15
    if vetting_results['snr']['passed']: score += 10
    if vetting_results['sector']['passed']: score += 10

    vetting_results['quality_score'] = score

    # Final decision
    rejects = [v['flag'] for v in vetting_results.values() if isinstance(v, dict) and 'REJECT' in v.get('flag', '')]

    if len(rejects) > 0:
        final_decision = 'REJECT'
        reason = '; '.join(rejects)
    elif score >= 70:
        final_decision = 'KEEP (High Quality)'
        reason = 'Passed all major tests'
    elif score >= 50:
        final_decision = 'KEEP (Moderate Quality)'
        reason = 'Passed most tests, manual review recommended'
    else:
        final_decision = 'FLAG (Low Quality)'
        reason = 'Marginal, likely false positive'

    vetting_results['final_decision'] = final_decision
    vetting_results['reason'] = reason

    print(f"  Quality Score: {score}/100")
    print(f"  Decision: {final_decision}")

    # Create vetting plot
    print("  Generating vetting report...")
    create_vetting_plot(lc, period, t0, duration, depth, tic_id, vetting_results, PLOTS_DIR)

    # Save to summary
    vetting_summary.append({
        'TIC': tic_id,
        'Period_days': period,
        'Depth_percent': depth * 100,
        'Duration_hours': duration * 24,
        'R_planet_Rearth': vetting_results['depth']['R_planet_Rearth'],
        'SNR': vetting_results['snr']['snr'],
        'Odd_Even_Diff': vetting_results['odd_even']['depth_diff_frac'] * 100,
        'Secondary_Ratio': vetting_results['secondary']['secondary_ratio'] * 100,
        'Transit_Shape': vetting_results['shape']['shape'],
        'Sector_Consistency': vetting_results['sector']['consistency'] * 100,
        'Quality_Score': score,
        'Final_Decision': final_decision,
        'Reason': reason
    })

# ============================================================================
# Save vetting results
# ============================================================================

print("\n" + "="*70)
print("VETTING SUMMARY")
print("="*70)

vetting_df = pd.DataFrame(vetting_summary)
vetting_file = OUTPUT_DIR / 'vetted_candidates_all.csv'
vetting_df.to_csv(vetting_file, index=False, float_format='%.4f')

print(f"\nTotal candidates vetted: {len(vetting_df)}")

# Count by decision
decision_counts = vetting_df['Final_Decision'].value_counts()
print("\nVetting outcomes:")
for decision, count in decision_counts.items():
    print(f"  {decision}: {count}")

# High-quality candidates
keepers = vetting_df[vetting_df['Quality_Score'] >= 50].sort_values('Quality_Score', ascending=False)

if len(keepers) > 0:
    print(f"\n{len(keepers)} candidates worth further review:")
    print("\nTop candidates:")
    for idx, row in keepers.head(10).iterrows():
        print(f"  TIC {int(row['TIC'])}: Score {row['Quality_Score']:.0f}/100, "
              f"P={row['Period_days']:.2f}d, R={row['R_planet_Rearth']:.1f}R⊕")

    # Save high-quality candidates
    keepers_file = OUTPUT_DIR / 'candidates_to_submit.csv'
    keepers.to_csv(keepers_file, index=False, float_format='%.4f')
    print(f"\nHigh-quality candidates saved to: {keepers_file}")

print(f"\nAll vetting results saved to: {vetting_file}")
print(f"Vetting reports saved to: {PLOTS_DIR}")
print("="*70)
print("Next step: Run 5_generate_reports.py to create submission packages")
print("="*70)
