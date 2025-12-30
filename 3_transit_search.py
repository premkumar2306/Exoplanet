#!/usr/bin/env python3
"""
TESS Transit Search using Box Least Squares (BLS)
==================================================

Performs period search on processed TESS light curves to identify
potential planetary transit signals.

Method: Box Least Squares (BLS)
- Searches for periodic box-shaped dips in light curve
- Returns strongest periodicities and their properties

Input: data/processed/*.fits
Output: data/bls_results/*.csv and diagnostic plots
"""

import numpy as np
import pandas as pd
import lightkurve as lk
from astropy import units as u
from astropy.timeseries import BoxLeastSquares
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

print("="*70)
print("TESS TRANSIT SEARCH (BOX LEAST SQUARES)")
print("="*70)

# ============================================================================
# Configuration
# ============================================================================

DATA_DIR = Path('data/processed')
OUTPUT_DIR = Path('data/bls_results')
PLOTS_DIR = Path('data/plots/bls')

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

# BLS parameters
MIN_PERIOD = 0.5       # days - below this is hot Jupiter territory (pipeline finds these)
MAX_PERIOD_FACTOR = 0.4  # maximum period as fraction of baseline (need 2.5+ transits)
MIN_DURATION = 0.02    # days (0.5 hours) - shorter likely grazing or artifact
MAX_DURATION = 0.33    # days (8 hours) - longer likely blend or stellar activity
DURATION_GRID = 20     # number of transit durations to test

# Detection threshold
MIN_SDE = 4.0          # Signal Detection Efficiency threshold (pipeline uses 7.1)

# ============================================================================
# Helper functions
# ============================================================================

def run_bls(lc, min_period=0.5, max_period=None, durations=None):
    """
    Run Box Least Squares period search on light curve.

    Parameters
    ----------
    lc : lightkurve.LightCurve
        Input light curve (should be normalized and detrended)
    min_period : float
        Minimum period to search (days)
    max_period : float
        Maximum period to search (days)
    durations : array-like
        Grid of transit durations to test (days)

    Returns
    -------
    result : astropy.timeseries.BoxLeastSquaresResults
        BLS results object containing periodogram and best-fit parameters
    """
    # Convert to astropy TimeSeries format
    time = lc.time.value * u.day
    flux = lc.flux.value
    flux_err = lc.flux_err.value if hasattr(lc, 'flux_err') else None

    # Set maximum period if not specified
    if max_period is None:
        baseline = time.max() - time.min()
        max_period = (MAX_PERIOD_FACTOR * baseline.value) * u.day

    # Set duration grid if not specified
    if durations is None:
        durations = np.linspace(MIN_DURATION, MAX_DURATION, DURATION_GRID) * u.day

    # Initialize BLS
    bls = BoxLeastSquares(time, flux, dy=flux_err)

    # Run period search
    # Frequency resolution: finer than 1/T^2 to avoid missing peaks
    baseline_days = (time.max() - time.min()).value
    frequency_factor = 0.1 / baseline_days**2

    result = bls.autopower(
        durations=durations,
        minimum_period=min_period * u.day,
        maximum_period=max_period,
        frequency_factor=frequency_factor
    )

    return result


def calculate_sde(bls_result):
    """
    Calculate Signal Detection Efficiency (SDE).

    SDE = (peak_power - median_power) / std_power

    This is the standard metric for BLS signal significance.

    Parameters
    ----------
    bls_result : BoxLeastSquaresResults
        BLS results object

    Returns
    -------
    sde : float
        Signal Detection Efficiency of strongest peak
    """
    power = bls_result.power
    median_power = np.median(power)
    std_power = np.std(power)

    peak_power = power.max()
    sde = (peak_power - median_power) / std_power

    return sde


def get_top_peaks(bls_result, n_peaks=5):
    """
    Identify top N peaks in BLS periodogram.

    Useful for checking harmonics and alternative periods.

    Parameters
    ----------
    bls_result : BoxLeastSquaresResults
        BLS results object
    n_peaks : int
        Number of peaks to return

    Returns
    -------
    peaks : list of dict
        Top peaks with their properties
    """
    power = bls_result.power
    period = bls_result.period.value

    # Find peaks (local maxima)
    from scipy.signal import find_peaks
    peak_indices, properties = find_peaks(power, height=np.median(power))

    # Sort by power
    sorted_indices = peak_indices[np.argsort(power[peak_indices])[::-1]]

    # Get top N
    top_indices = sorted_indices[:n_peaks]

    peaks = []
    for idx in top_indices:
        peaks.append({
            'period': period[idx],
            'power': power[idx],
            'index': idx
        })

    return peaks


def fold_lightcurve(lc, period, epoch_time):
    """
    Fold light curve at specified period.

    Parameters
    ----------
    lc : lightkurve.LightCurve
        Light curve to fold
    period : float
        Period in days
    epoch_time : float
        Reference time for phase=0 (typically first transit)

    Returns
    -------
    lc_folded : lightkurve.FoldedLightCurve
        Phase-folded light curve
    """
    lc_folded = lc.fold(period=period, epoch_time=epoch_time)
    return lc_folded


def plot_bls_results(lc, bls_result, tic_id, output_dir):
    """
    Create comprehensive BLS diagnostic plot.

    4-panel figure:
    1. Full light curve with best-fit transits marked
    2. BLS periodogram
    3. Folded light curve at best period
    4. Zoomed view of folded transit

    Parameters
    ----------
    lc : lightkurve.LightCurve
        Original light curve
    bls_result : BoxLeastSquaresResults
        BLS results
    tic_id : int
        TIC identifier
    output_dir : Path
        Directory to save plot
    """
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)

    # Get best-fit parameters
    best_period = bls_result.period[np.argmax(bls_result.power)].value
    best_t0 = bls_result.transit_time[np.argmax(bls_result.power)].value
    best_duration = bls_result.duration[np.argmax(bls_result.power)].value
    best_depth = bls_result.depth[np.argmax(bls_result.power)]

    sde = calculate_sde(bls_result)

    # Panel 1: Full light curve with transits marked
    ax1 = fig.add_subplot(gs[0, :])
    ax1.plot(lc.time.value, lc.flux.value, 'k.', markersize=1, alpha=0.5)

    # Mark predicted transit times
    time_min = lc.time.value.min()
    time_max = lc.time.value.max()
    transit_times = np.arange(best_t0, time_max, best_period)
    transit_times = transit_times[(transit_times >= time_min) & (transit_times <= time_max)]

    for t in transit_times:
        ax1.axvline(t, color='red', alpha=0.3, linewidth=1)

    ax1.set_xlabel('Time (BJD - 2457000)', fontsize=11)
    ax1.set_ylabel('Normalized Flux', fontsize=11)
    ax1.set_title(f'TIC {tic_id} - Full Light Curve (Transits Marked)', fontsize=13, fontweight='bold')
    ax1.grid(alpha=0.3)

    # Panel 2: BLS Periodogram
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.plot(bls_result.period, bls_result.power, 'k-', linewidth=0.5)
    ax2.axvline(best_period, color='red', linestyle='--', linewidth=2, label=f'Best: {best_period:.3f} d')
    ax2.axhline(np.median(bls_result.power), color='gray', linestyle=':', alpha=0.5)
    ax2.set_xlabel('Period (days)', fontsize=11)
    ax2.set_ylabel('BLS Power', fontsize=11)
    ax2.set_title(f'BLS Periodogram (SDE = {sde:.2f})', fontsize=12, fontweight='bold')
    ax2.legend(loc='upper right', fontsize=9)
    ax2.grid(alpha=0.3)

    # Panel 3: Folded light curve
    ax3 = fig.add_subplot(gs[1, 1])
    lc_folded = fold_lightcurve(lc, best_period, best_t0)
    ax3.plot(lc_folded.phase.value, lc_folded.flux.value, 'b.', markersize=2, alpha=0.5)

    # Bin for visibility
    lc_binned = lc_folded.bin(time_bin_size=0.01)
    ax3.plot(lc_binned.phase.value, lc_binned.flux.value, 'r-', linewidth=2, alpha=0.8)

    ax3.set_xlabel('Phase', fontsize=11)
    ax3.set_ylabel('Normalized Flux', fontsize=11)
    ax3.set_title(f'Folded at P = {best_period:.4f} days', fontsize=12, fontweight='bold')
    ax3.axvline(0, color='gray', linestyle=':', alpha=0.5)
    ax3.grid(alpha=0.3)

    # Panel 4: Zoomed folded transit
    ax4 = fig.add_subplot(gs[2, :])
    # Zoom to ±2 transit durations
    phase_width = 2 * best_duration / best_period
    mask = np.abs(lc_folded.phase.value) < phase_width

    ax4.plot(lc_folded.phase.value[mask], lc_folded.flux.value[mask], 'b.', markersize=3, alpha=0.6)

    # Binned version
    lc_binned_zoom = lc_folded[mask].bin(time_bin_size=0.005)
    ax4.plot(lc_binned_zoom.phase.value, lc_binned_zoom.flux.value, 'r-', linewidth=2.5, alpha=0.9)

    ax4.set_xlabel('Phase', fontsize=11)
    ax4.set_ylabel('Normalized Flux', fontsize=11)
    ax4.set_title(f'Transit Zoom (Depth = {best_depth:.4f}, Duration = {best_duration*24:.2f} hrs)',
                  fontsize=12, fontweight='bold')
    ax4.axvline(0, color='gray', linestyle=':', alpha=0.5)
    ax4.axhline(1, color='gray', linestyle=':', alpha=0.5)
    ax4.grid(alpha=0.3)

    plt.suptitle(f'TIC {tic_id} - BLS Transit Search', fontsize=15, fontweight='bold', y=0.995)

    plot_file = output_dir / f'TIC_{tic_id}_bls.png'
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    plt.close()


# ============================================================================
# Main processing loop
# ============================================================================

print("\nSearching for transit signals...")
print("-" * 70)

# Find all processed light curves
lc_files = sorted(DATA_DIR.glob('TIC_*_processed.fits'))
print(f"Found {len(lc_files)} processed light curves")

if len(lc_files) == 0:
    print("ERROR: No processed light curves found. Run 2_data_ingestion.py first.")
    exit(1)

results_summary = []

for idx, lc_file in enumerate(lc_files):
    # Extract TIC ID from filename
    tic_id = int(lc_file.stem.split('_')[1])

    print(f"\n[{idx+1}/{len(lc_files)}] TIC {tic_id}")

    # Load light curve
    try:
        lc = lk.read(lc_file)
    except Exception as e:
        print(f"  Error loading light curve: {e}")
        continue

    # Calculate maximum period
    baseline = (lc.time.value.max() - lc.time.value.min())
    max_period = MAX_PERIOD_FACTOR * baseline

    print(f"  Baseline: {baseline:.1f} days")
    print(f"  Period range: {MIN_PERIOD:.2f} - {max_period:.2f} days")

    # Run BLS
    print(f"  Running BLS...")
    try:
        bls_result = run_bls(lc, min_period=MIN_PERIOD, max_period=max_period)
    except Exception as e:
        print(f"  BLS failed: {e}")
        continue

    # Get best period
    best_idx = np.argmax(bls_result.power)
    best_period = bls_result.period[best_idx].value
    best_t0 = bls_result.transit_time[best_idx].value
    best_duration = bls_result.duration[best_idx].value
    best_depth = bls_result.depth[best_idx]
    best_power = bls_result.power[best_idx]

    # Calculate SDE
    sde = calculate_sde(bls_result)

    print(f"  Best period: {best_period:.4f} days")
    print(f"  Transit depth: {best_depth:.5f} ({best_depth*100:.3f}%)")
    print(f"  Duration: {best_duration*24:.2f} hours")
    print(f"  SDE: {sde:.2f}")

    # Calculate number of transits
    n_transits = int(baseline / best_period)

    # Estimate planet radius (assumes stellar radius from TIC)
    # R_planet = R_star * sqrt(depth)
    # For now, assume typical M dwarf R = 0.5 R_sun
    R_star = 0.5  # Solar radii (placeholder - should load from target list)
    R_planet_Rjup = R_star * np.sqrt(best_depth)
    R_planet_Rearth = R_planet_Rjup * 11.2  # 1 R_jup = 11.2 R_earth

    print(f"  Estimated planet radius: {R_planet_Rearth:.2f} R⊕")

    # Detection status
    if sde >= MIN_SDE:
        status = "DETECTED"
        print(f"  Status: {status} (SDE > {MIN_SDE})")

        # Create diagnostic plot
        print(f"  Generating diagnostic plots...")
        plot_bls_results(lc, bls_result, tic_id, PLOTS_DIR)

    else:
        status = "NOT DETECTED"
        print(f"  Status: {status} (SDE < {MIN_SDE})")

    # Get top 5 peaks for harmonic check
    top_peaks = get_top_peaks(bls_result, n_peaks=5)

    # Save results
    results_summary.append({
        'TIC': tic_id,
        'Period_days': best_period,
        'Epoch_BJD': best_t0,
        'Duration_hours': best_duration * 24,
        'Depth_ppt': best_depth * 1000,  # parts per thousand
        'Depth_percent': best_depth * 100,
        'BLS_Power': best_power,
        'SDE': sde,
        'N_transits': n_transits,
        'R_planet_Rearth': R_planet_Rearth,
        'Status': status,
        'Second_peak_period': top_peaks[1]['period'] if len(top_peaks) > 1 else np.nan,
        'Third_peak_period': top_peaks[2]['period'] if len(top_peaks) > 2 else np.nan
    })

# ============================================================================
# Save results summary
# ============================================================================

print("\n" + "="*70)
print("TRANSIT SEARCH SUMMARY")
print("="*70)

results_df = pd.DataFrame(results_summary)
results_file = OUTPUT_DIR / 'bls_results_all.csv'
results_df.to_csv(results_file, index=False, float_format='%.6f')

print(f"\nTotal stars analyzed: {len(results_df)}")

detections = results_df[results_df['Status'] == 'DETECTED']
print(f"Signals detected (SDE ≥ {MIN_SDE}): {len(detections)}")

if len(detections) > 0:
    print(f"\nDetection summary:")
    print(f"  Period range: {detections['Period_days'].min():.2f} - {detections['Period_days'].max():.2f} days")
    print(f"  Depth range: {detections['Depth_percent'].min():.3f} - {detections['Depth_percent'].max():.3f}%")
    print(f"  SDE range: {detections['SDE'].min():.2f} - {detections['SDE'].max():.2f}")
    print(f"  Planet size range: {detections['R_planet_Rearth'].min():.2f} - {detections['R_planet_Rearth'].max():.2f} R⊕")

    # Save detections to separate file
    detections_file = OUTPUT_DIR / 'bls_detections.csv'
    detections.to_csv(detections_file, index=False, float_format='%.6f')
    print(f"\nDetections saved to: {detections_file}")
    print(f"Diagnostic plots saved to: {PLOTS_DIR}")

print(f"\nFull results saved to: {results_file}")
print("="*70)
print("Next step: Run 4_vetting.py to eliminate false positives")
print("="*70)
