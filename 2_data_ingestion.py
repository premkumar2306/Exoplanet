#!/usr/bin/env python3
"""
TESS Data Ingestion and Preprocessing
======================================

Downloads TESS light curves for shortlisted targets and performs preprocessing:
- Download all sectors for each TIC ID
- Stitch sectors together
- Normalize flux
- Remove systematic trends
- Sigma-clip outliers
- Save processed light curves

Input: tess_shortlist_100.csv
Output: data/processed/ directory with individual light curve files
"""

import numpy as np
import pandas as pd
import lightkurve as lk
from astropy.io import fits
from astropy.time import Time
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

print("="*70)
print("TESS DATA INGESTION AND PREPROCESSING")
print("="*70)

# ============================================================================
# Configuration
# ============================================================================

INPUT_FILE = 'tess_shortlist_100.csv'
OUTPUT_DIR = Path('data/processed')
RAW_DIR = Path('data/raw')
PLOTS_DIR = Path('data/plots/ingestion')

# Create directories
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
RAW_DIR.mkdir(parents=True, exist_ok=True)
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

# Processing parameters
FLATTEN_WINDOW = 0.75  # days - removes trends longer than this
SIGMA_CLIP = 5.0       # standard deviations for outlier removal
MIN_DATA_POINTS = 1000 # minimum points after cleaning

# ============================================================================
# STEP 1: Load target list
# ============================================================================

print("\nSTEP 1: Loading target list...")
print("-" * 70)

try:
    targets = pd.read_csv(INPUT_FILE)
    print(f"Loaded {len(targets)} targets from {INPUT_FILE}")
except FileNotFoundError:
    print(f"ERROR: {INPUT_FILE} not found. Run 1_target_selection.py first.")
    exit(1)

# ============================================================================
# STEP 2: Download and process each target
# ============================================================================

print("\nSTEP 2: Downloading and processing light curves...")
print("-" * 70)

def download_tess_lightcurve(tic_id):
    """
    Download all available TESS light curves for a given TIC ID.

    Returns
    -------
    lc_collection : lightkurve.LightCurveCollection or None
        Collection of light curves from all sectors
    """
    try:
        # Search for all TESS observations
        search_result = lk.search_lightcurve(
            f'TIC {tic_id}',
            author='SPOC',        # Use SPOC pipeline (official NASA)
            exptime=120           # 2-minute cadence (120 seconds)
        )

        if len(search_result) == 0:
            print(f"  No data found for TIC {tic_id}")
            return None

        # Download all sectors
        lc_collection = search_result.download_all()

        if lc_collection is None:
            return None

        print(f"  Downloaded {len(lc_collection)} sectors")
        return lc_collection

    except Exception as e:
        print(f"  Download failed: {e}")
        return None


def stitch_and_normalize(lc_collection):
    """
    Stitch multiple TESS sectors together and normalize.

    Process:
    1. Remove NaN values from each sector
    2. Normalize each sector to median flux = 1
    3. Concatenate all sectors
    4. Return combined light curve

    Returns
    -------
    lc_stitched : lightkurve.LightCurve
        Combined, normalized light curve
    """
    cleaned_lcs = []

    for lc in lc_collection:
        # Remove NaNs
        lc = lc.remove_nans()

        if len(lc) < 100:  # Skip sectors with too little data
            continue

        # Normalize to median flux = 1
        lc = lc.normalize()

        cleaned_lcs.append(lc)

    if len(cleaned_lcs) == 0:
        return None

    # Stitch all sectors together
    lc_stitched = cleaned_lcs[0]
    for lc in cleaned_lcs[1:]:
        lc_stitched = lc_stitched.append(lc)

    return lc_stitched


def remove_trends(lc, window_length=0.75):
    """
    Remove long-term trends while preserving short-term signals (transits).

    Uses Savitzky-Golay filter to model and remove trends longer than window_length.

    Parameters
    ----------
    lc : lightkurve.LightCurve
        Input light curve
    window_length : float
        Window length in days for trend removal

    Returns
    -------
    lc_flat : lightkurve.LightCurve
        Flattened light curve with trends removed
    """
    try:
        # Lightkurve's flatten method uses Savitzky-Golay filter
        # Window length in days - removes trends longer than this
        lc_flat = lc.flatten(window_length=window_length)
        return lc_flat
    except Exception as e:
        print(f"  Flattening failed: {e}")
        return lc


def sigma_clip_outliers(lc, sigma=5.0):
    """
    Remove outliers using sigma clipping.

    Conservative threshold (5σ) to avoid clipping real transits.

    Parameters
    ----------
    lc : lightkurve.LightCurve
        Input light curve
    sigma : float
        Number of standard deviations for clipping

    Returns
    -------
    lc_clipped : lightkurve.LightCurve
        Light curve with outliers removed
    """
    # Use lightkurve's remove_outliers method
    lc_clipped = lc.remove_outliers(sigma=sigma)
    return lc_clipped


def save_processed_lightcurve(lc, tic_id, output_dir):
    """
    Save processed light curve to file.

    Saves as both FITS (standard format) and CSV (easy to read).

    Parameters
    ----------
    lc : lightkurve.LightCurve
        Processed light curve
    tic_id : int
        TIC identifier
    output_dir : Path
        Output directory
    """
    # Save as FITS
    fits_file = output_dir / f'TIC_{tic_id}_processed.fits'
    lc.to_fits(fits_file, overwrite=True)

    # Save as CSV for easy inspection
    csv_file = output_dir / f'TIC_{tic_id}_processed.csv'
    df = pd.DataFrame({
        'time': lc.time.value,
        'flux': lc.flux.value,
        'flux_err': lc.flux_err.value if hasattr(lc, 'flux_err') else np.nan
    })
    df.to_csv(csv_file, index=False)

    return fits_file


def plot_processing_steps(lc_raw, lc_processed, tic_id, plots_dir):
    """
    Create diagnostic plot showing processing steps.

    Parameters
    ----------
    lc_raw : lightkurve.LightCurve
        Raw stitched light curve
    lc_processed : lightkurve.LightCurve
        Fully processed light curve
    tic_id : int
        TIC identifier
    plots_dir : Path
        Directory for plots
    """
    fig, axes = plt.subplots(2, 1, figsize=(14, 8))

    # Raw light curve
    axes[0].plot(lc_raw.time.value, lc_raw.flux.value, 'k.', markersize=1, alpha=0.5)
    axes[0].set_ylabel('Normalized Flux', fontsize=12)
    axes[0].set_title(f'TIC {tic_id} - Raw Stitched Light Curve', fontsize=14, fontweight='bold')
    axes[0].grid(alpha=0.3)

    # Processed light curve
    axes[1].plot(lc_processed.time.value, lc_processed.flux.value, 'b.', markersize=1, alpha=0.7)
    axes[1].set_xlabel('Time (BJD - 2457000)', fontsize=12)
    axes[1].set_ylabel('Normalized Flux', fontsize=12)
    axes[1].set_title('Processed (Detrended + Sigma Clipped)', fontsize=14, fontweight='bold')
    axes[1].grid(alpha=0.3)

    plt.tight_layout()
    plot_file = plots_dir / f'TIC_{tic_id}_processing.png'
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    plt.close()


# ============================================================================
# Main processing loop
# ============================================================================

processing_log = []

print(f"\nProcessing {len(targets)} targets...")
print("(This will take 30-90 minutes depending on server speed)")
print("-" * 70)

for idx, row in targets.iterrows():
    tic_id = int(row['TIC'])
    print(f"\n[{idx+1}/{len(targets)}] TIC {tic_id}")

    # Check if already processed
    output_file = OUTPUT_DIR / f'TIC_{tic_id}_processed.fits'
    if output_file.exists():
        print(f"  Already processed - skipping")
        processing_log.append({
            'TIC': tic_id,
            'Status': 'Already processed',
            'N_points': 'N/A',
            'Time_span_days': 'N/A'
        })
        continue

    # Download
    lc_collection = download_tess_lightcurve(tic_id)
    if lc_collection is None:
        processing_log.append({
            'TIC': tic_id,
            'Status': 'Download failed',
            'N_points': 0,
            'Time_span_days': 0
        })
        continue

    # Stitch and normalize
    print(f"  Stitching sectors...")
    lc_raw = stitch_and_normalize(lc_collection)
    if lc_raw is None:
        processing_log.append({
            'TIC': tic_id,
            'Status': 'Stitching failed',
            'N_points': 0,
            'Time_span_days': 0
        })
        continue

    # Remove trends
    print(f"  Removing trends (window={FLATTEN_WINDOW} days)...")
    lc_flat = remove_trends(lc_raw, window_length=FLATTEN_WINDOW)

    # Sigma clip
    print(f"  Sigma clipping outliers ({SIGMA_CLIP}σ)...")
    lc_processed = sigma_clip_outliers(lc_flat, sigma=SIGMA_CLIP)

    # Quality check
    n_points = len(lc_processed)
    if n_points < MIN_DATA_POINTS:
        print(f"  WARNING: Only {n_points} points after processing (min={MIN_DATA_POINTS})")
        processing_log.append({
            'TIC': tic_id,
            'Status': 'Insufficient data',
            'N_points': n_points,
            'Time_span_days': 0
        })
        continue

    # Calculate time span
    time_span = (lc_processed.time.value.max() - lc_processed.time.value.min())

    # Save
    print(f"  Saving processed light curve ({n_points} points, {time_span:.1f} days)...")
    save_processed_lightcurve(lc_processed, tic_id, OUTPUT_DIR)

    # Create diagnostic plot
    print(f"  Creating diagnostic plot...")
    plot_processing_steps(lc_raw, lc_processed, tic_id, PLOTS_DIR)

    processing_log.append({
        'TIC': tic_id,
        'Status': 'Success',
        'N_points': n_points,
        'Time_span_days': time_span
    })

    print(f"  ✓ Complete")

# ============================================================================
# STEP 3: Save processing log
# ============================================================================

print("\n" + "="*70)
print("PROCESSING SUMMARY")
print("="*70)

log_df = pd.DataFrame(processing_log)
log_file = 'data/processing_log.csv'
log_df.to_csv(log_file, index=False)

# Summary statistics
status_counts = log_df['Status'].value_counts()
print("\nProcessing results:")
for status, count in status_counts.items():
    print(f"  {status}: {count}")

successful = log_df[log_df['Status'] == 'Success']
if len(successful) > 0:
    print(f"\nSuccessfully processed: {len(successful)} targets")
    print(f"Average data points: {successful['N_points'].mean():.0f}")
    print(f"Average time span: {successful['Time_span_days'].mean():.1f} days")
    print(f"\nProcessed light curves saved to: {OUTPUT_DIR}")
    print(f"Diagnostic plots saved to: {PLOTS_DIR}")
    print(f"\n{'='*70}")
    print("Next step: Run 3_transit_search.py")
    print("="*70)
else:
    print("\nWARNING: No targets successfully processed.")
    print("Check MAST server status or target list quality.")
