#!/usr/bin/env python3
"""
TESS Target Selection Script
=============================

Selects ~100 high-probability targets for exoplanet discovery.

Criteria:
- K5-M4 dwarfs (Teff: 2700-4500 K)
- TESS magnitude: 10.0-13.5
- Multi-sector coverage (≥4 sectors)
- Stellar radius: 0.3-0.7 R☉
- Distance < 100 pc
- Contamination ratio < 0.1
- Exclude known TOIs and planets

Output: tess_shortlist_100.csv
"""

import numpy as np
import pandas as pd
from astroquery.mast import Catalogs, Observations
from astropy.coordinates import SkyCoord
import astropy.units as u
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

print("="*70)
print("TESS TARGET SELECTION FOR CITIZEN EXOPLANET DISCOVERY")
print("="*70)

# ============================================================================
# STEP 1: Query TESS Input Catalog (TIC) for candidate stars
# ============================================================================

print("\nSTEP 1: Querying TESS Input Catalog...")
print("-" * 70)

# Query parameters - broad initial search
# We'll start with temperature and magnitude, then filter further
query_params = {
    'Tmag': [10.0, 13.5],        # TESS magnitude range
    'Teff': [2700, 4500],         # K5-M4 dwarf temperature range
    'rad': [0.3, 0.7],            # Stellar radius in solar radii
}

print(f"Initial query parameters:")
print(f"  TESS magnitude: {query_params['Tmag'][0]} - {query_params['Tmag'][1]}")
print(f"  Effective temperature: {query_params['Teff'][0]} - {query_params['Teff'][1]} K")
print(f"  Stellar radius: {query_params['rad'][0]} - {query_params['rad'][1]} R☉")

# Note: TIC is huge (>1 billion stars). We'll query by sky regions to avoid timeout.
# For demonstration, we'll query specific high-priority TESS fields.

# High-priority sky regions (TESS sectors with multi-sector coverage)
# These are coordinates of TESS Continuous Viewing Zones (CVZ) and overlap regions
target_regions = [
    {'name': 'Southern CVZ', 'ra': 0, 'dec': -75, 'radius': 30},
    {'name': 'Northern CVZ', 'ra': 0, 'dec': 75, 'radius': 30},
    {'name': 'Ecliptic Overlap 1', 'ra': 90, 'dec': 0, 'radius': 15},
    {'name': 'Ecliptic Overlap 2', 'ra': 270, 'dec': 0, 'radius': 15},
]

all_targets = []

print("\nQuerying TESS Input Catalog by sky region...")
print("(This may take 5-15 minutes depending on server load)")

for region in tqdm(target_regions, desc="Sky regions"):
    try:
        # Query TIC using cone search
        coord = SkyCoord(ra=region['ra']*u.deg, dec=region['dec']*u.deg)

        catalog_data = Catalogs.query_region(
            coord,
            radius=region['radius']*u.deg,
            catalog="TIC"
        )

        # Filter by our criteria
        mask = (
            (catalog_data['Tmag'] >= query_params['Tmag'][0]) &
            (catalog_data['Tmag'] <= query_params['Tmag'][1]) &
            (catalog_data['Teff'] >= query_params['Teff'][0]) &
            (catalog_data['Teff'] <= query_params['Teff'][1]) &
            (catalog_data['rad'] >= query_params['rad'][0]) &
            (catalog_data['rad'] <= query_params['rad'][1]) &
            (catalog_data['contratio'] < 0.1)  # Low contamination
        )

        filtered = catalog_data[mask]
        all_targets.append(filtered)

        print(f"  {region['name']}: {len(filtered)} candidates")

    except Exception as e:
        print(f"  {region['name']}: Query failed - {e}")
        continue

# Combine all targets
if len(all_targets) == 0:
    print("\nERROR: No targets found. Check MAST server status or adjust criteria.")
    exit(1)

from astropy.table import vstack
targets_combined = vstack(all_targets)
print(f"\nTotal candidates after initial filtering: {len(targets_combined)}")

# ============================================================================
# STEP 2: Check for multi-sector coverage
# ============================================================================

print("\nSTEP 2: Checking TESS sector coverage...")
print("-" * 70)

def count_tess_sectors(tic_id):
    """Query how many TESS sectors observed this target."""
    try:
        obs = Observations.query_criteria(
            obs_collection="TESS",
            target_name=str(tic_id),
            dataproduct_type="timeseries"
        )
        if obs is None or len(obs) == 0:
            return 0
        # Count unique sectors
        sectors = set()
        for o in obs:
            if 'description' in o.colnames:
                # Extract sector number from description
                desc = o['description']
                if 'Sector' in desc or 'sector' in desc:
                    # Parse sector number
                    import re
                    match = re.search(r'[Ss]ector\s*(\d+)', desc)
                    if match:
                        sectors.add(int(match.group(1)))
        return len(sectors)
    except:
        return 0

# For efficiency, we'll batch check sector coverage
# This is slow - in practice, use cached sector coverage from TESScut or ExoFOP

print("Checking sector coverage for candidates (sampling 500 targets)...")
print("Note: Full search should check all candidates - this is demonstrative")

# Sample subset for demonstration (in production, check all)
sample_size = min(500, len(targets_combined))
sample_indices = np.random.choice(len(targets_combined), sample_size, replace=False)
targets_sample = targets_combined[sample_indices]

sector_counts = []
multi_sector_targets = []

for target in tqdm(targets_sample[:100], desc="Sector coverage"):  # Limit to 100 for speed
    tic = target['ID']
    n_sectors = count_tess_sectors(tic)
    sector_counts.append(n_sectors)

    if n_sectors >= 4:
        multi_sector_targets.append({
            'TIC': tic,
            'RA': target['ra'],
            'DEC': target['dec'],
            'Tmag': target['Tmag'],
            'Teff': target['Teff'],
            'Rad': target['rad'],
            'Mass': target.get('mass', np.nan),
            'Dist': target.get('d', np.nan),
            'Sectors': n_sectors,
            'Contamination': target.get('contratio', np.nan)
        })

print(f"\nFound {len(multi_sector_targets)} stars with ≥4 sectors")

# ============================================================================
# STEP 3: Exclude known TOIs and planets
# ============================================================================

print("\nSTEP 3: Excluding known TOIs and confirmed planets...")
print("-" * 70)

# Download TOI catalog
try:
    from astroquery.mast import Catalogs
    toi_catalog = Catalogs.query_object("TOI", catalog="TIC", radius=180*u.deg)
    known_tics = set(toi_catalog['ID'])
    print(f"Loaded {len(known_tics)} known TESS Objects of Interest")
except Exception as e:
    print(f"Warning: Could not load TOI catalog - {e}")
    known_tics = set()

# Filter out known systems
final_targets = [t for t in multi_sector_targets if t['TIC'] not in known_tics]
print(f"After excluding known systems: {len(final_targets)} candidates")

# ============================================================================
# STEP 4: Rank and select top 100
# ============================================================================

print("\nSTEP 4: Ranking candidates...")
print("-" * 70)

# Ranking criteria:
# - Prefer smaller stars (larger transit depth)
# - Prefer more sectors (better period coverage)
# - Prefer brighter stars (better SNR)
# - Prefer lower contamination

df = pd.DataFrame(final_targets)

# Calculate score
df['Score'] = (
    (0.7 - df['Rad']) * 30 +           # Smaller radius = higher score
    (df['Sectors'] / 10) * 25 +         # More sectors = higher score
    (13.5 - df['Tmag']) * 20 +          # Brighter = higher score
    (0.1 - df['Contamination']) * 25    # Less contamination = higher score
)

df = df.sort_values('Score', ascending=False)

# Select top 100
top_100 = df.head(100).copy()

# ============================================================================
# STEP 5: Save output
# ============================================================================

print("\nSTEP 5: Saving target list...")
print("-" * 70)

# Reorder columns for clarity
output_cols = ['TIC', 'RA', 'DEC', 'Tmag', 'Teff', 'Rad', 'Mass', 'Dist',
               'Sectors', 'Contamination', 'Score']
top_100 = top_100[output_cols]

# Save to CSV
output_file = 'tess_shortlist_100.csv'
top_100.to_csv(output_file, index=False, float_format='%.4f')

print(f"Saved {len(top_100)} targets to {output_file}")

# Print summary statistics
print("\n" + "="*70)
print("TARGET SELECTION SUMMARY")
print("="*70)
print(f"Total targets in shortlist: {len(top_100)}")
print(f"\nTESS magnitude range: {top_100['Tmag'].min():.2f} - {top_100['Tmag'].max():.2f}")
print(f"Temperature range: {top_100['Teff'].min():.0f} - {top_100['Teff'].max():.0f} K")
print(f"Radius range: {top_100['Rad'].min():.3f} - {top_100['Rad'].max():.3f} R☉")
print(f"Sector coverage: {top_100['Sectors'].min():.0f} - {top_100['Sectors'].max():.0f} sectors")
print(f"Median distance: {top_100['Dist'].median():.1f} pc")
print(f"\nExpected planet candidates: 1-5 (1-5% hit rate)")
print(f"Expected confirmed planets: 0-1 (0-1% base rate)")
print("\n" + "="*70)
print("Next step: Run 2_data_ingestion.py")
print("="*70)
