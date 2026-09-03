import numpy as np
import rasterio
from scipy.ndimage import generic_filter
from scipy.linalg import fractional_matrix_power

CLASSES = [1, 2, 3, 4, 5]  # Built-up, Vegetation, Wetland, Water, Bareland
N_CLASSES = len(CLASSES)

BASE_YEAR = 2025
TARGET_YEARS = {
    2030: 2030 - BASE_YEAR,   # 5 years
    2050: 2050 - BASE_YEAR,   # 25 years
    2100: 2100 - BASE_YEAR,   # 75 years
}

# CA step size in years (smaller = smoother/more realistic spatial growth,
# but slower to run). 5-year steps is a good balance.
CA_STEP_YEARS = 5


# ---------------------------------------------------------------
# 1. LOAD RASTERS
# ---------------------------------------------------------------
def load_raster(path):
    with rasterio.open(path) as src:
        arr = src.read(1)
        profile = src.profile
    return arr, profile


def save_raster(arr, profile, out_path):
    profile.update(dtype=rasterio.uint8, count=1, compress='lzw', nodata=0)
    with rasterio.open(out_path, 'w', **profile) as dst:
        dst.write(arr.astype(np.uint8), 1)
    print(f"Saved: {out_path}")


# ---------------------------------------------------------------
# 2. MARKOV TRANSITION MATRIX (from two historical maps)
# ---------------------------------------------------------------
def build_transition_matrix(from_arr, to_arr, classes=CLASSES):
    """Cross-tabulate pixels: rows = class at time1, cols = class at time2."""
    n = len(classes)
    matrix = np.zeros((n, n), dtype=np.float64)
    class_index = {c: i for i, c in enumerate(classes)}

    mask = (from_arr > 0) & (to_arr > 0)  # ignore nodata
    from_valid = from_arr[mask]
    to_valid = to_arr[mask]

    for i in range(len(from_valid)):
        r = class_index.get(int(from_valid[i]))
        c = class_index.get(int(to_valid[i]))
        if r is not None and c is not None:
            matrix[r, c] += 1

    row_sums = matrix.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    prob_matrix = matrix / row_sums
    return prob_matrix


def annualize_matrix(prob_matrix, interval_years):
    """Convert an N-year transition matrix into a 1-year transition matrix."""
    annual = fractional_matrix_power(prob_matrix, 1.0 / interval_years)
    annual = np.real(annual)
    annual = np.clip(annual, 0, 1)
    row_sums = annual.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    return annual / row_sums


def project_class_proportions(current_counts, annual_matrix, n_years):
    """Project class pixel counts n_years into the future."""
    step_matrix = np.linalg.matrix_power(annual_matrix, n_years) \
        if isinstance(n_years, (int, np.integer)) else \
        np.real(fractional_matrix_power(annual_matrix, n_years))
    step_matrix = np.clip(step_matrix, 0, 1)
    row_sums = step_matrix.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    step_matrix = step_matrix / row_sums

    proportions = current_counts / current_counts.sum()
    projected = proportions @ step_matrix
    return projected  # proportions summing to 1, one per class


# ---------------------------------------------------------------
# 3. CELLULAR AUTOMATA — SPATIAL ALLOCATION
# ---------------------------------------------------------------
def neighborhood_density(arr, target_class, size=5):
    """Fraction of neighboring pixels (size x size window) equal to target_class."""
    def frac(window):
        return np.mean(window == target_class)
    return generic_filter(arr, frac, size=size, mode='nearest')


def compute_suitability(arr, target_class):
    """
    Suitability score for converting a pixel to target_class.
    Base rule: higher score if the pixel is already near that class
    (growth spreads from existing patches). Add extra terms here
    (e.g. distance-to-road, slope) if you have that data.
    """
    density = neighborhood_density(arr, target_class, size=5)
    noise = np.random.default_rng(42).random(arr.shape) * 0.05  # small stochastic term
    return density + noise


def ca_step(arr, current_counts, target_counts):
    """
    One CA iteration: reallocate pixels so that class totals move from
    current_counts toward target_counts, using neighborhood suitability
    to decide WHICH pixels change.
    """
    new_arr = arr.copy()
    diff = target_counts - current_counts  # +ve = class needs to grow

    # Classes that need to shrink become "donor" pool; classes that grow are "receivers"
    growing = [CLASSES[i] for i in range(N_CLASSES) if diff[i] > 0]
    shrinking = [CLASSES[i] for i in range(N_CLASSES) if diff[i] < 0]

    for g_idx, g_class in enumerate(growing):
        n_needed = int(round(diff[CLASSES.index(g_class)]))
        if n_needed <= 0:
            continue

        suitability = compute_suitability(new_arr, g_class)
        # Only allow conversion FROM classes that are shrinking
        eligible_mask = np.isin(new_arr, shrinking) & (new_arr != g_class)
        if not eligible_mask.any() or n_needed <= 0:
            continue

        flat_suit = suitability[eligible_mask]
        flat_idx = np.argwhere(eligible_mask)

        # Take the top-N most suitable eligible pixels
        n_take = min(n_needed, len(flat_suit))
        top_idx_positions = np.argpartition(-flat_suit, n_take - 1)[:n_take]
        chosen_pixels = flat_idx[top_idx_positions]

        for (row, col) in chosen_pixels:
            new_arr[row, col] = g_class

    return new_arr


def run_ca_markov(base_arr, annual_matrix, base_year, target_year, step_years=CA_STEP_YEARS):
    """Iteratively step the CA-Markov model from base_year to target_year."""
    arr = base_arr.copy()
    years_remaining = target_year - base_year
    current_year = base_year

    while years_remaining > 0:
        step = min(step_years, years_remaining)

        counts = np.array([np.sum(arr == c) for c in CLASSES], dtype=np.float64)
        proportions = project_class_proportions(counts, annual_matrix, step)
        target_counts = proportions * counts.sum()

        arr = ca_step(arr, counts, target_counts)

        years_remaining -= step
        current_year += step
        print(f"  CA step complete -> simulated year {current_year}")

    return arr


# ---------------------------------------------------------------
# 4. MAIN WORKFLOW
# ---------------------------------------------------------------
def main():
    # --- Load historical classified rasters (from GEE exports) ---
    arr_2006, _ = load_raster('LULC_2006.tif')
    arr_2025, profile_2025 = load_raster('LULC_2025.tif')

    # --- Build & annualize the Markov transition matrix ---
    print("Building transition matrix from 2006 -> 2025 ...")
    trans_matrix = build_transition_matrix(arr_2006, arr_2025)
    annual_matrix = annualize_matrix(trans_matrix, interval_years=2025 - 2006)
    print("Annual transition matrix:\n", np.round(annual_matrix, 4))

    # --- Project & export for each target year ---
    for year, n_years in TARGET_YEARS.items():
        print(f"\nProjecting LULC for {year} ({n_years} years from {BASE_YEAR}) ...")
        projected_arr = run_ca_markov(
            base_arr=arr_2025,
            annual_matrix=annual_matrix,
            base_year=BASE_YEAR,
            target_year=year,
        )
        save_raster(projected_arr, profile_2025.copy(), f'LULC_{year}.tif')

    print("\nDone. Outputs: LULC_2030.tif, LULC_2050.tif, LULC_2100.tif")


if __name__ == '__main__':
    main()
