import math
from typing import Optional
import numpy as np
from numpy.typing import ArrayLike


def _hist_bins(data: ArrayLike, bin_size: float = 30., fov_size: float = 3000.) -> np.ndarray:
    """Convenience function to calculate bins needed for 2D histogram.

    Parameters
    ----------
    data
        1D array of coordinates (either x or y).
    bin_size
        Size of bins in same units as coordinates specified by data.
    fov_size
        Size of field of view (histogram range) in same units as coordinates specified by data. Automatically set using
        range of data when set to None.

    Returns
    -------
    numpy.ndarray
        1D array containing bin edges.
    """
    bin_start = np.min(data)
    if fov_size is None:
        bin_steps = int(np.ceil((np.max(data) - np.min(data)) / bin_size)) + 1
        bin_end = bin_start + bin_steps * bin_size
    else:
        bin_end = bin_start + fov_size + bin_size

    return np.arange(bin_start, bin_end, bin_size)


def _hist1d(coordinates: np.ndarray, bins: tuple[np.ndarray, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    """Calculates a 1D histogram for the number of detections within a binned grid.

    Parameters
    ----------
    coordinates
        2D array containing xy coordinates of detections.
    bins
        2D histogram bins for each dimension, in same units as coordinates specified by data.
    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Bins and corresponding histogram for detection counts per pixel of binned image.
    """
    # Create binned image and flatten to 1D array.
    hist_img = np.histogram2d(coordinates[:, 0], coordinates[:, 1], bins)[0].flatten().astype(int)
    # Remove pixels with zero detections.
    hist_img = hist_img[hist_img != 0]
    # Length of 1D histogram equal to maximum number of detections across all pixels.
    hist_n_bins = np.max(hist_img)
    # Bins for 1D histogram have spacing of 1 detection and centred on integer values.
    return np.histogram(hist_img, bins=np.arange(0.5, stop=hist_n_bins + 1, step=1))


def freq_hist(coordinates: np.ndarray, bin_size: float = 30., fov_size: Optional[float] = 3000.,
              thinning_density: Optional[float] = None, thinning_repeats: int = 10) -> np.ndarray:
    """Calculates a 1D histogram for the number of detections within a binned grid with optional data thinning.

    First calculates a 2D image/histogram from a 2D point cloud. Each pixel in this image specifies the number of
    detections. A 1D histogram of detection counts per pixel is calculated from the image excluding pixels with zero
    detections. Before calculating 2D histograms the point cloud can be thinned (by randomly removing detections) to a
    specified density.

    Parameters
    ----------
    coordinates
        2D array containing xy coordinates of detections.
    bin_size
        Size of bins in same units as coordinates specified by data.
    fov_size
        Size of field of view (histogram range) in same units as coordinates specified by data. Automatically set using
        range of data when set to None.
    thinning_density
        Optional float specifying detection density to randomly thin data to.
    thinning_repeats
        Number of times the data thinning is repeated.
    Returns
    -------
    numpy.ndarray
        Histogram for detection counts per pixel of binned image.
    """

    # Bins for  image / 2D histogram
    bins = (_hist_bins(coordinates[:, 0], bin_size, fov_size),
            _hist_bins(coordinates[:, 1], bin_size, fov_size))

    if thinning_density:
        if not fov_size:
            # Area of smallest rectangle which contains all detections.
            detection_coverage = (np.max(coordinates[:, 0]) - np.min(coordinates[:, 0])) * (
                    np.max(coordinates[:, 1]) - np.min(coordinates[:, 1]))
        else:
            detection_coverage = fov_size * fov_size
        # Current density of detections.
        detection_density = coordinates.shape[0] / detection_coverage
        # Calculate the number of detections needed to get the specified thinning density.
        n_detections_keep = int(round(thinning_density * detection_coverage))
        # Check the current detection density isn't less than the thinning density.
        assert detection_density >= thinning_density, (f"Detection density: {detection_density:.03g}, less than "
                                                       f"specified thinning density: {thinning_density:.03g}")

        # Create the specified number of 1D histograms
        # using randomly selected detections from the original point cloud.
        histograms = [
            _hist1d(coordinates[np.random.choice(coordinates.shape[0], n_detections_keep, replace=False), :], bins)[0]
            for i in range(thinning_repeats)]
        # Pad the end of the histograms with zeros such that they are all the same length.
        histograms_max_length = np.max([histograms[i].size for i in range(thinning_repeats)])
        histograms = [np.pad(hist, (0, histograms_max_length - hist.size), "constant") for hist in histograms]
        # Return the average histogram
        return np.mean(np.stack(histograms, axis=1), axis=1)

    else:
        # Calculate and return the 1D histogram using all the detections from the original point cloud.
        return _hist1d(coordinates, bins)[0]


def ks_similarity(hists: tuple[np.ndarray, np.ndarray], norm: bool = False, alpha: float = 0.05) -> float:
    """Calculates the Kolmogorov-Smirnov (KS) similarity metric.

    Compares two normalised cumulative frequency histograms to calculate the KS metric. Values greater than one reject
    the null hypothesis at level alpha. Optional normalisation mode to account for offsets in total image intensity
    (total number of detections). In this case the brighter histogram is rebuilt with a bin size larger than 1
    detection before computing the KS metric.

    Parameters
    ----------
    hists
        Tuple containing two histograms for comparison.
    norm
        Should normalisation mode to account for offsets in intensity be used.
    alpha
        Significance level for rejecting the null hypothesis.
    Returns
    -------
    float
        The Kolmogorov-Smirnov (KS) dissimilarity metric. Values greater than one reject the null hypothesis.
    """
    assert len(hists) == 2

    if norm:
        hists = list(hists)
        # Expand histogram back to flattened 2D binned histogram.
        hist_expanded = [np.concatenate([(j + 1) * np.ones(hists[i][j]) for j in range(len(hists[i]))]) for i in
                         range(2)]
        # Total histogram intensity / number of detections is the sum of the expanded histogram.
        hist_intensity = [np.sum(hist_expanded[i]) for i in range(2)]
        # Make sure the brighter histogram is the second one, if not reverse the order.
        if hist_intensity[0] > hist_intensity[1]:
            hists.reverse()
            hist_expanded.reverse()
            hist_intensity.reverse()
        # The ratio between the two intensities defined the bin size for the normalised histogram.
        norm_fac = hist_intensity[1] / hist_intensity[0]
        # Compute bins for the normalised histogram.
        bin_steps = int(np.ceil(len(hists[1]) / norm_fac)) + 1
        bin_end = 1 + bin_steps * norm_fac
        bins_norm = np.arange(1, bin_end, norm_fac)
        # Recalculate histogram with new bins.
        hists[1] = np.histogram(hist_expanded[1], bins=bins_norm)[0]

    # Calculate cumulative histograms
    hist_cum = [np.cumsum(hists[i]) for i in range(2)]
    # Sample size is last value in cumulative histogram (total number of pixels excluding zero counts).
    sample_size = [hist_cum[i][-1] for i in range(2)]
    # Extend cumulative histograms by padding with last value such that they are the same length.
    hist_len = max(hist_cum[0].size, hist_cum[1].size)
    hist_cum = [np.pad(hist_cum[i], (0, max(hist_len - hist_cum[i].size, 0)), "edge") for i in range(2)]

    # Calculate KS metric using normalised cumulative histograms
    mf = np.max(np.abs(hist_cum[0] / sample_size[0] - hist_cum[1] / sample_size[1]))
    J = (sample_size[0] + sample_size[1]) / (sample_size[0] * sample_size[1])
    c = math.sqrt(-math.log(alpha / 2) / 2)
    return mf / math.sqrt(J) / c
