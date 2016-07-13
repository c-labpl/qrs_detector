# -*- coding: utf-8 -*-

import numpy as np
from scipy.signal import butter, lfilter


# TODO: Find a way (create a module) for keeping both online and offline algo version synced.
class QRSDetector(object):
    """QRS complex detector."""

    def __init__(self, file):
        """Variables initialization."""
        # TODO: refactor file to path
        ## Data file path.
        self.file = file

        ## Signal processing variables.
        self.raw_signal = np.array([])
        self.filtered = np.array([])
        self.differentiated_signal = np.array([])
        self.squared_signal = np.array([])
        self.integrated_signal = np.array([])

        ## Peak detection variables.
        self.fiducial_mark_val_i = np.array([])
        self.fiducial_mark_idx = np.array([])

        ## Peak thresholding variables.
        self.qrs_peak_i = np.array([])
        self.noise_peak_i = np.array([])

        ## Integrated signal detection and thresholding params.
        self.spk_i = 0.0
        self.npk_i = 0.0
        self.threshold_i_1 = 0.0
        self.threshold_i_2 = 0.0

        ## Params.
        self.signal_freq = 255  # signal frequency
        self.filter_lowcut = 0.0  # band pass filter low cut value
        self.filter_highcut = 15.0  # band pass filter high cut value
        self.integration_window = 15  # signal integration window length in samples


    ## Tool methods.

    def butter_bandpass_filter(self, data, lowcut, highcut, signal_freq, order):
        """Constructs signal filter and uses it to given dataset."""
        nyquist_freq = 0.5 * signal_freq
        low = lowcut / nyquist_freq
        high = highcut / nyquist_freq
        b, a = butter(order, [low, high], btype="band")
        y = lfilter(b, a, data)
        return y

    # Janko SLavic peak detection algorithm and implementation.
    # https://github.com/jankoslavic/py-tools/tree/master/findpeaks
    def findpeaks(self, data, spacing=1, limit=None):
        """Finds peaks in `data` which are of `spacing` width and >=`limit`.
        :param data: values
        :param spacing: minimum spacing to the next peak (should be 1 or more)
        :param limit: peaks should have value greater or equal
        :return:
        """
        len = data.size
        x = np.zeros(len + 2 * spacing)
        x[:spacing] = data[0] - 1.e-6
        x[-spacing:] = data[-1] - 1.e-6
        x[spacing:spacing + len] = data
        peak_candidate = np.zeros(len)
        peak_candidate[:] = True
        for s in range(spacing):
            start = spacing - s - 1
            h_b = x[start: start + len]  # before
            start = spacing
            h_c = x[start: start + len]  # central
            start = spacing + s + 1
            h_a = x[start: start + len]  # after
            peak_candidate = np.logical_and(peak_candidate, np.logical_and(h_c > h_b, h_c > h_a))

        ind = np.argwhere(peak_candidate)
        ind = ind.reshape(ind.size)
        if limit is not None:
            ind = ind[data[ind] > limit]
        return ind


    ## Data processing methods.

    def load_data(self):
        """Loads and cleans data."""
        with open(self.file) as f:
            content = f.readlines()
            for line in content:
                self.raw_signal = np.append(self.raw_signal, float(line.rstrip().split(';')[1]))

    def process_data(self):
        """Proceses received data."""
        # Signal filtering - band pass 0-15 Hz.
        self.filtered_signal = self.butter_bandpass_filter(self.raw_signal, lowcut=self.filter_lowcut,
                                                           highcut=self.filter_highcut, signal_freq=self.signal_freq,
                                                           order=1)

        # Derivative - provides QRS slope info.
        self.differentiated_signal = np.ediff1d(self.filtered_signal)

        # Squaring.
        self.squared_signal = self.differentiated_signal ** 2

        # Moving-window integration.
        self.integrated_signal = np.convolve(self.squared_signal, np.ones(self.integration_window))

        # Fiducial mark - peak detection - integrated signal.
        self.peaks_indices = self.findpeaks(self.integrated_signal, limit=0.30, spacing=50)

        for peak_index in self.peaks_indices:
            self.fiducial_mark_idx = np.append(self.fiducial_mark_idx, peak_index)
            self.fiducial_mark_val_i = np.append(self.fiducial_mark_val_i, self.integrated_signal[peak_index])

    def threshold_peaks(self):
        """Thresholding detected peaks - integrated - signal."""
        for peak_idx, peak_val_i in zip(self.fiducial_mark_idx, self.fiducial_mark_val_i):
            if peak_val_i > self.threshold_i_1:
                # TODO: Move these filter params to params section.
                self.spk_i = 0.125 * peak_val_i + 0.875 * self.spk_i
                self.qrs_peak_i = np.append(self.qrs_peak_i, peak_idx)
            else:
                self.npk_i = 0.125 * peak_val_i + 0.875 * self.npk_i
                self.noise_peak_i = np.append(self.noise_peak_i, peak_idx)

            self.threshold_i_1 = self.npk_i + 0.25 * (self.spk_i - self.npk_i)
            # TODO: Check whether this threshold is used.
            self.threshold_i_2 = 0.5 * self.threshold_i_1

if __name__ == "__main__":
    qrs_detector = QRSDetector("data250hz/pulse5.csv")
    qrs_detector.load_data()
    qrs_detector.process_data()
    qrs_detector.threshold_peaks()