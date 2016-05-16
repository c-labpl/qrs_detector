import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter, find_peaks_cwt

class QRSDetector():
    """QRS complex detector."""

    def __init__(self, file):
        """Variables initialization."""
        ## Data file path.
        self.file = file
        
        ## Signal processing variables.
        self.raw_signal = np.array([])
        self.filtered = np.array([])
        self.differentiated_signal = np.array([])
        self.squared_signal = np.array([])
        self.integrated_signal = np.array([])

        # Peak detection variables.
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

        ## Moving-window integration params.
        self.integration_window = 3
        
    ## Tool methods.
    def butter_bandpass_filter(self, data, lowcut, highcut, signal_freq, order):
        """Constructs signal filter."""
        nyquist_freq = 0.5 * signal_freq
        low = lowcut / nyquist_freq
        high = highcut / nyquist_freq
        b, a = butter(order, [low, high], btype='band')
        y = lfilter(b, a, data)
        return y
    
    ## Data processing methods.
    def load_data(self):
        """Loads and cleans data."""
        with open(self.file) as f:
            content = f.readlines()
        for line in content:
            self.raw_signal = np.append(self.raw_signal, float(line.rstrip().split(';')[1]))

        # for data in reading:
        #     self.process_data()
        #     self.detector.threshold_peaks()
    
    def process_data(self):
        """Process received data."""
        ## Signal filtering - pass band 0 - 12 Hz.
        self.filtered_signal = self.butter_bandpass_filter(self.raw_signal, lowcut=0.0, highcut=12.0, signal_freq=84.0, order=1)
        
        ## Derivative - provides QRS slope info.
        self.differentiated_signal = np.diff(self.filtered_signal)

        ## Squaring.
        self.squared_signal = np.power(self.differentiated_signal, 2)

        ## Moving-window integration.
        self.integrated_signal = np.convolve(self.squared_signal, np.ones((self.integration_window,)) / self.integration_window)

        ## Fiducial mark - peak detection - integrated signal
        self.peaks_indices = find_peaks_cwt(self.integrated_signal[:-1], np.arange(10, 15), noise_perc=0.1)

        for peak_index in self.peaks_indices:
            self.fiducial_mark_idx = np.append(self.fiducial_mark_idx, peak_index)
            self.fiducial_mark_val_i = np.append(self.fiducial_mark_val_i, self.integrated_signal[peak_index])
    
    def threshold_peaks(self):
        """Thresholding detect peaks - integrated signal."""
        for peak_idx, peak_val_i in zip(self.fiducial_mark_idx, self.fiducial_mark_val_i):
            if peak_val_i > self.threshold_i_1:
                self.spk_i = 0.125 * peak_val_i + 0.875 * self.spk_i
                self.qrs_peak_i = np.append(self.qrs_peak_i, peak_idx)

            else:
                self.npk_i = 0.125 * peak_val_i + 0.875 * self.npk_i
                self.noise_peak_i = np.append(self.noise_peak_i, peak_idx)

            self.threshold_i_1 = self.npk_i + 0.25 * (self.spk_i - self.npk_i)
            self.threshold_i_2 = 0.5 * self.threshold_i_1 
        # print self.qrs_peak_i

if __name__ == "__main__":
    qrs_detector = QRSDetector("data/pulse1.csv")
    qrs_detector.load_data()
    # qrs_detector.process_data()
    # qrs_detector.threshold_peaks()
