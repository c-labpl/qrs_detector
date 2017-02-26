import serial
import numpy as np
from scipy.signal import butter, lfilter
from collections import deque

class QRSDetector(object):
    """QRS complex detector."""

    def __init__(self, port, baud_rate):
        """Variables initialization and start reading ECG measurements."""

        # Configuration parameters.
        self.number_of_most_recent_samples = 200  # samples
        self.signal_freq = 255  # signal frequency
        self.filter_lowcut = 0.0  # band pass filter low cut value
        self.filter_highcut = 15.0  # band pass filter high cut value
        self.filter_order = 1
        # TODO: This value should be dynamic when dynamic freq will be implemented.
        self.integration_window = 15  # signal integration window length in samples
        self.findpeaks_limit = 0.40
        # TODO: This value should be dynamic when dynamic freq will be implemented.
        self.findpeaks_spacing = 50 # samples
        # TODO: This value should be dynamic when dynamic freq will be implemented.
        self.refractory_period = 120  # samples
        # TODO: This value should be dynamic when dynamic freq will be implemented.
        self.buffer_detection_window = 40  # samples
        # TODO: Rename these to more meaningful names.
        self.spk_i_measurement_weight = 0.125 # detection and thresholding params
        self.npk_i_measurement_weight = 0.125 # detection and thresholding params
        self.threshold_i_diff_weight = 0.25 # detection and thresholding params

        # Measured and calculated values.
        self.most_recent_measurements = deque([], self.number_of_most_recent_samples)  # most recent measurements array
        # TODO: Could this be time based instead of sample number based? If sample based - need to be dynamic when dynamic freq will be implemented.
        self.qrs_interval = 0  # samples
        self.spk_i = 0.0
        self.npk_i = 0.0
        self.threshold_i = 0.0

        # Run the detector.
        self.connect_to_ecg(port=port, baud_rate=baud_rate)

    # ECG interfacing methods.
    def connect_to_ecg(self, port, baud_rate):
        try:
            serial_port = serial.Serial(port, baud_rate)
            print("Connected! Starting reading ECG measurements.")
        except:
            print("Cannot connect to provided port!")
            return

        while True:
            raw_measurement = serial_port.readline()
            self.process_measurement(raw_measurement=raw_measurement)


    # Data processing methods.
    def process_measurement(self, raw_measurement):
        """Parsing raw data line."""

        # TODO: Add encoding for Python 3.
        raw_measurement_split = raw_measurement.rstrip().split(';')

        # TODO: Change it to == 2 - need to test it with real data.
        if len(raw_measurement_split) < 2:
            return
        try:
            timestamp = float(raw_measurement_split[0])
            measurement = float(raw_measurement_split[1])
        except Exception:
            return

        # Not physiologically possible ECG error measurements rejection.
        if measurement > 10:
            return

        self.most_recent_measurements.append(measurement)

        # TODO: Why?
        if len(self.most_recent_measurements) == 1:
            return

        # TODO: Maybe pass here most_recent_measurements as argument instead of calling self in method? More consistent with others.
        self.extract_peaks()


    def extract_peaks(self):
        """Proceses received data."""

        # Signal filtering - band pass 0-15 Hz.
        filtered_signal = self.bandpass_filter(self.most_recent_measurements, lowcut=self.filter_lowcut,
                                                    highcut=self.filter_highcut, signal_freq=self.signal_freq,
                                                    filter_order=self.filter_order)

        # Derivative - provides QRS slope info.
        differentiated_signal = np.ediff1d(filtered_signal)

        # Squaring.
        squared_signal = differentiated_signal ** 2

        # Moving-window integration.
        integrated_signal = np.convolve(squared_signal, np.ones(self.integration_window))

        # Fiducial mark - peak detection - integrated signal.
        detected_peaks_indices = self.findpeaks(integrated_signal, limit=self.findpeaks_limit, spacing=self.findpeaks_spacing)

        # TODO: Check whether detected_peaks_indices == fiducial_mark_idx_i - if yes use fiducial_mark_idx_i to create np.array()
        # TODO: Refactor fiducial_mark_idx_i to detected_peaks_indices and detected_peaks_values.
        # TODO: Refactor it to one 2D array with idx and vals in one array.
        fiducial_mark_idx_i = np.array([])
        fiducial_mark_val_i = np.array([])
        for detected_peak_index in detected_peaks_indices:
            fiducial_mark_idx_i = np.append(fiducial_mark_idx_i, detected_peak_index)
            fiducial_mark_val_i = np.append(fiducial_mark_val_i, integrated_signal[detected_peak_index])

        self.detect_qrs(detected_peaks_indices=fiducial_mark_idx_i, detected_peaks_values=fiducial_mark_val_i)

    # Detection methods.
    def detect_qrs(self, detected_peaks_indices, detected_peaks_values):
        """Thresholding detected peaks - integrated - signal."""
        # Check whether refractory period has passed.
        # After a valid QRS complex detection, there is a 200 ms refractory period before the next one can be detected.

        # TODO: Check whether it can be done differently.
        # TODO: Verify whether it does not break things being here.
        self.qrs_interval += 1

        # Add comment: refactory period - no detection if we are to close physiologically to last detected pulse.
        if self.qrs_interval > self.refractory_period:

            # Check whether any peak was detected in analysed samples window.
            if len(detected_peaks_indices) > 0:

                # Take the last one detected in analysed samples window as the most recent.
                # TODO: Refactor these names.
                peak_idx_i, peak_val_i = detected_peaks_indices[-1], detected_peaks_values[-1]

                # Check whether detected peak occured within defined window from end of samples buffer - making sure that the same peak will not be detected twice.
                # TODO: This is strange. Something can be done here to make it more obvious.
                # TODO: Also - why do we need this? Analyse this.
                if peak_idx_i > self.number_of_most_recent_samples - self.buffer_detection_window:

                    # Peak must be classified as a noise peak or a signal peak. To be a signal peak it must exceed threshold_i_1.
                    if peak_val_i > self.threshold_i:
                        self.qrs_interval = 0
                        self.handle_detection()
                        self.spk_i = self.spk_i_measurement_weight * peak_val_i + (1 - self.spk_i_measurement_weight) * self.spk_i
                    else:
                        self.npk_i = self.npk_i_measurement_weight * peak_val_i + (1 - self.npk_i_measurement_weight) * self.npk_i

                    self.threshold_i = self.npk_i + self.threshold_i_diff_weight * (self.spk_i - self.npk_i)


    def handle_detection(self):
        print("Pulse")

    # Tools methods.
    def bandpass_filter(self, data, lowcut, highcut, signal_freq, filter_order):
        """Constructs signal filter and uses it to given dataset."""
        nyquist_freq = 0.5 * signal_freq
        low = lowcut / nyquist_freq
        high = highcut / nyquist_freq
        b, a = butter(filter_order, [low, high], btype="band")
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


if __name__ == "__main__":
    qrs_detector = QRSDetector(port="/dev/cu.usbmodem1411", baud_rate="115200")
    qrs_detector.connect_to_ecg()
