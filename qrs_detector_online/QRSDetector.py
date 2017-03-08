import serial
import numpy as np
from scipy.signal import butter, lfilter
from collections import deque
from time import gmtime, strftime

LOG_DIR = "logs/"

# TODO: Class documentation
class QRSDetector(object):
    """QRS complex detector."""
    # TODO: Write something about which device it is dedicated to and with what params.

    def __init__(self, port, baud_rate):
        """
        QRSDetector class initialisation method.
        :param str port: port to which ECG device is connected
        :param str baud_rate: baud rate of data received from ECG device
        """
        # Configuration parameters.
        self.signal_frequency = 250                 # Set ECG device frequency in samples per second.

        self.number_of_samples_stored = 200         # Change proportionally when adjusting frequency (in samples).
        self.possible_measurement_upper_limit = 10  # ECG device physiologically possible upper measurement limit.

        self.filter_lowcut = 0.0
        self.filter_highcut = 15.0
        self.filter_order = 1

        self.integration_window = 15                # Change proportionally when adjusting frequency (in samples).

        self.findpeaks_limit = 0.30
        self.findpeaks_spacing = 50                 # Change proportionally when adjusting frequency (in samples).
        self.detection_window = 40                  # Change proportionally when adjusting frequency (in samples).

        self.refractory_period = 120                # Change proportionally when adjusting frequency (in samples).
        self.signal_peak_filtering_factor = 0.125
        self.noise_peak_filtering_factor = 0.125
        self.signal_noise_diff_weight = 0.25

        # Measurements and calculated values.
        self.timestamp = 0
        self.measurement = 0
        self.detected_qrs = 0
        self.most_recent_measurements = deque([0], self.number_of_samples_stored)  # most recent measurements array
        self.samples_since_last_detected_qrs = 0
        self.signal_peak_value = 0.0
        self.noise_peak_value = 0.0
        self.threshold_value = 0.0

        # Data logging.
        self.log_path = "{:s}QRS_detector_log_{:s}.csv".format(LOG_DIR, strftime("%Y_%m_%d_%H_%M_%S", gmtime()))
        self.log_data(self.log_path, "timestamp,measurement,qrs_detected\n")

        # Connect to ECG device and start the detector.
        self.connect_to_ecg(port=port, baud_rate=baud_rate)

    """Setting connection to ECG device methods."""

    def connect_to_ecg(self, port, baud_rate):
        """
        Method responsible for connecting to ECG device and starting reading ECG measurements.
        :param str port: port to which ECG device is connected
        :param str baud_rate: baud rate of data received from ECG device
        """
        try:
            serial_port = serial.Serial(port, baud_rate)
            print("Connected! Starting reading ECG measurements.")
        except:
            print("Cannot connect to provided port!")
            return

        while True:
            raw_measurement = serial_port.readline()
            self.process_measurement(raw_measurement=raw_measurement)

            self.log_data(self.log_path, "{:d},{:.10f},{:d}\n".format(int(self.timestamp), self.measurement, self.detected_qrs))

    """Measured data processing methods."""

    def process_measurement(self, raw_measurement):
        """
        Method responsible for parsing and initial processing of ECG measured data sample.
        :param str raw_measurement: ECG most recent raw measurement in "timestamp;measurement" format
        """
        raw_measurement_split = raw_measurement.decode().rstrip().split(';')

        # Parsing raw ECG data - modify this part in accordance to your device data format.
        if len(raw_measurement_split) != 2:
            return
        try:
            self.detected_qrs = 0
            self.timestamp = float(raw_measurement_split[0])
            self.measurement = float(raw_measurement_split[1])
        except Exception:
            return

        # Not physiologically possible ECG measurements rejection - filtering out device measurements errors.
        if self.measurement > self.possible_measurement_upper_limit:
            return

        # Appending measurements to deque used for rotating most recent samples for further analysis and detection.
        self.most_recent_measurements.append(self.measurement)

        self.extract_peaks(self.most_recent_measurements)

    def extract_peaks(self, most_recent_measurements):
        """
        Method responsible for extracting peaks from recently received ECG measurements data through signal processing.
        :param deque most_recent_measurements: most recent ECG measurements array
        """
        # Signal filtering - 0-15 Hz band pass filter.
        filtered_signal = self.bandpass_filter(most_recent_measurements, lowcut=self.filter_lowcut,
                                               highcut=self.filter_highcut, signal_freq=self.signal_frequency,
                                               filter_order=self.filter_order)

        # Derivative - provides QRS slope information.
        differentiated_signal = np.ediff1d(filtered_signal)

        # Squaring - intensifies values received in derivative.
        squared_signal = differentiated_signal**2

        # Moving-window integration.
        integrated_signal = np.convolve(squared_signal, np.ones(self.integration_window))

        # Fiducial mark - peak detection on integrated signal.
        detected_peaks_indices = self.findpeaks(integrated_signal, limit=self.findpeaks_limit, spacing=self.findpeaks_spacing)
        detected_peaks_indices = detected_peaks_indices[detected_peaks_indices > self.number_of_samples_stored - self.detection_window]
        detected_peaks_values = integrated_signal[detected_peaks_indices]

        self.detect_qrs(detected_peaks_indices=detected_peaks_indices, detected_peaks_values=detected_peaks_values)

    """Detection methods."""

    def detect_qrs(self, detected_peaks_indices, detected_peaks_values):
        """
        Method responsible for classifying detected ECG signal peaks either as noise or as QRS complex (heart beat).
        :param array detected_peaks_indices: detected peaks indices array
        :param array detected_peaks_values: detected peaks values array
        """
        self.samples_since_last_detected_qrs += 1

        # After a valid QRS complex detection, there is a 200 ms refractory period before next one can be detected.
        if self.samples_since_last_detected_qrs > self.refractory_period:

            # Check whether any peak was detected in analysed samples window.
            if len(detected_peaks_indices) > 0:

                # Take the last one detected in analysed samples window as the most recent.
                most_recent_peak_idx, most_recent_peak_value = detected_peaks_indices[-1], detected_peaks_values[-1]

                # Peak must be classified either as a noise peak or a signal peak.
                # To be classified as a signal peak (QRS peak) it must exceed dynamically set threshold value.
                if most_recent_peak_value > self.threshold_value:
                    self.handle_detection()
                    self.detected_qrs = 1
                    self.samples_since_last_detected_qrs = 0

                    # Adjust signal peak value used later for setting QRS-noise threshold.
                    self.signal_peak_value = self.signal_peak_filtering_factor * most_recent_peak_value + (1 - self.signal_peak_filtering_factor) * self.signal_peak_value
                else:
                    # Adjust noise peak value used later for setting QRS-noise threshold.
                    self.noise_peak_value = self.noise_peak_filtering_factor * most_recent_peak_value + (1 - self.noise_peak_filtering_factor) * self.noise_peak_value

                # Adjust QRS-noise threshold value based on previously detected QRS or noise peaks value.
                self.threshold_value = self.noise_peak_value + self.signal_noise_diff_weight * (self.signal_peak_value - self.noise_peak_value)

    def handle_detection(self):
        """
        Method responsible for generating any kind of response for detected QRS complex.
        """
        print("Pulse")

    """Tools methods."""

    def log_data(self, path, data):
        """
        Method responsible for logging measured ECG and detection results to a log file.
        :param str path: path to a log file
        :param str data: data line to log
        """
        with open(path, "a") as fin:
            fin.write(data)

    def bandpass_filter(self, data, lowcut, highcut, signal_freq, filter_order):
        """
        Method responsible for creating and applying Butterworth digital filter for received ECG signal.
        :param array data: raw data
        :param int lowcut: filter lowcut frequency value
        :param int highcut: filter highcut frequency value
        :param int signal_freq: signal frequency in samples per second (Hz)
        :param int filter_order: filter order
        :return array: filtered data
        """
        """Constructs signal filter and uses it to given dataset."""
        nyquist_freq = 0.5 * signal_freq
        low = lowcut / nyquist_freq
        high = highcut / nyquist_freq
        b, a = butter(filter_order, [low, high], btype="band")
        y = lfilter(b, a, data)
        return y

    def findpeaks(self, data, spacing=1, limit=None):
        """
        Janko Slavic peak detection algorithm and implementation.
        https://github.com/jankoslavic/py-tools/tree/master/findpeaks
        Finds peaks in `data` which are of `spacing` width and >=`limit`.
        :param array data: data
        :param int spacing: minimum spacing to the next peak (should be 1 or more)
        :param int limit: peaks should have value greater or equal
        :return array: detected peaks indexes array
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
    qrs_detector = QRSDetector(port="/dev/cu.usbmodem14311", baud_rate="115200")
