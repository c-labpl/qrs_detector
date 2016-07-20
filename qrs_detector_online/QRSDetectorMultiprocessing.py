# -*- coding: utf-8 -*-
import numpy as np
from scipy.signal import butter, lfilter
import serial
from collections import deque
from Logger import Logger
from AudioPlayer import AudioPlayer


class QRSDetectorMultiprocessing(object):
    """QRS complex detector."""

    def __init__(self, port, baud_rate, play_sound):
        """Variables initialization."""

        # General params.
        self.signal_freq = 255  # signal frequency
        self.filter_lowcut = 0.0  # band pass filter low cut value
        self.filter_highcut = 15.0  # band pass filter high cut value
        self.integration_window = 15  # signal integration window length in samples

        # Realtime params.
        self.cycling_window = 100  # samples
        self.r_interval = 0  # samples
        self.refractory_period = 50  # samples
        self.buffer_detection_window = 10  # samples

        # Connection details.
        self.port = port
        self.baud_rate = baud_rate
        self.serial = None
        self.update_data = True

        # Received data.
        self.data_line = None
        self.timestamp = None
        self.measurement = None

        # Signal processing variables.
        self.raw_signal = deque([], self.cycling_window)
        self.filtered = np.array([])
        self.differentiated_signal = np.array([])
        self.squared_signal = np.array([])
        self.integrated_signal = np.array([])

        # Integrated signal detection and thresholding params.
        self.spk_i = 0.0
        self.npk_i = 0.0
        self.threshold_i = 0.0
        self.spk_i_measurement_weight = 0.125
        self.npk_i_measurement_weight = 0.125
        self.threshold_i_diff_weight = 0.25

        # Peak thresholding variables.
        self.qrs_peak_i = np.array([])
        self.noise_peak_i = np.array([])
        self.detected_beat_indicator = 0

        # Data logger set up.
        self.logger = Logger("QRS", " ", "timestamp", "ecg", "beat_detected")

        # Audio player set up.
        self.play_sound = play_sound
        self.player = AudioPlayer(file_path="audio/beep.wav")


    # Lifecycle handling methods - public interface.

    def connect_to_arduino(self):
        self.serial = serial.Serial(self.port, self.baud_rate)
        print "Connected!"

    def start_updating_data(self):
        print "Detecting!"
        self.update_data = True
        while self.update_data:
            self.data_line = self.serial.readline()
            self.process_line()

    def stop_updating_data(self):
        self.update_data = False

    def disconnect_arduino(self):
        self.serial.close()


    # Data processing methods.

    def process_line(self):
        """Parsing raw data line."""

        self.detected_beat_indicator = 0

        update = self.data_line.rstrip().split(';')

        if len(update) < 2:
            return
        try:
            self.timestamp = float(update[0])
            self.measurement = float(update[1])
        except Exception:
            return

        self.r_interval += 1
        self.raw_signal.append(self.measurement)

        if len(self.raw_signal) == 1:
            return

        self.process_data()

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

        ## Peak detection variables.
        self.fiducial_mark_idx_i = np.array([])
        self.fiducial_mark_val_i = np.array([])

        for peak_index in self.peaks_indices:
            self.fiducial_mark_idx_i = np.append(self.fiducial_mark_idx_i, peak_index)
            self.fiducial_mark_val_i = np.append(self.fiducial_mark_val_i, self.integrated_signal[peak_index])

        self.threshold_peaks()

    def threshold_peaks(self):
        """Thresholding detected peaks - integrated - signal."""
        # Check whether refractory period has passed.
        # After a valid QRS complex detection, there is a 200 ms refractory period before the next one can be detected.
        if self.r_interval > self.refractory_period:

            # Check whether any peak was detected in analysed samples window.
            if len(self.fiducial_mark_idx_i) > 0:

                # Take the last one detected in analysed samples window.
                peak_idx_i, peak_val_i = self.fiducial_mark_idx_i[-1], self.fiducial_mark_val_i[-1]

                # Check whether detected peak occured within defined window from end of samples buffer - making sure that the same peak will not be detected twice.
                if peak_idx_i > self.cycling_window - self.buffer_detection_window:

                    # Peak must be classified as a noise peak or a signal peak. To be a signal peak it must exceed threshold_i_1.
                    if peak_val_i > self.threshold_i:
                        print "QRS detected!"
                        self.handle_detection()
                        self.detected_beat_indicator = 1
                        self.r_interval = 0

                        self.spk_i = self.spk_i_measurement_weight * peak_val_i + (
                                                                                      1 - self.spk_i_measurement_weight) * self.spk_i
                        self.qrs_peak_i = np.append(self.qrs_peak_i, peak_idx_i)
                    else:
                        print "NOISE detected!"
                        self.npk_i = self.npk_i_measurement_weight * peak_val_i + (
                                                                                      1 - self.npk_i_measurement_weight) * self.npk_i
                        self.noise_peak_i = np.append(self.noise_peak_i, peak_idx_i)

                    self.threshold_i = self.npk_i + self.threshold_i_diff_weight * (self.spk_i - self.npk_i)

        self.logger.log(str(self.timestamp), str(self.measurement), str(self.detected_beat_indicator))


    # Tool methods.

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

    def handle_detection(self):
        if self.play_sound:
            self.player.play()


if __name__ == "__main__":
    qrs_detector = QRSDetectorMultiprocessing(port="/dev/cu.usbmodem1411", baud_rate="115200", play_sound=True)
    qrs_detector.connect_to_arduino()
    qrs_detector.start_updating_data()
