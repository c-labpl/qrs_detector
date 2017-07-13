# Python Online and Offline ECG QRS Detector based on the Pan-Tomkins algorithm 

## Intro

The modules published in this repository are Python implementations of online and offline QRS complex detectors in ECG signal, based on the Pan-Tomkins algorithm (Pan J., Tompkins W. J., _A real-time QRS detection algorithm,_ IEEE Transactions on Biomedical Engineering, Vol. BME-32, No. 3, March 1985, pp. 230-236).

The QRS complex corresponds to the depolarization of the right and left ventricles of the human heart. It is the most visually obvious part of the ECG signal. QRS complex detection is essential for time-domain ECG signal analyses, namely heart rate variability. It makes it possible to compute inter-beat interval (RR interval) values that correspond to the time between two consecutive R peaks. Thus, a QRS complex detector is an ECG-based heart contraction detector.

You can find out more about cardiac cycle and QRS complex [here](https://en.wikipedia.org/wiki/Cardiac_cycle) and [here](https://en.wikipedia.org/wiki/QRS_complex).

This repository contains two versions of the Pan-Tomkins QRS detection algorithm implementation:
* __QRSDetectorOnline__ - Online version detects QRS complexes in a real-time acquired ECG signal. Therefore, it requires an ECG device to be plugged in and receiving a signal in real-time.
* __QRSDetectorOffline__ - Offline version detects QRS complexes in a pre-recorded ECG signal dataset (e.g. stored in _.csv_ format).

__This implementation of a QRS Complex Detector is by no means a certified medical tool and should not be used in health monitoring. It was created and used for experimental purposes in psychophysiology and psychology.__

## Algorithm

The published QRS Detector module is an implementation of the QRS detection algorithm known as the Pan-Tomkins QRS detection algorithm, first described in a paper by Jiapu Pan and Willis J. Tomkins titled _"A Real-Time QRS Detection Algorithm"_ (1985). The full version of the paper is accessible [here](http://www.robots.ox.ac.uk/~gari/teaching/cdt/A3/readings/ECG/Pan+Tompkins.pdf).

The direct input to the algorithm is a raw ECG signal. The detector processes the signal in two stages: filtering and thresholding.

First, in the filtering stage each raw ECG measurement is filtered using a cascade of low-pass and high-pass filters that together form a band-pass filter. This filtering mechanism ensures that only parts of the signal related to heart activity can pass through. The filters eliminate most of the measurement noise that could cause false positive detection. The band-pass filtered signal is then differentiated to identify signal segments with high signal change values. These changes are then squared and integrated to make them more distinct. In the last step of the processing stage, the integrated signal is screened by a peak detection algorithm to identify potential QRS complexes within the integrated signal.

In the next stage, the identified QRS complex candidates are classified by means of dynamically set thresholds, either as QRS complexes or as noise peaks. The thresholds are real-time adjusted: a threshold in a given moment is based on the signal value of the previously detected QRS and noise peaks. The dynamic thresholding accounts for changes in the noise level. The dynamic thresholding and complex filtering ensure sufficient detection sensitivity with relatively few false positive QRS complex detections.

Importantly, not all of the features presented in the original Pan and Tomkins paper were implemented in this module. Specifically, we decided not to implement supplementary mechanisms proposed in the paper that are not core elements of QRS detection. Therefore, we did not implement the following features: fiducial mark on filtered data, use of another set of thresholds based on the filtered ECG, irregular heart rate detection, and missing QRS complex detection search-back mechanism. Despite the lack of these supplementary features, implementation of the core features proposed by Pan and Tompkins allowed us to achieve a sufficient level of QRS detection. 

## Dependencies

Modules published here consist of the following dependencies:
* jupyter
* matplotlib
* numpy
* pyserial
* scipy

All the dependencies are in the _requirements.py_ file. 

The modules are implemented for use with Python 3.x. However, they are relatively easy to convert to work with Python 2.x: 
- import division module with 
```
from __future__ import division
```
- import print function with
```
from __future__ import print_function
```
- remove the _decode()_ function call when reading data from an ECG device in the online version of the QRS Detector, i.e. use 
```
raw_measurement.rstrip().split(',') 
```
instead of 
```
raw_measurement.decode().rstrip().split(',')
```

## Repository directory structure 
```
├── LICENSE
│
├── README.md          				 <- The top-level README for developers using this project.
│
├── arduino_ecg_sketch 				 <- E-health ECG device Arduino sketch source code and library.
│
├── ecg_data           				 <- Pre-recorded ECG datasets in .csv format.
│
├── logs               				 <- Data logged by Online and Offline QRS Detector in .csv format.
│
├── plots          	   			 <- Plots generated by Offline QRS Detector.
│
├── qrs_detector_offline_example.ipynb  	 <- Jupyter notebook with Offline QRS Detector usage.
│
├── QRSDetectorOffline.py   			 <- Offline QRS Detector module.
│
├── QRSDetectorOnline.py    			 <- Online QRS Detector module.
│
└── requirements.txt  	 			 <- The requirements file containing module dependencies.
```

## Installation and usage

The QRS Detector module was implemented in two separate versions: Online and Offline. Each has a different application and method of use.

### Online QRS Detector

The Online version is designed to work with a directly connected ECG device. As an input it uses an ECG signal received in real-time, detects QRS complexes, and outputs them to be used by other scripts in order to trigger external events. For example, the Online QRS Detector can be used to trigger visual, auditory, or tactile stimuli. It has already been successfully implemented in PsychoPy (Peirce, J. W. (2009). Generating stimuli for neuroscience using PsychoPy (Peirce, J. (2009). _Generating stimuli for neuroscience using PsychoPy._ Frontiers in Neuroinformatics, 2 (January), 1–8. [http://doi.org/10.3389/neuro.11.010.2008](http://doi.org/10.3389/neuro.11.010.2008)) and tested to study cardioceptive (interoceptive) abilities, namely, in Schandry’s heartbeat tracking task (Schandry, R. (1981). _Heart beat perception and emotional experience._ Psychophysiology, 18(4), 483–488. 
[http://doi.org/10.1111/j.1469-8986.1981.tb02486.x](http://doi.org/10.1111/j.1469-8986.1981.tb02486.x)) and heartbeat detection training based on that proposed by Schandry and Weitkunat (Schandry, R., & Weitkunat, R. (1990). _Enhancement of heartbeat-related brain potentials through cardiac awareness training._ The International Journal of Neuroscience, 53(2-4), 243–53. [http://dx.doi.org/10.3109/00207459008986611](http://dx.doi.org/10.3109/00207459008986611)).

The Online version of the QRS Detector module has been implemented to work with the Arduino-based e-Health Sensor Platform V2.0 ECG device. You will find more about this device [here](https://www.cooking-hacks.com/documentation/tutorials/ehealth-biometric-sensor-platform-arduino-raspberry-pi-medical#step4_2).  

An Arduino e-Health ECG device sketch is also provided in this repository. The sampling rate of the ECG signal acquisition is set to 250 (Hz) samples per second. Measurements are sent in a real-time in a string format, _"timestamp,measurement"_, and have to be parsed by the QRS Detector module.

To use the Online QRS Detector module, the ECG device with the loaded Arduino sketch has to be connected to a USB port. Then QRS Complex Detector object is initialized with the port name where the ECG device is connected and the measurement baud rate is set. No further calibration or configuration is needed. The Online QRS Detector starts detection immediately after initialization.

Below is example code of how to run the Online QRS Detector:

```
from QRSDetectorOnline import QRSDetectorOnline

qrs_detector = QRSDetectorOnline(port="/dev/cu.usbmodem14311", baud_rate="115200")
```

If you want to use the Online QRS Detector in the background with other processes running on the layer above it (e.g. display a visual stimulus or play a tone triggered by the detected QRS complexes), we suggest using a Python multiprocessing mechanism. Multiprocessing offers both local and remote concurrency, effectively side-stepping the Global Interpreter Lock by using sub-processes instead of threads. Here is example code of the many ways to achieve this:

```
from QRSDetectorOnline import QRSDetectorOnline

qrs_detector_process = subprocess.Popen(["python", "QRSDetectorOnline.py", "/dev/cu.usbmodem14311"], shell=False)
```

Even though the Online QRS Detector was implemented to be used with an Arduino-based e-Health Sensor Platform ECG device with 250 Hz sampling rate and a specified data format, it can be easily modified to be used with any other ECG device, sampling rate or data format. For more information, check the "Customization" section. 


### Offline QRS Detector

The Offline version of the detector works with ECG measurement datasets stored in _.csv_ format. These can be loaded into the detector to perform QRS detection.

The Offline QRS Detector loads the data, analyses it, and detects QRS complexes in the same way as the online version, but it uses an entire existing dataset instead of real-time measurements directly from an ECG device.

This module is intended for offline QRS detection in ECG measurements when real-time QRS-based event triggering is not crucial, but when more complex data analysis, visualisation and detection information are needed. It can also be used simply as a debugging tool to verify whether the online version works as intended, or to check the behaviour of QRS Detector intermediate data processing stages. 

The offline version of the QRS Detector module was implemented to work with any kind of ECG measurements data that can be loaded from a file. Unlike the online version, it was not designed to work with some specific device. The Offline QRS Detector expects _"timestamp,measurement"_ data format stored in _.csv_ format and is tuned for measurement data acquired with 250 Hz sampling rate. Both the format of the data and the sampling rate can be easily modified. For more information, check the "Customization" section. 

The Offline QRS Detector requires initialization with a path to the ECG measurements file. The QRS Detector will load the dataset, analyse measurements, and detect QRS complexes. It outputs a detection log file with marked detected QRS complexes. In the file, the detected QRS complexes are marked with a '1' flag in the 'qrs_detected' log data column. Additionally, the Offline QRS Detector stores detection results internally as an ecg_data_detected attribute of an Offline QRS Detector object. Optionally, it produces plots with all intermediate signal-processing steps and saves it to a *.csv* file.

Below is example code showing how to run the offline version of the QRS Detector module:

```
from QRSDetectorOnline import QRSDetectorOnline

qrs_detector = QRSDetectorOffline(ecg_data_path="ecg_data/ecg_data_1.csv", verbose=True, log_data=True, plot_data=True, show_plot=False)
```

Check _qrs_detector_offline_example.ipynb_ Jupyter notebook for an example usage of the Offline QRS Detector with generated plots and logs.

## Customization

### QRSDetectorOnline 
The Online QRS Detector module can be easily modified to work with other ECG devices, different sampling rates, and different data formats:

- QRSDetectorOnline works on a real-time received signal from an ECG device. Data is sent to the module in a _"timestamp,measurement"_ format string. If a different ECG data format is expected, the *process_measurement()* function requires some changes in order to enable correct data parsing. The Online QRS Detector works even if only measurement values without corresponding timestamps are available.

- QRSDetectorOnline is tuned for 250 Hz sampling rate by default; however, this can be customized by changing seven configuration attributes (marked in the code) according to the desired signal_frequency.  For example, to change the signal sampling rate from 250 to 125 samples per second, simply divide all the parameters by 2: set *signal_frequency* to 125, *number_of_samples_stored * to 100 samples, *integration_window* to 8 samples, *findpeaks_spacing* to 25 samples, *detection_window* to 20 samples and * refractory_period* to 60 samples.

### QRSDetectorOffline
The Offline QRS Detector is hardware independent because it uses an ECG measurements dataset instead of real-time acquired ECG signal. Therefore, the only parameters that need to be adjusted are sampling rate and data format:

- QRSDetectorOffline is by default tuned for a sampling rate of 250 samples per second. It can be customized by changing 4 configuration attributes (marked in the code) according to the desired *signal_frequency*. For example, to change signal sampling rate from 250 to 125 samples per second, divide all parameters by 2: set *signal_frequency* value to 125, *integration_window* to 8 samples, *findpeaks_spacing* to 25 samples and *refractory_period* to 60 samples.  

- QRSDetectorOffline uses a loaded ECG measurements dataset. Data is expected to be in csv format with each line in *"timestamp,measurement"* format. If a different ECG data format is expected, changes needs to be made in the *load_ecg_data()* function, which loads the dataset, or in the *detect_peaks()* function, which processes measurement values in:
```
ecg_measurements = ecg_data[:, 1] 
```
The algorithm will work fine even if only measurement values without timestamps are available.

## Authors
* Michał Sznajder (Jagiellonian University) - technical contact (msznajder@gmail.com)
* Marta Łukowska (Jagiellonian University)
 

## Citation information
If you use these modules in a research project, please consider citing it:

[![DOI](https://zenodo.org/badge/55516257.svg)](https://zenodo.org/badge/latestdoi/55516257)

If you use these modules in any other project, please refer to MIT open-source license.

## Acknowledgements
The following modules and repository were created as a part of the project “Relationship between interoceptive awareness and metacognitive abilities in near-threshold visual perception” supported by the National Science Centre Poland, PRELUDIUM 7, grant no. 2014/13/N/HS6/02963. Special thanks for Michael Timberlake for proofreading. 

