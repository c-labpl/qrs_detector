# import sys, os, random, cfg
# cwd = os.getcwd()
# os.chdir("../base")
# sys.path += ["../base"]
# from media import get_counter
# os.chdir(cwd)
import serial
import numpy as np
from scipy.signal import butter, lfilter, find_peaks_cwt
from pygame import mixer
#from matplotlib import pyplot as plt
from multiprocessing import Process, freeze_support
from time import gmtime, strftime

from libs.findpeaks import findpeaks # janko_slavic
from libs import detect_peaks # detect_peaks tony berta

class QRSDetectorMultiprocessing(Process):

	'''
	Initialization methods.
	'''
	def __init__(self, name, port):
         
         super(QRSDetectorMultiprocessing, self).__init__()

         self.name = name
         self.port = port

         # Connection.
         self.arduino = None
                  
         # Receiving data.
         self.updateData = False
         self.timestamp = None
         self.measurement = None   

         # Signal processing phases.
         self.rawSignal = []
         self.bandPassSignal = []
         self.differentiatedSignal = []
         self.squaredSignal = []
         self.integratedSignal = []
         
         # Detection.
         self.SPKI = 0.4
         self.NPKI = 0.1
         self.PEAKI = 0.0
         self.THRESHOLDI1 = 0.06
         self.THRESHOLDI2 = 0
         self.QRSPeak = []
         self.NOISEPeak = []
         self.RRInterval = 37 # previously 17
         self.RRCurrent = 0
         self.QRSInterval = 12
         self.recentBPM = []

         # Data logger set up.
         self.timestamp = 0.0
         self.measurement = 0.0
         self.detectedPulse = 0
         self.interbeatInterval = 0.0

         # Audio set up.
         self.playSound = 0
         self.s = None

         # Files set up.
         self.f = None
         
	def run(self):

         # Audio initialization.
         mixer.init()
         self.s = mixer.Sound('beep_ok.wav')
         self.audio_control_file = open("audio.txt")

         # Log file initialization.
         filename = "data/qrs_%s_%s_datalog.txt" % (strftime("%Y-%m-%d_%H-%M-%S", gmtime()), self.name)
         self.f = open(filename, 'w')
         self.f.truncate()

         self.f.write("%s %s %s %s %s \n" % ("sysTime", "timestamp", "ecg", "beat", "ibi"))

         # Arduino init.
         self.connectToArduino()
         self.startUpdatingData()

	'''
	Arduino connection methods.
	'''
	def connectToArduino(self):
		print "Initializing serial port."
#		 /dev/cu.usbmodem1411
#            COM5
		self.arduino = serial.Serial(self.port, 115200)
		self.updateData = True
		print "Ready."

	def disconnectFromArduino(self):
		self.updateData = False
  
   	'''
	Tools methods.
	'''
	# Data update method.
	def cycleList(self, N, collection, updateValue):
	    collection.append(updateValue)    
	    if len(collection) > N:
	        collection.pop(0)

	# Bandpass filter methods.
	def butterBandpass(self, lowcut, highcut, fs, order=5):
	    nyq = 0.5 * fs
	    low = lowcut / nyq
	    high = highcut / nyq
	    b, a = butter(order, [low, high], btype='band')
	    return b, a

	def butterBandpassFilter(self, data, lowcut, highcut, fs, order=5):
	    b, a = self.butterBandpass(lowcut, highcut, fs, order=order)
	    y = lfilter(b, a, data)
	    return y

	'''
	Loading data.
	'''
     	def startUpdatingData(self):
           while self.updateData:
            with open("audio.txt") as aud:
             line = aud.readline()
             if not(line == '0' or line == '1'):
                line = '0'
             self.playSound = int(line) 
         
            update = self.arduino.readline().rstrip().split(';')
         
            if len(update) < 2:
                continue
            try:
                self.timestamp = float(update[0])
                self.measurement = float(update[1])
            except Exception:
                continue

            self.cycleList(100, self.rawSignal, self.measurement)
            self.RRCurrent = self.RRCurrent + 1

            if (len(self.rawSignal) == 1):
		        continue;

            ## Bandpass filter - pass band 5-15 Hz.
            fs = 84.0
            lowcut = 5.0
            highcut = 15.0
            self.bandPassSignal = self.butterBandpassFilter(self.rawSignal, lowcut, highcut, fs, order=1)

            ## Derivative - provide QRS slope info. Five point derivative. Delay 1 sample.
            self.differentiatedSignal = np.diff(self.bandPassSignal)
            self.differentiatedSignal = self.differentiatedSignal / max(self.differentiatedSignal)

            ## Squaring - signal is squared point by point. Non-linear apflification of derivative.
            self.squaredSignal = self.differentiatedSignal * self.differentiatedSignal
            self.squaredSignal = self.squaredSignal / max(self.squaredSignal)

            # Moving-window integration - to obtain waveform feature information in addition to the slope.
            # Time window should be of the length of longest possible QRS complex - 30 samples (150 ms) for 200 samples/s.
            N = 12
            self.integratedSignal = np.convolve(self.squaredSignal, np.ones((N,)) / N)

            '''
            QRS detection phase.
            '''
            # Fiducial mark - peak detection.
            fiducialMark = []
            peaksIndices = find_peaks_cwt(self.integratedSignal[:-1], np.arange(10,15), noise_perc=0.1)
        
            # 2
            # N = 30
            # peaksIndices = findpeaks(self.integratedSignal, limit=0.06, spacing=1)

            # 3
            # N = 30
            # peaksIndices = detect_peaks.detect_peaks(self.integratedSignal, mph=0.06, mpd=50)

            for peakIndex in peaksIndices:
		        fiducialMark.append((peakIndex, self.integratedSignal[peakIndex]))

            # Thresholding detected peaks.
            if (self.RRCurrent > self.RRInterval):
		        if len(fiducialMark) > 0 and fiducialMark[-1][0]  > len(self.integratedSignal) - self.QRSInterval:
		            if fiducialMark[-1][1] > self.THRESHOLDI1:
				     self.detectedPulse = 1
				     self.QRSPeak.append(fiducialMark[-1][0])
				     self.SPKI = 0.125 * fiducialMark[-1][1] + 0.875 * self.SPKI
				     print "PULSE", fiducialMark[-1][1]
				     if self.playSound == 1:
						self.s.play()
				     currentBPM = 60.0 / ((1.0 / fs) * self.RRCurrent)
				     self.interbeatInterval = ((1.0 / fs) * self.RRCurrent)
				     self.cycleList(5, self.recentBPM, currentBPM)
				     print "BPM", sum(self.recentBPM) / len(self.recentBPM)
				     self.RRCurrent = 0
		            else:
		                self.NOISEPeak.append(fiducialMark[-1][0])
		                self.NPKI = 0.125 * fiducialMark[-1][1] + 0.875 * self.NPKI
		                
		            self.THRESHOLDI1 = self.NPKI + 0.25 * (self.SPKI - self.NPKI)
		            self.THRESHOLDI2 = 0.5 * self.THRESHOLDI1

            # self.f.write("%d %f %f %d %f \n" % (get_counter(), self.timestamp, self.measurement, self.detectedPulse, self.interbeatInterval))
            self.f.write("%d %f %f %d %f \n" % (0, self.timestamp, self.measurement, self.detectedPulse, self.interbeatInterval))
            self.detectedPulse = 0
            self.interbeatInterval = 0.0
            self.playSound = 0
              
if __name__ == "__main__":
    freeze_support() 
    QRSDetector = QRSDetectorMultiprocessing("test", "/dev/cu.usbmodem1411")
    QRSDetector.playSound = 1
    QRSDetector.daemon = True
    QRSDetector.start()
    QRSDetector.join() 
