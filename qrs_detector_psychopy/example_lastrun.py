#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy2 Experiment Builder (v1.84.0rc5),
    on Wed Sep 21 23:08:02 2016
If you publish work using this script please cite the PsychoPy publications:
    Peirce, JW (2007) PsychoPy - Psychophysics software in Python.
        Journal of Neuroscience Methods, 162(1-2), 8-13.
    Peirce, JW (2009) Generating stimuli for neuroscience using PsychoPy.
        Frontiers in Neuroinformatics, 2:10. doi: 10.3389/neuro.11.010.2008
"""

from __future__ import absolute_import, division

import psychopy
psychopy.useVersion('latest')

from psychopy import locale_setup, visual, core, data, event, logging, sound, gui
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED,
                                STOPPED, FINISHED, PRESSED, RELEASED, FOREVER)
import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import (sin, cos, tan, log, log10, pi, average,
                   sqrt, std, deg2rad, rad2deg, linspace, asarray)
from numpy.random import random, randint, normal, shuffle
import os  # handy system and path functions
import sys  # to get file system encoding

# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__)).decode(sys.getfilesystemencoding())
os.chdir(_thisDir)

# Store info about the experiment session
expName = 'example'  # from the Builder filename that created this script
expInfo = {u'session': u'001', u'participant': u''}
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName

# Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
filename = _thisDir + os.sep + u'data/%s_%s_%s' % (expInfo['participant'], expName, expInfo['date'])

# An ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(name=expName, version='',
    extraInfo=expInfo, runtimeInfo=None,
    originPath=u'/Users/michalsznajder/Dropbox (Personal)/github/qrsdetector/qrs_detector_psychopy/example.psyexp',
    savePickle=True, saveWideText=True,
    dataFileName=filename)
# save a log file for detail verbose info
logFile = logging.LogFile(filename+'.log', level=logging.EXP)
logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file

endExpNow = False  # flag for 'escape' or other condition => quit the exp

# Start Code - component code to be run before the window creation

# Setup the Window
win = visual.Window(
    size=[1440, 900], fullscr=False, screen=0,
    allowGUI=True, allowStencil=False,
    monitor='testMonitor', color=[0,0,0], colorSpace='rgb',
    blendMode='avg', useFBO=True)
# store frame rate of monitor if we can measure it
expInfo['frameRate'] = win.getActualFrameRate()
if expInfo['frameRate'] != None:
    frameDur = 1.0 / round(expInfo['frameRate'])
else:
    frameDur = 1.0 / 60.0  # could not measure, so guess

# Initialize components for Routine "trial"
trialClock = core.Clock()
ISI = core.StaticPeriod(win=win, screenHz=expInfo['frameRate'], name='ISI')
image = visual.ImageStim(
    win=win, name='image',
    image='resources/bird.jpg', mask=None,
    ori=0, pos=(0, 0), size=(0.5, 0.5),
    color=[1,1,1], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-1.0)
import subprocess
import glob

# Start QRSDetector as subprocess. Change port to correct one.
p = subprocess.Popen(["python", "QRSDetector.py", "/dev/cu.usbmodem1411"], shell=False)

# Initialize audio feedback setup.
# sound.Sound()
#audio = psychopy.sound.SoundPyo(value='C', secs=0.2, volume=5.0)
audio = sound.Sound(value='500', secs=0.2)
 
# Initialize QRSDetector log file.
qrs_log_file = max(glob.iglob('data/*.csv'), key=os.path.getctime)
print "QRSLog:", qrs_log_file
fin = open(qrs_log_file, "r")


# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.CountdownTimer()  # to track time remaining of each (non-slip) routine 

# ------Prepare to start Routine "trial"-------
t = 0
trialClock.reset()  # clock
frameN = -1
continueRoutine = True
routineTimer.add(8.000000)
# update component parameters for each repeat

# keep track of which components have finished
trialComponents = [ISI, image]
for thisComponent in trialComponents:
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED

# -------Start Routine "trial"-------
while continueRoutine and routineTimer.getTime() > 0:
    # get current time
    t = trialClock.getTime()
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *image* updates
    if t >= 0.0 and image.status == NOT_STARTED:
        # keep track of start time/frame for later
        image.tStart = t
        image.frameNStart = frameN  # exact frame index
        image.setAutoDraw(True)
    frameRemains = 0.0 + 8- win.monitorFramePeriod * 0.75  # most of one frame period left
    if image.status == STARTED and t >= frameRemains:
        image.setAutoDraw(False)
    # Log last qrs measurement in every frame.
    lines = fin.readlines()
    fin.seek(0)
    split_last_line = lines[-1].strip().split(" ")
    if len(split_last_line) == 4:
        string = "LOG %s %s %s %s" % (split_last_line[0], split_last_line[1], split_last_line[2], split_last_line[3])  
        logging.log(level=logging.DATA, msg=string)
    elif len(lines) >= 2 and lines[-2].strip().split(" ") == 4:
        split_lastbutone_line = lines[-2].strip().split(" ")
        string = "LOG %s %s %s %s" % (split_lastbutone_line[0], split_lastbutone_line[1], split_lastbutone_line[2], split_lastbutone_line[3])  
        logging.log(level=logging.DATA, msg=string)
    else:
        logging.log(level=logging.DATA, msg="No data")
    
    
    # Check if QRS peak was detected and trigger events if yes.
    peak_detected = os.path.exists("flag.txt")
    
    if peak_detected:
        logging.log(level=logging.DATA, msg="test")
        # Log detection frame.
        with open("flag.txt") as flag:
            line = flag.readline()
            logging.log(level=logging.DATA, msg="QRS " + line)
        
        # Play audio.
        # audio.play()
        audio.play()
        
        # Start listening for keys.
        psychopy.event.getKeys(keyList=["z"], timeStamped=True)
        
        # Do whatever else.
        print "Peak detected!"
        
        # Do not remove this.
        os.remove("flag.txt")
    
    # *ISI* period
    if t >= 0.0 and ISI.status == NOT_STARTED:
        # keep track of start time/frame for later
        ISI.tStart = t
        ISI.frameNStart = frameN  # exact frame index
        ISI.start(0.5)
    elif ISI.status == STARTED:  # one frame should pass before updating params and completing
        ISI.complete()  # finish the static period
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in trialComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # check for quit (the Esc key)
    if endExpNow or event.getKeys(keyList=["escape"]):
        core.quit()
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "trial"-------
for thisComponent in trialComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
peak_detected = os.path.exists("flag.txt")
if peak_detected:
    os.remove("flag.txt")

p.kill()
fin.close()
# these shouldn't be strictly necessary (should auto-save)
thisExp.saveAsWideText(filename+'.csv')
thisExp.saveAsPickle(filename)
logging.flush()
# make sure everything is closed down
thisExp.abort()  # or data files will save again on exit
win.close()
core.quit()
