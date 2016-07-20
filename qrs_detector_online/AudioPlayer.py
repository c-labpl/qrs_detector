# -*- coding: utf-8 -*-

from pygame import mixer


class AudioPlayer(object):

    def __init__(self, file_path):
        mixer.init()
        self.mix = mixer.Sound(file_path)

    def play(self):
        self.mix.play()

if __name__ == "__main__":
    player = AudioPlayer(file_path="audio/beep.wav")
    player.play()
