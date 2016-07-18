from pygame import mixer


class AudioPlayer(object):

    def __init__(self, filepath):
        mixer.init()
        self.mix = mixer.Sound(filepath)

    def play(self):
        self.mix.play()

if __name__ == "__main__":
    player = AudioPlayer(filepath="audio/beep.wav")
    player.play()

