# -*- coding: utf-8 -*-

from time import gmtime, strftime


class Logger(object):

    def __init__(self, name, separator, *columns):
        self.path = "data/%s_log_%s.csv" % (name, strftime("%Y_%m_%d_%H_%M_%S", gmtime()))
        self.fout = open(self.path, "w")

        self.separator = separator
        self.fout.write(self.separator.join(columns) + "\n")

    def log(self, *args):
        self.fout.write(self.separator.join(args) + "\n")

    def close(self):
        self.fout.close()

if __name__ == "__main__":
    logger = Logger("test", ",", "test1", "tes2")
    logger.log("1", "2")
