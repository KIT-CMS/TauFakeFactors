
import sys


class Logger(object):
    def __init__(self, log_file_name="logfile.log"):
        self.terminal = sys.stdout
        self.log = open(log_file_name, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        pass