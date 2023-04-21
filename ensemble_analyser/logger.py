

import logging
import os, sys

LOG_FORMAT = "%(message)s"

DEBUG = os.getenv('DEBUG')

ordinal = lambda n: "%d-%s" % (n,"tsnrhtdd"[(n//10%10!=1)*(n%10<4)*n%10::4])

def save_snapshot(output, confs, log):

    log.debug('Saving snapshot of the ensemble')
    with open(output, 'w') as f:
        for i in confs:
            f.write(f'{i.write_xyz()}\n')

    return None


class StreamToLogger():
    """
    Fake file-like stream object that redirects writes to a logger instance.
    """
    def __init__(self, logger, level):
       self.logger = logger
       self.level = level

    def write(self, buf):
        self.logger.log(self.level, ''.join([i for i in buf.splitlines() if i]))

    def flush(self):
        pass


def create_log(output):

    logging.basicConfig(
        filename=output,
        level=logging.DEBUG if DEBUG else logging.INFO,
        format=LOG_FORMAT,
        filemode='w'
    )

    log = logging.getLogger()
    sys.stdout = StreamToLogger(log,logging.INFO)
    sys.stderr = StreamToLogger(log,logging.ERROR)

    return logging.getLogger()



if __name__ == '__main__':

    log = create_log('tests/output_test.out')
    print('Test to standard out\n\n')
    raise Exception('Test to standard error')
