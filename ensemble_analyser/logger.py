import logging
import os

LOG_FORMAT = "%(message)s"

DEBUG = os.getenv('DEBUG')

logging.basicConfig(
    filename='output.out',
    level=logging.DEBUG if DEBUG else logging.INFO,
    format=LOG_FORMAT,
    filemode='w'
    )
log = logging.getLogger()


ordinal = lambda n: "%d-%s" % (n,"tsnrhtdd"[(n//10%10!=1)*(n%10<4)*n%10::4])

def save_snapshot(output, confs):

    log.debug('Saving snapshot of the ensemble')
    with open(output, 'w') as f:
        for i in confs:
            f.write(f'{i}\n')

    return None