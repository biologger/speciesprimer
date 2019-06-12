#!/usr/bin/python3
# -*- coding:utf-8 -*-

import os
import sys
import time
import atexit
from mydaemon import Daemon
import argparse

DEBUG = 1


def get_args():
    '''
    >>> get_args()
    ('start', 5)
    >>> get_args()
    ('stop', 4)
    '''

    try:
        parser = argparse.ArgumentParser()
        parser.add_argument('action', help='action',
                            choices=['start', 'stop', 'restart'])
        parser.add_argument('uniqid', help='Unique ID')

        args = parser.parse_args()

        result = (args.action, args.uniqid)
    except Exception as err:
        if DEBUG:
            raise
        else:
            sys.stderr.write('%s\n' % (err))

        result = (None, None)

    return result


class Pipeline(Daemon):

    def start_download(self):
        import getblastdb
        getblastdb.get_DB(mode="auto")

    def run(self):
        RUNFILE = self.runfile
        '''
        Process start to run here.
        '''
        self.start_download()
        # Delete the run file at exit. Maybe there will be no stop request.
        atexit.register(self.delrun)

        # Run while there is no stop request.
        n = 0
        while os.path.exists(RUNFILE):
            print('.', end='')
            n += 1

            if (n > 5):
                break

            time.sleep(1)


def main_daemon():
    DAEMON_NAME = 'BLAST DB background (id: #ID#)'
    DAEMON_STOP_TIMEOUT = 10
    PIDFILE = '/tmp/blastdb_#ID#.pid'
    RUNFILE = '/tmp/blastdb__#ID#.run'
    DEBUG = 1
    try:
        (action, uniqid) = get_args()
        # Create daemon object.
        DAEMON_NAME = DAEMON_NAME.replace('#ID#', uniqid)
        PIDFILE = PIDFILE.replace('#ID#', uniqid)
        RUNFILE = RUNFILE.replace('#ID#', uniqid)
        d = Pipeline(name=DAEMON_NAME, pidfile=PIDFILE, runfile=RUNFILE,
                     stoptimeout=DAEMON_STOP_TIMEOUT, debug=DEBUG)

        # Action requested.
        if action == 'start':
            d.start()
        elif action == 'stop':
            d.stop()
        elif action == 'restart':
            d.restart()
        else:
            raise NameError('Unknown action')

        sys.exit(0)
    except Exception as err:
        if DEBUG:
            raise
        else:
            sys.stderr.write('%s\n' % err)

    sys.exit(1)


# -----------------------------------------------------------------------------
if __name__ == '__main__':
    main_daemon()
