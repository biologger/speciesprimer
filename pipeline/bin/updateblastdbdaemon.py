#!/usr/bin/python3
# -*- coding:utf-8 -*-

import os
import sys
import time
import atexit
from mydaemon import Daemon
import argparse
import subprocess
import logging

DEBUG = 0


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
        today = time.strftime("%Y_%m_%d", time.localtime())
        logging.basicConfig(
            filename=os.path.join(
                "/", "home", "primerdesign",
                "speciesprimer_" + today + ".log"),
            level=logging.DEBUG, format="%(message)s")
        blastdb_dir = "/home/blastdb"
        os.chdir(blastdb_dir)
        command = [
            "update_blastdb.pl", "--passive", "--decompress",
            "--blastdb_version", "5", "nt_v5"]
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        def check_output():
            while True:
                output = process.stdout.readline().decode().strip()
                if output:
                    print(output)
                    logging.info(
                        time.strftime("%d %b %Y %H:%M:%S: ", time.localtime())
                        + str(output).strip())
                else:
                    break
        while process.poll() is None:
            check_output()

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
    DEBUG = 0
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
