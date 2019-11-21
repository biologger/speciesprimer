#!/usr/bin/python3
# -*- coding:utf-8 -*-

import sys
import argparse
from mydaemon import Daemon

DEBUG = 1


def get_args():

    try:
        parser = argparse.ArgumentParser()
        parser.add_argument('action', help='action',
                            choices=['start', 'stop', 'restart'])
        parser.add_argument('uniqid', help='Unique ID')
        parser.add_argument('dclass', help='daemonclass.py class',
                            choices=[
                                'SpeciesPrimer', 'Update_prokDB',
                                'Update_ntDB', 'GetBlastDB'])
        args = parser.parse_args()

        result = (args.action, args.uniqid, args.dclass)
    except Exception as err:
        if DEBUG:
            raise

        sys.stderr.write('%s\n' % (err))

        result = (None, None, None)

    return result


def main_daemon():
    DAEMON_NAME = 'SpeciesPrimer background (id: #ID#)'
    DAEMON_STOP_TIMEOUT = 10
    PIDFILE = '/programs/tmp/pipeline_#ID#.pid'
    RUNFILE = '/programs/tmp/pipeline_#ID#.run'
    DEBUG = 1
    try:
        (action, uniqid, dclass) = get_args()

        import importlib
        daemonclass = importlib.import_module('daemonclass')
        Runner = getattr(daemonclass, dclass)

        DAEMON_NAME = DAEMON_NAME.replace('#ID#', uniqid)
        PIDFILE = PIDFILE.replace('#ID#', uniqid)
        RUNFILE = RUNFILE.replace('#ID#', uniqid)
        d = Runner(name=DAEMON_NAME, pidfile=PIDFILE, runfile=RUNFILE,
                   stoptimeout=DAEMON_STOP_TIMEOUT, debug=DEBUG)

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

        sys.stderr.write('%s\n' % err)

    sys.exit(1)


if __name__ == '__main__':
    main_daemon()
