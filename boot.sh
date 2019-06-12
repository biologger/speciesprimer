#!/bin/sh
cd /home/pipeline/bin
exec gunicorn -b :5000 --access-logfile - --error-logfile - speciesprimergui:app

