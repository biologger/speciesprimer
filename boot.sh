#!/bin/sh
cd /pipeline/gui
exec gunicorn -b :5000 --access-logfile - --error-logfile - speciesprimergui:app

