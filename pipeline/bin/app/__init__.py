#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from flask import Flask
from guiconfig import Config

app = Flask(__name__)
app.config.from_object(Config)

from app import routes


