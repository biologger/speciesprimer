#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

BASE_PATH = os.path.dirname(os.path.abspath(__file__))
pipe_dir = os.path.join(BASE_PATH.split("tests")[0], "pipeline")
sys.path.append(pipe_dir)
