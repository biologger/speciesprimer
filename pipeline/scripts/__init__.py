#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys

script_dir = os.path.dirname(os.path.abspath(__file__))
pipe_dir, tail = os.path.split(script_dir)
sys.path.insert(0, pipe_dir)
