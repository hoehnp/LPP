#!/usr/bin/python

# Script:  clogview.py
# Purpose: plots of ChemCell log-file concentration data
# Syntax:  clogview.py gnu/matlab files ...
#          gnu/matlab = style of plots to create
#          files = one or more log files
# Example: clogview.py gnu log.*
# Author:  Steve Plimpton (Sandia)

# enable script to run from Python directly w/out Pizza.py

import sys
from clog import clog
from plotview import plotview
from gnu import gnu
from matlab import matlab
if "argv" not in globals(): argv = sys.argv

# main script

if len(argv) < 3:
  raise Exception("Syntax: clogview.py gnu/matlab files  ...")

style = argv[1]
files = ' '.join(argv[2:])

c = clog(files)
exec("plot = %s()" % style)
p = plotview(c,plot)
p.x = "Time"

