#!/usr/bin/python

# Script:  dview.py
# Purpose: launch vcr tool on LAMMPS dump files
# Syntax:  dview.py dump.1 dump.2 ...
#          files = one or more dump files
# Example: dview.py dump.*
# Author:  Steve Plimpton (Sandia)

# enable script to run from Python directly w/out Pizza.py

import sys
from dump import dump

# w/out Pizza.py these lines need to come before import of gl tool
import tkinter
tkroot = tkinter.Tk()
tkroot.withdraw()

from gl import gl
from vcr import vcr
if "argv" not in globals(): argv = sys.argv

# main script

if len(argv) < 2:
  raise Exception("Syntax: dview.py dump.1 ...")

files = ' '.join(argv[1:])

d = dump(files)
g = gl(d)
v = vcr(g)
