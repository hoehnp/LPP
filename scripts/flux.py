#!/usr/bin/python 

# Script:  flux.py
# Purpose: flux of atoms through a user-defined plane
# Syntax:  flux.py x/y/z plane outfile files ...
#          x/y/z = measure flux in x, y, or z direction
#          plane = plane where flux is measured, fraction of L from 0 to 1
#          outfile = file to write flux stats to
#          files = series of dump files
# Example: flux.py z 0.5 flux.out dump.*
# Author:  Paul Crozier (Sandia)

# enable script to run from Python directly w/out Pizza.py

import sys
from dump import dump
if "argv" not in globals(): argv = sys.argv

# main script

if len(argv) < 5:
  raise Exception("Syntax: flux.py x/y/z plane outfile files ...")

direction = argv[1]
scaled_plane = float(argv[2])
outfile = argv[3]
files = ' '.join(argv[4:])
d = dump(files)
d.unwrap()
tmp,ntypes = d.minmax("type")
ntypes = int(ntypes)
ntypes += 1

# loop over snapshots and compute net flux vs. first snapshot

f = open(outfile,"w")

jconfig = crossings = 0
flag = 0

while 1:
  which,time,flag = d.iterator(flag)
  if flag == -1: break
  
  if direction == "x":
    id,type,x = d.vecs(time,"id","type","x")
    lo = d.snaps[which].xlo
    hi = d.snaps[which].xhi
  elif direction == "y":
    id,type,x = d.vecs(time,"id","type","y")
    lo = d.snaps[which].ylo       
    hi = d.snaps[which].yhi       
  elif direction == "z":   
    id,type,x = d.vecs(time,"id","type","z")
    lo = d.snaps[which].zlo       
    hi = d.snaps[which].zhi       
  
  prd = hi - lo                
  plane = lo + scaled_plane*prd          
  
  print(time, end=' ')
  sys.stdout.flush()
    
  natoms = len(x)
  if jconfig == 0: x_initial = (natoms+1) * [0]
  jconfig += 1

  typeflux = ntypes * [0]
      
  for i in range(natoms):
    id[i] = int(id[i])
    type[i] = int(type[i])
    if jconfig == 1: x_initial[id[i]] = x[i]
    if x_initial[id[i]] < plane and x[i] > plane :
       crossings = int((x[i] - plane)/prd) + 1 
       typeflux[type[i]] += crossings
    elif x_initial[id[i]] > plane and x[i] < plane :
       crossings = int((plane - x[i])/prd) + 1
       typeflux[type[i]] -= crossings

  print(time, end=' ', file=f)
  for j in range(ntypes-1):
    print(typeflux[j+1], end=' ', file=f)
  print(file=f)
print()
  
f.close()
