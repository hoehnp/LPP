#!/usr/bin/python 

# Script:  bond_distribute.py 
# Purpose: binned bond length distributions by bond type
# Syntax:  bond_distribute.py datafile nbin rmin rmax outfile files ...
#          datafile = lammps data file
#          nbin = # of bins per bond type
#          rmin = min expected bond length
#          rmax = max expected bond length
#          outfile = file to write stats to
#          files = series of dump files
# Example: bond_distribute.py pore.data 1000 1.5 1.85 bonds.out pore.dump.1
# Author:  Paul Crozier (Sandia)

# enable script to run from Python directly w/out Pizza.py

import sys
from dump import dump
from math import sqrt
if "argv" not in globals(): argv = sys.argv

# main script

if len(argv) < 7:
  raise Exception("Syntax: bond_distribute.py datafile nbin rmin rmax outfile files ...")

dt = data(argv[1])	
nbins = int(argv[2])
rmin = float(argv[3])
rmax = float(argv[4])
outfile = argv[5]
files = ' '.join(argv[6:])

# get the bonds from the data file

bond = dt.get("Bonds")
nbonds = len(bond)
btype = nbonds * [0]
iatom = nbonds * [0]
jatom = nbonds * [0]
for i in range(nbonds):
  btype[i] = int(bond[i][1] - 1)
  iatom[i] = int(bond[i][2] - 1)
  jatom[i] = int(bond[i][3] - 1)

ntypes = 0
for i in range(nbonds): ntypes = max(bond[i][1],ntypes)
ntypes = int(ntypes)
ncount = ntypes * [0]
bin = nbins * [0]
for i in range(nbins): 
  bin[i] = ntypes * [0] 

# read snapshots one-at-a-time

d = dump(files,0)
d.map(1,"id",2,"type",3,"x",4,"y",5,"z")

while 1:
  time = next(d)
  if time == -1: break
   
  box = (d.snaps[-1].xlo,d.snaps[-1].ylo,d.snaps[-1].zlo,
         d.snaps[-1].xhi,d.snaps[-1].yhi,d.snaps[-1].zhi)
         
  xprd = box[3] - box[0]
  yprd = box[4] - box[1] 
  zprd = box[5] - box[2]
  
  d.unscale() 
  d.sort()
  x,y,z = d.vecs(time,"x","y","z")
   
  for i in range(nbonds):
  
    delx = x[jatom[i]] - x[iatom[i]]
    dely = y[jatom[i]] - y[iatom[i]]
    delz = z[jatom[i]] - z[iatom[i]]
        
    if abs(delx) > 0.5*xprd:
      if delx < 0.0:
        delx += xprd
      else:
        delx -= xprd
    if abs(dely) > 0.5*yprd:
      if dely < 0.0:
        dely += yprd
      else:
        dely -= yprd
    if abs(delz) > 0.5*zprd:
      if delz < 0.0:
        delz += zprd
      else:
        delz -= zprd
        
    r = sqrt(delx*delx + dely*dely + delz*delz)
    
    ibin = int(nbins*(r - rmin)/(rmax - rmin) + 0.5)
    if ((ibin >= 0) and (ibin <= nbins-1)): 
      bin[ibin][btype[i]] += nbins
      ncount[btype[i]] += 1
    else:
      print("Warning: bond distance outside specified range")
      print("Bond type:", btype[i]+1)
      print("Bond number:", i)
  print(time, end=' ')    
      
print()
print("Printing bond distance normalized distribution to",outfile)
    
fp = open(outfile,"w")
rrange = rmax - rmin
for i in range(nbins):
  print(rmin + rrange*float(i)/float(nbins), end=' ', file=fp) 
  for j in range(ntypes):
    if (ncount[j] > 0):
      print(float(bin[i][j])/float(ncount[j])/rrange, end=' ', file=fp)
    else:
      print(0.0, end=' ', file=fp)    
  print(file=fp) 
fp.close()
