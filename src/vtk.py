# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under 
# the GNU General Public License.

# vtk tool

oneline = "Convert LAMMPS snapshots to VTK format"

docstr = """
v = vtk(d)	        d = object containing atom coords (dump, data)

v.one()                 write all snapshots to tmp.vtk
v.one("new")            write all snapshots to new.vtk
v.many()                write snapshots to tmp0000.vtk, tmp0001.vtk, etc
v.many("new")           write snapshots to new0000.vtk, new0001.vtk, etc
v.single(N)             write snapshot for timestep N to tmp.vtk 
v.single(N,"file")      write snapshot for timestep N to file.vtk

  surfaces in snapshot will be written to SURF1.vtk, SURF2.vtk, etc
    where each surface (triangle type) is in a different file
"""

# History
#   8/05, Steve Plimpton (SNL): original version

# ToDo list

# Variables
#   data = data file to read from

# Imports and external programs

import sys, re

# Class definition

class vtk:
  
  # --------------------------------------------------------------------

  def __init__(self,data):
    self.data = data
  
  # --------------------------------------------------------------------

  def one(self,*args):
    if len(args) == 0: file = "tmp.vtk"
    elif args[0][-4:] == ".vtk": file = args[0]
    else: file = args[0] + ".vtk"

    n = flag = 0
    which,time,flag = self.data.iterator(flag)
    time,box,atoms,bonds,tris,lines = self.data.viz(which)
    print(time, end=' ')
    sys.stdout.flush()

    if len(tris): surface(tris)
    
    allatoms = []
    for atom in atoms:
      allatoms.append(atom)
    
    while 1:
      which,time,flag = self.data.iterator(flag)
      if flag == -1: break
      time,box,atoms,bonds,tris,lines = self.data.viz(which)
      
      for atom in atoms: allatoms.append(atom)
      print(time, end=' ')
      sys.stdout.flush()
      n += 1

    particle(file,allatoms)
    print("\nwrote %d snapshots to %s in VTK format" % (n,file))

  # --------------------------------------------------------------------

  def many(self,*args):
    if len(args) == 0: root = "tmp"
    else: root = args[0]

    surfflag = 0
    n = flag = 0
    while 1:
      which,time,flag = self.data.iterator(flag)
      if flag == -1: break
      time,box,atoms,bonds,tris,lines = self.data.viz(which)

      if surfflag == 0 and len(tris):
        surfflag = 1
        surface(tris)
        
      if n < 10:
        file = root + "000" + str(n) + ".vtk"
      elif n < 100:
        file = root + "00" + str(n) + ".vtk"
      elif n < 1000:
        file = root + "0" + str(n) + ".vtk"
      else:
        file = root + str(n) + ".vtk"

      particle(file,atoms)
      
      print(time, end=' ')
      sys.stdout.flush()
      n += 1
      
    print("\nwrote %s snapshots in VTK format" % n)
  # --------------------------------------------------------------------
  def manyGran(self,*args,**kwargs):
    
    # check whether to output or not
    outputfl = True
    if "output" in kwargs: outputfl = kwargs["output"]
    
    # read startIndex (offset for filename due to parallel processing)
    startIndex = 0
    fileNos = []
    if "fileNos" in kwargs:
      fileNos = kwargs["fileNos"]
    else:
      fileNos = range(len(self.data.snaps))
    
    # output name
    if len(args) == 0: root = "tmp"
    else: root = args[0]

    surfflag = 0
    n = flag = 0
    
    # iterate over snaps
    while 1:
      which,time,flag = self.data.iterator(flag)
      if flag == -1: break
      time,box,atoms,bonds,tris,lines = self.data.viz(which)
      xlo=self.data.snaps[n].xlo
      xhi=self.data.snaps[n].xhi
      ylo=self.data.snaps[n].ylo
      yhi=self.data.snaps[n].yhi
      zlo=self.data.snaps[n].zlo
      zhi=self.data.snaps[n].zhi
      
      atoms=self.data.snaps[n].atoms
      names=self.data.names

      if surfflag == 0 and len(tris):
        surfflag = 1
        surface(tris)
        
      file, file_bb, file_walls = generateFilename(root,fileNos,n)

      boundingBox(file_bb,xlo,xhi,ylo,yhi,zlo,zhi)
      nvalues = 0
      try: nvalues = len(self.data.snaps[0].atoms[0])
      except: nvalues = 0
      particleGran(file,atoms,names,nvalues)
      
      if outputfl: print(time, end=' ')
      if outputfl: sys.stdout.flush()
      n += 1
      
    if outputfl: print("\nwrote %s granular snapshots in VTK format" % n)
  # --------------------------------------------------------------------

  def single(self,time,*args):
    if len(args) == 0: file = "tmp.vtk"
    elif args[0][-4:] == ".vtk": file = args[0]
    else: file = args[0] + ".vtk"

    which = self.data.findtime(time)
    time,box,atoms,bonds,tris,lines = self.data.viz(which)
    if len(tris): surface(tris)
    particle(file,atoms)

# ----------------------------------------------------------------------------
# generates the filename of the output-vtk-files from
# - a root string,
# - the a string of numbers (timestamps) and
# - the index, i.e. which of those numbers is going to be used
# ----------------------------------------------------------------------------
def generateFilename(root,fileNos,n):
  if fileNos[n] < 10:
    file = root + "000" + str(fileNos[n]) + ".vtk"
    file_bb= root + "000" + str(fileNos[n]) + "_boundingBox.vtk"
    file_walls= root + "000" + str(fileNos[n]) + "_walls.vtk"
  elif fileNos[n] < 100:
    file = root + "00" + str(fileNos[n]) + ".vtk"
    file_bb= root + "00" + str(fileNos[n]) + "_boundingBox.vtk"
    file_walls= root + "00" + str(fileNos[n]) + "_walls.vtk"
  elif fileNos[n] < 1000:
    file = root + "0" + str(fileNos[n]) + ".vtk"
    file_bb= root + "0" + str(fileNos[n]) + "_boundingBox.vtk"
    file_walls= root + "0" + str(fileNos[n]) + "_walls.vtk"
  else:
    file = root + str(fileNos[n]) + ".vtk"
    file_bb= root + str(fileNos[n]) + "_boundingBox.vtk"
    file_walls= root + str(fileNos[n]) + "_walls.vtk"

  return (file, file_bb, file_walls)

# --------------------------------------------------------------------
# write list of triangles into VTK surface files: SURF1.vtk, SURF2.vtk, ...
# all triangles of one type constitute 1 surface = 1 file
# create list of unique vertices (via dictionary) from triangle list

def surface(tris):
  ntypes = tris[-1][1]
  
  for i in range(ntypes):
    itype = i+1
    v = {}
    nvert = ntri = 0
    for tri in tris:
      if tri[1] == itype:
        ntri += 1
        vert = (tri[2],tri[3],tri[4])
        if vert not in v:
          v[vert] = nvert
          nvert += 1
        vert = (tri[5],tri[6],tri[7])
        if vert not in v:
          v[vert] = nvert
          nvert += 1
        vert = (tri[8],tri[9],tri[10])
        if vert not in v:
          v[vert] = nvert
          nvert += 1

    keys = list(v.keys())
    vinverse = {}
    for key in keys:
      vinverse[v[key]] = key
      
    filename = "SURF" + str(itype) + ".vtk"
    f = open(filename,"w")
  
    print("# vtk DataFile Version 3.0", file=f)
    print("Generated by pizza.py", file=f)
    print("ASCII", file=f)
    print("DATASET POLYDATA", file=f)
    print("POINTS %d float" % nvert, file=f)
    for i in range(nvert):
      tuple = vinverse[i]
      print(tuple[0],tuple[1],tuple[2], file=f)
    print("POLYGONS",ntri,4*ntri, file=f)
    for tri in tris:
      if tri[1] == itype:
        vert = (tri[2],tri[3],tri[4])
        ivert1 = v[vert]
        vert = (tri[5],tri[6],tri[7])
        ivert2 = v[vert]
        vert = (tri[8],tri[9],tri[10])
        ivert3 = v[vert]
        print(3,ivert1,ivert2,ivert3, file=f)
    print(file=f)
    print("CELL_DATA",ntri, file=f)
    print("POINT_DATA",nvert, file=f)

    f.close()

# --------------------------------------------------------------------
# write atoms from one snapshot into file in VTK format

def particle(file,atoms):
  f = open(file,"w")
  
  print("# vtk DataFile Version 2.0", file=f)
  print("Generated by pizza.py", file=f)
  print("ASCII", file=f)
  print("DATASET POLYDATA", file=f)
  print("POINTS %d float" % len(atoms), file=f)
  for atom in atoms:
    print(atom[2],atom[3],atom[4], file=f)
  print("VERTICES",len(atoms),2*len(atoms), file=f)
  for i in range(len(atoms)):
    print(1,i, file=f)
  print("POINT_DATA",len(atoms), file=f)
  print("SCALARS atom_type int 1", file=f)
  print("LOOKUP_TABLE default", file=f)
  for atom in atoms:
    itype = int(atom[1])
    print(itype, end=' ', file=f)
  print(file=f)
  
  f.close()


def boundingBox(file,xlo,xhi,ylo,yhi,zlo,zhi):
  f = open(file,"w")
  
  print("# vtk DataFile Version 2.0", file=f)
  print("Generated by pizza.py", file=f)
  print("ASCII", file=f)
  print("DATASET RECTILINEAR_GRID", file=f)
  print("DIMENSIONS 2 2 2", file=f)
  print("X_COORDINATES 2 float", file=f)
  print(xlo,xhi, file=f)
  print("Y_COORDINATES 2 float", file=f)
  print(ylo,yhi, file=f)
  print("Z_COORDINATES 2 float", file=f)
  print(zlo,zhi, file=f)

def typestr(o):
  string = str(type(o))
  sp = string.split('\'')
  return sp[1]

def particleGran(file,atoms,names,n_values):
  f = open(file,"w")
  
  # if no atoms are present
  if atoms is None:
    atoms = []
  
  # find indices of scalars and vectors
  scalars, vectors = findScalarsAndVectors(names)
  
  # print head
  print("# vtk DataFile Version 2.0", file=f)
  print("Generated by lpp.py", file=f)
  print("ASCII", file=f)
  print("DATASET POLYDATA", file=f)
  print("POINTS %d float" % len(atoms), file=f)
  for atom in atoms:
    print(atom[vectors['x']], atom[vectors['x']+1], atom[vectors['x']+2] , file=f) #atom[3],atom[4],atom[5]  #write x,y,z  [atom[0]=id, atom[1]=type]
  print("VERTICES", len(atoms), 2*len(atoms), file=f)
  for i in range(len(atoms)):
    print(1,i, file=f)
  print("POINT_DATA",len(atoms), file=f)
  
  if len(atoms) == 0:
    print('', file=f)
    f.close()
    return

  # print VECTORS
  for key in vectors.keys():
    
    # don't print coodinates again
    if key == 'x':
      continue
    
    vectortype = 'float'
    if atoms != []:
      vectortype = typestr(atoms[0][vectors[key]])
      if 'float' in vectortype: vectortype = 'float'
      elif 'int' in vectortype: vectortype = 'int'
      else: vectortype = 'float'
    else: # if no atoms are present
      pass
    
    print("VECTORS",key,vectortype, file=f)
    for atom in atoms:
      print(atom[vectors[key]], atom[vectors[key]+1], atom[vectors[key]+2], file=f)
  
  # print SCALARS
  for key in scalars.keys():
    scalartype =''
    if atoms != []:
      scalartype = typestr(atoms[0][scalars[key]])
      if 'float' in scalartype: scalartype = 'float'
      elif 'int' in scalartype: scalartype = 'int'
      else: scalartype = 'int'
    else: # if no atoms are present
      pass
    
    print("SCALARS",key,scalartype,1, file=f)
    print("LOOKUP_TABLE default", file=f)
    for atom in atoms:
      print(atom[scalars[key]], file=f)
    
  print('', file=f)
  f.close()

def findScalarsAndVectors(names):
  
  vectors={}
  scalars={}
  
  # create reversed dictionary {position:name}
  indices = {}
  for name in names:
    indices[names[name]]=name
  
  # fill missing indices (occurrs e.g. if output is like vx vy vz fx fy fz vx vy vz)
  for i in range(max(indices)):
    if i not in indices:
      indices[i]=""
  
  # compile regexes to find vectors
  regvx = re.compile(".*x")
  regvy = re.compile(".*y")
  regvz = re.compile(".*z")
  regf = re.compile("f_.*\[[0-9]+\]")
  regc = re.compile("c_.*\[[0-9]+\]")
  regv = re.compile("v_.*\[[0-9]+\]")
  
  # loop over all indices and look if their names represent a vector (if not: it's a scalar)
  i = 0
  while i<= max(indices):
    
    if i+2 <= max(indices) and regvx.match(indices[i]) != None and regvy.match(indices[i+1]) != None and regvz.match(indices[i+2]) != None:
      
      newname=''
      if len(indices[i]) == 1:
        newname=indices[i]
      else:
        newname=indices[i][:-1]
      
      vectors[newname]=i
      i+=3
      continue
    
    if regf.match(indices[i]) != None or regc.match(indices[i]) != None or regv.match(indices[i]) != None:
      
      name = indices[i]
      number = int( name.split('[')[1].split(']')[0] )
      nextName = name.split('[')[0]+'['+str(number+1)+']'
      nextButOneName = name.split('[')[0]+'['+str(number+2)+']'
      
      newname = name[2:-(len(name.split('[')[1])+1)]
      
      if i+2 <= max(indices) and indices[i+1] == nextName and indices[i+2] == nextButOneName:
        vectors[newname]=i
        i+=3
        continue
      
      else:
        scalars[newname]=i
        i+=1
        continue
    
    # program only here if not a vector
    if indices[i] != '':
      newname = indices[i]
      scalars[newname]=i
    i+=1
  
  if 'x' not in vectors.keys():
    print("vector x y z has to be contained in dump file. please change liggghts input script accordingly.")
    exit()
  
  return scalars, vectors
      
     
  
