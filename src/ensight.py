# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under 
# the GNU General Public License.

# ensight tool

oneline = "Convert LAMMPS snapshots or meshes to Ensight format"

docstr = """
e = ensight(d)	     d = object with atoms or elements (dump,data,mdump)
e.change = 1         set to 1 if element nodal xyz change with time (def = 0)
e.maxtype = 10       max particle type, set if query to data will be bad

e.one()
e.one("new")
e.one("cns","Centro","eng","Energy")
e.one("new","cns","Centro","eng","Energy")
                     write all snapshots as an Ensight data set
                     Ensight header file = tmp.case (no 1st arg) or new.case
                     Ensight coord file = tmp.xyz or new.xyz
                     additional pairs of args create auxiliary files:
                       tmp.cns, tmp.eng or new.cns, new.eng
                     cns,eng = column name in dump file and file name suffix
                     Centro,Energy = Ensight name for the variable

e.increment()        same args as one(), but process dump out-of-core

e.many()             same args as one(), but create multiple Ensight files
                     tmp0000.xyz, tmp0001.xyz, etc
                     new0000.cns, new0001.cns, etc
                     new0000.eng, new0001.eng, etc

e.single(N)          same args as one() prepended by N, but write a single snap
"""

# History
#   10/06, Steve Plimpton (SNL): original version

# ToDo list
#   binary files
#   create vector or tensor variable files, not just scalar
#     via pair of args like ["vx","vy","vz"],"vel"

# Variables
#   data = data file to read from
#   which = 0 for particles, 1 for elements
#   change = 0 for unchanging mesh coords, 1 for changing mesh coords (def = 0)

# Imports and external programs

import sys, types

# Class definition

class ensight:

  # --------------------------------------------------------------------

  def __init__(self,data):
    self.change = 0
    self.maxtype = 0
    self.data = data
    if type(data) is types.InstanceType and ".dump" in str(data.__class__):
      self.which = 0
    elif type(data) is types.InstanceType and ".data" in str(data.__class__):
      self.which = 0
    elif type(data) is types.InstanceType and ".mdump" in str(data.__class__):
      self.which = 1
    elif type(data) is types.InstanceType and ".cdata" in str(data.__class__):
      self.which = 1
    else:
      raise Exception("unrecognized object passed to ensight")
    
  # --------------------------------------------------------------------

  def one(self,*args):
    if len(args) % 2 == 0: root = "tmp"
    else: 
      root = args[0]
      args = args[1:]

    pairs = []
    for i in range(0,len(args),2): pairs.append([args[i],args[i+1]])

    # max # of types for all steps in Ensight files

    if self.which == 0 and self.maxtype == 0:
      self.maxtype = self.data.maxtype()
    
    # write Ensight *.case header file

    f = open("%s.case" % root,"w")
    times = self.data.time()
    self.case_file(f,root,pairs,0,len(times),times)
    f.close()

    # open additional files

    f = open(root + ".xyz","w")
    vfiles = []
    for pair in pairs: vfiles.append(open(root + "." + pair[0],"w"))

    # loop over snapshots
    # write coords into xyz file, variables into their files

    first = 1
    n = flag = etype = 0
    while 1:
      which,time,flag = self.data.iterator(flag)
      if flag == -1: break

      if self.which == 0:
        print("BEGIN TIME STEP", file=f)
        time,box,atoms,bonds,tris,lines = self.data.viz(which)
        self.coord_file_atoms(f,box,atoms)
        print("END TIME STEP", file=f)
      elif self.change == 0 and first:
        print("BEGIN TIME STEP", file=f)
        time,box,nodes,elements,nvalues,evalues = self.data.mviz(which)
        self.coord_file_elements(f,box,nodes,elements)
        etype = len(elements[0])
        first = 0
        print("END TIME STEP", file=f)
      elif self.change:
        print("BEGIN TIME STEP", file=f)
        time,box,nodes,elements,nvalues,evalues = self.data.mviz(which)
        self.coord_file_elements(f,box,nodes,elements)
        etype = len(elements[0])
        print("END TIME STEP", file=f)

      for i in range(len(pairs)):
        print("BEGIN TIME STEP", file=vfiles[i])
        values = self.data.vecs(time,pairs[i][0])
        if self.which == 0:
          self.variable_file_atoms(vfiles[i],pairs[i][1],atoms,values)
        else:
          self.variable_file_elements(vfiles[i],pairs[i][1],etype,values)
        print("END TIME STEP", file=vfiles[i])
        
      print(time, end=' ')
      sys.stdout.flush()
      n += 1

    # close additional files

    f.close()
    for f in vfiles: f.close()

    print("\nwrote %s snapshots in Ensight format" % n)

  # --------------------------------------------------------------------

  def increment(self,*args):
    if len(args) % 2 == 0: root = "tmp"
    else: 
      root = args[0]
      args = args[1:]

    pairs = []
    for i in range(0,len(args),2): pairs.append([args[i],args[i+1]])

    # max # of types for all steps in Ensight files

    if self.which == 0 and self.maxtype == 0:
      self.maxtype = self.data.maxtype()

    # open additional files

    f = open(root + ".xyz","w")
    vfiles = []
    for pair in pairs: vfiles.append(open(root + "." + pair[0],"w"))

    # loop over snapshots
    # write coords into xyz file, variables into their files

    times = []
    first = 1
    n = etype = 0
    while 1:
      time = next(self.data)
      if time == -1: break
      times.append(time)
      self.data.tselect.one(time)
      self.data.delete()

      if self.which == 0:
        print("BEGIN TIME STEP", file=f)
        time,box,atoms,bonds,tris,lines = self.data.viz(0)
        self.coord_file_atoms(f,box,atoms)
        print("END TIME STEP", file=f)
      elif self.change == 0 and first:
        print("BEGIN TIME STEP", file=f)
        time,box,nodes,elements,nvalues,evalues = self.data.mviz(0)
        self.coord_file_elements(f,box,nodes,elements)
        etype = len(elements[0])
        first = 0
        print("END TIME STEP", file=f)
      elif self.change:
        print("BEGIN TIME STEP", file=f)
        time,box,nodes,elements,nvalues,evalues = self.data.mviz(0)
        self.coord_file_elements(f,box,nodes,elements)
        etype = len(elements[0])
        print("END TIME STEP", file=f)

      for i in range(len(pairs)):
        print("BEGIN TIME STEP", file=vfiles[i])
        values = self.data.vecs(time,pairs[i][0])
        if self.which == 0:
          self.variable_file_atoms(vfiles[i],pairs[i][1],atoms,values)
        else:
          self.variable_file_elements(vfiles[i],pairs[i][1],etype,values)
        print("END TIME STEP", file=vfiles[i])
        
      print(time, end=' ')
      sys.stdout.flush()
      n += 1

    # close additional files

    f.close()
    for f in vfiles: f.close()

    # write Ensight *.case header file now that know all timesteps

    f = open("%s.case" % root,"w")
    self.case_file(f,root,pairs,0,len(times),times)
    f.close()

    print("\nwrote %s snapshots in Ensight format" % n)

  # --------------------------------------------------------------------

  def many(self,*args):
    if len(args) % 2 == 0: root = "tmp"
    else: 
      root = args[0]
      args = args[1:]

    pairs = []
    for i in range(0,len(args),2): pairs.append([args[i],args[i+1]])

    # max # of types for all steps in Ensight files

    if self.which == 0 and self.maxtype == 0:
      self.maxtype = self.data.maxtype()

    # write Ensight *.case header file

    f = open("%s.case" % root,"w")
    times = self.data.time()
    self.case_file(f,root,pairs,1,len(times),times)
    f.close()

    # loop over snapshots
    # generate unique filenames
    # write coords into one xyz file per snapshot, variables into their files

    first = 1
    n = flag = etype = 0
    while 1:
      which,time,flag = self.data.iterator(flag)
      if flag == -1: break

      files = []
      if n < 10:
        file = root + "000" + str(n) + ".xyz"
	for pair in pairs:
	  files.append(root + "000" + str(n) + "." + pair[0])
      elif n < 100:
        file = root + "00" + str(n) + ".xyz"
	for pair in pairs:
	  files.append(root + "00" + str(n) + "." + pair[0])
      elif n < 1000:
        file = root + "0" + str(n) + ".xyz"
	for pair in pairs:
	  files.append(root + "0" + str(n) + "." + pair[0])
      else:
        file = root + str(n) + ".xyz"
	for pair in pairs:
	  files.append(root + str(n) + "." + pair[0])

      if self.which == 0:
        f = open(file,"w")
        time,box,atoms,bonds,tris,lines = self.data.viz(which)
        self.coord_file_atoms(f,box,atoms)
        f.close()
      elif self.change == 0 and first:
        f = open(root + ".xyz","w")
        time,box,nodes,elements,nvalues,evalues = self.data.mviz(which)
        self.coord_file_elements(f,box,nodes,elements)
        etype = len(elements[0])
        first = 0
        f.close()
      elif self.change:
        f = open(file,"w")
        time,box,nodes,elements,nvalues,evalues = self.data.mviz(which)
        self.coord_file_elements(f,box,nodes,elements)
        etype = len(elements[0])
        f.close()

      for i in range(len(pairs)):
        values = self.data.vecs(time,pairs[i][0])
        f = open(files[i],"w")
        if self.which == 0:
          self.variable_file_atoms(f,pairs[i][1],atoms,values)
        else:
          self.variable_file_elements(f,pairs[i][1],etype,values)
	f.close()

      print(time, end=' ')
      sys.stdout.flush()
      n += 1
  
    print("\nwrote %s snapshots in Ensight format" % n)

  # --------------------------------------------------------------------

  def single(self,time,*args):
    if len(args) % 2 == 0: root = "tmp"
    else: 
      root = args[0]
      args = args[1:]

    pairs = []
    for i in range(0,len(args),2): pairs.append([args[i],args[i+1]])

    # max # of types for all steps in Ensight files

    if self.which == 0 and self.maxtype == 0:
      self.maxtype = self.data.maxtype()

    # write Ensight *.case header file

    f = open("%s.case" % root,"w")
    self.case_file(f,root,pairs,0,1,[time])
    f.close()

    # write coords into xyz file, variables into their files

    which = self.data.findtime(time)
    etype = 0
    
    f = open(root + ".xyz","w")
    if self.which == 0:
      time,box,atoms,bonds,tris,lines = self.data.viz(which)
      self.coord_file_atoms(f,box,atoms)
    else:
      time,box,nodes,elements,nvalues,evalues = self.data.mviz(which)
      self.coord_file_elements(f,box,nodes,elements)
      etype = len(elements[0])
    f.close()
    
    for i in range(len(pairs)):
      values = self.data.vecs(time,pairs[i][0])
      f = open(root + "." + pairs[i][0],"w")
      if self.which == 0:
        self.variable_file_atoms(f,pairs[i][1],atoms,values)
      else:
        self.variable_file_elements(f,pairs[i][1],etype,values)
      f.close()
      
  # --------------------------------------------------------------------
  # write Ensight case file

  def case_file(self,f,root,pairs,multifile,nsnaps,times):
    print("# Ensight case file\n", file=f)
    print("FORMAT\ntype: ensight gold\n", file=f)
    
    if self.which == 0:
      if multifile:
#        print >>f,"GEOMETRY\nmodel: %s****.xyz change_coords_only\n" % root
        print("GEOMETRY\nmodel: %s****.xyz\n" % root, file=f)
      else:
#        print >>f,"GEOMETRY\nmodel: 1 1 %s.xyz change_coords_only\n" % root
        print("GEOMETRY\nmodel: 1 1 %s.xyz\n" % root, file=f)
    else:
      if self.change == 0:
        print("GEOMETRY\nmodel: %s.xyz\n" % root, file=f)
      elif multifile:
        print("GEOMETRY\nmodel: %s****.xyz\n" % root, file=f)
      else:
        print("GEOMETRY\nmodel: 1 1 %s.xyz\n" % root, file=f)

    if len(pairs):
      print("VARIABLE", file=f)
      for pair in pairs:
        if self.which == 0:
          if multifile:
            print("scalar per node: %s %s****.%s" % (pair[1],root,pair[0]), file=f)
          else:
            print("scalar per node: 1 1 %s %s.%s" % (pair[1],root,pair[0]), file=f)
        else:
          if multifile:
            print("scalar per element: %s %s****.%s" % (pair[1],root,pair[0]), file=f)
          else:
            print("scalar per element: 1 1 %s %s.%s" % (pair[1],root,pair[0]), file=f)
      print(file=f)

    print("TIME", file=f)
    print("time set: 1", file=f)
    print("number of steps:",nsnaps, file=f)
    print("filename start number: 0", file=f)
    print("filename increment: 1", file=f)
    print("time values:", file=f)
    for i in range(nsnaps):
      print(times[i], end=' ', file=f)
      if i % 10 == 9: print(file=f)
    print(file=f)
    print(file=f)
    
    if not multifile:
      print("FILE", file=f)
      print("file set: 1", file=f)
      print("number of steps:",nsnaps, file=f)

  # --------------------------------------------------------------------
  # write Ensight coordinates for atoms
  # partition into "parts"
  # one part = coords for all atoms of a single type

  def coord_file_atoms(self,f,box,atoms):
    print("Particle geometry\nfor a collection of atoms", file=f)
    print("node id given", file=f)
    print("element id off", file=f)
    print("extents", file=f)
    print("%12.5e%12.5e" % (box[0],box[3]), file=f)
    print("%12.5e%12.5e" % (box[1],box[4]), file=f)
    print("%12.5e%12.5e" % (box[2],box[5]), file=f)

    for type in range(1,self.maxtype+1):
      print("part", file=f)
      print("%10d" % type, file=f)
      print("type",type, file=f)
      print("coordinates", file=f)
      group = [atom for atom in atoms if int(atom[1]) == type]
      print("%10d" % len(group), file=f)
      for atom in group: print("%10d" % int(atom[0]), file=f)
      for atom in group: print("%12.5e" % atom[2], file=f)
      for atom in group: print("%12.5e" % atom[3], file=f)
      for atom in group: print("%12.5e" % atom[4], file=f)
      print("point", file=f)
      print("%10d" % len(group), file=f)
      for i in range(1,len(group)+1): print("%10d" % i, file=f)

  # --------------------------------------------------------------------
  # write Ensight coordinates for elements

  def coord_file_elements(self,f,box,nodes,elements):
    print("Element geometry\nfor a collection of elements", file=f)
    print("node id given", file=f)
    print("element id given", file=f)
    print("extents", file=f)
    print("%12.5e%12.5e" % (box[0],box[3]), file=f)
    print("%12.5e%12.5e" % (box[1],box[4]), file=f)
    print("%12.5e%12.5e" % (box[2],box[5]), file=f)

    print("part", file=f)
    print("%10d" % 1, file=f)
    print("all elements", file=f)
    print("coordinates", file=f)
    print("%10d" % len(nodes), file=f)
    for node in nodes: print("%10d" % int(node[0]), file=f)
    for node in nodes: print("%12.5e" % node[2], file=f)
    for node in nodes: print("%12.5e" % node[3], file=f)
    for node in nodes: print("%12.5e" % node[4], file=f)

    if len(elements[0]) == 5: print("tria3", file=f)
    elif len(elements[0]) == 6: print("tetra4", file=f)
    else: raise Exception("unrecognized element type")
    print("%10d" % len(elements), file=f)

    for element in elements: print("%10d" % int(element[0]), file=f)
    if len(elements[0]) == 5:
      for element in elements:
        print("%10d%10d%10d" % \
              (int(element[2]),int(element[3]),int(element[4])), file=f)
    elif len(elements[0]) == 6:
      for element in elements:
        print("%10d%10d%10d%10d" % \
              (int(element[2]),int(element[3]),int(element[4]),int(element[5])), file=f)

  # --------------------------------------------------------------------
  # write Ensight variable values for atoms
  # partition into "parts"
  # one part = values for all atoms of a single type

  def variable_file_atoms(self,f,name,atoms,values):
    print("Particle %s" % name, file=f)
    for type in range(1,self.maxtype+1):
      print("part", file=f)
      print("%10d" % type, file=f)
      print("coordinates", file=f)
      group = [values[i] for i in range(len(atoms))
               if int(atoms[i][1]) == type]
      for value in group: print("%12.5e" % value, file=f)

  # --------------------------------------------------------------------
  # write Ensight variable values for elements

  def variable_file_elements(self,f,name,etype,values):
    print("Element %s" % name, file=f)
    print("part", file=f)
    print("%10d" % 1, file=f)
    if etype == 5: print("tria3", file=f)
    elif etype == 6: print("tetra4", file=f)
    for value in values: print("%12.5e" % value, file=f)
