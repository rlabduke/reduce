from __future__ import division
import os, sys

#This script updates the het dictionary used by Reduce to add hydrogens to het groups
#It obsoletes and replaces update_het_dict.sh
#Outline:
#Download the lasted version of the het_dictionary from the PDB
#Reformat certain columns of the het dictionary
#Edit certain functional groups
#Regularize endlines for Windows
#Set filenames and clen up

#----- Download het dict from PDB -----
import urllib
import time

print >> sys.stderr, "Downloading most recent het dict"
#Download the lasted version of the het_dictionary from the PDB
urllib.urlretrieve('ftp://ftp.wwpdb.org/pub/pdb/data/monomers/het_dictionary.txt', 'reduce_wwPDB_het_dict_download.txt')
#Store time of download for versioning
fetchtime = time.strftime("%Y%m%d_%H%M%S")
#----- End download -----

#----- Adjust columns in het dict -----
def make_pdb3_format_atom(atomstring):
  #put atoms into a PDB3-like format
  #Note that this reformatting does not distinguish between calcium 'CA  ' and
  #  alpha carbons ' CA '.  Both will be rendered as ' CA '.
  #No issues are known to result, but the possibility remains.
  atomstring = atomstring.strip()
  l = len(atomstring)
  if   l == 4: return atomstring
  elif l == 3: return ' '+atomstring
  elif l == 2: return ' '+atomstring+' '
  elif l == 1: return ' '+atomstring+'  '
  else:
    sys.stdout.write("ERROR: Strange atom name |"+atomstring+"|")
    sys.exit()

infile = open("reduce_wwPDB_het_dict_download.txt")
outfile = open("reduce_wwPDB_het_dict_adjusted.txt","w")

print >> sys.stderr, "Adjusting atom columns"
for line in infile:
  if line.startswith("CONECT"):
    x = line.split()
    numconnections = '%5s' % x[2]
    connections = x[3:]
    atomlist = []
    for atom in connections:
      atomlist.append(make_pdb3_format_atom(atom))
    outlist = [x[0],"     ",make_pdb3_format_atom(x[1]),numconnections,' '.join(atomlist).rstrip()]
    print >> outfile, ''.join(outlist).strip()
  elif line.startswith("RESIDUE"):
    x = line.split()
    resname = x[1].rjust(3)
    members = x[2].strip().rjust(7)
    print >> outfile, x[0]+"   "+resname+ members
  elif line.startswith("END"):
    print >> outfile, line.strip(),"  "
  elif line.startswith("HETNAM") or line.startswith("HETSYN") or line.startswith("FORMUL"):
    outfile.write(line)
  else:
    print >> outfile, line.strip()

infile.close()
outfile.close()
#----- End adjustment -----

#----- Remove H from Phosphate OH -----
#The default het dict contains some hydrogens that are not present at
#  physiological conditions. These hydrogens are usually present to make the
#  ligands charge-neutral.
#The most prevalent of these is the H on phosphate OH.
#Here we comment out those hydrogens so Reduce will not add them.
import re

class record():
  def __init__(self):
    self.lines = []
    self.has_O = False
    self.has_P = False

  def addline(self, line):
    self.lines.append(line)

  def output_record(self,outfile):
    for line in self.lines:
      outfile.write(line)
    outfile.write("\n")

  def remove_OH_on_P(self):
    line_index = 0
    hydrogens_to_remove = []
    while line_index < len(self.lines):
      line = self.lines[line_index]
      if line.startswith("CONECT"):
        parent_atom = line[11:15].strip()
        number_atoms = line[16:20].strip()
        connected_atoms = line[20:].split()
        if parent_atom in hydrogens_to_remove:
          self.lines[line_index] = "xx"+line #assumes hydrogen always comes after oxygen
        elif parent_atom.startswith("O") and number_atoms=="2": #need 2 connected atoms
          atom_1 = connected_atoms[0].strip()
          atom_2 = connected_atoms[1].strip()
          if atom_1.startswith("P") and atom_2.startswith("H"):
            hydrogens_to_remove.append(atom_2)
            #matches pattern for a phosphate oxygen, need to comment out the attached hydrogen
      line_index += 1

  def process(self, formul_line, outfile):
    phosphate_pattern = re.compile('(FORMUL)(.........)(.*H[1-9])(.*O[1-9])(.*P[1-9])')
    if phosphate_pattern.match(formul_line):
      self.remove_OH_on_P()
    self.output_record(outfile)

infile = open("reduce_wwPDB_het_dict_adjusted.txt")
outfile = open("reduce_wwPDB_het_dict_no_poh.txt","w")

print >> sys.stderr, "Removing phosphate hydorgens"
residue_line_pattern = re.compile('(RESIDUE)(...)(?P<resname>...)')
formula_line_pattern = re.compile('FORMUL')
for line in infile:
  if residue_line_pattern.match(line):
    #resname = residue_line_pattern.match(line).group('resname')
    rec = record()
    rec.addline(line)
  elif formula_line_pattern.match(line):
    rec.addline(line)
    rec.process(line, outfile)
  else:
    rec.addline(line)

infile.close()
outfile.close()
#----- End remove H from POH -----

#----- Regularized endlines for Windows -----
#Swap newline characters for Windows-compatible ones to ensure usability on Win
infile = open("reduce_wwPDB_het_dict_no_poh.txt")
outfile = open("reduce_wwPDB_het_dict_endlines.txt","w")

print >> sys.stderr, "Switching to Windows endlines"
lines = infile.readlines()
for line in lines:
  line=line.replace("\n","")
  outfile.write("%s\r\n" % line)

infile.close()
outfile.close()
#----- End regularize endlines -----

#----- Rename and remove files -----
#The default path Reduce uses to find the het dict file can be checked by
#  invoking reduce -h and looking under the -DB flag.
#Generally, Reduce is looking for a file named reduce_wwPDB_het_dict.txt.
#Renaming the updated file is renamed to match that name allows it to be used
#  automatically in most circumstances.
#The previous file is renamed with a datestamp.
print >> sys.stderr, "File cleanup"
existing_files = os.listdir(".")
if "reduce_wwPDB_het_dict.txt" in existing_files:
  print >> sys.stderr, "  The previous het dict was renamed to reduce_wwPDB_het_dict_obsoleted_"+fetchtime+".txt"
  os.rename("reduce_wwPDB_het_dict.txt","reduce_wwPDB_het_dict_obsoleted_"+fetchtime+".txt")
os.rename("reduce_wwPDB_het_dict_endlines.txt","reduce_wwPDB_het_dict.txt")
os.remove("reduce_wwPDB_het_dict_download.txt")
os.remove("reduce_wwPDB_het_dict_adjusted.txt")
os.remove("reduce_wwPDB_het_dict_no_poh.txt")
#----- End rename and remove -----
