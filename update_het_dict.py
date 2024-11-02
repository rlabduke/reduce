from __future__ import division, print_function
import os, sys

#This script updates the het dictionary used by Reduce to add hydrogens to het groups
#It obsoletes and replaces update_het_dict.sh
#Outline:
#Download the lasted version of the het_dictionary from the PDB
#Reformat certain columns of the het dictionary
#Edit certain functional groups
#Regularize endlines for Windows
#Set filenames and clean up

#----- Download het dict from PDB -----
import urllib
import time

print("Downloading most recent het dict", file=sys.stderr)
#Download the lasted version of the het_dictionary from the PDB
urllib.urlretrieve('https://files.wwpdb.org/pub/pdb/data/monomers/het_dictionary.txt', 'reduce_wwPDB_het_dict_download.txt')
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

print("Adjusting atom columns", file=sys.stderr)
for line in infile:
  if line.startswith("CONECT"):
    x = line.split()
    numconnections = '%5s' % x[2]
    connections = x[3:]
    atomlist = []
    for atom in connections:
      atomlist.append(make_pdb3_format_atom(atom))
    outlist = [x[0],"     ",make_pdb3_format_atom(x[1]),numconnections,' '.join(atomlist).rstrip()]
    print(''.join(outlist).strip(), file=outfile)
  elif line.startswith("RESIDUE"):
    x = line.split()
    resname = x[1].rjust(3)
    members = x[2].strip().rjust(7)
    print(x[0]+"   "+resname+ members, file=outfile)
  elif line.startswith("END"):
    print(line.strip(),"  ", file=outfile)
  elif line.startswith("HETNAM") or line.startswith("HETSYN") or line.startswith("FORMUL"):
    outfile.write(line)
  else:
    print(line.strip(), file=outfile)

infile.close()
outfile.close()
#----- End adjustment -----

#----- Remove H from Phosphate OH -----
#The default het dict contains some hydrogens that are not present at
#  physiological conditions. These hydrogens are usually present to make the
#  ligands charge-neutral.
#The most prevalent of these is the H on phosphate OH.
#Here we comment out those hydrogens so Reduce will not add them.
#We also adjust SPD, SPM, and PAR nitrogens to have appropriate charge states
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

  def fix_spd_charge(self):
    line_index = 0
    while line_index < len(self.lines):
      line = self.lines[line_index]
      if line.startswith("CONECT      N1 "):
        self.lines[line_index] = "CONECT      N1     4 C2  HN11 HN12 HN13\n"
      elif line.startswith("CONECT      N6 "):
        self.lines[line_index] = "CONECT      N6     4 C5   C7  HN61 HN62\n"
      elif line.startswith("CONECT      N10 "):
        self.lines[line_index] = "CONECT      N10    4 C9  H101 H102 H103\n"
      elif line.startswith("CONECT     HN12"):
        self.lines[line_index] = "CONECT     HN12    1 N1\nCONECT     NH13    1 N1\n"
      elif line.startswith("CONECT      HN6 "):
        self.lines[line_index] = "CONECT     HN61    1 N6\nCONECT     HN62    1 N6\n"
      elif line.startswith("CONECT     H102    1 N10"):
        self.lines[line_index] = "CONECT     H102    1 N10\nCONECT     H103    1 N10\n"
      line_index += 1

  def fix_spm_charge(self):
    line_index = 0
    while line_index < len(self.lines):
      line = self.lines[line_index]
      if line.startswith("CONECT      N1 "):
        self.lines[line_index] = "CONECT      N1     4 C2  HN11 HN12 HN13\n"
      elif line.startswith("CONECT      N5 "):
        self.lines[line_index] = "CONECT      N5     4 C4   C6  HN51 HN52\n"
      elif line.startswith("CONECT      N10"):
        #HN0 -> H101
        self.lines[line_index] = "CONECT      N10    4 C9   C11 H101 H102\n"
      elif line.startswith("CONECT      N14"):
        self.lines[line_index] = "CONECT      N14    4 C13 HN41 HN42 HN43\n"

      elif line.startswith("CONECT     HN12"):
        self.lines[line_index] = "CONECT     HN12    1 N1\nCONECT     NH13    1 N1\n"
      elif line.startswith("CONECT      HN5"):
        self.lines[line_index] = "CONECT     HN51    1 N5\nCONECT     HN52    1 N5\n"
      elif line.startswith("CONECT      HN0"):
        #HN0 -> H101
        self.lines[line_index] = "CONECT     H101    1 N10\nCONECT     H102    1 N10\n"
      elif line.startswith("CONECT     HN42"):
        self.lines[line_index] = "CONECT     HN42    1 N14\nCONECT     HN43    1 N14\n"
      line_index += 1

  def fix_par_charge(self):
    line_index = 0
    while line_index < len(self.lines):
      line = self.lines[line_index]
      if line.startswith("CONECT      N21"):
        self.lines[line_index] = "CONECT      N21    4 C21 HN21 HN22 HN23\n"
      elif line.startswith("CONECT      N12"):
        self.lines[line_index] = "CONECT      N12    4 C12 H121 H122 H123\n"
      elif line.startswith("CONECT      N32"):
        self.lines[line_index] = "CONECT      N32    4 C32 H321 H322 H323\n"
      elif line.startswith("CONECT      N24"):
        self.lines[line_index] = "CONECT      N24    4 C24 H241 H242 H243\n"
      elif line.startswith("CONECT      N64"):
        self.lines[line_index] = "CONECT      N64    4 C64 HN61 HN62 HN63\n"

      elif line.startswith("CONECT     HN22"):
        self.lines[line_index] = "CONECT     HN22    1 N21\nCONECT     HN23    1 N21\n"
      elif line.startswith("CONECT     H122"):
        self.lines[line_index] = "CONECT     H122    1 N12\nCONECT     H123    1 N12\n"
      elif line.startswith("CONECT     H322"):
        self.lines[line_index] = "CONECT     H322    1 N32\nCONECT     H323    1 N32\n"
      elif line.startswith("CONECT     H242"):
        self.lines[line_index] = "CONECT     H242    1 N24\nCONECT     H243    1 N24\n"
      elif line.startswith("CONECT     HN62"):
        self.lines[line_index] = "CONECT     HN62    1 N64\nCONECT     HN63    1 N64\n"
      line_index += 1

  def process(self, formul_line, outfile):
    phosphate_pattern = re.compile('(FORMUL)(.........)(.*H[1-9])(.*O[1-9])(.*P[1-9])')
    if phosphate_pattern.match(formul_line):
      self.remove_OH_on_P()
    if self.lines[0].startswith("RESIDUE   SPD"):
      self.fix_spd_charge()
    elif self.lines[0].startswith("RESIDUE   SPM"):
      self.fix_spm_charge()
    elif self.lines[0].startswith("RESIDUE   PAR"):
      self.fix_par_charge()
    self.output_record(outfile)

infile = open("reduce_wwPDB_het_dict_adjusted.txt")
outfile = open("reduce_wwPDB_het_dict_no_poh.txt","w")

print("Adjusting charges", file=sys.stderr)
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

print("Switching to Windows endlines", file=sys.stderr)
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
print("File cleanup", file=sys.stderr)
existing_files = os.listdir(".")
if "reduce_wwPDB_het_dict.txt" in existing_files:
  print("  The previous het dict was renamed to reduce_wwPDB_het_dict_obsoleted_"+fetchtime+".txt", file=sys.stderr)
  os.rename("reduce_wwPDB_het_dict.txt","reduce_wwPDB_het_dict_obsoleted_"+fetchtime+".txt")
os.rename("reduce_wwPDB_het_dict_endlines.txt","reduce_wwPDB_het_dict.txt")
os.remove("reduce_wwPDB_het_dict_download.txt")
os.remove("reduce_wwPDB_het_dict_adjusted.txt")
os.remove("reduce_wwPDB_het_dict_no_poh.txt")
#----- End rename and remove -----
