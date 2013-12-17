#!/bin/bash

# update_het_dict.sh
# unix shell script which invokes the scripts to update the Reduce program's
# dictionary file.  Run script from a writable directory containing this script
# and the appropriate perl scripts (see below),
#   $>  ./update_het_dict.sh 
# Also required is:
#     * the capability to download the current "het_dictionary.txt".  Some 
#       simple modifications to the commands below allows for use of other than
#       the curl download machinery used here.
#     * perl version 5 (maybe 6).  Gotta have it or some other version of Perl
#       to run the perl scripts: adjust_het_dict.pl and comment_out_OH_on_P_v3.0.pl
#       for the typical run yielding a PDB V3 atom naming dictionary.
# 
# scripts authors: Robert Immormino, Gary Kapral
# date written: April 2011
# update history at bottom of this file

# *nix builtin to generate date used in file names
DATESTAMP=$(date +%Y%m%d)

# invokes curl to download current het_dict from www PDB
curl --output reduce_wwPDB_het_dict$DATESTAMP.unmod.txt --stderr curlerr.log \
     ftp://ftp.wwpdb.org/pub/pdb/data/monomers/het_dictionary.txt

# perl script to reformat dictionary
perl ./adjust_het_dict.pl reduce_wwPDB_het_dict$DATESTAMP.unmod.txt \
     reduce_het_dict$DATESTAMP.adjusted.txt

# perl script to remove phosphate hydroxyl hydrogen added for charge neutralization
perl ./comment_out_OH_on_P_v3.0.pl reduce_het_dict$DATESTAMP.adjusted.txt

# rename (and perhaps move) dictionary
mv reduce_het_dict$DATESTAMP.adjusted.txt_mod reduce_wwPDB_het_dict$DATESTAMP.txt
# Further renaming and movement to default het_dict location for the 
# reduce executable may be needed.  Run reduce with the -h flag for default 
# location and file name.
