#!/bin/bash
######################################
#Gary Kapral 04/20/2011
#update the het dictionary in reduce on the cluster; downloads the het dictionary from PDB, date stamps it
#makes it readable by reduce, and links the reduce_wwPDB_het_dict.txt used by reduce to the new het dictionary
#use: sign in as srcer, run from /home/srcer/src/reduce/reduce_trunk/
#other programs used: adjust_het_dict.pl, comment_out_OH_on_P_v3.0.pl
######################################
DATESTAMP=$(date +%Y%m%d)
curl --output reduce_wwPDB_het_dict$DATESTAMP.unmod.txt --stderr curlerr.log ftp://ftp.wwpdb.org/pub/pdb/data/monomers/het_dictionary.txt	#download new het_dict from PDB

perl ./adjust_het_dict.pl reduce_wwPDB_het_dict$DATESTAMP.unmod.txt reduce_het_dict$DATESTAMP.adjusted.txt					#get rid of unnecessary spaces
perl ./comment_out_OH_on_P_v3.0.pl reduce_het_dict$DATESTAMP.adjusted.txt									#comment out the H on phosphate that PDB puts in for neutral charge
mv reduce_het_dict$DATESTAMP.adjusted.txt_mod reduce_wwPDB_het_dict$DATESTAMP.txt								#rename the new file follow the standard name format
#rm /home/srcer/bin/reduce_wwPDB_het_dict.txt 													#remove the old symbolic link
#ln -s /home/srcer/src/reduce/reduce_trunk/reduce_wwPDB_het_dict$DATESTAMP.txt /home/srcer/bin/reduce_wwPDB_het_dict.txt				#add a new symbolic link that points to the latest het_dict

