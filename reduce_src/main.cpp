// Name: main.cpp
// Author: J. Michael Word
// Date Written: 7/15/97
// Purpose: add hydrogens to a Protein Data Bank file and
//          writing a new Protein Data Bank file to standard output
//
// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999-2016 J. Michael Word
// **************************************************************
//
//  reduceChanges now contains the CHANGELOG or history info
//

#if defined(_MSC_VER)
#pragma warning(disable:4786)
#pragma warning(disable:4305)
#pragma warning(disable:4800)
#endif

#include "reduce.h"
#include "DotSph.h"
#include <iostream>
using std::cout;
using std::cin;
using std::cerr;
using std::endl;
#include <fstream>
#include <sstream>

#ifndef NOSYSSTATS
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif

#define DIRECTORY_SEP_CHAR '/'

#ifndef HET_DICTIONARY
//#define HET_DICTIONARY "reduce_het_dict.txt"
#define HET_DICTIONARY "reduce_wwPDB_het_dict.txt"
#endif
#ifndef HET_DICTOLD
#define HET_DICTOLD "reduce_het_dict.txt"
#endif
std::string DBfilename(HET_DICTIONARY);

void reduceHelp(bool showAll) { /*help*/
   cerr << versionString << endl;
   cerr << shortVersion << endl;
   cerr << "arguments: [-flags] filename or -" << endl;
//   cerr << "040509 reduce.C ln 455: NO renumber, 1459: NO RXR msg" << endl;
//   cerr << "041113 rework main to do first and loop over other NMR models if model# not specified."<< endl;
//   cerr << "Adds hydrogens to a PDB format file and writes to standard output." << endl;
//   cerr << "(note: By default, HIS sidechain NH protons are not added. See -BUILD)" << endl;
   cerr << endl;
   cerr << "Suggested usage:" << endl;
   cerr << "reduce -FLIP myfile.pdb > myfileFH.pdb (do NQH-flips)" << endl;
   cerr << "reduce -NOFLIP myfile.pdb > myfileH.pdb (do NOT do NQH-flips)" << endl << endl;
   cerr << "Flags:" << endl;
   cerr << "-FLIP             add H and rotate and flip NQH groups" << endl;
   cerr << "-NOFLIP           add H and rotate groups with no NQH flips" << endl;
   cerr << "-Trim             remove (rather than add) hydrogens and skip all optimizations" << endl;

  if (showAll) {
   cerr << endl;
   cerr << "-NUClear          use nuclear X-H distances rather than default" << endl;
   cerr << "                  electron cloud distances" << endl;
   cerr << "-NOOH             remove hydrogens on OH and SH groups" << endl;
   cerr << "-OH               add hydrogens on OH and SH groups (default)" << endl;
   cerr << endl;
   cerr << "-HIS              create NH hydrogens on HIS rings" << endl;
//   cerr << "-FLIPs            allow complete ASN, GLN and HIS sidechains to flip" << endl; - removed 130724 JJH
   cerr << "                        (usually used with -HIS)" << endl;
   cerr << "-NOHETh           do not attempt to add NH proton on Het groups" << endl;
//   cerr << "-ADDNHATGAP            add \"amide\" hydrogen on chain breaks" <<endl;
//   cerr << "-GAPERROR#.#       sets the half width for allowed peptide bond lengths variations around 1.4 Angstroms: default 0.3" << endl;
   cerr << "-ROTNH3           allow lysine NH3 to rotate (default)" << endl;
   cerr << "-NOROTNH3         do not allow lysine NH3 to rotate" << endl;
   cerr << "-ROTEXist         allow existing rotatable groups (OH, SH, Met-CH3) to rotate" << endl;
   cerr << "-ROTEXOH          allow existing OH & SH groups to rotate" << endl;
//   cerr << "-ALLMETHYLS       allow all methyl groups to rotate" << endl;
   cerr << "-ALLALT           process adjustments for all conformations (default)" << endl;
   cerr << "-ONLYA            only adjust 'A' conformations" << endl;
   cerr << "-CHARGEs          output charge state for appropriate hydrogen records" << endl;
//   cerr << "-NOROTMET         do not rotate methionine methyl groups" << endl; //now default behavior -cjw 160602
   cerr << "-DOROTMET         allow methionine methyl groups to rotate (not recommended)" << endl;
   cerr << "-NOADJust         do not process any rot or flip adjustments" << endl;
  }
   cerr << endl;
   cerr << "-NOBUILD#.#       build with a given penalty often 200 or 999" << endl;
   cerr << "-BUILD            add H, including His sc NH, then rotate and flip groups" << endl;
   cerr << "                  (except for pre-existing methionine methyl hydrogens)" << endl;
   cerr << endl;
  if (showAll) {
   cerr << "                  (same as: -OH -ROTEXOH -HIS -FLIP)" << endl;
   cerr << "-Keep             keep bond lengths as found" << endl;
   cerr << "-MAXAromdih#      dihedral angle cutoff for aromatic ring planarity check (default="<< MaxAromRingDih <<")" << endl; // - Aram 07/24/12
   cerr << "-NBonds#          remove dots if cause within n bonds (default="<< NBondCutoff <<")" << endl;
   cerr << "-Nterm#           max number of nterm residue (default="<<MinNTermResNo<<")" << endl;
   cerr << "-DENSity#.#       dot density (in dots/A^2) for VDW calculations (Real, default="<<VdwDotDensity<<")" << endl;
   cerr << "-RADius#.#        probe radius (in A) for VDW calculations (Real, default="<<ProbeRadius<<")" << endl;
   cerr << "-OCCcutoff#.#     occupancy cutoff for adjustments (default="<<OccupancyCutoff<<")" << endl;
   cerr << "-H2OOCCcutoff#.#  occupancy cutoff for water atoms (default="<<WaterOCCcutoff<<")" << endl;
   cerr << "-H2OBcutoff#      B-factor  cutoff for water atoms (Integer, default="<<WaterBcutoff<<")" << endl;
   cerr << "-PENalty#.#       fraction of std. bias towards original orientation (default="<<PenaltyMagnitude<<")" << endl;
   cerr << "-HBREGcutoff#.#   over this gap regular HBonds bump (default="<<MinRegHBgap<<")" << endl;
   cerr << "-HBCHargedcut#.#  over this gap charged HBonds bump (default="<<MinChargedHBgap<<")" << endl;
   cerr << "-BADBumpcut#.#    at this gap a bump is 'bad' (default="<<BadBumpGapCut<<")" << endl;
   cerr << "-METALBump#.#     H 'bumps' metals at radius plus this (default="<<MetalBumpBias<<")" << endl;
   cerr << "-NONMETALBump#.#  'bumps' nonmetal at radius plus this (default="<<NonMetalBumpBias<<")" << endl;
   cerr << "-SEGIDmap \"seg,c...\"  assign chainID based on segment identifier field" << endl;
   cerr << "-Xplor            use Xplor conventions for naming polar hydrogens" << endl;
   cerr << "-OLDpdb 	      use the pre-remediation names for hydrogens" << endl;
   cerr << "-BBmodel	      expects a backbone only model and will build HA hydrogens on Calpha truncated residues" <<endl;
   cerr << "-NOCon            drop conect records" << endl;
   cerr << "-LIMIT#           max seconds to spend in exhaustive search (default="<< ExhaustiveLimit <<")" << endl;
   cerr << "-NOTICKs          do not display the set orientation ticker during processing" << endl;
   cerr << "-SHOWSCore        display scores for each orientation considered during processing" << endl;
   cerr << "-DROP_HYDROGENS_ON_ATOM_RECORDS drop hydrogens on incoming ATOM records before other processing" << endl;
   cerr << "-DROP_HYDROGENS_ON_OTHER_RECORDS drop hydrogens on incoming non-ATOM records before other processing" << endl;
   cerr << "-NO_ADD_WATER_HYDROGENS don't add hydrogens on incoming HOH records but do other processing" << endl;
   cerr << "                   (the default is to add water hydrogens even if hydrogens have been dropped)" << endl;
   cerr << "-NO_ADD_OTHER_HYDROGENS don't add hydrogens on incoming non-HOH records but do other processing" << endl;
   cerr << "                   (the default is to add other hydrogens even if hydrogens have been dropped)" << endl;
   cerr << "-NOOPT            do not perform optimizations, only drop/add hydrogens" << endl;
   cerr << "-FIX \"filename\"   if given, file specifies orientations for adjustable groups" << endl;
   cerr << "-DB \"filename\"    file to search for het info" << endl;
   cerr << "                        (default=\""<<DBfilename<<"\")" << endl;
   cerr << "note: can also redirect with unix environment variable: REDUCE_HET_DICT" << endl;
  }
   cerr << endl;
   cerr << "-STRING           pass reduce a string object from a script, must be quoted" << endl;
   cerr << "usage: from within a script, reduce -STRING \"_name_of_string_variable_\"" << endl << endl;
   cerr << "-Quiet            do not write extra info to the console" << endl;
   cerr << "-REFerence        display citation reference" << endl;
   cerr << "-Version          display the version of reduce" <<endl;
   cerr << "-Changes          display the change log" <<endl;
   cerr << "-Help             the more extensive description of command line arguments" << endl;
   exit(2);
}

void reduceChanges(bool showAll) { /*changes*/
   cerr  <<  "History: this new version of reduce has been completely re-written" << endl;
   cerr  <<  "  to include het groups and rotations" << endl ;
   cerr  << endl;
   cerr  << "10/ 1/97 - jmw - changed rotation score to dot score, modified geometry," << endl;
   cerr  << "                 of CH2 and mcNH angles to be compatible with ECEPP," << endl;
   cerr  << "                 and added option to permit all methyls to rotate." << endl;
   cerr  << "10/22/97 - jmw - fixed bug in rotation code which did not separate" << endl;
   cerr  << "                 alternate conformations" << endl;
   cerr  << "11/12/97 - jmw - support for flips of HIS and ASN/GLN residues," << endl;
   cerr  << "                 orientable waters and Hbonds to aromatics" << endl;
   cerr  << "                 (includes some analysis of pairing of flips and rots)" << endl;
   cerr  << "12/10/97 - jmw - fixed waters to have phantom H atoms," << endl;
   cerr  << "                 better rotational searching, metal covalent radii" << endl;
   cerr  << " 2/ 7/98 - jmw - broke out motions and re-wrote search," << endl;
   cerr  << "                 use only A conf, fix penalty" << endl;
   cerr  << " 2/18/98 - jmw - fixup ambiguous atom names, no creation of altB Hs," << endl;
   cerr  << "                 MTO waters, added -BUILD and other flags," << endl;
   cerr  << "                 better recognise S-S or C-O-C bonds, etc." << endl;
   cerr  << " 3/ 1/98 - jmw - fix orientations" << endl;
   cerr  << " 3/25/98 - jmw - updated metal radii, added HET dict env var," << endl;
   cerr  << "                 expanded fixed orientations, changed -OH default" << endl;
   cerr  << " 4/ 8/98 - jmw - fixed bug, occupancy of water phantomHs now > 0.0 !" << endl;
   cerr  << " 4/18/98 - jmw - bad bump check, short water h, new limit on HB, ..." << endl;
   cerr  << " 4/26/98 - jmw - new hires search strategy" << endl;
   cerr  << " 4/29/98 - jmw - fixed min rotation angle and penalty bugs" << endl;
   cerr  << " 5/13/98 - jmw - fixed another min rotation angle bug where the coarse score" << endl;
   cerr  << "                 is the best that can be obtained. Also updated output to" << endl;
   cerr  << "                 the header to document K/C/X/F categories." << endl;
   cerr  << " 5/15/98 - jmw - fine tuned cases where group is fixed by metal or modification" << endl;
   cerr  << " 7/ 3/98 - jmw - fixed bug: ASN/GLN sc C=O carbon radii was 1.75, now 1.65" << endl;
   cerr  << " 8/ 7/98 - jmw - stop putting notes in the segID, elem & charge fields" << endl;
   cerr  << "                 because PDB use of these fields are now being expanded" << endl;
   cerr  << "                 by naming convention differences with XPLOR and XPLORs" << endl;
   cerr  << "                 use of SEGID rather than CHAINID" << endl;
   cerr  << " 8/12/98 - jmw - changed how we figure num of permutations" << endl;
   cerr  << " 8/14/98 - jmw - now search for het dictionary in current directory" << endl;
   cerr  << " 9/ 1/98 - jmw - worked on portability by cleaning up g++ warnings" << endl;
   cerr  << " 9/13/98 - jmw - figured out a work-around to g++ template and" << endl;
   cerr  << "                 static const class initializer problems" << endl;
   cerr  << " 1/ 7/99 - jmw - consolidate Linux and Mac changes with SGI src" << endl;
   cerr  << " 3/16/99 - jmw - extended atom name parsing for wierd HETs with col1 ABCDEFGs" << endl;
   cerr  << " 4/ 6/99 - jmw - added OW to list of names for oxygen in water" << endl;
   cerr  << " 7/ 7/99 - jmw - Improved portability (near->nearr, List <T>::, Makefiles)," << endl;
   cerr  << "                 made penalty 1.0, added reference, changed basic io hooks" << endl;
   cerr  << " 9/ 2/99 - jmw - Updated Makefile list and Utility.h for DEC alpha support" << endl;
   cerr  << "10/ 5/99 - jmw - Updated a Makefiles and the main in reduce.C for sgi6.5 compiler" << endl;
   cerr  << " 8/ 4/00 - jmw - Added -segid flag to support segment identifiers" << endl;
   cerr  << "10/30/00 - jmw - Modified Seq mergesort const/non const stuff" << endl;
   cerr  << " 4/19/01 - jmw - Added support for left justified A/T/U/C/G for nucleic acids" << endl;
   cerr  << "                 which fixed a rare but nasty bug in the water recognition" << endl;
   cerr  << " 5/ 7/01 - jmw - Stopped output with -quiet while fixing orientation" << endl;
   cerr  << "                 and finally dealt with string literal conversion messages" << endl;
   cerr  << " 5/31/01 - jmw - Pass signal of abandoned clique search in rc," << endl;
   cerr  << "                 added new flag to allow mapping of segids to chains," << endl;
   cerr  << "                 changed properties of Nterminal fragment nitrogens to acceptor" << endl;
   cerr  << "10/ 4/01 - jmw - fiddled to get compiled on RH linux 7 (gcc 2.96)" << endl;
   cerr  << " 5/24/02 - jmw - added control over hbump" << endl;
   cerr  << " 1/23/03 - jmw -v2.16 - Changed global MinChargedHBgap: 0.4 => 0.8" << endl;
   cerr  << "                        Changed phosphorus properties to drop acceptor status" << endl;
   cerr  << " 3/06/03 - jmw -v2.17 - Cleaning declarations found by Leo and Andrew such as" << endl;
   cerr  << "                        adding compile flag -DOLD_CXX_DEFNS to select const longs" << endl;
   cerr  << "                        instead of std::_Ios_Fmtflags in AtomPositions.C" << endl;
   cerr  << " 3/31/03 - jmw -v2.18 - Fixed spelling mistake by renaming cuttoff to cutoff throughout." << endl;
   cerr  << "                        Edited help for -H2OBcutoff# to suggest an integer value" << endl;
   cerr  << "                        Made -H2OOCCcutoff#.# parse a real (was integer)" << endl;
   cerr  << " 4/ 3/03 - jmw -v2.19 - Moved -Help to the end of the command line parsing." << endl;
   cerr  << " 6/ 3/03 - jmw -v2.20 - updated to isoC++ style std includes: #include <cstring>," << endl;
   cerr  << "                        changed String class to Stringclass class," << endl;
   cerr  << "                        nice speedups from Jack Snoeyink" << endl;
   cerr  << " 6/ 4/03 - jmw -v2.21 - added out for OLD_STD_HDRS for sgi plus other sgi polishing" << endl;
   cerr  << "11/ 7/03 - jmw -v2.22 - fixed bug with creating hydrogens at the end of a triple bond" << endl;
   cerr  << "                        and supressed the warnings about connecting nucleic acid bases" << endl;
   cerr  << endl;
   cerr  << "040509dcr reduce.C ln 455: NO renumber, 1459: NO RXR msg" << endl;
   cerr  << ""<< endl;
   cerr  << "041113dcr reduce.C reconstructed main to loop over any NMR models in the file" << endl;
   cerr  << "          a specific model can still be specified." << endl;
   cerr  << endl;
   cerr  << "changes" << endl;
   cerr  << endl;
   cerr  << "050314 dcr incorporated Mike's version of 030703, to wit:" << endl;
   cerr  << endl;
   cerr  << " 6/15/06 - apl -v3.0  - incorporated decomposition of scoring function into interactions of" << endl;
   cerr  << "                        atom singles, atom pairs, atom tripples and all other" << endl;
   cerr  << "                        higher-order atom interactions.  incorporated" << endl;
   cerr  << "                        dynamic programming. incorporated changes from v2.21.mod_dcr" << endl;
   cerr  << "                        disabling hydrogen atom numbering and skipInfo. removing" << endl;
   cerr  << "                        PDBrecNAMEout as in v2.999. Disabling hydrogen bonds to his" << endl;
   cerr  << "                        aromatic carbon atoms." << endl;
   cerr  << " 6/19/06 - apl -v3.01- decomposing the scoring function in terms of which dots should" << endl;
   cerr  << "                       be scored with which hyperedges.  Incorporating S3 reduction rules" << endl;
   cerr  << "                       into dynamic programming.  Incorporating code to handle 4-way overlap" << endl;
   cerr  << "                       (but not five way overlap) though I have not observed any 4-way" << endl;
   cerr  << "                       overlap using the new decomposition scheme (it would show up in the previous scheme)." << endl;
   cerr  << " 6/24/06 - apl -v3.02- incorporating the additions to v2.23 that deal with multiple NMR models" << endl;
   cerr  << "                       in a single file.  Altering output from dp to be less intrusive and a" << endl;
   cerr  << "                       little more informative." << endl;
   cerr  << "                       Adding new reduction rule for vertices with exactly 1 state (i.e. no" << endl;
   cerr  << "                       real options) which speeds up dynamic programming for the second-round" << endl;
   cerr  << "                       of optimizations that calculate the optimal network states for sub-optimal" << endl;
   cerr  << "                       flip states.  GLN and ASN will have 1 state in these optimizations." << endl;
   cerr  << "7/03/06 - apl -      - Fixing USER MOD Set ordering in output PDB.  Fixing 'flip' records in columns" << endl;
   cerr  << "                       82 to 85 for the atoms on flipped residues.  Fixing atom placement plan bug" << endl;
   cerr  << "                       in genHydrogens() that manifested itself as aberant behavior when presented" << endl;
   cerr  << "                       with several (3 or more) alternate conformations." << endl;
   cerr  << "7/09/06 - apl -      - Adding #include <cassert> for gcc3.3.4 builds" << endl;
   cerr  << "7/11/06 - apl -      - Fixing 'node with 0 states' bug." << endl;
   cerr  << "10/19/06 - apl - v3.03- changing HIS carbons to regular carbons and not arromatic carbons" << endl;
   cerr  << "10/20/06 - apl -      - fixing bug in optimization code that failed to keep optimal network states" << endl;
   cerr  << "                       when a network was forced to incur a penalty." << endl;
   cerr  << "3/ 7/07 - apl -        Bug fix: do not add hydrogens to 'N' on non-amino acids" << endl;
   cerr  << "                       fixes 3H bug on SAC in 1b0b.pdb" << endl;
   cerr  << "3/ 7/07 - apl -        Bug fix: march ResBlk::_insertPtr backwards at the end of constructor" << endl;
   cerr  << "                       even when _insertPtr has reached the end of the rlst.  This ensures all" << endl;
   cerr  << "                       residues are protonated, instead of all but the last one.  Fixes Loren's" << endl;
   cerr  << "                       bug in trying to protonate nicotine when it was by itself in a .pdb." << endl;
   cerr  << "3/ 7/07 - apl -        inline distanceSquared in toolclasses/Point3d.h for a 10% speedup." << endl;
   cerr  << endl;
   cerr  << "7/ 7/07 - jjh & rmi -  incoporated new hydrogen names in StdResH.cpp for remediated pdb files" << endl;
   cerr  << "               v3.10   followed the example of the -Xplor flag and added a -OLDpdb flag to allow" << endl;
   cerr  << "                       output of new (default) or old (pre-remediation) hydrogen names" << endl;
   cerr  << "7/13/07 - jjh & rmi -  added -BBmodel flag allows addition of hydrogens on Calpha of truncated" << endl;
   cerr  << "                       amino acids" << endl;
   cerr  << "7/16/07 - jjh & rmi -  added fixAtomName() to ElementInfo.cpp to check for Hg, Ho, and Hf atoms" << endl;
   cerr  << "                       the corresponding warnings are no longer output" << endl;
   cerr  << "7/31/07 - rmi -        Bug fix: changed the logic in selecting Xplor vs. PDBv2.3 vs. PDBv3.0 atom" << endl;
   cerr  << "                       names in reduce.cpp and flipmemo.cpp. Updated the Version and shortVersion" << endl;
   cerr  << "7/31/07 - rmi -        Bug fix: changed the logic in selecting Xplor vs. PDBv2.3 vs. PDBv3.0 atom" << endl;
   cerr  << "                       names in reduce.cpp and flipmemo.cpp. Updated the Version and shortVersion" << endl;
   cerr  << "                       strings and changed the year to the four digit year in the versionString." << endl;
   cerr  << "8/01/07 - rmi          Added the -Version and -Changes flags. The -Changes flag uses a new function" << endl;
   cerr  << "                       reduceChanges and follows the format of reduceHelp" << endl;
   cerr  << "                       Also commented out non-RNA non-DNA 'nucleic acids' from StdResH.cpp to fix" << endl;
   cerr  << "                       doubling of backbone hydrogens." << endl;
   cerr  << "8/14/07 - rmi          Fixed bug in StdResH.cpp to add hydrogens on MSE. The change was to add SE" << endl;
   cerr  << "                       as well as SED as the name for the selenium atom" <<endl;
   cerr  << "8/18/07 - rwgk         Patched Elementinfo.cpp for compiler problems: (a)Visual C++ warning and (b)Tru64 error" << endl;
   cerr  << "svn rev 67, 68         (a)threw runtime error on fixAtomName() (b)added 'using std::sprintf'" << endl;
   cerr  << "8/29/07 - rmi          Modified the reduce het dict so that hydrogens are not built on carboxylates" << endl;
//   cerr  << "9/25/07 - rmi          Added a flag ADDNHATGAP which allows a single hydrogen to be built at the N-termini of chain breaks" << endl;
   cerr  << "                       added break-amide to StdResH to treat these amides as a special case" << endl;
   cerr  << "10/3/07 - rmi          Added support for Hybrid36 atom and residue numbers" << endl;
   cerr  << "11/1/07 - rmi          Added support for two character chainIds" << endl;
   cerr  << "11/7/07 - rmi          Reverted changes to pdb_sscanf.cpp and changed write format for chains to %-2s" <<endl;
   cerr  << "11/14/07- rmi          Several BUG fixes:  in PDBrec.h getAtomDescr() wants seqNum not serialNum" <<endl;
   cerr  << "                         in AtomPositions.cpp change format to %-2.2s for two character chains" <<endl;
   cerr  << "                         fixed format strings in read_format.i and write_format.i"  <<endl;
   cerr  << "                         explicitly added Hy36seqNum to pdb++.h" << endl;
   cerr  << "02/14/08 - vbc & rmi   Bob's fixes for dealing with windows line ending files (no fix for mac files though)." <<endl;
   cerr  << "             & jjh       I (vbc) attempted to fix some of the warnings for new hydrogen names so molprobity isn't" <<endl;
   cerr  << "                         quite as swamped." <<endl;
   cerr  << "02/20/08 - vbc & rmi   Fixed double H bug from Bob's previous correction" << endl;
   cerr  << "02/28/08 - jmw & jjh   Fixed altID bug for H(alpha) when connecting atoms only non-blank altID and only one total conformation" <<endl;
   cerr  << "04/11/08 - jjh         Added -STRING flag to allow scripts in Perl/Python to pass a string to reduce for processing.  Output still directed to standard out." << endl;
   cerr  << "04/28/08 - jjh          fixed 4 character Deuterium recognition w/ PDB 3.0 names" << endl;
   cerr  << "08/21/08 - jjh          added -CHARGEs flag to control charge state output - off by default" << endl;
   cerr  << "11/06/09 - jjh         added -FLIP and -NOFLIP flag" << endl;
   cerr  << "09/15/10 - wba         'reducer' versions are for evaluation of changes to H bond-distances" << endl;
   cerr  << "11/18/11 - jjh         Overhauled handling of alternate conformers" << endl;
   cerr  << "08/15/12 - aram        added -MAXAromdih cutoff to support rotating aromatic methyls" << endl;
   cerr  << "2012/08/23 - lnd & vbc New, shorter, H bond distances and van der Waals - version 3.17" << endl;
   cerr  << "2012/09/05 - gjk       New, shorter, H bond distances and van der Waals for nucleic acids" << endl;
   cerr  << "2013/01/10 - jjh       Remove N-H atom from N-terminal residue of amino acid chain" << endl;
   cerr  << "2013/01/16 - jjh       Added -NUCLEAR flag, which uses nuclear distances/vdW rather" << endl;
   cerr  << "                         the default electron cloud distances/vdW for H placement" << endl;
   cerr  << "2013/01/22 - jjh       fixed handling of group rotation for alternates" << endl;
   cerr  << "2013/02/19 - wba       updated version number and date" << endl;
   cerr  << "2013/03/26 - jjh       fixed bugs related to aromatic methyl rotations" << endl;
   cerr  << "2013/03/27 - jjh v3.22 fixed bug where number of brute force node searches" << endl;
   cerr  << "                        exceeded the system size of an int type" << endl;
   cerr  << "2013/05/09 - jjh v3.23 support for segid instead of chainid added" << endl;
   cerr  << "2013/05/20 - jjh       fixed -trim handling of H5'' atom in RNA" << endl;
   cerr  << "2013/07/05 - jjh v3.24 fixed handling of ANISOU records when SEGIDs in use," << endl;
   cerr  << "                        and fixed removal of redundant N-terminal H atoms" << endl;
   cerr  << "2013/07/11 - jjh       fixed calculation of neighbor atoms with invalidated records" << endl;
   cerr  << "2013/07/11 - jjh       fixed enforcement of time limit for clique search" << endl;
   cerr  << "2013/07/18 - jjh       fixed handling of multiple models when first model" << endl;
   cerr  << "                        is not model 1" << endl;
   cerr  << "2013/07/24 - jjh       cleaned up command line parsing function, removed redundant -FLIPs option" << endl;
   cerr  << "2015/08/19 - sj        fixed multimodel handling. Header now printed only once, and" << endl;
   cerr  << "                       everything else not within MODEL ENDMDL printed in the end" << endl;
   cerr  << "2015/??    - v3.3, sj & cjw, New flip method to replace nqh_minimize" << endl;
   cerr  << "2015/??    - bjh       test system" << endl;
   cerr  << "2016/06/02 - cjw       set default behavior *not* to rotate methionine methyls, added -DOROTMET flag" << endl;
   cerr  << "2020/10/07 - rmt       Added fine-grained command-line control over hydrogen drop/add and optimization" << endl;
   cerr  << "2020/10/08 - rmt       Adjusted vector operations to enable running in debug compile" << endl;
   cerr  << "2020/10/22 - sbstnk    Fix SIGUIJ output format" << endl;
   cerr  << "2020/10/22 - sbstnk    Return code '2' for critical errors to distinguish from timeouts" << endl;
   cerr  << "2020/10/29 - sbstnk    Avoid floating-point accuracy issues when doing some flips" << endl;
   cerr  << "2020/11/03 - cjw       Updated het dictionary" << endl;
   cerr  << "2020/11/19 - rmt       Verifies that atoms are valid before adjusting them" << endl;
   cerr  << "2020/12/03 - rmt       Version 3.8 adds hydrogens on all alternates even if not adjusting" << endl;
   cerr  << "2020/12/03 - rmt       Version 3.9 fixed bugs in abandoned-clique reporting" << endl;
   cerr  << "2020/12/18 - rmt       Version 3.10 wraps the library for use in Python" << endl;
   cerr  << endl;
   exit(2);
}

void establishHetDictionaryFileName(void) {
	int i =DBfilename.find_last_of(DIRECTORY_SEP_CHAR);
	std::string localfile = DBfilename.substr(i+1, DBfilename.length());
#ifdef NOSYSSTATS
	const int rc = 1; // force the name on the MAC
#else
	struct stat filestatbuf;
	const int rc = stat(localfile.c_str(), &filestatbuf);
#endif
	if (rc == 0) {
  		DBfilename = localfile;
	}
	else {
		const char *hetdbAlias = getenv("REDUCE_HET_DICT");
		if (hetdbAlias && (strlen(hetdbAlias) > 0)) {
			DBfilename = hetdbAlias;
		}
	}
}

#ifdef MACAPPL
char *getOptionsOnTheMac() {
   establishHetDictionaryFileName(); // this may need modifying

   BuildHisHydrogens  = TRUE; // this block is the same as the old -build
   SaveOHetcHydrogens = TRUE;
   RotExistingOH      = TRUE;
   DemandFlipAllHNQs  = TRUE;

   return "test.pdb";
}
#endif


char* parseCommandLine(int argc, char **argv) {
  char *pdbfile = NULL;
  int nfile = 0, n;

  establishHetDictionaryFileName();

  if (argc <= 1) {
    reduceHelp(FALSE);
  }
  for (int i = 1; i < argc; i++) {
    char *p = argv[i];
    if (p[0] == '-') {
      if (p[1] == '\0') {
        nfile = 1;
        pdbfile = NULL; // i.e. standard input
      }
      else if(compArgStr(p+1, "STRING", 6)){
        StringInput = TRUE;
      }
      else if(compArgStr(p+1, "FLIP", 4)){
        BuildHisHydrogens  = TRUE;
        SaveOHetcHydrogens = TRUE;
        RotExistingOH      = TRUE;
        DemandFlipAllHNQs  = TRUE;
      }
      else if(compArgStr(p+1, "NOFLIP", 6)){
        PenaltyMagnitude=9999;
        BuildHisHydrogens  = TRUE;
        SaveOHetcHydrogens = TRUE;
        RotExistingOH      = TRUE;
        DemandFlipAllHNQs  = TRUE;
      }
      else if(compArgStr(p+1, "BUILD", 5)){
        BuildHisHydrogens  = TRUE;
        SaveOHetcHydrogens = TRUE;
        RotExistingOH      = TRUE;
        DemandFlipAllHNQs  = TRUE;
      }
      else if((n = compArgStr(p+1, "NOBUILD", 7))){
        PenaltyMagnitude = parseReal(p, n+1, 10, PenaltyMagnitude);
        if(PenaltyMagnitude < 0){
          cerr << "!!ERROR!!" << endl;
          cerr << "Penalty Magnitude for -NOBUILD must be > 0" << endl << endl;
          reduceHelp(FALSE);
        }
        // PenaltyMagnitude = 200;      9999 in molprobity
        BuildHisHydrogens  = TRUE;
        SaveOHetcHydrogens = TRUE;
        RotExistingOH      = TRUE;  //  not used in molprobity
        DemandFlipAllHNQs  = TRUE;
      }
      else if ((n = compArgStr(p+1,"RENAMEFLIP",10))){ // SJ - 09/25/2015 added to set the RenameFlip flag to TRUE. See top of the file for intended behavior of the flag
          RenameFlip=TRUE;
      }
      else if((n = compArgStr(p+1, "Version", 1))){
        cerr << shortVersion << endl;
        exit(2);
      }
      else if((n = compArgStr(p+1, "Changes", 1))) {
        reduceChanges(TRUE);
      }
      else if((n = compArgStr(p+1, "Quiet", 1))){
        Verbose = FALSE;
      }
      else if((n = compArgStr(p+1, "NOTICKs", 6))){
        ShowCliqueTicks = FALSE;
      }
      else if((n = compArgStr(p+1, "SHOWSCore", 6))){
        ShowOrientScore = TRUE;
      }
      else if((n = compArgStr(p+1, "NOCon", 3))){
        KeepConnections = FALSE;
      }
      else if((n = compArgStr(p+1, "NOROTMET", 8))){
        OKProcessMetMe = FALSE;
      }
      else if((n = compArgStr(p+1, "DOROTMET", 8))){
        OKProcessMetMe = TRUE;
      }
      else if((n = compArgStr(p+1, "NOADJust", 5))){
        OKtoAdjust = FALSE;
      }
      else if((n = compArgStr(p+1, "HIS", 3))){
        BuildHisHydrogens = TRUE;
      }
      else if((n = compArgStr(p+1, "OH", 2))){
        SaveOHetcHydrogens = TRUE;
      }
      else if((n = compArgStr(p+1, "NOOH", 4))){
        SaveOHetcHydrogens = FALSE;
      }
      else if((n = compArgStr(p+1, "Xplor", 1)) && !UseOldNames){
        UseXplorNames = TRUE;
      }
      else if((n = compArgStr(p+1, "Xplor", 1)) && UseOldNames){
        cerr << "Cannot use both -Xplor and -OLDpdb flags" << endl;
        exit(2);
      }
      else if((n = compArgStr(p+1, "OLDpdb", 3)) && ! UseXplorNames){
        UseOldNames = TRUE;
        DBfilename = HET_DICTOLD;
      }
      else if((n = compArgStr(p+1, "OLDpdb", 3)) && UseXplorNames){
        cerr << "Cannot use both -Xplor and -OLDpdb flags" << endl;
        exit(2);
      }
      else if((n = compArgStr(p+1, "BBmodel", 2))){
        BackBoneModel = TRUE;
      }
      else if((n = compArgStr(p+1, "Trim", 1))){
        RemoveATOMHydrogens = TRUE;
	RemoveOtherHydrogens = TRUE;
	AddWaterHydrogens = FALSE;
	AddOtherHydrogens = FALSE;
	StopBeforeOptimizing = TRUE;
      }
      else if((n = compArgStr(p+1, "Keep", 1))){
        StandardizeRHBondLengths = FALSE;
      }
      else if((n = compArgStr(p+1, "ALLMETHYLS", 10))){
        DemandRotAllMethyls = TRUE;
      }
      else if((n = compArgStr(p+1, "ROTEXist", 5))){
        DemandRotExisting = TRUE;
      }
      else if((n = compArgStr(p+1, "ADDNHATGAP", 10))){
        NeutralTermini = TRUE;
      }
      else if((n = compArgStr(p+1, "ROTNH3", 6))){
        DemandRotNH3 = TRUE;
      }
      else if((n = compArgStr(p+1, "NOROTNH3", 8))){
        DemandRotNH3 = FALSE;
      }
      else if((n = compArgStr(p+1, "ROTEXOH", 7))){
        RotExistingOH = TRUE;
      }
      /*else if((n = compArgStr(p+1, "FLIPs", 4))){
        DemandFlipAllHNQs = TRUE;
      }*/ //removed redundant "FLIPs" flag - JJH 130724
      else if((n = compArgStr(p+1, "SEGIDmap", 5))){
        if (++i < argc) {
          UseSEGIDtoChainMap = TRUE;
          PDBrec::InstallMapOfSEGIDstoChains(argv[i]);
        }
        else {
          cerr << "no mapping info after -SEGIDmap flag" << endl;
        }
      }
      else if((n = compArgStr(p+1, "MAXAromdih", 1))){ // - Aram 07/24/12
        MaxAromRingDih = parseInteger(p, n+1, 10);
      }
      else if((n = compArgStr(p+1, "Nterm", 1))){
        MinNTermResNo = parseInteger(p, n+1, 10);
      }
      else if((n = compArgStr(p+1, "ONLTA", 5))){
        DoOnlyAltA = TRUE;
      }
      else if((n = compArgStr(p+1, "ALLALT", 6))){
        DoOnlyAltA = FALSE;
      }
      else if((n = compArgStr(p+1, "CHARGEs", 6))){
        ShowCharges = TRUE;
      }
      else if((n = compArgStr(p+1, "NOHETh", 5))){
        ProcessConnHydOnHets = FALSE;
      }
      else if((n = compArgStr(p+1, "DENSity", 4))){
        VdwDotDensity = parseReal(p, n+1, 10, VdwDotDensity);
      }
      else if((n = compArgStr(p+1, "PENalty", 3))){
        PenaltyMagnitude = parseReal(p, n+1, 10, PenaltyMagnitude);
      }
      else if((n = compArgStr(p+1, "RADius", 3))){
        ProbeRadius = parseReal(p, n+1, 10, ProbeRadius);
      }
      else if((n = compArgStr(p+1, "NBonds", 2))){
        NBondCutoff = parseInteger(p, n+1, 10);
      }
      else if((n = compArgStr(p+1, "OCCcutoff", 3))){
        OccupancyCutoff = parseReal(p, n+1, 10, OccupancyCutoff);
      }
      else if((n = compArgStr(p+1, "H2OBcutoff", 4))){
        WaterBcutoff = 1.0 * parseInteger(p, n+1, 10);
      }
      else if((n = compArgStr(p+1, "H2OOCCcutoff", 6))){
        WaterOCCcutoff = parseReal(p, n+1, 10, WaterOCCcutoff);
      }
      else if((n = compArgStr(p+1, "HBREGcutoff", 5))){
        MinRegHBgap = parseReal(p, n+1, 10, MinRegHBgap);
      }
      else if((n = compArgStr(p+1, "HBCHargedcutoff", 4))){
        MinChargedHBgap = parseReal(p, n+1, 10, MinChargedHBgap);
      }
      else if((n = compArgStr(p+1, "BADBumpcutoff", 4))){
        BadBumpGapCut = parseReal(p, n+1, 10, BadBumpGapCut);
      }
      else if((n = compArgStr(p+1, "NONMETALBump", 9))){
        NonMetalBumpBias = parseReal(p, n+1, 10, NonMetalBumpBias);
      }
      else if((n = compArgStr(p+1, "METALBump", 6))){
        MetalBumpBias = parseReal(p, n+1, 10, MetalBumpBias);
      }
      else if((n = compArgStr(p+1, "GAPERROR", 8))){
        GapWidth = parseReal(p, n+1, 10, GapWidth);
        if (GapWidth > 1.4) {
          cerr << "Max allowed HalfGapWidth is 1.4" << endl;
          exit(2);
        }
      }
      else if((n = compArgStr(p+1, "REFerence", 3))){
        cerr << "Please cite: " << referenceString << endl;
        cerr << "For more information see " << electronicReference << endl;
        exit(2);
      }
      else if((n = compArgStr(p+1, "FIX", 3))){
        if (++i < argc) {
          OFile = argv[i];
        }
        else {
          cerr << "no filename after -FIX flag" << endl;
        }
      }
      else if((n = compArgStr(p+1, "DB", 2))){
        if (++i < argc) {
          DBfilename = argv[i];
        }
        else {
          cerr << "no filename after -DB flag" << endl;
        }
      }
      else if((n = compArgStr(p+1, "LIMITsearch", 5))){
        ExhaustiveLimit = parseInteger(p, n+1, 10);
      }
      else if((n = compArgStr(p+1, "NUClear", 3))){
        UseNuclearDistances = TRUE;
      }
      else if(compArgStr(p+1, "Help", 1)){ // has to be after all the other -HXXXs
        reduceHelp(TRUE);
      }
      else if ((n = compArgStr(p + 1, "NOOPT", 4))) {
	StopBeforeOptimizing = TRUE;
      }
      else if ((n = compArgStr(p + 1, "NO_ADD_WATER_HYDROGENS", 8))) {
        AddWaterHydrogens = FALSE;
      }
      else if ((n = compArgStr(p + 1, "NO_ADD_OTHER_HYDROGENS", 8))) {
	AddOtherHydrogens = FALSE;
      }
      else if ((n = compArgStr(p + 1, "DROP_HYDROGENS_ON_ATOM_RECORDS:", 19))) {
        RemoveATOMHydrogens = TRUE;
      }
      else if ((n = compArgStr(p + 1, "DROP_HYDROGENS_ON_OTHER_RECORDS:", 19))) {
	RemoveOtherHydrogens = TRUE;
      }
      else {
        cerr << "unrecognized flag, \"" << p << "\", ignored." << endl;
      }
    }
    else if (StringInput == TRUE){
      pdbfile = p;
      nfile = 1;
    }
    else if (nfile <= 0) {
      pdbfile = p;
      nfile = 1;
    }
    else {
      cerr << "unrecognized parameter, \"" << p << "\", ignored." << endl;
    }
  }
  if (nfile != 1) { reduceHelp(FALSE); }
  return pdbfile;
}

int main(int argc, char **argv) {

    int ret = -1;
    
    char *pdbFile = parseCommandLine(argc, argv);
    
    // See whether we're supposed to read from a PDB file, (pdbFile is a string name and
    // StringInput is FALSE), from a string (StringInput is TRUE and pdbFile is set to
    // the actual string to read from), or from standard input (the - option tells us to
    // read from standard input by setting pdbfile to NULL and StringInput to FALSE).
    // In any case, read the input into a string.
    std::istream *ifPtr;
    if (StringInput == TRUE) {
      if (Verbose) cerr << "Processing input string" << endl;
      ifPtr = new std::istringstream(pdbFile);
    } else if (pdbFile) {
      if (Verbose) cerr << "Processing file: \"" << pdbFile << "\"" << endl;
      ifPtr = new std::ifstream(pdbFile);
    } else {
      if (Verbose) cerr << "Processing file: --standard input--" << endl;
      ifPtr = &cin;
    }
    std::istreambuf_iterator<char> endBuf;
    std::string s(std::istreambuf_iterator<char>(*ifPtr), endBuf);
    if (ifPtr != &cin) {
      delete ifPtr;
    }

    // Read all models from the file into a list of lists of records.  This gives
    // us one list of PDB records for each model in the file.
    std::vector< std::list< std::shared_ptr<PDBrec> > > models = inputModels(s);
    if (models.size() == 0) {
      cerr << "Error: no input records" << endl;
      return 100;
    }

    // Process each model.
    for (std::vector< std::list< std::shared_ptr<PDBrec> > >::iterator it = models.begin();
         it != models.end(); it++) {
      std::list< std::shared_ptr<PDBrec> > &m = *it;
      //=====================================================================
      // Removing hydrogens

      if (RemoveATOMHydrogens || RemoveOtherHydrogens) {
        dropHydrogens(m, RemoveATOMHydrogens, RemoveOtherHydrogens);
      }

      //=====================================================================
      // Adding hydrogens

      UseSEGIDasChain = checkSEGIDs(m);

      if (AddWaterHydrogens || AddOtherHydrogens) {
        if (Verbose) {
          if (BuildHisHydrogens) {
            cerr << "Building His ring NH Hydrogens." << endl;
          }
          if (SaveOHetcHydrogens) {
            cerr << "Building or keeping OH & SH Hydrogens." << endl;
            if (RotExistingOH && !DemandRotExisting) {
              cerr << "Rotating existing OH & SH Hydrogens" << endl;
            }
          }
          else {
            cerr << "Dropping existing OH & SH Hydrogens." << endl;
          }
          if (NeutralTermini) {
            cerr << "Adding \"amide\" Hydrogens to chain breaks." << endl;
          }
        }
      }
      /// @todo Only load this if we have one...
      CTab hetdatabase(DBfilename);

      DotSphManager dotBucket(VdwDotDensity);

      AtomPositions xyz(2000, DoOnlyAltA, UseXplorNames, UseOldNames, BackBoneModel,
        NBondCutoff, MinRegHBgap, MinChargedHBgap,
        BadBumpGapCut, dotBucket, ProbeRadius,
        PenaltyMagnitude, OccupancyCutoff,
        Verbose, ShowOrientScore, ShowCliqueTicks, cerr);

      //    NonConstListIter<PDBrec> infoPtr(records); // info on changes can be
                                                     // inserted before infoPtr
      std::list< std::shared_ptr<PDBrec> >::iterator infoPtr = m.begin();

      scanAndGroupRecords(m, xyz, infoPtr);

      // if some sidechain needs adjustment...
      std::vector<std::string> adjNotes;
      Tally._num_adj = 0;

      reduceList(hetdatabase, m, xyz, adjNotes);
      // SJ - All Hydrogens are generated at this time, but no flips and methyl rotations have happened.

      //=====================================================================
      // Performing optimizations

      if (!StopBeforeOptimizing) {
        ret = optimize(xyz, adjNotes);

        if (OKtoAdjust && xyz.numChanges() > 0) {
          xyz.describeChanges(m, infoPtr, adjNotes);
        }
      }

      //=====================================================================
      // Describing operations

      if (Verbose) {
        if (AddWaterHydrogens || AddOtherHydrogens) {
          cerr << "Found " << Tally._H_found << " hydrogens ("
            << Tally._H_HET_found << " hets)" << endl;
          cerr << "Standardized " << Tally._H_standardized << " hydrogens ("
            << Tally._H_HET_standardized << " hets)" << endl;
          cerr << "Added " << Tally._H_added << " hydrogens ("
            << Tally._H_HET_added << " hets)" << endl;
        }
        if (RemoveATOMHydrogens || RemoveOtherHydrogens) {
          cerr << "Removed " << Tally._H_removed << " hydrogens ("
            << Tally._H_HET_removed << " hets)" << endl;
        }
        if (!StopBeforeOptimizing) {
          if (Tally._num_adj > 0) {
            cerr << "Adjusted " << Tally._num_adj << " group(s)" << endl;
          }
          if (Tally._num_renamed > 0) {
            cerr << "Renamed and marked " << Tally._num_renamed
              << " ambiguous 'A' atom name(s)" << endl;
          }
        }
      }
    }

    // SJ This is where the outputrecords should be called for all records.
    //This function now prints all models
    outputRecords_all(models, cout);
    
    if (Verbose) { // copied over from processPDBfile to here, will be printed only once
     cerr << "If you publish work which uses reduce, please cite:"
     << endl << referenceString << endl;
     cerr << "For more information see " << electronicReference << endl;
     }
    
   return ret;
}
