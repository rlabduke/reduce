TOOLHOME    =../toolclasses
PDBHOME     =../libpdb
INSTALLHOME =..
DICT_HOME   =/usr/local/reduce_het_dict.txt
DICT_LOC    = -DHET_DICTIONARY='"$(DICT_HOME)"'
PROG_FLGS   = -DAROMATICS_ACCEPT_HBONDS \
              -DCHARFUNCMACROS \
	      -DBOOLPREDEFINED  -DBRACKETOPERPARMS -DLEFT_JUSTIFY_NUC_RES_OK \
              -DLISTFRIENDFIX -DINCTEMPLATEDEFNS
#		-DOLD_CXX_DEFNS 
#              -DDEBUGLOOKFORTHREEWAYOVERLAP 
#              -DDEBUGMOVABLENODE -DDEBUGDOTCOUNTS
#              -DMINPREDEFINED -DMAXPREDEFINED\

CFLAGS = $(OPT) $(DEBUG) $(PROG_FLGS) $(DICT_LOC)\
         -I$(PDBHOME) -I$(TOOLHOME) 
#         -Wall

LFLAGS = -L$(PDBHOME) -lpdb++ -L$(TOOLHOME) -ltoolclasses -lm

OBJLIST = GraphToHoldScores.o reduce.o CTab.o ElementInfo.o StdResH.o ResBlk.o AtomConn.o \
            AtomPositions.o DotSph.o Mover.o RotMethyl.o RotDonor.o \
            FlipMemo.o CliqueList.o AtomDescr.o PDBrec.o MoveableNode.o \
            

CXX		= g++

.C.o:
	$(CXX) -c $*.C $(CFLAGS)

OPT             = -O3
DEBUG           = $(CXXDEBUGFLAGS)

all:    reduce

reduce: $(OBJLIST)
	$(CXX) -o $@ $(OBJLIST) $(LFLAGS)

clean:
	@rm -rf *.o *.ckp ii_files

install: reduce
	mv reduce $(INSTALLHOME)
	gunzip -c -d ../reduce_het_dict.txt.gz > $(DICT_HOME)

# DO NOT DELETE THIS LINE -- make depend uses it
AtomConn.o: $(TOOLHOME)/Point3d.h $(TOOLHOME)/utility.C \
         $(TOOLHOME)/utility.h AtomConn.C \
	 AtomConn.h ElementInfo.h
AtomDescr.o: AtomDescr.C AtomDescr.h
AtomPositions.o: $(PDBHOME)/pdb++.h $(TOOLHOME)/DisjointSets.h \
	 $(TOOLHOME)/Point3d.h $(TOOLHOME)/utility.C \
	 $(TOOLHOME)/utility.h AtomConn.h AtomPositions.C \
	 AtomPositions.h BumperPoint.h CliqueList.h DotSph.h \
	 ElementInfo.h FlipMemo.h LocBlk.h Mover.h PDBrec.h ResBlk.h \
	 RotDonor.h RotMethyl.h StdResH.h neighbors.C neighbors.h \
	 AtomDescr.h MoveableNode.h GraphToHoldScores.h
CTab.o: $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.C $(TOOLHOME)/utility.h AtomConn.h \
	 CTab.C CTab.h ElementInfo.h
CliqueList.o: $(PDBHOME)/pdb++.h \
	 $(TOOLHOME)/Point3d.h $(TOOLHOME)/utility.C \
	 $(TOOLHOME)/utility.h BumperPoint.h CliqueList.C \
	 CliqueList.h DotSph.h ElementInfo.h LocBlk.h Mover.h PDBrec.h \
	 neighbors.C neighbors.h
DotSph.o: $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.C $(TOOLHOME)/utility.h DotSph.C \
	 DotSph.h
ElementInfo.o: $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.C $(TOOLHOME)/utility.h AtomConn.h \
	 ElementInfo.C ElementInfo.h StdResH.h
FlipMemo.o: $(PDBHOME)/pdb++.h $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.C $(TOOLHOME)/utility.h AtomConn.h \
	 AtomPositions.h BumperPoint.h CliqueList.h DotSph.h \
	 ElementInfo.h FlipMemo.C FlipMemo.h LocBlk.h Mover.h PDBrec.h \
	 ResBlk.h neighbors.C neighbors.h
GraphToHoldScores.o: GraphToHoldScores.h GraphToHoldScores.C AtomPositions.h \
    AtomDescr.h
Mover.o: $(PDBHOME)/pdb++.h $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.C $(TOOLHOME)/utility.h \
	 BumperPoint.h DotSph.h ElementInfo.h LocBlk.h Mover.C Mover.h \
	 PDBrec.h ResBlk.h neighbors.C neighbors.h
MoveableNode.o: MoveableNode.C MoveableNode.h AtomDescr.h \
    GraphToHoldScores.h AtomPositions.h
PDBrec.o: $(PDBHOME)/pdb++.h \
	 $(TOOLHOME)/Point3d.h $(TOOLHOME)/utility.C \
	 $(TOOLHOME)/utility.h ElementInfo.h AtomDescr.h PDBrec.C PDBrec.h
ResBlk.o: $(PDBHOME)/pdb++.h $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.C $(TOOLHOME)/utility.h \
	 ElementInfo.h PDBrec.h ResBlk.C ResBlk.h
RotDonor.o: $(PDBHOME)/pdb++.h $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.C $(TOOLHOME)/utility.h \
	 AtomPositions.h BumperPoint.h CliqueList.h DotSph.h \
	 ElementInfo.h FlipMemo.h LocBlk.h Mover.h PDBrec.h ResBlk.h \
	 RotDonor.C RotDonor.h neighbors.C neighbors.h
RotMethyl.o: $(PDBHOME)/pdb++.h $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.C $(TOOLHOME)/utility.h \
	 AtomPositions.h BumperPoint.h CliqueList.h DotSph.h \
	 ElementInfo.h LocBlk.h Mover.h PDBrec.h ResBlk.h RotMethyl.C \
	 RotMethyl.h neighbors.C neighbors.h
StdResH.o: $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.C $(TOOLHOME)/utility.h AtomConn.h \
	 ElementInfo.h StdResH.C StdResH.h
neighbors.o: $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.C $(TOOLHOME)/utility.h LocBlk.h \
	 neighbors.C neighbors.h
reduce.o: $(PDBHOME)/pdb++.h $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.C $(TOOLHOME)/utility.h AtomConn.h \
	 AtomPositions.h BumperPoint.h CTab.h CliqueList.h DotSph.h \
	 ElementInfo.h LocBlk.h Mover.h PDBrec.h ResBlk.h StdResH.h \
	 neighbors.C neighbors.h reduce.C
# DO NOT DELETE THIS 2nd LINE -- make depend uses it
