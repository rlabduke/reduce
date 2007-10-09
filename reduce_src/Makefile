TOOLHOME    =../toolclasses
PDBHOME     =../libpdb
INSTALLHOME =..
DICT_HOME   =/usr/local/reduce_wwPDB_het_dict.txt
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
            FlipMemo.o CliqueList.o AtomDescr.o PDBrec.o MoveableNode.o hybrid_36_c.o\
            

CXX		= g++

.cpp.o:
	$(CXX) -c $*.cpp $(CFLAGS)

OPT             = -O3
DEBUG           = $(CXXDEBUGFLAGS)

all:    reduce

reduce: $(OBJLIST)
	$(CXX) -o $@ $(OBJLIST) $(LFLAGS)

clean:
	@rm -rf *.o *.ckp ii_files

install: reduce
	mv reduce $(INSTALLHOME)
	gunzip -c -d ../reduce_wwPDB_het_dict.txt.gz > $(DICT_HOME)

# DO NOT DELETE THIS LINE -- make depend uses it
AtomConn.o: $(TOOLHOME)/Point3d.h $(TOOLHOME)/utility.cpp \
         $(TOOLHOME)/utility.h AtomConn.cpp \
	 AtomConn.h ElementInfo.h
AtomDescr.o: AtomDescr.cpp AtomDescr.h
AtomPositions.o: $(PDBHOME)/pdb++.h $(TOOLHOME)/DisjointSets.h \
	 $(TOOLHOME)/Point3d.h $(TOOLHOME)/utility.cpp \
	 $(TOOLHOME)/utility.h AtomConn.h AtomPositions.cpp \
	 AtomPositions.h BumperPoint.h CliqueList.h DotSph.h \
	 ElementInfo.h FlipMemo.h LocBlk.h Mover.h PDBrec.h ResBlk.h \
	 RotDonor.h RotMethyl.h StdResH.h neighbors.cpp neighbors.h \
	 AtomDescr.h MoveableNode.h GraphToHoldScores.h
CTab.o: $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.cpp $(TOOLHOME)/utility.h AtomConn.h \
	 CTab.cpp CTab.h ElementInfo.h
CliqueList.o: $(PDBHOME)/pdb++.h \
	 $(TOOLHOME)/Point3d.h $(TOOLHOME)/utility.cpp \
	 $(TOOLHOME)/utility.h BumperPoint.h CliqueList.cpp \
	 CliqueList.h DotSph.h ElementInfo.h LocBlk.h Mover.h PDBrec.h \
	 neighbors.cpp neighbors.h
DotSph.o: $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.cpp $(TOOLHOME)/utility.h DotSph.cpp \
	 DotSph.h
ElementInfo.o: $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.cpp $(TOOLHOME)/utility.h AtomConn.h \
	 ElementInfo.cpp ElementInfo.h StdResH.h
FlipMemo.o: $(PDBHOME)/pdb++.h $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.cpp $(TOOLHOME)/utility.h AtomConn.h \
	 AtomPositions.h BumperPoint.h CliqueList.h DotSph.h \
	 ElementInfo.h FlipMemo.cpp FlipMemo.h LocBlk.h Mover.h PDBrec.h \
	 ResBlk.h neighbors.cpp neighbors.h
GraphToHoldScores.o: GraphToHoldScores.h GraphToHoldScores.cpp AtomPositions.h \
    AtomDescr.h
hybrid_36_c.o: hybrid_36_c.c hybrid_36_c.h
Mover.o: $(PDBHOME)/pdb++.h $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.cpp $(TOOLHOME)/utility.h \
	 BumperPoint.h DotSph.h ElementInfo.h LocBlk.h Mover.cpp Mover.h \
	 PDBrec.h ResBlk.h neighbors.cpp neighbors.h
MoveableNode.o: MoveableNode.cpp MoveableNode.h AtomDescr.h \
    GraphToHoldScores.h AtomPositions.h
PDBrec.o: $(PDBHOME)/pdb++.h \
	 $(TOOLHOME)/Point3d.h $(TOOLHOME)/utility.cpp \
	 $(TOOLHOME)/utility.h ElementInfo.h AtomDescr.h PDBrec.cpp PDBrec.h hybrid_36_c.h
ResBlk.o: $(PDBHOME)/pdb++.h $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.cpp $(TOOLHOME)/utility.h \
	 ElementInfo.h PDBrec.h ResBlk.cpp ResBlk.h
RotDonor.o: $(PDBHOME)/pdb++.h $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.cpp $(TOOLHOME)/utility.h \
	 AtomPositions.h BumperPoint.h CliqueList.h DotSph.h \
	 ElementInfo.h FlipMemo.h LocBlk.h Mover.h PDBrec.h ResBlk.h \
	 RotDonor.cpp RotDonor.h neighbors.cpp neighbors.h
RotMethyl.o: $(PDBHOME)/pdb++.h $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.cpp $(TOOLHOME)/utility.h \
	 AtomPositions.h BumperPoint.h CliqueList.h DotSph.h \
	 ElementInfo.h LocBlk.h Mover.h PDBrec.h ResBlk.h RotMethyl.cpp \
	 RotMethyl.h neighbors.cpp neighbors.h
StdResH.o: $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.cpp $(TOOLHOME)/utility.h AtomConn.h \
	 ElementInfo.h StdResH.cpp StdResH.h
neighbors.o: $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.cpp $(TOOLHOME)/utility.h LocBlk.h \
	 neighbors.cpp neighbors.h
reduce.o: $(PDBHOME)/pdb++.h $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.cpp $(TOOLHOME)/utility.h AtomConn.h \
	 AtomPositions.h BumperPoint.h CTab.h CliqueList.h DotSph.h \
	 ElementInfo.h LocBlk.h Mover.h PDBrec.h ResBlk.h StdResH.h \
	 neighbors.cpp neighbors.h reduce.cpp
# DO NOT DELETE THIS 2nd LINE -- make depend uses it
