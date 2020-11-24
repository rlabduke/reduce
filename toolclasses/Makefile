#
#	On BSD machines, RANLIB should be 'ranlib'
#
#	On System V machines, RANLIB should be ':'
#
SHELL		= /bin/sh
RANLIB		= :

CXX		= g++

.SUFFIXES:	.cpp

.cpp.o:
	$(CXX) $(CFLAGS) -c $< -o $@

OPT		= -O
DEBUG		= $(CXXDEBUGFLAGS)
CFLAGS		= $(OPT) $(DEBUG) \
                   -DCHARFUNCMACROS -DINCTEMPLATEDEFNS \
		   -DBOOLPREDEFINED -DBRACKETOPERPARMS -DLISTFRIENDFIX

LIBRARY		= toolclasses

LIBARCH		= lib$(LIBRARY).a

OBJS		= UseCount.o Point3d.o Hash.o utility.o \
		  column_sprintf.o column_sscanf.o DisjointSets.o

all:		$(LIBARCH)

install:	$(LIBARCH)
		@echo ""

$(LIBARCH):     $(OBJS)
		@echo "Loading $(LIBARCH) ... "
		@ar cru $(LIBARCH) $(OBJS)
		@$(RANLIB) $(LIBARCH)
		@echo "done"

clean:;		@rm -f $(OBJS)

spotless:;	@rm -f $(OBJS) $(LIBARCH)

# DO NOT DELETE THIS LINE -- make depend uses it
Dict.o: Dict.cpp Dict.h Hash.h utility.h
DisjointSets.o: DisjointSets.cpp DisjointSets.h
Hash.o: Hash.cpp Hash.h
Hdl.o: Hdl.cpp Hdl.h UseCount.h
List.o: List.cpp List.h
ListIter.o: Iter.h List.h ListIter.cpp ListIter.h
MultiDict.o: Hash.h MultiDict.cpp MultiDict.h Seq.h utility.h
Point3d.o: Point3d.cpp Point3d.h
Seq.o: Seq.cpp Seq.h
UseCount.o: UseCount.cpp UseCount.h
Vector.o: Vector.cpp Vector.h utility.h
column_sprintf.o: column_sprintf.cpp
column_sscanf.o: column_sscanf.cpp
utility.o: utility.cpp
# DO NOT DELETE THIS 2nd LINE -- make depend uses it
