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
DisjointSets.o: DisjointSets.cpp DisjointSets.h
Hash.o: Hash.cpp Hash.h
Point3d.o: Point3d.cpp Point3d.h
UseCount.o: UseCount.cpp UseCount.h
column_sprintf.o: column_sprintf.cpp
column_sscanf.o: column_sscanf.cpp
utility.o: utility.cpp
# DO NOT DELETE THIS 2nd LINE -- make depend uses it
