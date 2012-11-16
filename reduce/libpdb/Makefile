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
CFLAGS		= $(OPT) $(DEBUG) -DBOOLPREDEFINED -DHASSSCANFEXTERN

LIBRARY		= pdb++

LIBARCH		= lib$(LIBRARY).a
OBJS		= pdb_read.o pdb_sprntf.o pdb_sscanf.o pdb_chars.o \
		  pdb_type.o pdb++.o pdbinput.o
SRCS		= pdb_read.cpp pdb_sprntf.cppc pdb_sscanf.cpp pdb_chars.cpp \
		  pdb_type.cpp pdb++.cppc pdbinput.cpp

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
pdb++.o: pdb++.cpp pdb++.h
pdb_chars.o: pdb++.h pdb_chars.cpp pdbrun5_write.i pdbrun6_write.i \
	 write_format.i
pdb_read.o: pdb++.h pdb_read.cpp pdbrun5_read.i pdbrun6_read.i \
	 read_format.i
pdb_sprntf.o: pdb++.h pdb_sprntf.cpp
pdb_sscanf.o: pdb++.h pdb_sscanf.cpp
pdb_type.o: pdb++.h pdb_type.cpp
pdbinput.o: pdb++.h pdbinput.cpp
# DO NOT DELETE THIS 2nd LINE -- make depend uses it
