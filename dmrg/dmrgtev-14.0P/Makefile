############################################################################
#	Makefile			(C) 2009 Fabio Ortolani	fbo 090801 #
############################################################################
#
#	Name and version.
#		It is recomended that $(NAME)-$(VERSION) is the name of 
#		actual directory.
#
NAME 		:= dmrgtev
VERSION		:= 10.1
#
############################################################################
#
#	Copyright string
#
COPYRIGHT 	:= (C) 2011 Fabio Ortolani
#
#	Host, operating system and machine types.
#		(These can be hardcoded inside the program including header
#		version.h or version.hh, see below)
#
HOSTTYPE 	:= $(shell echo $$HOSTTYPE)
OSTYPE 		:= $(shell echo $$OSTYPE)
MACHTYPE 	:= $(shell echo $$MACHTYPE)
#
#	Executable and archive name
#
EXECUTABLE	:= $(NAME)
ARCHIVE		:= $(NAME)-$(VERSION).tbz
#
#	Sources
#
SOURCES 	:= $(wildcard *.h *.hh *.c *.cc)
SOURCES 	:= $(filter-out version.h version.hh, $(SOURCES))
#
#	Reordering of sources (first include files, then implementation files)
#
SOURCES 	:= $(patsubst %.c,  %.z,  $(SOURCES))
SOURCES 	:= $(patsubst %.cc, %.zz, $(SOURCES))
SOURCES 	:= $(sort $(SOURCES))
SOURCES 	:= $(patsubst %.z,  %.c,  $(SOURCES))
SOURCES 	:= $(patsubst %.zz, %.cc, $(SOURCES))
#
#	Completing list of sources
#
SOURCES 	:= Makefile $(SOURCES)
SOURCES 	:= $(SOURCES) $(wildcard *.inp)
SOURCES 	:= $(SOURCES) $(wildcard *.sh)
#
#	C, C++ units
#
CSOURCES 	:= $(filter %.c,  $(SOURCES))
CXXSOURCES	:= $(filter %.cc, $(SOURCES))
#
#	Compilation objects
#
COBJECTS	:= $(patsubst %.c,  %.o,   $(CSOURCES))       
CXXOBJECTS 	:= $(patsubst %.cc, %.o, $(CXXSOURCES))
OBJECTS 	:= $(sort $(COBJECTS) $(CXXOBJECTS))
#
#	Add path to sources (for archive creation)
#
PATHSOURCES	:= $(patsubst %, $(shell basename `pwd`)/%, $(SOURCES)) 
#
#	Targets, not sources files and intermediate files
#
TARGETS  	 = $(wildcard $(EXECUTABLE) *.out)
NOTSOURCES	 = $(filter-out $(SOURCES), $(wildcard *))
GARBAGE  	 = $(filter-out $(TARGETS) $(ARCHIVE), $(NOTSOURCES))
#
#	Actual date
#
NOW		 = $(shell date '+%y%m%d')
#
############################################################################
#
#	Compilation control variables 
#
CPP		:= gcc
CC 		:= gcc
LD 		:= gcc
#
#CPP		:= /Users/naldesi/Downloads/usr/local/bin/gcc
#CC 		:= /Users/naldesi/Downloads/usr/local/bin/gcc
#LD 		:= /Users/naldesi/Downloads/usr/local/bin/gcc
#
#	Options
#
DFLAGS		:= -E -MM
CPPFLAGS	:=
CFLAGS 		:= -O3 -Wall -Winline -mfpmath=sse -msse3
CXXFLAGS 	:= -O3 -Wall -Winline -mfpmath=sse -msse3
LDFLAGS		:= 
LOADLIBES	:= -lstdc++ -lpthread -lm 
#

#CFLAGS 	:= -g -faddress-sanitizer
#CXXFLAGS 	:= -g -faddress-sanitizer
#CFLAGS 		:= -g -fsanitize=address
#CXXFLAGS 	:= -g -fsanitize=address

ifeq  ($(strip $(OSTYPE)),darwin9.0)
#
#	Mac version of gcc often defaults to 32 bit architecture
#
CFLAGS		:= -m64 $(CFLAGS)
CXXFLAGS	:= -m64 $(CXXFLAGS)
LDFLAGS		:= -m64 $(LDFLAGS)
#
endif
#
ifeq  ($(strip $(HOSTTYPE)),powerpc)
#
#	on Mac powerpc we must add explicitly this library (I think this is 
#	a bug of Mac version of gcc)
#
LOADLIBES 	:= -lgcc_s.10.4 $(LOADLIBES)
endif
#
############################################################################
#
#	Rules
#
.SUFFIXES:
.SUFFIXES:	.c .cc .o .inp .out .tgz .tbz .pdf
#
#	Compilation rules
#
%.o:	%.c
	$(CC) -c $(CPPFLAGS) $(CFLAGS)   $< -o $@
%.o:  	%.cc
	$(CC) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@
#
#	Execution rules
#		Generation of output files (file.out or standard output) from 
#		input files (file.inp)
#
%.out:	%.inp $(EXECUTABLE)
	./$(EXECUTABLE) < $< >& $@ &
%.pdf:  %.inp $(EXECUTABLE)
	@echo "Generating $@ from $<" 
	./$(EXECUTABLE) <$< 2>&1 | a2ps -M a4 -o - 2>/dev/null | \
	ps2pdf -sPAPERSIZE=a4 - $@  &
%:	%.inp $(EXECUTABLE)
	./$(EXECUTABLE) < $< 2>&1 
#
#	Archiving rules (tar gzipped or tar bzipped)
#	Archives are created in upper directory. 
#
%.tgz:	$(SOURCES)
	@echo "Creating archive $@ containing:"
	@tar -cvf - -C ../ $(PATHSOURCES) | gzip -9 > $@
%.tbz:	$(SOURCES)
	@echo "Creating archive $@ containing:"
	@tar -cvf - -C ../ $(PATHSOURCES) | bzip2   > $@
#
############################################################################
#
#	Phony targets
#
.PHONY:		all clean cleanall distclean save restore list
#
#	Default target (the executable)
#
all:		$(EXECUTABLE)
#
#	Clean all intermediate files
#
clean:
	rm -f $(GARBAGE)
#
#	Full cleaning (only source files are preserved)
#
cleanall:
	rm -f $(NOTSOURCES)
distclean:	cleanall
#
#	Save and restore (to and from upper directory)
#
save:		../$(ARCHIVE)
restore:	../$(ARCHIVE)
	tar xjpvf ../$(ARCHIVE) -C ../
#
#	Create a listing of all sources
#
list:		$(NAME).pdf
#
##############################################################################
#
#	Explicite dependencies and rules.
#
#	Listing creation
#
$(NAME).pdf:	$(SOURCES)
#	@a2ps -M a4 -o - $^ 2>/dev/null | ps2pdf -sPAPERSIZE=a4 - $@
	@a2ps --pro=color -C -2 -M a4 -g --pretty-print='C++' -o - $^ 2>/dev/null | ps2pdf -sPAPERSIZE=a4 - $@
	@printf "Generated sources listing (pdf format): $@"
	@echo "$(shell wc -l $^ | tail -n 1 | sed -e 's/total/lines/')"
#
#	Executable creation
#
$(EXECUTABLE):	$(OBJECTS)
	$(LD) -o $@ $(LDFLAGS) $^ $(LOADLIBES)
#
############################################################################
#
#	System dependent include files (can be included to hardcode
#	system informations).
#
#		C version include file
#
version.h:	Makefile
	@printf "%-8s%-32s%-24s%-10s\n" 			\
		"/*" "version.h" "$(COPYRIGHT)" "fbo $(NOW)"	>  $@
	@printf "        ============================="		>> $@
	@printf "=====================================  */\n"	>> $@
	@echo "#ifndef VERSION_H"				>> $@
	@echo "#define VERSION_H"				>> $@
	@echo "/*  Automatic generated file. Don't edit!  */"	>> $@
	@echo "#define NAME      \"$(NAME)\""			>> $@
	@echo "#define VERSION   \"$(VERSION)\""		>> $@
	@echo "#define COPYRIGHT \"$(COPYRIGHT) fbo $(NOW)\""	>> $@
	@echo "#define HOSTTYPE  \"$(HOSTTYPE)\""		>> $@
	@echo "#define OSTYPE    \"$(OSTYPE)\""			>> $@
	@echo "#define MACHTYPE  \"$(MACHTYPE)\""		>> $@
	@echo "#endif /* VERSION_H */"				>> $@
	@printf "/* ==================================="	>> $@
	@printf "====================================  */\n"	>> $@
#
#		C++ version include file 
#
version.hh:	Makefile
	@printf "%-8s%-32s%-24s%-10s\n" 			\
		"//" "version.hh" "$(COPYRIGHT)" "fbo $(NOW)"	>  $@
	@printf "//      ============================="		>> $@
	@printf "=====================================\n"	>> $@
	@echo "#ifndef VERSION_HH"				>> $@
	@echo "#define VERSION_HH"				>> $@
	@echo "//      Automatic generated file. Don't edit!"	>> $@
	@echo "#define NAME      \"$(NAME)\""			>> $@
	@echo "#define VERSION   \"$(VERSION)\""		>> $@
	@echo "#define COPYRIGHT \"$(COPYRIGHT) fbo $(NOW)\""	>> $@
	@echo "#define HOSTTYPE  \"$(HOSTTYPE)\""		>> $@
	@echo "#define OSTYPE    \"$(OSTYPE)\""			>> $@
	@echo "#define MACHTYPE  \"$(MACHTYPE)\""		>> $@
	@echo "#endif // VERSION_HH"				>> $@
	@printf "//==================================="		>> $@
	@printf "=====================================\n"	>> $@
#
############################################################################
#
#	Dependency list creation and inclusion
#
depend:	version.h version.hh $(CSOURCES) $(CXXSOURCES)	
	@$(CPP) $(DFLAGS) $(CPPFLAGS) $(CSOURCES) $(CXXSOURCES) > $@
-include	depend
#
############################################################################
