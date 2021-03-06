#----------------------------------------------------------------------
#
#
# 980303: Bottom-level makefile for a subproject.
#	  If the amount of comments is exasperating, do 'make undocumented'.
#	  User maintainable parts are placed near the top. Readability
#	  decreases towards the end.
#
#	  Hardcoded stuff is limited to:
#
#	    libraries, paths for libraries and header files,
#	    compiler options, directory names.
#
#	  Make invocation:
#
#	    make [MODE=opt|MODE=nopt]    \
#                [LIBMODE=true]          \
#		 [PROFILING=true]        \
#                target
#
#	    (Defaults are: MODE=nopt, LIBMODE=true, PROFILING=false.)
#	    (Note that a user changing MODE becomes responsible for
#	     making [real]clean.)
#	    (LIBMODE=true ensures that executables are linked with
#	     libraries, not object files.)
#
#	  Some of the explicit targets:
#
#	    lib:		Collects all object code into a library.
#	    all:		Makes applications and libs.
#	    clean:		Removes objects and libs for specific platform.
#	    realclean:		Removes all object code, libs and emacs
#				backups.
#	    depend:		Updates all dependency files.
#	    app, apps		Forces a remaking of all applications. Useful
#				when changing platform and one does not want
#				to make clean.
#
#	  Bugs and features:
#
#	    - Will not work unless there is at least one source file.
#	      Still true?
#	    - Will not work with applications using the same name
#	      and source in different languages.
#	    - When switching platform/optimization level, may have to
#	      rerun make, due to a small problem with make reevaluating
#	      and updating the makefiles.
#	    - As of now, applications are not separated according to
#	      $(HOSTTYPE), even if the object files and libraries are.
#	    - Only one library for each sub-project directory is supported.
#	    - 980402: Added support for GLOptimizer stuff.
#	    - 980515: Mesa tested and found to work.
#	    - 980521: Added '-p' to mkdir to make '-j 2' possible.
#	    - 9837:   Initiated removal of Cosmo/Optimizer from Makefiles.
#	    - 980914: Trying to merge yvh/hke's PC-changes into master.
#           - 20050222: So many changes have been made to the system
#                       that the comments above may be out of date... 
#
#
#
#
#
#----------------------------------------------------------------------


#------------------------------------------------------------
#
#	afr 11/01 2000: I do not understand what controls
#	the HOSTTYPE variable on Linux, so instead of using
#	it we check for linux and set it to 'linux'.
#
#------------------------------------------------------------

UNAME=$(shell uname)
ifeq "$(UNAME)" "Linux"
  PLATFORM=linux
else
  ifeq "$(UNAME)" "HP-UX"
    PLATFORM=hp-pa
  else
    ifeq "$(UNAME)" "IRIX64"
      PLATFORM=sgi
    else
      PLATFORM=$(HOSTTYPE)
    endif
  endif
endif



#------------------------------------------------------------
#
#	MODE will be overridden by the MODE input parameter.
#	Same for LIBMODE which is explained at the top.
#       ... and for PROFILING.
#
#------------------------------------------------------------

MODE=nopt
PROFILING=false
LIBMODE=true



#------------------------------------------------------------
#
#	C specifics.
#
#------------------------------------------------------------

CC			="A_C_compiler_for_use_on_$(PLATFORM)"



ifeq "$(PLATFORM)" "winnt"

  CC			=cl
  LD			=link

  ifeq "$(MODE)" "opt"
    CFLAGS		=-GX -G6 -GA -Gs -Gf -Gy -Ox -Ob2 -nologo -MD
    CDEFS		=-DWIN32 -DMICROSOFT -DNDEBUG
    LDFLAGS		=-nologo -opt:ref

  else
    CFLAGS		=-GX -nologo -Yd -Z7 -MDd
    CDEFS		=-DWIN32 -DMICROSOFT
    LDFLAGS		=-nologo

  endif
endif

ifeq "$(PLATFORM)" "borland"

  CC			=bcc32
  LD			=ilink32

  ifeq "$(MODE)" "opt"
    CFLAGS		=-6 -O2
    CDEFS		=-DWIN32 -DNDEBUG
    LDFLAGS		=

  else
    CFLAGS		=-6 -v
    CDEFS		=-DWIN32
    LDFLAGS		=

  endif
endif


ifeq "$(PLATFORM)" "hp-pa"

  CC			=/usr/bin/cc
  LD			=$(CC)
#	LD           = /SI/Programs/purify/purify-5.1-hpux/purify CC

  ifeq "$(MODE)" "opt"
    CFLAGS		=+DAportable +Oall
    CDEFS		=-DHPUX -DNDEBUG
  else
    CFLAGS		=+DAportable -g
    CDEFS		=-DHPUX
  endif

  LDFLAGS		=+DAportable
endif



ifeq "$(PLATFORM)" "sgi"

  CC			=cc
  LD			=$(CC)

  ifeq "$(MODE)" "opt"
    CFLAGS		=-O3 -n32
    CDEFS		=-DSGI -DNDEBUG
  else
    CFLAGS		=-g -n32
    CDEFS		=-DSGI
  endif
endif


ifeq "$(PLATFORM)" "linux"

  CC		=gcc
  LD		=$(CC)

  ifeq "$(MODE)" "opt"
    CFLAGS	=-O2 -Wall
    CDEFS	=-DNDEBUG
  else
    CFLAGS	=-g -ggdb -Wall -pthread -fPIC
    CDEFS	=
  endif

endif



#------------------------------------------------------------
#
#	C++ specifics.
#
#------------------------------------------------------------

CXX			="A_C++_compiler_for_use_on_$(PLATFORM)"



ifeq "$(PLATFORM)" "winnt"

  CXX			=cl
  LDXX			=link

  ifeq "$(MODE)" "opt"
    CXXFLAGS		=-TP -GX -GR -G6 -GA -Gs -Gf -Gy -Ox -Ob2 -nologo -MD
    CXXDEFS		=-DWIN32 -DMICROSOFT -DNDEBUG
    LDXXFLAGS		=-nologo -opt:ref

  else
    CXXFLAGS		=-TP -GX -GR -nologo -Yd -Z7 -MDd
    CXXDEFS		=-DWIN32 -DMICROSOFT
    LDXXFLAGS		=-nologo -DEBUGTYPE:BOTH 

  endif
endif


ifeq "$(PLATFORM)" "borland"

  CXX			=bcc32
  LDXX			=ilink32

  ifeq "$(MODE)" "opt"
    CXXFLAGS		=-6 -P -O2
    CXXDEFS		=-DWIN32 -DNDEBUG
    LDXXFLAGS		=

  else
    CXXFLAGS		=-6 -P -v
    CXXDEFS		=-DWIN32
    LDXXFLAGS		=

  endif
endif


ifeq "$(PLATFORM)" "hp-pa"

  CXX			=/opt/aCC/bin/aCC
  LDXX		=$(CXX)
#  LDXX		=/SI/Programs/purify/purify-5.1-hpux64/purify $(CXX)

  ifeq "$(MODE)" "opt"
    CXXFLAGS		=-O -AA
    CXXDEFS		=-DHPUX -DHP_ACC_Cplusplus -DHP_Cplusplus -DNDEBUG
  else
    CXXFLAGS		=-g -AA
    CXXDEFS		=-DHPUX -DHP_ACC_Cplusplus -DHP_Cplusplus
  endif

  LDXXFLAGS		=-g -z -G +DAportable -AA
endif


ifeq "$(PLATFORM)" "sgi"

  CXX		=CC
  LDXX		=$(CXX)

  ifeq "$(MODE)" "opt"
    CXXFLAGS	=-O3 -n32 -LANG:std -ptused
    CXXDEFS	=-DSGI -DNDEBUG
  else
    CXXFLAGS	=-g -n32 -LANG:std -ptused
    CXXDEFS	=-DSGI
  endif

  LDXXFLAGS		=-LANG:std -ptused
endif


ifeq "$(PLATFORM)" "linux"

    CXX		=g++
    LDXX	=$(CXX)

  ifeq "$(MODE)" "opt"
    CXXFLAGS	=-O2 -Wall
    CXXDEFS	=-DNDEBUG
  else
    ifeq "$(MODE)" "maxopt"
      CXXFLAGS	=-O3 -ffast-math -Wall
      CXXDEFS	=-DNDEBUG
    else
      CXXFLAGS	=-g -ggdb -Wall -pthread -fPIC
      CXXDEFS	=
    endif
  endif
endif



ifeq "$(PROFILING)" "true"
  CFLAGS +=-pg
  CXXFLAGS +=-pg
  LDFLAGS +=-pg
  LDXXFLAGS +=-pg
endif



#------------------------------------------------------------
#
#	Libraries to include by default.
#
#------------------------------------------------------------

# afr, 10/5 1999:
# Changed this section. Now, DEPLIBSLOCAL should contain
# a list of local libraries to link (local meaning that it
# is a subproject on the same level as the current project).
# Examples include GoGeometry, utils, brep etc.
# DEPLIBSGLOBAL should contain the non-local non-system libraries, such
# as siscat and basictools.
# DEPLIBSSYSTEM should be system and OpenGL libraries.
# All three variables will be read from Makefile.local, if it exist.
# Any variables not set will be set in this file instead, to a
# default value.

DEPLIBSSYSTEM     = m pthread X11 cln

SYSTEMINCLUDEDIRS = /usr/X11R6/include

# SYSTEMLIBDIRS     = /usr/X11R6/lib64 
SYSTEMLIBDIRS     = /usr/X11R6/lib

# Here follows default values. These should be made platform dependant.

ifndef DEPLIBSLOCAL
DEPLIBSLOCAL	=
endif

ifndef DEPLIBSSYSTEM
DEPLIBSSYSTEM =
endif

ifndef SYSTEMINCLUDEDIRS
SYSTEMINCLUDEDIRS =
endif

ifndef SYSTEMLIBDIRS
SYSTEMLIBDIRS     =
endif


CLIBS		= $(foreach lib, $(DEPLIBSLOCAL),-l$(lib)) $(foreach lib, $(DEPLIBSSYSTEM),-l$(lib))



CXXLIBS		= $(foreach lib, $(DEPLIBSLOCAL),-l$(lib)) $(foreach lib, $(DEPLIBSSYSTEM),-l$(lib))



#------------------------------------------------------------
#
#	Include file locations.
#	Note: '-I.' is not the src-directory, but the one
#	      above, i.e. the one containing the makefile.
#
#------------------------------------------------------------

CINCLUDES	+=-Iinclude

CXXINCLUDES +=-Iinclude $(foreach inc, $(DEPLIBSLOCAL),-I../$(inc)/include) $(foreach inc, $(SYSTEMINCLUDEDIRS),-I$(inc))


#------------------------------------------------------------
#
#	Library locations.
#	980629: Hmm... Some winnt stuff seems to be
#	      	hardcoded further down... Look for oglsdk.
#
#------------------------------------------------------------

# For some reason, we use hp-pa as PLATFORM, but siscat et al use hp-aCC, so
# the PLATFORM for non-local, Sintef libraries (="global") is set especially
# for hp.
#
# Same applies to sgi. BasicTools etc. use IRIX for -n32, sgi for -o32 and
# IRIX64 for -64 compiled objects.
#

ifeq "$(PLATFORM)" "hp-pa"
  GLOBALPLATFORM = hp-aCC
else
  ifeq "$(PLATFORM)" "sgi"
    GLOBALPLATFORM = IRIX
  else
    GLOBALPLATFORM = $(PLATFORM)
  endif
endif


ifeq "$(PLATFORM)" "winnt"
  CXXLIBPATH	=-LIBPATH:lib/$(PLATFORM)/$(MODE) $(foreach lib, $(DEPLIBSLOCAL),-LIBPATH:../$(lib)/lib/$(PLATFORM)/$(MODE)) $(foreach lib, $(SYSTEMLIBDIRS),-LIBPATH:$(lib))
else
  CXXLIBPATH	=-Llib/$(PLATFORM)/$(MODE) $(foreach lib, $(DEPLIBSLOCAL),-L../$(lib)/lib/$(PLATFORM)/$(MODE)) $(foreach lib, $(SYSTEMLIBDIRS),-L$(lib))
endif

ifeq "$(PLATFORM)" "winnt"
  CLIBPATH	=-LIBPATH:lib/$(PLATFORM)/$(MODE) $(foreach lib, $(DEPLIBSLOCAL),-LIBPATH:../$(lib)/lib/$(PLATFORM)/$(MODE)) $(foreach lib, $(SYSTEMLIBDIRS),-LIBPATH:$(lib))
else
  CLIBPATH	=-Llib/$(PLATFORM)/$(MODE) $(foreach lib, $(DEPLIBSLOCAL),-L../$(lib)/lib/$(PLATFORM)/$(MODE)) $(foreach lib, $(SYSTEMLIBDIRS),-L$(lib))
endif


#------------------------------------------------------------
#
#	Auxiliary variables.
#	Concatenation of options.
#	Wildcard expanded and derived variables.
#	(Note to myself: Why did I use simply expanded variables here?)
#
#------------------------------------------------------------

COPTS		=$(CFLAGS) $(CINCLUDES) $(CDEFS)
CXXOPTS		=$(CXXFLAGS) $(CXXINCLUDES) $(CXXDEFS)
LDOPTS		=$(LDFLAGS) $(CLIBPATH)
LDXXOPTS	=$(LDXXFLAGS) $(CXXLIBPATH)



CPROGS		=$(basename $(notdir $(wildcard app/*.c)))
#CXXPROGS	=$(basename $(notdir $(wildcard app/*.C)))
CXXPROGS	=$(basename $(notdir $(wildcard app/*.C) $(wildcard app/*.cpp)))
CXXPROGSRCS	=$(notdir $(wildcard app/*.C) $(wildcard app/*.cpp))

LIBNAME		=$(notdir $(shell pwd))

CSRCS		:=$(wildcard src/*.c) $(addsuffix .c, $(CPROGS))
#CXXSRCS		:=$(wildcard src/*.C) $(addsuffix .C, $(CXXPROGS))
CXXSRCS		:=$(wildcard src/*.C) $(wildcard src/*.cpp) $(CXXPROGSRCS)
CXXSRCSBASE	:=$(notdir $(basename $(CXXSRCS)))
QFORMBASE	=$(basename $(notdir $(wildcard include/*.ui)))
QNOFORMBASE	=$(basename $(notdir $(shell grep -ls Q_OBJECT include/*.h | grep -v _form)))

COBJS		=$(addprefix lib/$(PLATFORM)/$(MODE)/, \
			     $(notdir $(CSRCS:.c=.o)))
#CXXOBJS		=$(addprefix lib/$(PLATFORM)/$(MODE)/, \
#			     $(notdir $(CXXSRCS:.C=.o)))
CXXOBJS		=$(addprefix lib/$(PLATFORM)/$(MODE)/, \
			     $(addsuffix .o, $(CXXSRCSBASE)))
QFORMOBJS	=$(addprefix lib/$(PLATFORM)/$(MODE)/, \
			     $(addsuffix .o, $(QFORMBASE))) 
QNOFORMOBJS	=$(addprefix lib/$(PLATFORM)/$(MODE)/, \
			     $(addsuffix .o, $(QFORMBASE))) 
QFORMMOCOBJS	=$(addprefix lib/$(PLATFORM)/$(MODE)/moc_, \
			     $(addsuffix .o, $(QFORMBASE))) 
QNOFORMMOCOBJS	=$(addprefix lib/$(PLATFORM)/$(MODE)/moc_, \
			     $(addsuffix .o, $(QNOFORMBASE))) 
ALLQOBJS	=$(QFORMOBJS) $(QNOFORMOBJS) $(QFORMMOCOBJS) $(QNOFORMMOCOBJS)


CDEPFILES	= $(addprefix dep/, $(addsuffix .d, $(notdir $(CSRCS:.c=))))
#CXXDEPFILES	= $(addprefix dep/, $(addsuffix .d, $(notdir $(CXXSRCS:.C=))))
CXXDEPFILES	= $(addprefix dep/, $(addsuffix .d, $(CXXSRCSBASE)))

CPROGOBJS	=$(addprefix lib/$(PLATFORM)/$(MODE)/, \
			     $(addsuffix .o, $(CPROGS)))
CXXPROGOBJS	=$(addprefix lib/$(PLATFORM)/$(MODE)/, \
			     $(addsuffix .o, $(CXXPROGS)))

NOAPPCOBJS	=$(filter-out $(CPROGOBJS), $(COBJS))
NOAPPCXXOBJS	=$(sort $(filter-out $(CXXPROGOBJS), $(CXXOBJS)) $(ALLQOBJS))
NOAPPNOQCXXOBJS	=$(filter-out $(ALLQOBJS), $(NOAPPCXXOBJS))


HOSTMODELIB	=$(PLATFORM)/$(MODE)/lib$(LIBNAME)




#------------------------------------------------------------
#
#	Targets. Let 'default' be the first one.
#
#------------------------------------------------------------

.PHONY: 	all clean realclean depend varcheck Makefile lib appreport \
		libreport update_from_master update_the_master apps


default:	
		@echo No default target supported. Please make something.


#
# This is for reporting application and lib names to the makefile on the
# level above, it needs to be aware of them for copying stuff around...
#
appreport:
		@echo $(CXXPROGS) $(CPROGS)


libreport:
		@echo lib/$(HOSTMODELIB).a


undocumented:	Makefile
		cat Makefile | egrep -v '^#'


all:		$(addprefix app/, $(CPROGS) $(CXXPROGS)) \
		lib/$(HOSTMODELIB).a


varcheck:	
		@echo ""
		@echo "USER=" $(USER)
		@echo ""
		@echo "PWD=" $(PWD)
		@echo "(notdir (PWD))=" $(notdir $(PWD))
		@echo ""
		@echo "CXXDEPFILES="$(CXXDEPFILES)
		@echo "DEPLIBS="$(DEPLIBS)
		@echo "DEPLIBOBJS="$(DEPLIBOBJS)
		@echo ""
		@echo "QQQ="$(QQQ)
		@echo "QQQ="$(word 1, $(QQQ))
		@echo ""
		@echo "HOST="$(HOST)
		@echo "HOSTNAME="$(HOSTNAME)
		@echo "PLATFORM="$(PLATFORM)
		@echo "LIBMODE="$(LIBMODE)
ifeq "$(LIBMODE)" "true"
		@echo "LIBMODE=true er true"
else
		@echo "LIBMODE=true er false"
endif
		@echo ""
		@echo "BONES_MODE="$(BONES_MODE)
		@echo ""
		@echo "CLIBPATH="$(CLIBPATH)
		@echo "CXXLIBPATH="$(CXXLIBPATH)
		@echo ""
		@echo "CXXLIBS="$(CXXLIBS)
		@echo "CXXDEFS="$(CXXDEFS)
		@echo "CXXOPTS="$(CXXOPTS)
		@echo ""
		@echo "LIBNAME="$(LIBNAME)
		@echo ""
		@echo "MAKE=" $(MAKE)
		@echo "MFLAGS=" $(MFLAGS)
		@echo "MAKEFLAGS=" $(MAKEFLAGS)
		@echo "UNAME="$(UNAME)
		@echo "PLATFORM="$(PLATFORM)
		@echo "HOSTTYPE="$(HOSTTYPE)
		@echo ""
		@echo "NOAPPCXXOBJS="$(NOAPPCXXOBJS)
		@echo "CXXOBJS="$(CXXOBJS)
		@echo "CXXPROGOBJS="$(CXXPROGOBJS)
		@echo "CXXSRCSBASE="$(CXXSRCSBASE)
		@echo "QFORMBASE="$(QFORMBASE)
		@echo "QNOFORMBASE="$(QNOFORMBASE)
		@echo "ALLQOBJS="$(ALLQOBJS)
		@echo "NOAPPNOQCXXOBJS="$(NOAPPNOQCXXOBJS)

clean:
		@echo ""
		@echo "Remember that MODE=[opt|nopt] *does* make a difference"
		@echo "for the 'clean' target."
		@echo ""
		@if [ ! -d lib ]; then mkdir -p lib; fi
		@if [ ! -d lib/$(PLATFORM) ]; then mkdir -p lib/$(PLATFORM); fi
		@if [ ! -d lib/$(PLATFORM)/$(MODE) ]; then \
		  mkdir -p lib/$(PLATFORM)/$(MODE); fi
		rm -f lib/$(PLATFORM)/$(MODE)/*.o lib/$(PLATFORM)/$(MODE)/*.a
		rm -rf lib/$(PLATFORM)/$(MODE)/ii_files core
		rm -f $(addprefix app/, $(CXXPROGS) $(CPROGS) core)
		rm -f lib/$(PLATFORM)/$(MODE)/ii_files
		rm -f src/moc_*.C
		rm -f src/*_form.C
		rm -f include/*_form.h


realclean:
		@echo ""
		@echo "This might remake dependencies. Nothing to worry about."
		@echo "Yes. It is unnecessary. Yes. It is annoying."
		@echo "No. I don't know how to avoid it, yet."
		@echo ""
		rm -rf $(addprefix app/, $(CXXPROGS) \
					 $(CPROGS) core)
		rm -f dep/*.d
		rm -rf lib core
		find . -name '*~' -exec rm -f {} \;
		mkdir -p dep
		mkdir -p lib
		rm -f lib/$(PLATFORM)/$(MODE)/ii_files
		rm -f src/moc_*.C
		rm -f src/*_form.C
		rm -f include/*_form.h


ifeq "$(LIBMODE)" "true"

  ifeq "$(PLATFORM)" "winnt"
    lib:	lib/$(HOSTMODELIB).lib
  else
    lib:	lib/$(HOSTMODELIB).a
  endif

else

  lib:		$(NOAPPCOBJS) $(NOAPPCXXOBJS)

endif


lib/$(HOSTMODELIB).a:	$(NOAPPCOBJS) $(NOAPPCXXOBJS)
		ar -rv lib/$(HOSTMODELIB).a $^

lib/$(HOSTMODELIB).lib:	$(NOAPPCOBJS) $(NOAPPCXXOBJS)
		link -lib $(LDXXOPTS) \
		/out:"lib/$(HOSTMODELIB).lib" $^


#
# 980304: Dependencies. To do it perfect is perhaps non-trivial?
#
#	  1) Executables. The object file of an executable is easy, but
#	     the executable itself might depend on libraries or
#	     object files not known to
#	     the compiler. (They might be known through headerfiles
#	     containing prototypes, but libraries might have been rebuilt
#	     without changing the headerfiles. The libraries depend on the
#	     headerfiles, not the other way around, therefore; a problem.)
#	     Solution: Make executables dependent on all libraries and all
#	     objects. This is not perfect.
#	  2) As of now, there is one set of dependency files for all
#	     platforms, so 'make depend' should be issued when moving to
#	     another platform, I think. Perhaps. Better: A set of dependency
#	     files for each platform? Not very important. Problems might
#	     appear, though. Ex.: #include statements wrapped in platform-
#	     dependent #ifdefs.
#	  3) If compilation speed is of concern, one might want to remove
#	     all dependencies on the dependency files. This will stop
#	     the dependency files from being (re)made, but may
#	     cause some hard-to-find problems. If wanted, comment out
#	     the 'Makefile' rule.
#	  4) Note that 'make depend' will always say
#	     "make: Nothing to be done for `depend'.", since the makefile
#	     itself depends on 'depend'. In fact, this makes it completely
#	     unnecessary to do 'make depend' explicitly.
#
# 980309: By not using dependecy on the file dep/hope*.dep, but rather
#	  including the rule to make it here, we forces this to be remade
#	  *always*, which is necessary since we might have changed opt-mode
#	  from the last time it was produced.
#	  Also, we avoid infinite loops, since 'depend' will be up to date.
# 980309: Aha. Not sure if this comment fits in here.
#	  Hmm. Seems that make does not detect that hopefully* was
#	  remade as a makefile, and therefore does not detect that Makefile
#	  should be "reconsidered". Possible to solve this by forcing
#	  make to regard hopefully* as outdated from the start?
# 980310: This is relevant to the problem of having to "make" twice, so I'll
#	  let the comment remain for a while.
#
depend:		$(CXXDEPFILES) $(CDEPFILES)
		@if [ ! -d dep ]; then mkdir -p dep; fi
		@rm -f dep/hopefully_a_unique_name.dep
		@( $(foreach app, hopefully_another_unique_name \
				  $(CXXPROGS) $(CPROGS), \
		     echo app/$(app): lib/$(PLATFORM)/$(MODE)/$(app).o; \
		     echo ""; ) ) > \
		   dep/hopefully_a_unique_name.dep


Makefile:	depend


#
# 980305: Wow! This works. Means that pattern matching is done not only on
#	  target, but also on dependencies! (Ok, '(make.info)Pattern Intro'
#	  explains this.)
# 980310: Cannot see that it is possible to combine these rules into one,
#	  in a neat way.
#	  We have two (three) languages and two possible directories for
#	  sources.
#
dep/%.d:	src/%.C
		@if [ ! -d dep ]; then mkdir -p dep; fi
		g++ -MM -MG $(CXXINCLUDES) $(CXXDEFS) -USGI -UMICROSOFT $< | \
		sed 's/$*\.o/lib\/$(PLATFORM)\/opt\/& \
			     lib\/$(PLATFORM)\/nopt\/& dep\/$(@F)/' | \
		sed 's/include\/\(\w*_form\.h\)/\1/' | \
		sed 's/\(\w*_form\.h\)/include\/\1/' > $@

dep/%.d:	src/%.cpp
		@if [ ! -d dep ]; then mkdir -p dep; fi
		g++ -MM -MG $(CXXINCLUDES) $(CXXDEFS) -USGI -UMICROSOFT $< | \
		sed 's/$*\.o/lib\/$(PLATFORM)\/opt\/& \
			     lib\/$(PLATFORM)\/nopt\/& dep\/$(@F)/' > $@

dep/%.d:	app/%.C
		@if [ ! -d dep ]; then mkdir -p dep; fi
		g++ -MM -MG $(CXXINCLUDES) $(CXXDEFS) -USGI -UMICROSOFT $< | \
		sed 's/$*\.o/lib\/$(PLATFORM)\/opt\/& \
			     lib\/$(PLATFORM)\/nopt\/& dep\/$(@F)/' > $@

dep/%.d:	app/%.cpp
		@if [ ! -d dep ]; then mkdir -p dep; fi
		g++ -MM -MG $(CXXINCLUDES) $(CXXDEFS) -USGI -UMICROSOFT $< | \
		sed 's/$*\.o/lib\/$(PLATFORM)\/opt\/& \
			     lib\/$(PLATFORM)\/nopt\/& dep\/$(@F)/' > $@

dep/%.d:	src/%.c
		@if [ ! -d dep ]; then mkdir -p dep; fi
		gcc -MM -MG $(CINCLUDES) $(CDEFS) -UMICROSOFT $< | \
		sed 's/$*\.o/lib\/$(PLATFORM)\/opt\/& \
			     lib\/$(PLATFORM)\/nopt\/& dep\/$(@F)/' > $@

dep/%.d:	app/%.c
		@if [ ! -d dep ]; then mkdir -p dep; fi
		gcc -MM -MG $(CINCLUDES) $(CDEFS) -UMICROSOFT $< | \
		sed 's/$*\.o/lib\/$(PLATFORM)\/opt\/& \
			     lib\/$(PLATFORM)\/nopt\/& dep\/$(@F)/' > $@


#
# 980306: Ending .dep is chosen so as not to confuse with the other .d files.
# 980309: This is never used now, see comments to 'depend', but we can
#	  leave it here for later use. Mighty strange. This rule shouldn't
#	  be used, still get into (even more) problems when commented out.
#	  Ok, seems that the include-statement triggers a search for the
#	  rule. But still has to force a remake of it in 'depend'.
#
dep/hopefully_a_unique_name.dep:
		@rm -f dep/hopefully_a_unique_name.dep
		@( $(foreach app, hopefully_another_unique_name \
				  $(CXXPROGS) $(CPROGS), \
		     echo $(app): lib/$(PLATFORM)/$(MODE)/$(app).o; \
		     echo ""; ) ) > \
		   dep/hopefully_a_unique_name.dep


#------------------------------------------------------------
#
#	A Borland trick
#
#------------------------------------------------------------

ifeq "$(PLATFORM)" "borland"
  SPACEIFNOBORLAND= 
else
  SPACEIFNOBORLAND=
endif

#
# 980304: These are the (pattern) rules for making executables.
#	  Note! The targets have to be indented by spaces, not tabs!
# 980402: Note that for LIBMODE=false, we should probably make things
#	  dependent on the NOAPPOBJS... This was most likely a bug.
# 980609: In case of LIBMODE=true, dependency on libraries other than the
#	  current library would be nice!
#	  (Perhaps not desirable under all circumstances?)
#
ifeq "$(LIBMODE)" "true"

  ifeq "$(PLATFORM)" "winnt"

    $(addprefix app/, $(CXXPROGS)):	lib/$(HOSTMODELIB).lib \
					$(DEPLIBS)
	$(LDXX) $(LDXXOPTS) -out:$@ lib/$(PLATFORM)/$(MODE)/$(@F).o \
	        lib/$(HOSTMODELIB).lib \
		$(addsuffix .lib, $(addprefix lib, $(DEPLIBSLOCAL))) \
		$(addsuffix .lib, $(DEPLIBSGLOBAL)) \
		$(DEPLIBSSYSTEM)
	@cp -p $@ $@.exe


  else

    $(addprefix app/, $(CXXPROGS)):	lib/$(HOSTMODELIB).a $(DEPLIBS)
	$(LDXX) $(LDXXOPTS) -o$(SPACEIFNOBORLAND)$@ lib/$(PLATFORM)/$(MODE)/$(@F).o \
	        -l$(LIBNAME) $(CXXLIBS)

    $(addprefix app/, $(CPROGS)):	lib/$(HOSTMODELIB).a $(DEPLIBS)
	$(LD) $(LDOPTS) -o$(SPACEIFNOBORLAND)$@ lib/$(PLATFORM)/$(MODE)/$(@F).o \
	      -l$(LIBNAME) $(CLIBS)

  endif

else

  ifeq "$(LIBMODEX)" "true"

    ifeq "$(PLATFORM)" "winnt"

      $(addprefix app/, $(CXXPROGS)):	$(NOAPPCXXOB)
	      @echo "$(LDXXOPTS) /OUT:$@ \
		     lib/$(PLATFORM)/$(MODE)/$(@F).o \
		     $(NOAPPCXXOB) $(DEPLIBOBJS) \
		     gdi32.lib user32.lib glut32.lib opengl32.lib \
		    " > linkcmd
	      $(LDXX) @linkcmd
	      @cp -p $@ $@.exe

      $(addprefix app/, $(CPROGS)):;		@echo Not implemented.

    else
      $(addprefix app/, $(CXXPROGS)):;		@echo Not implemented.
      $(addprefix app/, $(CPROGS)):;		@echo Not implemented.
    endif

  else

    $(addprefix app/, $(CXXPROGS)):	$(NOAPPCXXOBJS)
	$(LDXX) $(LDXXOPTS) -o$(SPACEIFNOBORLAND)$@ lib/$(PLATFORM)/$(MODE)/$(@F).o \
	        $(NOAPPCXXOBJS) $(CXXLIBS)

    $(addprefix app/, $(CPROGS)):	$(NOAPPCOBJS)
	$(LD) $(LDOPTS) -o$(SPACEIFNOBORLAND)$@ lib/$(PLATFORM)/$(MODE)/$(@F).o \
	      $(NOAPPCOBJS) $(CLIBS)

  endif

endif


#
# 980310: The reason for the seemingly stupid rm is that we might have
#	  switched platform or opt-mode... This is undetectable. (?)
#
apps app:
	rm -f $(addprefix app/, $(CXXPROGS) $(CPROGS))
	$(MAKE) $(addprefix app/, $(CXXPROGS) $(CPROGS))

#
# 980304: The (important) dependencies are added at the end of the makefile.
#	  We need this one (main source file) dependency to decide
#	  which compiler to use... This is another case where combining the
#	  rules seems impossible.
# 980617: VC++ ignores the '-o' option.
#

# Experimenting with QT-stuff

# QFORMINC	=$(addsuffix .h, $(addprefix include/, $(QFORMBASE)))
# QFORMSRC	=$(addsuffix .C, $(addprefix src/, $(QFORMBASE)))
# $(QFORMINC) : include/%.h : include/%.ui
# 		@echo "This is in the uic rule. Target is " $@
# 		cd include; uic -o $*.h $*.ui
# $(QFORMSRC) : src/%.C : include/%.h include/%.ui
# 		cd include; uic -i $*.h -o ../src/$*.C $*.ui
# QALLMOC		=$(addsuffix .o, $(addprefix lib/$(PLATFORM)/$(MODE)/, $(QNOFORMMOCOBJS) $(QFORMMOCOBJS)))
# $(QALLMOC): lib/$(PLATFORM)/$(MODE)/%.o: include/$(removeprefix moc_ ,%).h
# 		moc -o $@ $<
# 		@if [ ! -d lib ]; then mkdir -p lib; fi
# 		@if [ ! -d lib/$(PLATFORM) ]; then mkdir -p lib/$(PLATFORM); fi
# 		@if [ ! -d lib/$(PLATFORM)/$(MODE) ]; then \
# 		  mkdir -p lib/$(PLATFORM)/$(MODE); fi
# 		$(CXX) $(CXXOPTS) -o$(SPACEIFNOBORLAND)$@ -c src/$*.C
# ifeq "$(PLATFORM)" "winnt"
# 		@mv $*.obj $@
# endif
# 		@if [ ! -d lib ]; then mkdir -p lib; fi
# 		@if [ ! -d lib/$(PLATFORM) ]; then mkdir -p lib/$(PLATFORM); fi
# 		@if [ ! -d lib/$(PLATFORM)/$(MODE) ]; then \
# 		  mkdir -p lib/$(PLATFORM)/$(MODE); fi
# 		$(CXX) $(CXXOPTS) -o$(SPACEIFNOBORLAND)$@ -c src/$*.C
# ifeq "$(PLATFORM)" "winnt"
# 		@mv $*.obj $@
# endif

include/%.h: include/%.ui
		cd include; uic $*.ui -o $*.h; cd ..

src/%.C: include/%.ui
		cd include; uic -impl $*.h $*.ui -o ../src/$*.C; cd ..

src/moc_%.C: include/%.h
		moc $< -o $@

lib/$(PLATFORM)/$(MODE)/moc_%.o:	src/moc_%.C
		@if [ ! -d lib ]; then mkdir -p lib; fi
		@if [ ! -d lib/$(PLATFORM) ]; then mkdir -p lib/$(PLATFORM); fi
		@if [ ! -d lib/$(PLATFORM)/$(MODE) ]; then \
		  mkdir -p lib/$(PLATFORM)/$(MODE); fi
		$(CXX) $(CXXOPTS) -o$(SPACEIFNOBORLAND)$@ -c src/moc_$*.C
ifeq "$(PLATFORM)" "winnt"
		@mv moc_$*.obj $@
endif

lib/$(PLATFORM)/$(MODE)/%.o:	src/%.C
		@if [ ! -d lib ]; then mkdir -p lib; fi
		@if [ ! -d lib/$(PLATFORM) ]; then mkdir -p lib/$(PLATFORM); fi
		@if [ ! -d lib/$(PLATFORM)/$(MODE) ]; then \
		  mkdir -p lib/$(PLATFORM)/$(MODE); fi
		$(CXX) $(CXXOPTS) -o$(SPACEIFNOBORLAND)$@ -c src/$*.C
ifeq "$(PLATFORM)" "winnt"
		@mv $*.obj $@
endif

lib/$(PLATFORM)/$(MODE)/%.o:	src/%.cpp
		@if [ ! -d lib ]; then mkdir -p lib; fi
		@if [ ! -d lib/$(PLATFORM) ]; then mkdir -p lib/$(PLATFORM); fi
		@if [ ! -d lib/$(PLATFORM)/$(MODE) ]; then \
		  mkdir -p lib/$(PLATFORM)/$(MODE); fi
		$(CXX) $(CXXOPTS) -o$(SPACEIFNOBORLAND)$@ -c src/$*.cpp
ifeq "$(PLATFORM)" "winnt"
		@mv $*.obj $@
endif

lib/$(PLATFORM)/$(MODE)/%.o:	app/%.C
		@if [ ! -d lib ]; then mkdir -p lib; fi
		@if [ ! -d lib/$(PLATFORM) ]; then mkdir -p lib/$(PLATFORM); fi
		@if [ ! -d lib/$(PLATFORM)/$(MODE) ]; then \
		  mkdir -p lib/$(PLATFORM)/$(MODE); fi
		$(CXX) $(CXXOPTS) -o$(SPACEIFNOBORLAND)$@ -c app/$*.C
ifeq "$(PLATFORM)" "winnt"
		@mv $*.obj $@
endif

lib/$(PLATFORM)/$(MODE)/%.o:	app/%.cpp
		@if [ ! -d lib ]; then mkdir -p lib; fi
		@if [ ! -d lib/$(PLATFORM) ]; then mkdir -p lib/$(PLATFORM); fi
		@if [ ! -d lib/$(PLATFORM)/$(MODE) ]; then \
		  mkdir -p lib/$(PLATFORM)/$(MODE); fi
		$(CXX) $(CXXOPTS) -o$(SPACEIFNOBORLAND)$@ -c app/$*.cpp
ifeq "$(PLATFORM)" "winnt"
		@mv $*.obj $@
endif

lib/$(PLATFORM)/$(MODE)/%.o:	src/%.c
		@if [ ! -d lib ]; then mkdir -p lib; fi
		@if [ ! -d lib/$(PLATFORM) ]; then mkdir -p lib/$(PLATFORM); fi
		@if [ ! -d lib/$(PLATFORM)/$(MODE) ]; then \
		  mkdir -p lib/$(PLATFORM)/$(MODE); fi
		$(CC) $(COPTS) -o$(SPACEIFNOBORLAND)$@ -c src/$*.c
ifeq "$(PLATFORM)" "winnt"
		@mv $*.obj $@
endif

lib/$(PLATFORM)/$(MODE)/%.o:	app/%.c
		@if [ ! -d lib ]; then mkdir -p lib; fi
		@if [ ! -d lib/$(PLATFORM) ]; then mkdir -p lib/$(PLATFORM); fi
		@if [ ! -d lib/$(PLATFORM)/$(MODE) ]; then \
		  mkdir -p lib/$(PLATFORM)/$(MODE); fi
		$(CC) $(COPTS) -o$(SPACEIFNOBORLAND)$@ -c app/$*.c
ifeq "$(PLATFORM)" "winnt"
		@mv $*.obj $@
endif


#
# 980609: This rule is handy for making libraries for other subprojects
#	  than the one we are currently working on. Hmm... How do we transfer
#	  this into a rule? There is a problem with trailing whitespaces here,
#	  at least... The ugly hack will have to do for the time being.
#

QQQ=$(foreach lib, $(DEPLIBS0), ../$(lib)/lib/$(PLATFORM)/$(MODE)/lib$(lib).a:\;\\t\(cd\\t../$(lib)\;make\\tlib\)\\n\\n)


#------------------------------------------------------------
#
#	Appending the dependencies here. Why is the wildcard
#	function used here? To ensure that the files actually
#	exist, probably. That piece is stolen.
#
# 980306: Tried to include 'hopefully*' the same way, but
#	  it didn't work! Had to 'make' twice to get it to
#	  read the file, seemingly...
#	  A simple 'include' works, but then we have a problem,
#	  if, for instance, we remove the file... Ok, not
#	  really a problem, only an annoying error message
#	  from the first make invocation. This is 'solved' by
#	  using the '-include' statement. The problem of
#	  understanding what happened in the first case,
#	  remains, though...
#	  Somebody knows the answer? Please tell me!
#
#------------------------------------------------------------

ifneq "$(foreach file, $(CXXDEPFILES) $(CDEPFILES), $(wildcard $(file)))" ""
  include $(foreach file, $(CXXDEPFILES) $(CDEPFILES), $(wildcard $(file)))
endif

-include dep/hopefully_a_unique_name.dep
