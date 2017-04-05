FILES = digorg stat crossover crossover2

CC    = gcc
CXX   = g++
# Add -DNDEBUG to skip asserts and stuff.
override CFLAGS := -c -std=c++11 -pedantic -Wall -Wextra $(CFLAGS)
# -Wno-switch -Wno-unused-parameter
GFXLIBS = -lcairo -lSDL
LIBS = -lgsl -lgslcblas

# make, make mode=parallel, make mode=debug, make mode=profile.
ifeq ($(mode),debug)
   FLAGS = $(CFLAGS) -O0 -g3
   LDFLAGS = $(LIBS)
else ifeq ($(mode),profile)
   FLAGS = $(CFLAGS) -O2 -pg -fprofile-arcs -ftest-coverage
   LDFLAGS = $(LIBS) -pg -fprofile-arcs -ftest-coverage
else ifeq ($(mode),parallel)
   FLAGS = $(CFLAGS) -O3 -fopenmp
   LDFLAGS = $(LIBS) -lgomp
else
   FLAGS = $(CFLAGS) -O2
   LDFLAGS = $(LIBS)
# -L/usr/local/lib
endif

SRC_D := $(filter-out crossover%.cc test% stat%, $(wildcard *.cc))
SRC_S := $(filter-out crossover%.cc digorg.cc test%, $(wildcard *.cc))
SRC_C := crossover.cc alignment.cc common.cc
SRC_C2 := crossover2.cc alignment.cc common.cc


SOURCE = $(wildcard *.cc)
OBJ_D = $(SRC_D:.cc=.o)
OBJ_S = $(SRC_S:.cc=.o)
OBJ_C = $(SRC_C:.cc=.o)
OBJ_C2 = $(SRC_C2:.cc=.o)
TESTS := $(patsubst %.cc, %, $(wildcard test_*.cc))

all: $(FILES)

tests: $(TESTS)

.cc.o:
	$(CXX) $(FLAGS) -o $@ $<

.c.o:
	$(CC) $(FLAGS) -o $@ $<

clean:
	-rm -f $(FILES) $(TESTS) *.o *~ *.rpo *.gcov *.gcno gmon.out *.gcda a.out \
		test_bitstring
	touch -d20020101 Makedep.rule
	make depend

depend:
	$(CXX) $(FLAGS) -MM $(SOURCE) >Makedep.rule

Makedep.rule: Makefile
	@echo "The makefile has changed since the last 'make depend'."

#For Emacs-lovers
tags:
	etags *.[cc,h]

tar:
	tar -cvzf evolution.tar.gz Makefile $(SOURCE) harvardforest *.h

# Dependencies for executables
digorg: $(OBJ_D)
	$(CXX) $(LDFLAGS) -o $@ $+

test_%: test_%.o
	$(CXX) $(LDFLAGS) -o $@ $+

stat: $(OBJ_S)
	$(CXX) $(LDFLAGS) -o $@ $+

crossover: $(OBJ_C)
	$(CXX) $(LDFLAGS) -o $@ $+

crossover2: $(OBJ_C2)
	$(CXX) $(LDFLAGS) -o $@ $+

doc: README.org
	emacs -q --batch --visit=README.org --execute='(org-export-as-ascii 3 1 nil nil)'
#	emacs --batch --execute='(setq vc-follow-symlinks nil)' --visit=README.org --execute='(org-export-as-html-and-open nil)' --visit=README.org --execute='(org-export-as-utf8)'

.PHONY: all clean depend tests tags tar doc

# Dependencies for object files
include Makedep.rule
