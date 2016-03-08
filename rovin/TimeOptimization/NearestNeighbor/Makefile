# Author: HJ Noh
# Start date: 2008/06/30
# Nearest Neighbor

#LIBPATH = /usr/lib/sineTest

WCPPCC = g++
WCPPFLAGS = -Wall -I.

#OBJS = RRTPlanner.o
#TARGET = RRTPlanner.cpp
#LIBTARGET = CArm
#APPTARGET = app

all:
	$(WCPPCC) -g -c $(WCPPFLAGS) *.cpp
	cp *.o ../
.cpp.o :
	$(WCPPCC) -g -c $(WCPPFLAGS) $<

install:
#	cp -rf lib$(LIBTARGET).so $(LIBPATH)
#	cp -rf $(LIBTARGET).h $(INSTALL_DIR)/include

clean: 
	rm -drf *.bak *.so *.o

