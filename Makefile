## ChIA-PET2
## Copyleft 2015 Tsinghua University
## Author: Guipeng Li
## Contact: guipeng.lee@gmail.com
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.


## DO NOT EDIT THE REST OF THIS FILE!!

MK_PATH = $(shell dirname $(abspath $(lastword $(MAKEFILE_LIST))))
## tmp = $(shell chmod +x $(MA_PATH)/bin/ChIA-PET2)
VNUM = $(shell $(MK_PATH)/bin/ChIA-PET2 --version | cut -d " " -f 3)

CC=g++
CFLAGS=-Wall -O2 -Wno-char-subscripts -fopenmp
LDFLAGS= -lz
#CFLAGS=-Wall -Wno-char-subscripts -g -lz -fpermissive
SOURCES=$(MK_PATH)/src
BIN=$(MK_PATH)/bin

all : trimLinker buildBedpe removeDup buildTagAlign tag2Depth bedpe2Interaction bedpe2Phased bedpe2Matrix
install : cp


######################################
## Config file
##
######################################
ifndef PREFIX
PREFIX = ${HOME}/bin
endif


######################################
## Compile
##
######################################

## Build C++ code
trimLinker: $(SOURCES)/trimLinker.cpp
	$(CC) $(CFLAGS) -o ${BIN}/$@ ${SOURCES}/trimLinker.cpp ${LDFLAGS}

buildBedpe: $(SOURCES)/buildBedpe.cpp
	$(CC) $(CFLAGS) -o ${BIN}/$@ ${SOURCES}/buildBedpe.cpp ${LDFLAGS}

removeDup: $(SOURCES)/removeDup.cpp
	$(CC) $(CFLAGS) -o ${BIN}/$@ ${SOURCES}/removeDup.cpp ${LDFLAGS}

buildTagAlign: $(SOURCES)/buildTagAlign.cpp
	$(CC) $(CFLAGS) -o ${BIN}/$@ ${SOURCES}/buildTagAlign.cpp ${LDFLAGS}

tag2Depth: $(SOURCES)/tag2Depth.cpp
	$(CC) $(CFLAGS) -o ${BIN}/$@ ${SOURCES}/tag2Depth.cpp ${LDFLAGS}

bedpe2Interaction: $(SOURCES)/bedpe2Interaction.cpp
	$(CC) $(CFLAGS) -o ${BIN}/$@ ${SOURCES}/bedpe2Interaction.cpp ${LDFLAGS}

bedpe2Phased: $(SOURCES)/bedpe2Phased.cpp
	$(CC) $(CFLAGS) -o ${BIN}/$@ ${SOURCES}/bedpe2Phased.cpp ${LDFLAGS}

bedpe2Matrix: $(SOURCES)/bedpe2Matrix.cpp
	$(CC) -Wall -O2 -std=c++0x -o ${BIN}/$@ ${SOURCES}/bedpe2Matrix.cpp ${LDFLAGS}

######################################
## Create installed version
##
######################################

cp:
ifneq ($(realpath $(MK_PATH)), $(realpath $(PREFIX))/ChIA-PET2_$(VNUM))
	if [ ! -d ${PREFIX} ]; then mkdir ${PREFIX}; fi
	chmod +x ${BIN}/ChIA-PET2
	cp -ri $(MK_PATH) $(PREFIX)/ChIA-PET2_$(VNUM)
endif
	@echo "Install ChIA-PET2 in $(realpath $(PREFIX))/ChIA-PET2_$(VNUM) ..."

clean:
	rm -rf $(PREFIX)/ChIA-PET2_$(VNUM)
