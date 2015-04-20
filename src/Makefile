# Main makefile

OBJ_CORE = Distance.o Util.o Readconf.o Output.o Cells.o Init.o WriteCfg.o Structure.o Interactions.o Thermo.o Energy.o Move.o Main.o
MKDIR_P = mkdir -p
BINDIR = ../bin
MODPATH = ../modules

BIN = ../bin/gpMC.exe

OBJ_ALL = $(OBJ_CORE) set_precision.o Definitions.o

FCOPTS = -fast -module $(MODPATH)

%.o : %.mod

.SUFFIXES : .f90

.f90.o:
	ifort -c $(FCOPTS)   $<

$(BIN): ${MODPATH} ${BINDIR} $(OBJ_CORE)
	ifort -o $(BIN)  $(OBJ_ALL)

$(OBJ_CORE) : set_precision.o Definitions.o

Definitions.o : set_precision.o

clean:
	rm -f  $(OBJ_ALL)  $(MODPATH)/*.mod

all:
	$(MAKE)  $(BIN)

debug:
	$(MAKE) FCOPTS="-debug -O0 -CB -module $(MODPATH)"

directories: ${BINDIR} ${MODPATH}

${BINDIR}:
	${MKDIR_P} ${BINDIR}

${MODPATH}:
	${MKDIR_P} ${MODPATH}
	