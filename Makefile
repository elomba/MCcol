# Main makefile for program gpMC (project MCcol)

RESPATH = results/
VPATH = src/
OBJDIR = obj
OBJ_CORE = $(addprefix $(OBJDIR)/, Distance.o Util.o Readconf.o Output.o Cells.o Init.o\
 WriteCfg.o Structure.o Interactions.o Dump.o Thermo.o \
 Energy.o Move.o VolumeMove.o Main.o)
MKDIR_P = mkdir -p
BINDIR = bin
MODPATH = modules

BIN = bin/gpMC.exe

OBJ_ALL = $(OBJ_CORE) $(OBJDIR)/set_precision.o $(OBJDIR)/Definitions.o  

FCOPTS = -fast -module $(MODPATH) 

%.o : %.mod

.SUFFIXES : .f90
	
$(OBJDIR)/%.o: %.f90 ${OBJDIR}
	ifort -c $(FCOPTS) $< -o $@

$(BIN): ${MODPATH} ${BINDIR} $(OBJ_CORE) $(OBJDIR)/libutil.o
	ifort $(LKOPTS) -o $(BIN)  $(OBJ_ALL) $(OBJDIR)/libutil.o ; ${MKDIR_P} ${RESPATH}

$(OBJ_CORE) : $(OBJDIR)/set_precision.o $(OBJDIR)/Definitions.o

$(OBJDIR)/Definitions.o : $(OBJDIR)/set_precision.o ${OBJDIR}

$(OBJDIR)/libutil.o: $(VPATH)/libutil.c
	$(CC) -c -o $(OBJDIR)/libutil.o $(VPATH)/libutil.c
	
clean:
	rm -f  $(OBJ_ALL)  $(MODPATH)/*.mod

all:
	$(MAKE)  $(BIN)

debug:
	$(MAKE) FCOPTS="-debug -O0 -CB -module $(MODPATH)"

debug-mac:
	$(MAKE) FCOPTS="-debug -O0 -CB -module $(MODPATH)" LKOPTS="-Wl,-no_pie"
	
directories: ${BINDIR} ${MODPATH} ${OBJDIR}

${BINDIR}:
	${MKDIR_P} ${BINDIR}

${MODPATH}:
	${MKDIR_P} ${MODPATH}
	
${OBJDIR}:
	${MKDIR_P} ${OBJDIR}
