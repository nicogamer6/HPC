# -------------- #
# -- Makefile -- #
# -------------- #


# -- Lile list ----------
FILE = main.c mouvement.c nrutil.c test_mouvement.c vnrutil.c morpho.c test_morpho.c matriceROC.c mouvement_SSE2.c test_mouvement_SSE2.c morpho_SSE2.C test_morpho_SSE2.C

# -- Paths ----------
SRC_PATH = src
OBJ_PATH = obj
EXE_PATH = exe
INC_PATH = include

IMG_PATH_FD = testFD/* testFDmorphoF/* testFDmorphoFO/* testFDmorphoO/* testFDmorphoOF/* testFD_SSE/* testFD_SSEmorphoF/* testFD_SSEmorphoFO/* testFD_SSEmorphoO/* testFD_SSEmorphoOF/* testFD_SSE_OMP/*
IMG_PATH_SD = testSD/* testSDmorphoF/* testSDmorphoFO/* testSDmorphoO/* testSDmorphoOF/* testSD_SSE/* testSD_SSEmorphoF/* testSD_SSEmorphoFO/* testSD_SSEmorphoO/* testSD_SSEmorphoOF/* testSD_SSE_OMP/*
IMG_PATH_morpho = testmorpho/* testmorphoSSE/*

# -- OS ----------
OS = MACH_OSX
#OS = LINUX

# -- Config ----------
# if CONFIG = CLI  (Command Line Interface, no Apple Framework)
CONFIG = CLI

# -- Macros ----------
CC = gcc
AR = ar -rc

# -- Flags ----------
C_DEBUG_FLAGS = -O0
C_CC_FLAGS = -std=c99 -DNOALIAS -DALIGNED -mssse3 -fopenmp 
C_SSE_FLAGS = -mfpmath=sse -mmmx -msse -msse2 -msse3
C_OPTIMISATION_FLAGS = -O3 -fstrict-aliasing

#C_ARCH_FLAGS = -xSSE4.2
C_ARCH_FLAGS =

C_OS_FLAGS = -D$(OS)
C_CONFIG_FLAGS = -D$(CONFIG)
C_INC_FLAGS = -I$(INC_PATH)

CFLAGS =  $(C_CC_FLAGS) $(C_DEBUG_FLAGS)        $(C_ARCH_FLAGS) $(C_OS_FLAGS) $(C_CONFIG_FLAGS) $(C_INC_FLAGS) $(LIB_INC_PATH) 
CFLAGS = $(C_CC_FLAGS) $(C_OPTIMISATION_FLAGS) $(C_ARCH_FLAGS) $(C_OS_FLAGS) $(C_CONFIG_FLAGS) $(C_INC_FLAGS) $(LIB_INC_PATH)
LDFLAGS = $(C_CC_FLAGS) $(C_OPTIMISATION_FLAGS) $(C_ARCH_FLAGS) $(C_OS_FLAGS) $(C_CONFIG_FLAGS) $(C_INC_FLAGS) $(LIB_LIB_PATH) 

# -- Final product ----------
PRODUCT   = hpc

# -- src and obj List ----------
SRC = $(addprefix ${SRC_PATH}/, $(FILE))
OBJ = $(addprefix ${OBJ_PATH}/, $(addsuffix .o, $(basename $(FILE))))

# -- Base rules ----------
$(OBJ_PATH)/%.o : $(SRC_PATH)/%.c
	$(CC) $(CFLAGS) -c $< -o $@
   
#-----Main rule ----------
$(EXE_PATH)/$(PRODUCT): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS) $(OPTFLAGS) $(CFG) $(INC) $(LIB) -lm

# -- Other stuff ----------
depend:
	makedepend $(CFLAGS) -Y $(SRC)

clean:
	rm -f $(OBJ)
	rm -f ${EXE_PATH}/${PRODUCT}
	rm -f $(IMG_PATH_FD) $(IMG_PATH_SD) $(IMG_PATH_morpho)

