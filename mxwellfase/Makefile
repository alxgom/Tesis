#
# This Makefile was automatically generated by Code::Blocks IDE.
#

SRCS_f90d1 = \
module.f90 \
legendre.f90 

SRCS_f95d1 = \
main.f95 

OBJS_f90d1 = \
module.o \
legendre.o 

OBJS_f95d1 = \
main.o 

SRC_DIR_f90d1 = 

SRC_DIR_f95d1 = 
OBJS_DIR = obj/Debug/
EXE_DIR = bin/Debug/

EXE = mxwellfase.exe
FC = mingw32-gfortran.exe
IDIR = 
CFLAGS = -Wall -g  -J$(OBJS_DIR) $(IDIR)
LFLAGS = 
LIBS = 

VPATH = $(SRC_DIR_f90d1):$(OBJS_DIR):$(SRC_DIR_f95d1):$(OBJS_DIR)
OBJS = $(addprefix $(OBJS_DIR), $(OBJS_f90d1) $(OBJS_f95d1))

all : $(EXE)

$(EXE) : $(OBJS_f90d1) $(OBJS_f95d1)
	@mkdir -p $(EXE_DIR)
	$(FC) -o $(EXE_DIR)$(EXE) $(OBJS) $(LFLAGS) $(LIBS)

$(OBJS_f90d1):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f90d1)$(@:.o=.f90) -o $(OBJS_DIR)$@

$(OBJS_f95d1):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f95d1)$(@:.o=.f95) -o $(OBJS_DIR)$@

clean :
	rm -f $(OBJS_DIR)*.*
	rm -f $(EXE_DIR)$(EXE)

# Dependencies of files
module.o: \
    module.f90
legendre.o: \
    legendre.f90 \
    module.o
main.o: \
    main.f95 \
    module.o

