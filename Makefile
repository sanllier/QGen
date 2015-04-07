all: qgen

#-----------------------------------------------------------------------------

GPU =
CURAND = 
MTRAND =
DEBUG =

DFLAG =

#-----------------------------------------------------------------------------

ifneq ($(GPU), )
	DFLAG += -DGPU
endif

ifneq ($(CURAND), )
	DFLAG += -DCURAND
else 
	ifneq ($(MTRAND), )
		DFLAG += -DMTRAND
	endif
endif

#-----------------------------------------------------------------------------

OBJDIR = obj/
INCDIR = include/
SRCDIR = src/
BINDIR = lib/

#-----------------------------------------------------------------------------

CC     = mpicxx
CUDACC = g++
CFLAG  = -I$(INCDIR)

ifeq ($(DEBUG), )
	CFLAG += -O3
endif

#-----------------------------------------------------------------------------

BINFILE =

ifneq ($(GPU), )
	BINFILE = qgen_gpu
else
	BINFILE = qgen_cpu
endif

ifneq ($(DEBUG), )
	BINFILE := $(BINFILE)_d
endif

BINFILE := $(BINFILE).a

#-----------------------------------------------------------------------------

FILES = mpicheck pugixml qgen \
        qindivid_base qindivid_cpu \
        qobservestate qrotoperator sparams \
        randomizer

ifneq ($(GPU), )
	FILES += qindivid_gpu
endif

ifneq ($(CURAND), )
	FILES += random_curand
endif

ifneq ($(MTRAND), )
	FILES += mtrand
endif

#-----------------------------------------------------------------------------

OBJECTS = $(addprefix $(OBJDIR), $(addsuffix .o, $(FILES)))
CFILES =  $(addprefix $(SRCDIR), $(addsuffix .cpp, $(FILES)))

ifneq ($(GPU), )
	OBJECTS += $(OBJDIR)kernels.o
	CFILES += $(SRCDIR)kernels.cu
endif

#-----------------------------------------------------------------------------

qgen: $(OBJECTS)
	@mkdir -p lib
	@echo "\033[30;1;41m "lib" dir was created \033[0m"
	@ar rs $(BINDIR)$(BINFILE) $(OBJECTS)
	@echo "\033[30;1;41m QGen builded successfully! \033[0m"
	@echo "\033[30;1;41m --> $(BINDIR)$(BINFILE) \033[0m"

$(OBJDIR)%.o: $(SRCDIR)%.cpp
	@mkdir -p $(OBJDIR)
	@$(CC) -c $(DFLAG) $(CFLAG) $^ -o $@
	@echo "\033[30;1;46m $@ - done \033[0m\n"

$(OBJDIR)%.o: $(SRCDIR)%.cu
	@mkdir -p $(OBJDIR)
	@nvcc -ccbin=$(CUDACC) -arch=sm_20 -maxrregcount=21 $(CFLAG) $< -c -o $@	
	@echo "\033[30;1;46m $@ - done \033[0m\n"

clean:
	rm -r -f lib
	rm -r -f obj

#-----------------------------------------------------------------------------
