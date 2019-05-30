# Makefile is in main folder
CC:= g++ 

VPATH = src/:serie_problemas/:lib/

CFLAGS:= -O3 -m64 -std=c++11 -g

ROOTINC:= $(shell root-config --cflags)

ROOTLIB:= $(shell root-config --libs)

MAIN:= $(notdir $(wildcard serie_problemas/*.cpp))

SRC:= $(notdir $(wildcard src/*.cpp)) $(notdir $(wildcard src/*.C))

LIB:= $(notdir $(wildcard lib/*.cpp)) $(notdir $(wildcard lib/*.C))

SRC_cpp := $(filter %.cpp, $(SRC))
SRC_C:=$(filter %.C, $(SRC))

LIB_cpp := $(filter %.cpp, $(LIB))
LIB_C :=$(filter %.C, $(LIB))

SRC_OBJ := $(addprefix bin/, $(SRC_cpp:.cpp=.o))
SRC_OBJ += $(addprefix bin/, $(SRC_C:.C=.o))

LIB_OBJ := $(addprefix bin/, $(LIB_cpp:.cpp=.o))
LIB_OBJ += $(addprefix bin/, $(LIB_C:.C=.o)) 

OBJS:= $(SRC_OBJ) $(LIB_OBJ)

EXE:= $(addprefix bin/, $(MAIN:.cpp=.x))

INCLUDE := -I src -I lib

all : $(EXE)

lib: $(LIB_OBJ)

src: $(SRC_OBJ)

obj: $(OBJS)

bin/%.x : bin/%.o $(LIB_OBJ)
	@echo "Making executable $< [$@]"
	@-$(CC) $(CFLAGS) -o $@  $^ $(INCLUDE) $(ROOTINC) $(ROOTLIB) 
	@echo
#@echo "Running executable"
#$@

bin/%.exe : bin/%.o $(OBJS)
	@echo "Making executable $< [$@]"
	@-$(CC) $(CFLAGS) -o $@  $^ $(INCLUDE) $(ROOTINC) $(ROOTLIB) 
	@echo 
#@echo "Running executable"
#$@

#necessary for extension .cpp
bin/%.o : %.cpp
	@echo "Compiling $< [$@]"
	@$(CC) $(CFLAGS) -c -o $@ $< $(INCLUDE) $(ROOTINC)

#needed for extension .C
bin/%.o : %.C
	@echo "Compiling $< [$@]"
	@$(CC) $(CFLAGS) -c -o $@ $< $(INCLUDE) $(ROOTINC)

clean:
	@rm -fv $(LIB_OBJ) $(SRC_OBJ) $(EXE)

clnlib:
	@rm -fv $(LIB_OBJ)

clnsrc:
	@rm -fv $(SRC_OBJ)

clnexe:
	@rm -fv $(EXE)

shwroot:
	@echo "ROOT INCLUDE"
	@echo

shwsrc:
	@echo "SOURCE .cpp"
	@echo $(SRC)
	@echo

shwlib:
	@echo "LIBRARY .cpp"
	@echo $(LIB)
	@echo

shwall:
	@echo "ROOT INCLUDE"
	@echo $(ROOTINC)
	@echo 

	@echo "SOURCE .cpp"
	@echo $(SRC)
	@echo  

	@echo "LIB .cpp"
	@echo $(LIB)
	@echo

