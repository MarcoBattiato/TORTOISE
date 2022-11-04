###############################################################################
#                                                                             #
#                  Makefile for projects utilising TORTOISE                   #
#                                                                             #
##################################### Usage ###################################
#                                                                             #
# 1) Place this makefile in the same folder as the main.cpp file              #
# 2) Specify the location of the dependency files                             #
# 3) Edit executable name, compiler flags and the like (optional)             #
# 4) Run make all                                                             #
#                                                                             #
# PS: Do not forget to specify the path to GNUPlot in the Configuration file  #
#                                                                             #
############################## Editable settings ##############################

# Location of the library TORTOISE
TORT_DIR := ./TORTOISE/
# Location of the library Eigen
EIGE_DIR := /usr/local/Cellar/eigen/3.4.0/
# Location of the project main file and subdirectories
USER_DIR := .
# Name of the main file
EXE      := main
# Compiler
CXX      := g++
# Compilation flags
FLAGS    := -std=c++20 -O3 -flto

###############################################################################

####################### Do not edit below this point ##########################

# sources, objects and headers from the base directory
SRCS_TORT := $(shell find $(TORT_DIR)src -type f \( -name "*.cpp" \) -print)
SRCS_USR := $(shell find $(USER_DIR) -not \( -path "$(TORT_DIR)*" -prune \) -type f \( -name "*.cpp" \) -print)
OBJS := $(SRCS_USR:%.cpp=%.o) $(SRCS_TORT:%.cpp=%.o)
DEPS := $(OBJS:.o=.d)

# folders to include with -I
INC_DIRS := $(shell find $(USER_DIR) -not \( -path "$(TORT_DIR)*" -prune \) -not \( -path "./TORTOISE" -prune \) -type d \( -not -name "." \) -print) \
			$(TORT_DIR)include  $(EIGE_DIR)
INC_FLAGS := $(addprefix -I ,$(INC_DIRS))

# c++ flags, -I and automatic dependency generator
CPPFLAGS := $(FLAGS) $(INC_FLAGS) -MMD -MP

########################### Makefile rules start here ########################

default :
	@echo 'Typing "make" with no arguments displays this message.'
	@echo 'Specific targets for compiling are:'
	@echo '  all       compiles the whole project'
	@echo '  clean     removes executable, objects, archive and dependencies'

.PHONY: clean

all : $(EXE)

$(EXE) : $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

-include $(DEPS)

clean :
	@rm -f $(OBJS) $(EXE) $(DEPS)

