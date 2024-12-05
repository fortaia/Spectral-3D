#
# adapted from: https://fortran-lang.org/learn/building_programs/project_make
#
# Disable the default rules
.SUFFIXES:
MAKEFLAGS += --no-builtin-rules --no-builtin-variables
SHELL=/bin/bash

# Project name
NAME := run3d

TARGET := $(NAME)
INPUT_FILE := THI_PJET

ROOT_DIR := .
SRCS_DIR := $(ROOT_DIR)/src
EXE_DIR := $(ROOT_DIR)/bin
LIBS_DIR := $(ROOT_DIR)/dependencies
LIBS :=
INCS :=

DEFINES :=

EXE := $(EXE_DIR)/$(TARGET)

# Configuration settings
FC := mpifort
FFLAGS :=
AR := ar rcs
LD := $(FC)
RM := rm -f
GD := $(SRCS_DIR)/gen-deps.awk
CPP := -cpp

# edit build.conf file desired
include $(ROOT_DIR)/build.conf

# List of all source files
SRCS := $(wildcard $(SRCS_DIR)/*.f90)

# Add source directory to search paths
vpath % .:$(SRCS_DIR)
vpath % $(patsubst -I%,%,$(filter -I%,$(INCS)))

# Define a map from each file name to its object file
obj = $(src).o
$(foreach src, $(SRCS), $(eval $(src) := $(obj)))

# Create lists of the build artefacts in this project
OBJS := $(addsuffix .o, $(SRCS))
LIB := $(patsubst %, lib%.a, $(NAME))
DEPS := $(SRCS_DIR)/.depend.mk

# Declare all public targets
.PHONY: all clean allclean libs libsclean 
all: $(EXE)

$(EXE): $(OBJS)
	$(FC) $(FFLAGS) $^ $(LIBS) $(INCS) -o $(EXE)
	@cp $(SRCS_DIR)/$(INPUT_FILE) $(EXE_DIR)
	@printf "\nDefault input file $(INPUT_FILE) copied to run folder $(EXE_DIR)\n"

# Create object files from Fortran source
$(OBJS): %.o: %
	$(FC) $(FFLAGS) $(CPP) $(DEFINES) $(INCS) $(FFLAGS_MOD_DIR) $(SRCS_DIR) -c -o $@ $<

# Process the Fortran source for module dependencies
$(DEPS):
	@echo '# This file contains the module dependencies' > $(DEPS)
	@$(foreach file, $(SRCS), $(GD) $(file) >> $(DEPS))

# Define all module interdependencies
-include $(DEPS)
$(foreach dep, $(OBJS), $(eval $(dep): $($(dep))))

# Cleanup, filter to avoid removing source code by accident
clean:
	$(RM) $(SRCS_DIR)/*.{i,mod,smod,d,o} $(EXE) $(DEPS) $(INPUT_FILE)

allclean:
	@make libsclean
	@make clean
#
# rules for building the external libraries (compile with 'make libs'):
#
include $(LIBS_DIR)/external_libs.mk
