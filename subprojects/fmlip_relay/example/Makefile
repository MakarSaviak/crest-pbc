# fmlip_relay/example/Makefile
# ─────────────────────────────────────────────────────────────────────────────

FC      = gfortran
CC      = gcc
FCFLAGS = -O2 -Wall -fcheck=all -fopenmp
CCFLAGS = -O2 -Wall
LDFLAGS = -fopenmp

SRCDIR  = ../src
BDIR    = build

C_SRC   = $(SRCDIR)/c/socket_utils.c
F90_SRC = $(SRCDIR)/fortran/mlip_client.f90
EX_SRC  = example_usage.f90

C_OBJ   = $(BDIR)/socket_utils.o
F90_OBJ = $(BDIR)/mlip_client.o
EX_OBJ  = $(BDIR)/example_usage.o
EXE     = $(BDIR)/run_example

.PHONY: all test-dummy clean

all: $(BDIR) $(EXE)

$(BDIR):
	mkdir -p $(BDIR)

$(C_OBJ): $(C_SRC)
	$(CC) $(CCFLAGS) -c $< -o $@

$(F90_OBJ): $(F90_SRC) $(C_OBJ)
	$(FC) $(FCFLAGS) -c $< -o $@ -J$(BDIR)

$(EX_OBJ): $(EX_SRC) $(F90_OBJ)
	$(FC) $(FCFLAGS) -I$(BDIR) -c $< -o $@

$(EXE): $(EX_OBJ) $(F90_OBJ) $(C_OBJ)
	$(FC) $(LDFLAGS) $^ -o $@

# Quick smoke-test with the dummy backend (no model needed)
test-dummy: $(EXE)
	./$(EXE) dummy

clean:
	rm -rf $(BDIR)
