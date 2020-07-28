CFLAGS = -Wuninitialized -g


ODIR=/quac/QuaC/obj
SRCDIR=/quac/QuaC/src
#TESTS=$(basename $(notdir $(wildcard $(TESTDIR)/*test*.c)))
CFLAGS += -isystem $(SRCDIR)

include ${SLEPC_DIR}/lib/slepc/conf/slepc_variables

_DEPS = quantum_gates.h dm_utilities.h operators.h solver.h operators_p.h quac.h quac_p.h kron_p.h qasm_parser.h error_correction.h qsystem.h qsystem_p.h qvec_utilities.h quantum_circuits.h
DEPS  = $(patsubst %,$(SRCDIR)/%,$(_DEPS))

_OBJ  = quac.o operators.o solver.o kron.o dm_utilities.o quantum_gates.o error_correction.o qasm_parser.o qsystem.o qvec_utilities.o quantum_circuits.o
 OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

_TEST_OBJ  = unity.o t_helpers.o
TEST_OBJ = $(patsubst %,$(ODIR)/%,$(_TEST_OBJ))

_TEST_DEPS = t_helpers.h
TEST_DEPS  = $(patsubst %,$(TESTDIR)/%,$(_TEST_DEPS))

obj/%.o: %.c $(DEPS)
	@mkdir -p $(@D)
	${PETSC_COMPILE} -c -o $@ $< $(CFLAGS) ${PETSC_KSP_LIB} ${PETSC_CC_INCLUDES} ${SLEPC_EPS_LIB}


.phony: clean_test test count_fails

count_fails:
	@echo "All failures listed below"
	@grep FAIL test_results || true

clean_test:
	rm -f $(TEST_OBJ)
	@rm -f test_results

simple_circuit: obj/simple_circuit.o $(OBJ)
	${CLINKER} -o $@ $^ $(CFLAGS) ${PETSC_KSP_LIB} ${SLEPC_EPS_LIB}

UChicagoCavity: obj/UChicagoCavity.o $(OBJ)
	${CLINKER} -o $@ $^ $(CFLAGS) ${PETSC_KSP_LIB} ${SLEPC_EPS_LIB}

.PHONY: clean

clean:
	rm -f simple_circuit
	rm -f obj/*
