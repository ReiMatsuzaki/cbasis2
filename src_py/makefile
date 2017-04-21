include local.mk
#VPATH=../src_cpp ../utils ../math

# ==== Common settings ====
CXXFLAGS=-Wall -fPIC

#ifeq (\$(OS),Linux)
PY_FLAGS=-I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -shared -fPIC -lboost_python 
#endif

ifeq (\$(OS), Darwin)
#PY_FLAGS=-I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -lboost_python -DPIC -bundle -fPIC -framework Python 
PY_FLAGS= -lboost_python -DPIC -bundle -fPIC -framework Python 
endif


BIN_ORIG=bin
ARCH=debug
BINDIR=./bin/${ARCH}

$(shell mkdir -p ${BINDIR})

ifeq (${ARCH},fast)
  # calculation speed
  CXXFLAGS=${INC_PATH} -Wall -O3 -DARG_NO_CHECK -DEIGEN_NO_DEBUG  -fPIC
endif
ifeq (${ARCH},debug)
  # debug mode
  CXXFLAGS=${INC_PATH} -Wall -g  -fPIC
endif
ifeq (${ARCH},mem)
  # for memory check. (very slow)
  CXXFLAGS=${INC_PATH} -Wall -pg -g -fno-inline -DARG_NO_CHECK -DEIGEN_NO_DEBUG  -fPIC
endif


SRCS:=$(wildcard *.cpp)
DEPS:=$(SRCS:%.cpp=${BINDIR}/%.d)
-include ${DEPS}

# -- basic relation --
${BINDIR}/%.o: ../src_cpp/%.cpp
	${CXX} -c -o $@ -MMD ${CXXFLAGS} $<

${BINDIR}/%.o: ../utils/%.cpp
	${CXX} -c -o $@ -MMD ${CXXFLAGS} $<

${BINDIR}/%.o: ../math/%.cpp
	${CXX} -c -o $@ -MMD ${CXXFLAGS} $<

# ==== L2(R1) ====
L2_OBJS=fact.o erfc.o lgamma.o math_utils.o exp_func.o cut_exp.o cip_exp.o delta.o
${BINDIR}/l2func_bind.so: l2func_bind.cpp $(foreach o, ${L2_OBJS}, ${BINDIR}/$o)
	${CXX} -o $@ $^ ${CXXFLAGS} ${PY_FLAGS}
.PHONY: check_l2
check_l2: ${BINDIR}/l2func_bind.so
	ln -s ${BINDIR}/l2func_bind.so .
	python test.py
	rm l2func_bind.so
install_l2: ${BINDIR}/l2func_bind.so
	cp $< ${PY_DIR}

# ==== SymMolInt ====
SYMMOLINT_OBJS=symmolint.o one_int.o fact.o cfunc.o mol_func.o bmatset.o angmoment.o eigen_plus.o symgroup.o molecule.o int_exp.o erfc.o
${BINDIR}/symmolint_bind.so: symmolint_bind.cpp $(foreach o, ${SYMMOLINT_OBJS}, ${BINDIR}/$o)
	${CXX} -o $@ $^ ${CXXFLAGS} ${PY_FLAGS} ${LIBS} -lgsl -lgslcblas
check_symmolint: ${BINDIR}/symmolint_bind.so
	rm -f symmolint_bind.so
	ln -s ${BINDIR}/symmolint_bind.so .
#	python test_solve_pi.py
	python test_symmolint.py
#	python test_part_wave.py
#	rm symmolint_bind.so

install_symmolint: ${BINDIR}/symmolint_bind.so
	cp $< ${PY_DIR}
	cp symmolint_plus.py ${PY_DIR}
	cp symmolint.py ${PY_DIR}
	cp part_wave.py ${PY_DIR}
	cp solve_pi.py ${PY_DIR}

# ==== R1GTOs ====
R1GTOINT_OBJS=fact.o erfc.o lgamma.o exp_int.o eigen_plus.o opt_alpha.o r1gtoint.o l_algebra.o op_driv.o
${BINDIR}/r1gtoint_bind.so: r1gtoint_bind.cpp $(foreach o, ${R1GTOINT_OBJS}, ${BINDIR}/$o)
	${CXX} -o $@ $^ ${CXXFLAGS} ${PY_FLAGS} ${LDFLAGS}
check_r1gtoint: ${BINDIR}/r1gtoint_bind.so
	rm -f ./r1gtoint_bind.so
	ln -s ${BINDIR}/r1gtoint_bind.so .
	python test_r1gtoint.py
	rm r1gtoint_bind.so
install_r1gtoint: ${BINDIR}/r1gtoint_bind.so
	cp ${BINDIR}/r1gtoint_bind.so ~/local/Python/
	cp r1gtoint.py ~/local/Python/

# ==== operation ====
.PHONY: clean
clean:
	rm -r -f ${BINDIR}

# ==== old ====
check_gto3d: l2func_bind.so
	python test_gto3d.py

l2func_molint_bind.o: l2func_molint_bind.cpp
	${CXX} -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -bundle -fPIC -o $@ -c $< ${CXXFLAGS} -framework Python

l2func_molint_bind.so: l2func_molint_bind.cpp ${OBJMOLINT} ../eigen_plus.o
	${CXX} -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -bundle -fPIC -o $@ $^  ${CXXFLAGS} -lboost_python -lboost_numpy -lgsl -lblas -framework Python

check_molint_bind: l2func_molint_bind.so
	python test_molint.py
