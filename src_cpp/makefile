include local.mk
VPATH=../utils ../math
ARCH=debug
VALGRIND=valgrind --error-limit=no --tool=memcheck --leak-check=full --show-reachable=no 

ifeq (${ARCH},fast)
  # calculation speed
  CXXFLAGS=${INC_PATH} -Wall -O3 -DARG_NO_CHECK -DEIGEN_NO_DEBUG 
  RUN=
endif
ifeq (${ARCH},debug)
  # debug mode
  CXXFLAGS=${INC_PATH} -Wall
  RUN=
#VALGRIND
endif
ifeq (${ARCH},mem)
  # for memory check. (very slow)
  CXXFLAGS=${INC_PATH} -Wall -pg -g -fno-inline -DARG_NO_CHECK -DEIGEN_NO_DEBUG 
  RUN=VALGRIND
endif
ifeq (${ARCH},time)
  # for time profile
  CXXFLAGS=${INC_PATH} -Wall -pg -g -DARG_NO_CHECK -DEIGEN_NO_DEBUG 
  RUN=iprofiler -timeprofiler -T 30
endif
BINDIR=bin/${ARCH}
$(shell mkdir -p ${BINDIR})


SRCS:=$(wildcard *.cpp)
DEPS:=$(SRCS:%.cpp=${BINDIR}/%.d)
-include ${DEPS}

# -- google test --
# read README in googletest
CPPFLAGS += -isystem ${GTEST_DIR}/include
CXXFLAGS += -pthread
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h

GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)
${BINDIR}/gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest-all.cc -o $@
${BINDIR}/gtest.a : ${BINDIR}/gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

# -- basic relation --
${BINDIR}/%.o: %.cpp
	${CXX} -c -o $@ -MMD ${CPPFLAGS} ${CXXFLAGS} $<

# -- symmolint --
OBJS=symmolint.o molecule.o one_int.o two_int.o symgroup.o bmatset.o angmoment.o eigen_plus.o cfunc.o mol_func.o b2eint.o fact.o erfc.o int_exp.o
${BINDIR}/symmolint.a: $(foreach o, ${OBJS}, ${BINDIR}/$o)
	ar rcs $@ $^ 
symmolint.a: ${BINDIR}/symmolint.a

# -- test math --
MATH_OBJS = \
	test_math.o cfunc.o lgamma.o fact.o mol_func.o \
	angmoment.o eigen_plus.o bmatset.o gtest.a
${BINDIR}/test_math: $(foreach o, ${MATH_OBJS}, ${BINDIR}/$o) 
	${CXX} -o $@ $^  ${CXXFLAGS} ${GTEST} ${LIBS} -lgsl -lgslcblas
.PHONY: check_math
check_math: ${BINDIR}/test_math
	${RUN} ./$<

# -- test symgroup --
SYMGROUP_OBJS=test_symgroup.o symgroup.o fact.o gtest.a
${BINDIR}/test_symgroup: $(foreach o, ${SYMGROUP_OBJS}, ${BINDIR}/$o)
	${CXX} -o $@ $^ ${CXXFLAGS} ${GTEST}
.PHONY: check_symgroup
check_symgroup: ${BINDIR}/test_symgroup
	${RUN} ./$<

# -- test ivo --
OBJS = test_ivo.o ivo.o symmolint.o symgroup.o molecule.o one_int.o two_int.o mo.o b2eint.o bmatset.o angmoment.o eigen_plus.o cfunc.o mol_func.o int_exp.o fact.o erfc.o gtest.a
${BINDIR}/test_ivo: $(foreach o, ${OBJS}, ${BINDIR}/$o)
	${CXX} -o $@ $^ ${CXXFLAGS} ${GTEST} -lgsl -lgslcblas
.PHONYS: check_ivo
check_ivo: ${BINDIR}/test_ivo
	${RUN} ./$<

# -- test symmolint --
SYMMOLINT_OBJS = \
	test_symmolint.o symmolint.o molecule.o one_int.o two_int.o symgroup.o bmatset.o angmoment.o eigen_plus.o cfunc.o mol_func.o b2eint.o fact.o erfc.o int_exp.o gtest.a 
#${BINDIR}/test_symmolint: $(foreach o, ${SYMMOLINT_OBJS}, ${BINDIR}/$o)
#	${CXX} -o $@ $^ ${CXXFLAGS} ${GTEST} ${LIBS} -lgsl -lgslcblas
${BINDIR}/test_symmolint: ${BINDIR}/test_symmolint.o ${BINDIR}/symmolint.a ${BINDIR}/gtest.a 
	${CXX} -o $@ $^ ${CXXFLAGS} ${GTEST} ${LIBS} -lgsl -lgslcblas
.PHONY: check_symmolint
check_symmolint: ${BINDIR}/test_symmolint
	${RUN} ./$<


# -- test 2e int --
TWOINT_OBJS = test_2eint.o b2eint.o symmolint.o symgroup.o one_int.o two_int.o bmatset.o angmoment.o  \
	eigen_plus.o molecule.o erfc.o int_exp.o fact.o cfunc.o mol_func.o timer.o gtest.a
${BINDIR}/test_2eint:  $(foreach o, ${TWOINT_OBJS}, ${BINDIR}/$o) 
	${CXX} -o $@ $^ ${CXXFLAGS} ${GTEST} ${LIBS} -lgsl -lgslcblas
.PHONY: check_2eint
check_2eint: ${BINDIR}/test_2eint
	${RUN} ./$<
time_2eint: ${BINDIR}/test_2eint
	iprofiler -timeprofiler ./$<

# -- test hf --
HF_OBJS = test_hf.o trans_eri.o mo.o b2eint.o symmolint.o one_int.o two_int.o \
	symgroup.o bmatset.o angmoment.o eigen_plus.o molecule.o int_exp.o erfc.o\
	fact.o cfunc.o mol_func.o timer.o gtest.a
${BINDIR}/test_hf: $(foreach o, ${HF_OBJS}, ${BINDIR}/$o)
	${CXX} -o $@ $^ ${CXXFLAGS} ${GTEST} ${LIBS} -lgsl -lgslcblas
.PHONY: check_hf
check_hf: ${BINDIR}/test_hf
	${RUNE} ./$<

clean:
	rm -f ${BINDIR}/*.o
	rm -f ${BINDIR}/*.a
	rm -f ${BINDIR}/*.d

# ==== old ====
# -- test L2(R1) --
#L2OBJS=fact.o erfc.o lgamma.o math_utils.o exp_func.o cut_exp.o exp_int.o delta.o
#${BINDIR}/l2.a: $(foreach o, ${L2OBJS}, ${BINDIR}/$o)
#	ar r $@ $^
#${BINDIR}/test_l2: ${BINDIR}/test_l2.o ${BINDIR}/l2.a
#	${CXX} -o $@ ${CXXFLAGS} $^ -lgtest
#.PHONY: /check_l2
#check_l2: ${BINDIR}/test_l2
#	${VALGRIND} ./$<


# -- not used now --
gto3dset.o: gto3dset.cpp gto3dset.hpp

molint.o: molint.cpp molint.hpp math_utils.hpp macros.hpp

# -- not used now --

test_gto3d.o: test_gto3d.cpp gto3d.hpp 
#test_gto3d: test_gto3d.o cints.o angmoment.o gto3dset.o math_utils.o molint.o cip_exp.o exp_func.o eigen_plus.o
test_gto3d: test_gto3d.o cints.o angmoment.o gto3dset.o molint.o spec_func.o eigen_plus.o ${OBJS} 
	${CXX} -o $@ $^  ${CXXFLAGS} ${LIBGTEST} -lgsl
.PHONY: check_gto3d
check_gto3d: test_gto3d
	./test_gto3d

test_gto3d_time.o: test_gto3d_time.cpp
test_gto3d_time:  test_gto3d_time.o cints.o angmoment.o gto3dset.o molint.o spec_func.o eigen_plus.o ${OBJS} timer.o
	${CXX} -o $@ ${CXXFLAGS} ${LIBGTEST} $^ -lgsl
.PHONY: check_gto3d_time
check_gto3d_time: test_gto3d_time
	./test_gto3d_time
profile_gto3d_time: test_gto3d_time
	iprofiler -timeprofiler ./test_gto3d_time

check_gto3d_time_10: test_gto3d_time
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 

inc_gamma_grid.o: inc_gamma_grid.cpp
inc_gamma_grid: inc_gamma_grid.o molint.o
	${CXX} -o $@ ${CXXFLAGS} $^ -lgsl
.PHONY: run_inc_gamma_grid
run_inc_gamma_grid: inc_gamma_grid
	./inc_gamma_grid 

l2func_bind.so: wrapper.cpp ${OBJS}
	${CXX} -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -bundle -fPIC -o $@ wrapper.cpp ${OBJS} ${CXXFLAGS} -lboost_python  -framework Python

utest_py: l2func_bind.so utest.py
	python utest.py


