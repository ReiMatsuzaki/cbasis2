include local.mk
F90=ifort
CXX=icc
VPATH=../utils ../math ../src_cpp
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
${BINDIR}/%.o: %.f90
	${F90} -c -o $@ $<

# -- test --
OBJS = \
	read_molint.o test_ccol.o cfunc.o lgamma.o fact.o mol_func.o \
	read_aoint_wrapper.o\
	angmoment.o eigen_plus.o bmatset.o  molecule.o symgroup.o symmolint.o gtest.a
#OBJS = \
#	test_ccol.o  read_molint.o read_aoint_wrapper.o  symmolint.o gtest.a
#${BINDIR}/symmolint.a:
#	echo "compile symmolint.a"
#	cd ../src_cpp; make ARCH=${ARCH} BINDIR=${BINDIR}
#symmolint.a: ${BINDIR}/symmolint.a

${BINDIR}/test: $(foreach o, ${OBJS}, ${BINDIR}/$o) 
	${CXX} -o $@ $^  ${CXXFLAGS} ${GTEST} ${LIBS} -lgsl -lgslcblas -lifport -lifcore

.PHONY: check
check: ${BINDIR}/test
	cd out/; ~/src/prog/cCOLUMBUS/molint<../molint_h2plus.in #> molint_h2plus.out
	${RUN} ./$<
#	cat molint_h2plus_cpp.in

.PHONY: clean
clean:
	rm -f ${BINDIR}/test
	rm -f ${BINDIR}/*.o
	rm -f ${BINDIR}/*.a
	rm -f ${BINDIR}/*.d
