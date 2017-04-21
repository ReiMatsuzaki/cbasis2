VPATH=../src_cpp ../utils ../math
CPPFLAGS+=-Wall -fPIC
#ifeq (\$(OS),Linux)
#PY_FLAGS=-I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -shared -fPIC -lboost_python 

## ==== build options ====
ARCH=fast
ifeq (${ARCH}, fast)
  CXXFLAGS+= -O3 -DARG_NO_CHECK -DEIGEN_NO_DEBUG
endif

## ==== build directory ====
BINORIG=bin
BINDIR=${BINORIG}/${ARCH}
$(shell mkdir -p ${BINDIR})

## ==== dependencies ====
SRCS:=$(wildcard *.cpp)
DEPS:=$(SRCS:%.cpp=${BINDIR}/%.d)
-include ${DEPS}

## ==== basic relation ====
${BINDIR}/%.o: %.cpp
	${CXX} -c -o $@ -MMD ${CPPFLAGS} ${CXXFLAGS} $<

## ==== MAIN ====
## -- one_eig --
OBJS=one_driv.o read_json.o symmolint.o molecule.o one_int.o two_int.o symgroup.o bmatset.o angmoment.o eigen_plus.o cfunc.o mol_func.o b2eint.o fact.o erfc.o int_exp.o timestamp.o
${BINDIR}/one_eig: $(foreach o, ${OBJS}, ${BINDIR}/$o)
	${CXX} -o $@ $^ ${CXXFLAGS} ${LIBS} -lgsl -lgslcblas
check: ${BINDIR}/one_eig
	../one_eig/bin/fast/one_eig hatom0.json
	./$^ hatom.json
	cat cs.csv
	@echo reference value at w=1 at length are 0.954
#	./$^ h2plus.json; echo "reference is -1.1026342144949"
#	./$^ h2plus_small.json
check_h2plus: ${BINDIR}/one_eig
	../one_eig/bin/fast/one_eig h2plus0.json
	./$^ h2plus.json;
	cat cs.csv
	@echo "ref : E=0.1,  piu=0.690953, sigu=0.007387, beta=1.34948"

## ==== Basic command ====
clean:
	rm -fr ${BINDIR}

