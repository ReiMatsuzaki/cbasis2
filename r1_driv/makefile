VPATH=../r1basis ../utils ../math
CPPFLAGS+=-Wall -fPIC

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
OBJS=r1_driv.o r1_lc.o r1basis.o read_json.o eigen_plus.o erfc.o int_exp.o fact.o timestamp.o
#R1BASIS_OBJS=r1_lc.o opt_green.o r1basis.o erfc.o eigen_plus.o int_exp.o fact.o
${BINDIR}/r1_driv: $(foreach o, ${OBJS}, ${BINDIR}/$o)
	${CXX} -o $@ $^ ${CXXFLAGS} ${LIBS} -lgsl -lgslcblas
check: ${BINDIR}/r1_driv
	$< h1skpl.json

## ==== Basic command ====
clean:
	rm -fr ${BINDIR}


