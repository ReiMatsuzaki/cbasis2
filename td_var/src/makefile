TARGET=td_var
VPATH=../../r1basis ../../utils ../../math
CPPFLAGS+=-Wall -fPIC

## ==== build options ====
ARCH=fast
#ARCH=debug
ifeq (${ARCH}, fast)
  CXXFLAGS+= -O3 -DARG_NO_CHECK -DEIGEN_NO_DEBUG
endif
ifeq (${ARCH}, debug)
  CXXFLAGS+= 
endif

## ==== build directory ====
TESTORIG=$(realpath ../test)
BINORIG=$(abspath ../bin)
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
OBJS=${TARGET}.o r1_lc.o r1basis.o nderiv.o read_json.o eigen_plus.o erfc.o int_exp.o fact.o timestamp.o

${BINDIR}/${TARGET}: $(foreach o, ${OBJS}, ${BINDIR}/$o)
	${CXX} -o $@ $^ ${CXXFLAGS} ${LIBS} -lgsl -lgslcblas
check: ${BINDIR}/${TARGET}
	cd ${TESTORIG}/one_sto; $< one_sto.json; head czs.csv
#	cd ${TESTORIG}/3sto; $< h1skpl.json


## ==== Basic command ====
clean:
	rm -fr ${BINDIR}


