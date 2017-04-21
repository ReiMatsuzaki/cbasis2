TARGET=opt_r1basis
VPATH=../r1basis ../utils ../math
CPPFLAGS+=-Wall -fPIC

## ==== build options ====
#ARCH=fast
ARCH=debug
ifeq (${ARCH}, fast)
  CXXFLAGS+= -O3 -DARG_NO_CHECK -DEIGEN_NO_DEBUG
endif
ifeq (${ARCH}, debug)
  CXXFLAGS+= 
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
OBJS=${TARGET}.o r1_lc.o r1basis.o nderiv.o read_json.o eigen_plus.o erfc.o int_exp.o fact.o timestamp.o

${BINDIR}/${TARGET}: $(foreach o, ${OBJS}, ${BINDIR}/$o)
	${CXX} -o $@ $^ ${CXXFLAGS} ${LIBS} -lgsl -lgslcblas
check: ${BINDIR}/${TARGET}
	$< test_grad_hess.json

check2: ${BINDIR}/${TARGET}
	./$^ h1skpl.json
	@echo  --------Convergence-----------
	@echo  energy:   0.500000000000000     
	@echo zeta_1:   0.9827984729  -0.0127337071
	@echo zeta_2:   0.9180784353  -0.6594777057
	@echo zeta_3:   0.3344314317  -1.0972055028
## ==== Basic command ====
clean:
	rm -fr ${BINDIR}


