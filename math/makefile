include ../local.mk
include ../common.mk

## ==== Main build ====
OBJS=test_math.o nderiv.o fact.o erfc.o int_exp.o 
${BINDIR}/test_math: $(foreach o, ${OBJS}, ${BINDIR}/$o) ${BINDIR}/gtest.a
	${CXX} -o $@ $^ ${CPPFLAGS} ${CXXFLAGS}

## ==== launch unit test ====
check: ${BINDIR}/test_math
	./$^

