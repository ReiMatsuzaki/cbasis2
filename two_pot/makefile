include ../local.mk
VPATH=../src_cpp ../utils ../math
CPPFLAGS+=-Wall -fPIC
#ifeq (\$(OS),Linux)
#PY_FLAGS=-I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -shared -fPIC -lboost_python 

## ==== build options ====
ARCH=fast
ifeq (${ARCH}, fast)
  CXXFLAGS+=-I${INC_PATH} -O3 -DARG_NO_CHECK -DEIGEN_NO_DEBUG
endif
ifeq (${ARCH},debug)
  CXXFLAGS=-I${INC_PATH} -Wall
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
OBJS=two_pot.o read_json.o symmol_read_json.o symmolint.o molecule.o one_int.o two_int.o symgroup.o bmatset.o angmoment.o eigen_plus.o cfunc.o mol_func.o b2eint.o fact.o erfc.o int_exp.o timestamp.o two_int.o mo.o ivo.o
${BINDIR}/two_pot: $(foreach o, ${OBJS}, ${BINDIR}/$o)
	${CXX} -o $@ $^ ${CXXFLAGS} ${LIBS} -lgsl -lgslcblas
check_h2plus: ${BINDIR}/two_pot
	../one_eig/bin/fast/one_eig h2plus0.json
	./$^ h2plus.json
	@echo "ref : E=0.1,  piu=0.690953, sigu=0.007387, beta=1.34948"
	cat cs.csv

check_h2plus_eig: ${BINDIR}/two_pot
	../one_eig/bin/fast/one_eig h2plus0.json
	./$^ h2plus_eig.json
	@echo "ref : E=0.1,  piu=0.690953, sigu=0.007387, beta=1.34948"
	cat cs.csv

# below is 2017/4/18
#Cs(total,alpha):        -0.6908653097,        -0.6949812508
#Cs(sigu,alpha):         -0.0072038318,        -0.0069877573
#Cs(piu,alpha):          -0.6836614779,        -0.6879934935
#coefficient A_lm
#1, -1,        -0.3862386383,         0.0408794718
#1, 0,         0.0004736234,         0.0000956235
#1, 1,        -0.3862386383,         0.0408794718
#3, -1,        -0.0038349705,        -0.0447900564
#3, 0,        -0.0057271232,        -0.0564848179
#3, 1,        -0.0038349705,        -0.0447900564
#Cs(total):         0.6908562458,         0.6949827261
#Cs(sig_u):         0.0072079351,         0.0069865847
#Cs(pi_u) :         0.6836483107,         0.6879961414
#beta     :         1.3878035287,         1.3672135054
#coef:         0.5338084163
#Coef(M=0) = 
#0:         0.0017207678,         0.0016679242
#1:         0.0000000000,         0.0000000000
#2:         0.0022774275,         0.0014104359
#3:         0.0000000000,         0.0000000000
#4:         0.0027929144,         0.0015255229
#5:         0.0000000000,         0.0000000000
#6:         0.0052140701,         0.0049112876
#Cs(M=0,fixed):         0.0072079351,         0.0069865847
#Coef(M=1) = 
#0:        -0.0816045060,        -0.0821234901
#1:         0.0000000000,         0.0000000000
#2:         0.0800458279,         0.0806520444
#3:         0.0000000000,         0.0000000000
#4:        -0.0008930265,        -0.0008852961
#5:         0.0000000000,         0.0000000000
#6:         0.0024517046,         0.0023567417
#Cs(M=1,fixed):        -0.6836483107,        -0.6879961414

EXAM=h2_stex_small
check: ${BINDIR}/two_pot
	../one_eig/bin/fast/one_eig ${EXAM}.oneeig.json
	../rhf/bin/fast/rhf ${EXAM}.rhf.json	
	./$^ ${EXAM}.twopot.json
#w_eV: 32.7252
#w_au: 1.20263
#E_eV: (16.6555,0)
#E_au: (0.612077,0)
#k_au: (1.10641,0)
#[  CalcDriv] Tue Apr 18 16:23:03 2017
#[CalcDriv_eigen_value] Tue Apr 18 16:23:03 2017
#calculate psi1
#calculate psi0. L = 1
#convert V, HV. L = 1
#calculate psi1
#calculate psi0. L = 1
#convert V, HV. L = 1
#calculate psi1
#calculate psi0. L = 1
#convert V, HV. L = 1
#[calc_braket] Tue Apr 18 16:23:03 2017
#[calc_alpha] Tue Apr 18 16:23:03 2017
#Cs(total,alpha):        -2.4812808959,        -1.6083991073
#Cs(sigu,alpha):         -0.5973002912,        -0.2034588730
#Cs(piu,alpha):          -1.8839806048,        -1.4049402343
#coefficient A_lm
#1, -1,         0.3611375502,         0.5391116956
#1, 0,         0.0536986019,         0.5136438775
#1, 1,         0.3611375502,         0.5391116956
#Cs(total):         2.4793724102,         1.6086774500
#Cs(sig_u):         0.5963746066,         0.2035607649
#Cs(pi_u) :         1.8829978037,         1.4051166851
#beta     :         1.8970085121,         1.8229532590
#coef:         0.5338084163
#Coef(M=0) = 
#0:         0.1423739499,         0.0485965529
#1:         0.0000000000,         0.0000000000
#2:         0.2847478997,         0.0971931059
#Cs(M=0,fixed):         0.5963746066,         0.2035607649
#Coef(M=1) = 
#0:        -0.2247663062,        -0.1677234495
#1:         0.0000000000,         0.0000000000
#2:         0.2247663062,         0.1677234495
#Cs(M=1,fixed):        -1.8829978037,        -1.4051166851

check_hatom: ${BINDIR}/two_pot
	../one_eig/bin/fast/one_eig hatom0.json
	./$^ hatom.json

calc_hf:
	../one_eig/bin/fast/one_eig h2_oneeig.json
	../rhf/bin/fast/rhf h2_rhf.json

calc_hf_small:
	../one_eig/bin/fast/one_eig h2_oneeig.json
	../rhf/bin/fast/rhf h2_rhf_small.json

check_stex: ${BINDIR}/two_pot
	./$^ h2.json

check_rpa: ${BINDIR}/two_pot
	@echo RPA
	./$^ h2_rpa_small.json
#oscillator strength
#(0.99796,0.00137469), (0.934428,0.000519194)
#(L,m)11
#(L,m)1-1
#(L,m)10
#[calc_alpha] Wed Apr 26 08:43:53 2017
#Cs(total,alpha):        -1.9655444218,        -1.7838145039
#Cs(sigu,alpha):         -0.3696543025,        -0.2373796170
#Cs(piu,alpha):          -1.5958901193,        -1.5464348869
#1-1(-1.24925,0)(-0.614421,0)(-0.12756,0.112096)(0.7194,0)(0.0597706,-0.132814)
#10(-1.24925,0)(-0.763776,0)(0.188447,-0.150676)(0.496396,0)(-0.106317,0.145656)
#11(-1.24925,0)(-0.614421,0)(-0.12756,0.112096)(0.7194,0)(0.0597706,-0.132814)
#coefficient A_lm
#1, -1,         0.3518934984,         0.5396341329
#1, 0,         0.0872689310,         0.5030735464
#1, 1,         0.3518934984,         0.5396341329
#Cs(total):         2.4389673450,         1.6213292625
#Cs(sig_u):         0.5829256387,         0.1975596425
#Cs(pi_u) :         1.8560417063,         1.4237696200
#beta     :         1.9214906466,         1.8069212818
#coef:         0.5338084163
#Coef(M=0) = 
#0:         0.1391632453,         0.0471638905
#1:         0.0000000000,         0.0000000000
#2:         0.2783264906,         0.0943277810
#Cs(M=0,fixed):         0.5829256387,         0.1975596425
#Coef(M=1) = 
#0:        -0.2215486591,        -0.1699499796
#1:         0.0000000000,         0.0000000000
#2:         0.2215486591,         0.1699499796
#Cs(M=1,fixed):        -1.8560417063,        -1.4237696200


check_propagator: ${BINDIR}/two_pot
	@echo propagator
	./$^ h2_propagator.json

check_propagator_ao: ${BINDIR}/two_pot
	@echo propagator_ao
	./$^ h2_propagator_ao.json

## ==== Basic command ====
clean:
	rm -fr ${BINDIR}

