CXXFLAGS=-I.. -Wall -O0

timer.o: timer.hpp timer.cpp

calc.o: calc.cpp 

calc_cpp: calc.o timer.o
	${CXX} -o $@ calc.o ../l2.a timer.o ${CXXFLAGS} 

test: calc_cpp
	echo '----------------------'
	./calc_cpp
	echo '----------------------'
	python calc.py
	echo '----------------------'
