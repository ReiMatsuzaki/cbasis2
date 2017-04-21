

#ifeq (\$(OS),Linux)
PY_FLAGS=-I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -shared -fPIC -lboost_python 
CXXFLAGS=-Wall -fPIC
#endif

ifeq (\$(OS), Darwin)
#PY_FLAGS=-I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -lboost_python -DPIC -bundle -fPIC -framework Python 
PY_FLAGS= -lboost_python -DPIC -bundle -fPIC -framework Python 
CXXFLAGS=-Wall
endif

coulomb_bind.so: coulomb.cpp 
	${CXX} ${PY_FLAGS} -o $@ coulomb.cpp -lgsl -lgslcblas ${CXXFLAGS} -lboost_python 

test: coulomb_bind.so
	python test.py
install: coulomb_bind.so coulomb.py
	cp coulomb.py ~/local/Python/
	cp coulomb_bind.so ~/local/Python/
