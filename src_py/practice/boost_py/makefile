include local.mk
basic.so: basic.cpp
	${CXX} -o $@ $^ ${CXXFLAGS} ${PY_FLAGS}

check:basic.so
	python -c "from basic import *; print add(1,1)"
	rm basic.so
