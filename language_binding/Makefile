PYTHON_INC=`python-config --includes`
#JAVA_INC=-I/usr/lib/j2se/1.4/include/ -I/usr/lib/j2se/1.4/include/linux/
#TARGET=-python
OPTDIR=../
OS=LINUX
CXX=g++
python:mcmc_py_wrap.cpp
	$(CXX) -shared $< -I../ $(PYTHON_INC) -o _mcmc.so -fPIC -g

mcmc_py_wrap.cpp:mcmc.i
	swig -o $@ -c++ -python -Wall -I../ mcmc.i 


clean:
	rm -f mcmc_py_wrap.*
	rm -f *~
#	rm -f *.java
#	rm -f *.class
#	rm -f *.cs
#	rm -f *.py


