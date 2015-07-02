PYTHON_INC=`python-config --includes`
JAVA_HOME=/usr/lib/jvm/java-7-openjdk
JAVA_INC=-I$(JAVA_HOME)/include/ -I$(JAVA_HOME)/include/linux/
#TARGET=-python
OS=LINUX
CXX=g++
_mcmc.so:mcmc_py_wrap.cpp
	$(CXX) -shared $< -I../../ $(PYTHON_INC) -o _mcmc.so -fPIC -g

mcmc_py_wrap.cpp:mcmc.i
	swig -o $@ -c++ -python -Wall -I../../ mcmc.i 


clean:
	rm -f mcmc_py_wrap.*
	rm -f *~
#	rm -f *.java
#	rm -f *.class
#	rm -f *.cs
#	rm -f *.py


