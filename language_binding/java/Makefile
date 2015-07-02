PYTHON_INC=`python-config --includes`
JAVA_HOME=/usr/lib/jvm/java-7-openjdk
JAVA_INC=-I$(JAVA_HOME)/include/ -I$(JAVA_HOME)/include/linux/
#TARGET=-python
OS=LINUX
CXX=g++
all:libmcmc.so classes

libmcmc.so:mcmc_java_wrap.cpp
	$(CXX) -shared $< -I../../ $(JAVA_INC) -o libmcmc.so -fPIC -g

mcmc_java_wrap.cpp:mcmc.i
	mkdir -p mcmc && swig -o $@ -c++ -java -Wall -package mcmc -outdir mcmc/ -I../../ mcmc.i 


classes:
	cd mcmc/; javac *.java

clean:
	rm -f mcmc_java_wrap.*
	rm -f *~
	rm -f *.java
	rm -f *.class
#	rm -f *.cs
#	rm -f *.py


