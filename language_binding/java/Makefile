targets=libmcmc.so

all:$(targets)


libmcmc.so:mcmc_Sampler.cpp mcmc_Sampler.h
	g++ -shared $< -I ../../ ${JAVA_INC} -o $@ -fPIC -g -O2

clean:
	rm -f $(targets)
	rm -f `find . -iname '*.class'`
