targets=libmcmc.so

all:$(targets)


libmcmc.so:mcmc_Sampler.cpp mcmc_Sampler.h
	g++ -shared $< -I ../../ ${JAVA_INC} -o $@ -fPIC -g

clean:
	rm -f $(targets)
