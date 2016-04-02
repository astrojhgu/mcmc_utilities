lib/liblogger.a:logger/logger.o
	mkdir -p lib && ar rv $@ $^

logger/logger.o:logger/logger.cpp logger/logger.hpp
	$(CXX) -fPIC $< -o $@ -c -std=c++11 -O3

clean:
	rm -fv logger/logger.o
