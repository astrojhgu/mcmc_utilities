all:split_jags_data split_jags_output

split_jags_data:split_jags_data.o
	$(CXX) -o $@ $< -std=c++11 -O2 -static

split_jags_data.o:split_jags_data.cpp jags_data_reader.hpp
	$(CXX) -o $@ -c $< -std=c++11 -O2

split_jags_output:split_jags_output.o
	$(CXX) -o $@ $< -std=c++11 -O2 -static

split_jags_output.o:split_jags_output.cpp jags_data_reader.hpp
	$(CXX) -o $@ -c $< -std=c++11 -O2

clean:
	rm -rf *.o
