CXXFLAGS = -Wall -Werror -g -std=c++14
ROOT = $(shell pwd)

linalg/matrix_test:
	@mkdir -p bin
	g++ $(CXXFLAGS) linalg/matrix_test.cc -I $(ROOT) -o bin/matrix_test
	bin/matrix_test

linalg/cholesky_test:
	@mkdir -p bin
	g++ $(CXXFLAGS) linalg/cholesky_test.cc -I $(ROOT) -o bin/cholesky_test
	bin/cholesky_test

liegroups/so3.o: liegroups/so3.cc
	@mkdir -p build/liegroups
	g++ $(CXXFLAGS) -c liegroups/so3.cc -I $(ROOT) -o build/liegroups/so3.o

liegroups/so3_test: liegroups/so3.o
	@mkdir -p bin
	g++ $(CXXFLAGS) build/liegroups/so3.o liegroups/so3_test.cc -I $(ROOT) -o bin/so3_test
	bin/so3_test

test: linalg/matrix_test linalg/cholesky_test liegroups/so3_test

.PHONY: clean
clean:
	$(RM) -r build bin
