CXXFLAGS = -Wall -Werror -g -std=c++14

#build/manifolds/so3.o: src/manifolds/so3.cc
#	@mkdir -p build/manifolds
#	g++ $(CXXFLAGS) -c src/manifolds/so3.cc -I include -o build/manifolds/so3.o
#
matrix_test:
	@mkdir -p bin
	g++ $(CXXFLAGS) linalg/matrix_test.cc -I linalg -o bin/matrix_test
	bin/matrix_test

cholesky_test:
	@mkdir -p bin
	g++ $(CXXFLAGS) linalg/cholesky_test.cc -I linalg -o bin/cholesky_test
	bin/cholesky_test

#build/liegroups/so3.o: liegroups/so3.cc
#    @mkdir -p build/liegroups
#    g++ $(CXXFLAGS) -c liegroups/so3.cc -I liegroups/so3.hh -o build/liegroups/so3.o
#
#so3_test: build/liegroups/so3.o
#	@mkdir -p bin
#	g++ $(CXXFLAGS) build/manifolds/so3.o tests/manifolds_test.cc -I include -o bin/manifolds_test
#	bin/manifolds_test
#
test: matrix_test cholesky_test

.PHONY: clean
clean:
	$(RM) -r build bin
