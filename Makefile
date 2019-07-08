DEPS := $(shell find . -name *.d)

INC_DIRS := $(shell pwd)
INC_FLAGS := $(addprefix, -I,$(INC_DIRS))

CXXFLAGS = -I$(INC_DIRS) -Wall -Werror -g -std=c++14 -MD $(INC_FLAGS)

build/linalg/matrix_test.cc.o: linalg/matrix_test.cc
	@mkdir -p $(shell dirname $@)
	g++ $(CXXFLAGS) -c $< -o $@

bin/matrix_test: build/linalg/matrix_test.cc.o
	@mkdir -p $(shell dirname $@)
	g++ $(CXXFLAGS) $< -o  $@

.PHONY: linalg/matrix_test
linalg/matrix_test: bin/matrix_test
	bin/matrix_test

build/linalg/cholesky_test.cc.o: linalg/cholesky_test.cc
	@mkdir -p $(shell dirname $@)
	g++ $(CXXFLAGS) -c $< -o $@

bin/cholesky_test: build/linalg/cholesky_test.cc.o
	@mkdir -p $(shell dirname $@)
	g++ $(CXXFLAGS) $< -o  $@

.PHONY: linalg/cholesky_test
linalg/cholesky_test: bin/cholesky_test
	bin/cholesky_test

LIEGROUPS_SRCS = $(shell find liegroups -name *.cc ! -name *_test.cc)
LIEGROUPS_OBJS = $(addprefix build/,$(LIEGROUPS_SRCS:%.cc=%.cc.o))

build/liegroups/%.cc.o: liegroups/%.cc
	@mkdir -p $(shell dirname $@)
	g++ $(CXXFLAGS) -c $< -o $@

build/liegroups/so3_test.cc.o: liegroups/so3_test.cc
	@mkdir -p $(shell dirname $@)
	g++ $(CXXFLAGS) -c $< -o $@

bin/so3_test: build/liegroups/so3_test.cc.o $(LIEGROUPS_OBJS)
	@mkdir -p $(shell dirname $@)
	g++ $(CXXFLAGS) $^ -o $@

.PHONY: liegroups/so3_test
liegroups/so3_test: bin/so3_test
	bin/so3_test

build/liegroups/se3_test.cc.o: liegroups/se3_test.cc
	@mkdir -p $(shell dirname $@)
	g++ $(CXXFLAGS) -c $< -o $@

bin/se3_test: build/liegroups/se3_test.cc.o $(LIEGROUPS_OBJS)
	@mkdir -p $(shell dirname $@)
	g++ $(CXXFLAGS) $^ -o $@

.PHONY: liegroups/se3_test
liegroups/se3_test: bin/se3_test
	bin/se3_test

IMAGE_SRCS = $(shell find image -name *.cc ! -name *_test.cc)
IMAGE_OBJS = $(addprefix build/,$(IMAGE_SRCS:%.cc=%.cc.o))

build/image/%.cc.o: image/%.cc
	@mkdir -p $(shell dirname $@)
	g++ $(CXXFLAGS) -c $< -o $@

build/image/load_image_from_png_test.cc.o: image/load_image_from_png_test.cc
	@mkdir -p $(shell dirname $@)
	g++ $(CXXFLAGS) -c $< -o $@

bin/load_image_from_png_test: build/image/load_image_from_png_test.cc.o $(IMAGE_OBJS)
	@mkdir -p $(shell dirname $@)
	g++ $(CXXFLAGS) $^ -o $@ -lpng

.PHONY: image/load_image_from_png_test
image/load_image_from_png_test: bin/load_image_from_png_test
	bin/load_image_from_png_test data/pngkey.com-cat-silhouette-png-1527338.png

.PHONY: test
test: linalg/matrix_test linalg/cholesky_test liegroups/so3_test liegroups/se3_test image/load_image_from_png_test

.PHONY: clean
clean:
	$(RM) -r bin build

.PHONY: install_deps
install_deps:
	sudo apt install libpng-dev

-include $(DEPS)
