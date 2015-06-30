OBJS = equations.o common.o symbolics.o containers.o combinatorics.o matrix.o linear.o quadratic.o dimensions.o filtering.o invariants.o
CC = g++
FLAGS = -std=c++0x -g -O0

vpath %.cpp src
vpath %.h src

RGIsearch: $(addprefix bin/, $(OBJS))
	@echo Linking
	@$(CC) $(FLAGS) $^ -o $@ 

bin/%.o : %.cpp
	@mkdir -p bin
	@echo Building $@
	@${CC} ${FLAGS} -c -o $@ $<


bin/common.o: common.cpp common.h symbolics.h
bin/symbolics.o: symbolics.cpp symbolics.h common.h
bin/containers.o: containers.cpp containers.h common.h
bin/combinatorics.o: combinatorics.cpp combinatorics.h common.h containers.h
bin/matrix.o: matrix.cpp matrix.h common.h settings.h
bin/linear.o: linear.cpp linear.h common.h
bin/quadratic.o: quadratic.cpp quadratic.h common.h
bin/dimensions.o: dimensions.cpp dimensions.h common.h containers.h linear.h settings.h
bin/filtering.o: filtering.cpp filtering.h common.h containers.h linear.h settings.h
bin/invariants.o: invariants.cpp common.h symbolics.h containers.h settings.h
bin/equations.o: equations.cpp invariants.h common.h symbolics.h containers.h combinatorics.h matrix.h quadratic.h dimensions.h filtering.h 

clean: 
	@-rm -f bin/*.o 
	@-rmdir bin 
	@-rm -f RGIsearch
