
all: lib/Diff.a bin/TestDiff bin/TestDiffToCCode

lib/Diff.a: src/*.cpp include/*.h
	g++ --std=c++11 -O2 -I ./include -mavx src/*.cpp -c
	mkdir -p lib
	ar -crv lib/Diff.a *.o
bin/TestDiffToCCode: bin/TestDiff lib/Diff.a
	g++ --std=c++11 -O2 -I ./include -I. -L=. -mavx -D TEST_CCODE -o bin/TestDiffToCCode test/TestDiff.cpp lib/Diff.a

bin/TestDiff: lib/Diff.a 
	mkdir -p bin
	g++ --std=c++11 -O2 -I ./include -L=.                     -o bin/TestDiff        test/TestDiff.cpp lib/Diff.a

