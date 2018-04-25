
all: bin/TestDiff

bin/TestDiff: TestDiff.cpp Diff.cpp Diff.h
	mkdir -p bin
	g++ --std=c++11 -o bin/TestDiff TestDiff.cpp Diff.cpp

