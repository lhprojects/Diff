
all: bin/TestDiff bin/TestDiffToCCode


bin/TestDiffToCCode: bin/TestDiff
	g++ --std=c++11 -O2 -L=. -mavx -D TEST_CCODE -o bin/TestDiffToCCode TestDiff.cpp Diff.cpp Quad.cpp ReplaceVariable.cpp ToCCode.cpp

bin/TestDiff: TestDiff.cpp Diff.cpp Quad.cpp Diff.h Quad.h Func.h ToCCode.cpp ReplaceVariable.cpp
	mkdir -p bin
	g++ --std=c++11 -O2 -o bin/TestDiff TestDiff.cpp Diff.cpp Quad.cpp ToCCode.cpp ReplaceVariable.cpp

