
all: bin/TestDiff bin/TestDiffToCCode


X_Ws_D0.h X_Ws_D1.h X_Ws_D2.h X_Ws_D3.h X_Ws_D4.h X_Ws_D5.h X_Ws_D6.h X_Ws_D0_avx.h X_Ws_D1_avx.h X_Ws_D2_avx.h X_Ws_D3_avx.h X_Ws_D4_avx.h X_Ws_D5_avx.h X_Ws_D6_avx.h: bin/TestDiff
	bin/TestDiff

bin/TestDiffToCCode: bin/TestDiff X_Ws_D0.h X_Ws_D1.h X_Ws_D2.h X_Ws_D3.h X_Ws_D4.h X_Ws_D5.h X_Ws_D6.h X_Ws_D0_avx.h X_Ws_D1_avx.h X_Ws_D2_avx.h X_Ws_D3_avx.h X_Ws_D4_avx.h X_Ws_D4_avx.h X_Ws_D6_avx.h
	g++ --std=c++11 -L=. -mavx -D TEST_CCODE -o bin/TestDiffToCCode TestDiff.cpp Diff.cpp

bin/TestDiff: TestDiff.cpp Diff.cpp Diff.h
	mkdir -p bin
	g++ --std=c++11 -o bin/TestDiff TestDiff.cpp Diff.cpp

