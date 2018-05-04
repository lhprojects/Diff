#pragma once

#include <stdint.h>
#include <exception>
#include <stdexcept>
#include <vector>
#include <stdlib.h>

namespace Diff {


	template<class T, size_t N>
	struct SmallVector
	{
	

		SmallVector() : fPtr((T*)fLocal), fN(0), fRes(N) {
		}

		void push_back(T const &v) {
			prepareForPush();
			new(fPtr + fN - 1) T(v);
		}

		void push_back(T &&v) {
			prepareForPush();
			new(fPtr + fN - 1) T(std::move(v));
		}

		void clear() {
			// never shink
			for (size_t i = fN; i != 0; --i) {
				fPtr[i - 1].~T();
			}
			fN = 0;			
		}

		~SmallVector() {
			clear();
			if (fPtr != (T*)fLocal) {
				free((void*)fPtr);
			}
			fPtr = (T*)fLocal;
			fN = 0;
			fRes = N;
		}

		T const & operator[](size_t i) const {
			return fPtr[i];
		}
		T & operator[](size_t i) {
			return fPtr[i];
		}
		T const & at(size_t i) const {
			if (i >= fN) throw std::out_of_range("SmallVector::at()");
			return fPtr[i];
		}
		T & at(size_t i) {
			if (i >= fN) throw std::out_of_range("SmallVector::at()");
			return fPtr[i];
		}
		T const* begin() const { return fPtr; }
		T const* end() const { return fPtr + fN; }
		T* begin() { return fPtr; }
		T* end() { return fPtr + fN; }
		size_t size() const { return fN; }
	private:
		void prepareForPush() {
			if (fN + 1 > fRes) {
				if (fPtr == (T*)fLocal) {
					fRes = N >= 7 ? N : 7;
				} else {
					fRes = fN * 2 - 1;
				}
				void *ptr = malloc(fRes * sizeof(T));
				for (size_t i = 0; i < fN; ++i) {
					new ((char*)ptr + sizeof(T)*i) T(std::move(fPtr[i]));
				}
				clear();
				fPtr = (T*)ptr;
			}
			fN += 1;
		}

		size_t fRes;
		size_t fN;
		T *fPtr;
		alignas(T) char fLocal[sizeof(T)*N];
	};
}
