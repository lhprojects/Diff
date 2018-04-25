#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <string>
#include <set>
#include <vector>

struct DExprImpl;
struct DVariable;
struct DVariableImpl;

// reference count information
// key: the pointer to Expression object
// value: how many references of this object
extern std::set<DExprImpl*> DCount;

// a handler of Expression object
struct DExpr
{
	// bind nothing
	DExpr() : fImpl(nullptr) { }
	DExpr(DExpr const &s);
	DExpr(DExpr &&s);
	DExpr const &operator=(DExpr const &s);
	DExpr const &operator=(DExpr &&s);
	~DExpr();

	// evaluate the value
	double V() const;
	double operator()() const { return V(); }
	// return the value of last call of V()
	// you must ensure that any relevant variable not changed
	double W() const;
	// make a Expression object reprent the differential respect` s`
	// `s` must be Variable object
	DExpr D(DExpr const &s) const;
	DExpr D(DVariable const &s) const;

	// number of nodes
	// a hits for compuation time
	size_t Nodes() const;

	// find all variables
	// not optimized, very slow
	// use it for debug
	std::vector<DVariable> GetVariablesList() const;

	// make a variable to constant
	DExpr FixVariable(DVariable const &s) const;

	// return string for debug
	// all sub expression will expanded
	// maybe it cause out of memory
	std::string ToString() const;

//private:
	DExpr(DExprImpl const*impl);
	DExprImpl const *fImpl;
};

// a handler of Constant object
struct DConst : DExpr
{
	// bind a constant
	DConst(double v);

	DConst(DConst const &s) : DExpr(s) { }
	DConst(DConst &&s) : DExpr(std::move(s)) { }
	DConst const &operator=(DConst const &s) = delete;
	DConst const &operator=(DConst &&s) = delete;
};


// a handler of Variable object
struct DVariable : DExpr
{
	// bind a variable with value of `v` name fo `name`
	// `name` is only used for refining the output of `ToString`
	DVariable(double v);
	DVariable(std::string const &name, double v);
	DVariable(DVariableImpl const *p);
	DVariable(DVariable const &s) : DExpr(s) { }
	DVariable(DVariable &&s) : DExpr(std::move(s)) { }
	DVariable const &operator=(DVariable const &s) = delete;
	DVariable const &operator=(DVariable &&s) = delete;

	// get name
	std::string const &GetName() const;
	// reset the value
	void SetV(double v);
};

DExpr operator*(DExpr const &s1, DExpr const &s2);
DExpr operator+(DExpr const &s1, DExpr const &s2);
DExpr operator-(DExpr const &s1, DExpr const &s2);
DExpr operator/(DExpr const &s1, DExpr const &s2);
DExpr sqrt(DExpr const &s);
DExpr log(DExpr const &s);
DExpr pow(DExpr const &s, double n);
DExpr exp(DExpr const &s);
DExpr sin(DExpr const &s);
DExpr cos(DExpr const &s);
DExpr sinh(DExpr const &s);
DExpr cosh(DExpr const &s);
inline DExpr POW2(DExpr const &s) { return pow(s, 2); }
inline DExpr POW4(DExpr const &s) { return pow(s, 4); }
inline DExpr pow2(DExpr const &s) { return pow(s, 2); }
inline DExpr pow4(DExpr const &s) { return pow(s, 4); }

inline DExpr operator+(DExpr const &s1, double s2) { 
	return s1 + DConst(s2);
}
inline DExpr operator+(double s1, DExpr const &s2) { 
	return DConst(s1) + s2;
}
inline DExpr operator-(DExpr const &s1, double s2) { 
	return s1 - DConst(s2);
}
inline DExpr operator-(double s1, DExpr const &s2) { 
	return DConst(s1) - s2;
}
inline DExpr operator*(DExpr const &s1, double s2) { 
	return s1 * DConst(s2);
}
inline DExpr operator*(double s1, DExpr const &s2) { 
	return DConst(s1) * s2;
}
inline DExpr operator/(DExpr const &s1, double s2) { 
	return s1 / DConst(s2);
}
inline DExpr operator/(double s1, DExpr const &s2) { 
	return DConst(s1) / s2;
}

inline DExpr operator-(DExpr const &s2) {
	return 0.0 - s2;
}

