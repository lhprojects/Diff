#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <string>
#include <set>
#include <vector>

namespace Diff {

	struct DExprImpl;
	struct DConstant;
	struct DVariableImpl;
	struct Var;

	// reference count information
	// key: the pointer to Expression object
	// value: how many references of this object
	extern std::set<DExprImpl*> DCount;

	// a handler of Expression object
	struct Expr
	{
		// bind nothing
		Expr() : fImpl(nullptr) { }
		Expr(Expr const &s);
		Expr(Expr &&s);
		Expr const &operator=(Expr const &s);
		Expr const &operator=(Expr &&s);
		~Expr();

		// evaluate the value
		double V() const;
		// make a Expression object reprent the differential respect` s`
		// `s` must be handler of Variable object
		Expr D(Expr const &s) const;


		// same as V()
		double operator()() const { return V(); }

		// return the value of last call of V()
		// please call V() at once
		double W() const;

		// same as D(Expr const &)
		Expr D(Var const &s) const;

		// number of nodes
		// a hits for compuation time
		size_t Nodes() const;

		// find all variables
		// not optimized, very slow
		// and some variables may be optimized out, and is not in the list
		// only use it for debug
		std::vector<Var> GetVariablesList() const;

		// make a variable to constant
		Expr FixVariable(Var const &s) const;
		// make a variable to constant
		Expr ReplaceVariable(Var const &s, Expr const &expr) const;

		// return string for debug
		// all sub expression will expanded
		// maybe it cause out of memory
		std::string ToString() const;

		//private:
		Expr(DExprImpl const &impl);
		DExprImpl const *fImpl;
	};

	// a handler of Constant object
	struct Const : Expr
	{
		// bind a constant
		Const(double v);

		Const(Const const &s) : Expr(s) { }
		Const(Const &&s) : Expr(std::move(s)) { }
		Const &operator=(Const const &s) = delete;
		Const &operator=(Const &&s) = delete;
	private:
		friend Const CastToConst(Expr const &);
		friend Const const &Zero();
		friend Const const &One();
		friend Const const &Two();
		Const(DConstant const &p);
	};

	// a handler of Variable object
	struct Var : Expr
	{
		// bind a variable with value of `v` name fo `name`
		// `name` is only used for refining the output of `ToString`
		Var(double v);
		Var(std::string const &name, double v);
		// bind same variable of s
		Var(Var const &s) : Expr(s) { }
		Var(Var &&s) : Expr(std::move(s)) { }
		// rebind is not allowed
		Var &operator=(Var const &s) = delete;
		Var &operator=(Var &&s) = delete;

		// get name
		std::string const &GetName() const;
		// reset the value
		void SetV(double v) const;
	private:
		friend Var CastToVar(Expr const &);
		friend struct Expr;
		Var(DVariableImpl const &p);
	};


	// reflection functions
	bool IsVar(Expr const &);
	bool IsConst(Expr const &);

	Var CastToVar(Expr const &);
	Const CastToConst(Expr const &);

	template<class T>
	T Cast(Expr const &);

	template< > inline Expr Cast<Expr>(Expr const &expr) { return expr; }
	template< > inline Var Cast<Var>(Expr const &expr) { return CastToVar(expr); }
	template< > inline Const Cast<Const>(Expr const &expr) { return CastToConst(expr); }

	Expr operator*(Expr const &s1, Expr const &s2);
	Expr operator+(Expr const &s1, Expr const &s2);
	Expr operator-(Expr const &s1, Expr const &s2);
	Expr operator/(Expr const &s1, Expr const &s2);
	Expr sqrt(Expr const &s);
	Expr log(Expr const &s);
	Expr pow(Expr const &s, double n);
	Expr exp(Expr const &s);
	Expr sin(Expr const &s);
	Expr cos(Expr const &s);
	Expr sinh(Expr const &s);
	Expr cosh(Expr const &s);
	inline Expr POW2(Expr const &s) { return pow(s, 2); }
	inline Expr POW4(Expr const &s) { return pow(s, 4); }
	inline Expr pow2(Expr const &s) { return pow(s, 2); }
	inline Expr pow4(Expr const &s) { return pow(s, 4); }
	// int_from^to dx y
	Expr Integrate(Expr const &x, Expr const &from, Expr const &to, Expr const &y);
	Expr Integrate(Var const &x, Expr const &from, Expr const &to, Expr const &y);

	inline Expr operator+(Expr const &s1, double s2) {
		return s1 + Const(s2);
	}
	inline Expr operator+(double s1, Expr const &s2) {
		return Const(s1) + s2;
	}
	inline Expr operator-(Expr const &s1, double s2) {
		return s1 - Const(s2);
	}
	inline Expr operator-(double s1, Expr const &s2) {
		return Const(s1) - s2;
	}
	inline Expr operator*(Expr const &s1, double s2) {
		return s1 * Const(s2);
	}
	inline Expr operator*(double s1, Expr const &s2) {
		return Const(s1) * s2;
	}
	inline Expr operator/(Expr const &s1, double s2) {
		return s1 / Const(s2);
	}
	inline Expr operator/(double s1, Expr const &s2) {
		return Const(s1) / s2;
	}

	inline Expr operator-(Expr const &s2) {
		return 0.0 - s2;
	}

}
