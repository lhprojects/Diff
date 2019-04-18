#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include "Num.h"
#include "Constants.h"
#include "SamllVector.h"

#include <string>
#include <set>
#include <vector>
#include <map>
#include <stdint.h>
#include <tuple>
#include <utility>

namespace Diff {

	struct Expr;
	struct DExprImpl;
	struct DConstant;
	struct DVariableImpl;
	struct Var;
	struct CCode;
	typedef SmallVector<Expr, 2> SubExpressionVector;
	typedef SmallVector<double, 2> ParameterVector;


	// a handler of Expression object
	struct Expr
	{
		// bind nothing, it is empty
		Expr() : fImpl(nullptr) { }
		Expr(Expr const &s);
		Expr(Expr &&s);
		// rebind is not allowed
		Expr &operator=(Expr const &s) = delete;
		Expr &operator=(Expr &&) = delete;
		~Expr();

		// if nothing is binded
		bool Empty() const { return fImpl == nullptr; }
		// evaluate the value
		double V() const;

		Expr D(Expr const &var) const;

		// number of nodes
		// a hits for compuation time
		size_t Nodes() const;

		// find all variables
		// not optimized, very slow
		// and some variables may be optimized out, and is not in the list
		// only use it for debug
		std::vector<Var> GetVariablesList() const;


		std::string const &GetTypeName() const;
		void GetSubExpressions(SubExpressionVector &expr) const;
		void GetParameters(ParameterVector &expr) const;

		// return string for debug
		// all sub expression will expanded
		// maybe it cause out of memory
		std::string ToString() const;

		// Unique id of an expression
		// use it as key of unordered_map/_set or map/set
		uint64_t Uid() const;
		//private:

		Expr(DExprImpl const &impl);
		DExprImpl const *fImpl;
	};


	struct RebindableExpr : Expr {
		// bind nothing
		RebindableExpr();
		RebindableExpr(Expr const &s) : Expr(s) { }
		RebindableExpr(Expr &&s) : Expr(std::move(s)) { }
		// rebind is not allowed
		RebindableExpr &operator=(Expr const &s);
		RebindableExpr &operator=(Expr &&);

	};

	// a handler of Constant object
	struct Const : Expr
	{
		// bind a constant
		Const(double v);

		Const(Const const &s) : Expr(s) { }
		Const(Const &&s) : Expr(std::move(s)) { }
		// rebind is not allowed
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


	struct ExprOrDouble : Expr
	{
		ExprOrDouble(Expr const &s) : Expr(s) { }
		ExprOrDouble(Expr &&s) : Expr(std::move(s)) { }
		ExprOrDouble(double x);
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

	Expr operator*(ExprOrDouble const &s1, ExprOrDouble const &s2);
	Expr operator+(ExprOrDouble const &s1, ExprOrDouble const &s2);
	Expr operator-(ExprOrDouble const &s1, ExprOrDouble const &s2);
	Expr operator/(ExprOrDouble const &s1, ExprOrDouble const &s2);
	inline Expr operator-(Expr const &s2) { return 0.0 - s2; }
	Expr sqrt(Expr const &s);
	Expr log(Expr const &s);
	Expr pow(Expr const &s, double n);
	Expr exp(Expr const &s);
	Expr sin(Expr const &s);
	Expr cos(Expr const &s);
	Expr tan(Expr const &s);
	Expr sinh(Expr const &s);
	Expr cosh(Expr const &s);
	inline Expr POW2(Expr const &s) { return pow(s, 2); }
	inline Expr POW4(Expr const &s) { return pow(s, 4); }
	inline Expr pow2(Expr const &s) { return pow(s, 2); }
	inline Expr pow4(Expr const &s) { return pow(s, 4); }

	// make a Expression object reprent the differential respect` s`
	// `var` must be handler of Variable object
	Expr D(Expr const &expr, Expr const &var);
	Expr D(Expr const &expr, Expr const &var, int i);
	Expr D(Expr const &expr, std::pair<Expr, int> const &pair);

	Expr Integrate(ExprOrDouble const &y, Expr const &x, ExprOrDouble const &from, ExprOrDouble const &to);
	Expr GaussLegendre64PointsIntegrate(ExprOrDouble const &y, Expr const &x, ExprOrDouble const &from, ExprOrDouble const &to);

	Expr Integrate(ExprOrDouble const &y, std::tuple<Expr, ExprOrDouble, ExprOrDouble> const &x);
	Expr GaussLegendre64PointsIntegrate(ExprOrDouble const &y, std::tuple<Expr, ExprOrDouble, ExprOrDouble> const &x);


	struct SumSecondArg
	{
		Expr fExpr;
		double fFist;
		double fSecond;
		double fInc;
		SumSecondArg(ExprOrDouble const &a, double b, double c, double d = 1) : fExpr(a), fFist(b), fSecond(c), fInc(d) { }
	};

	Expr Sum(Expr const &expr, Expr const &var, double first, double last, double inc = 1);
	Expr Sum(Expr const &expr, SumSecondArg const &arg);

	// reference count information
	// key: the pointer to Expression object
	// value: how many references of this object
	extern std::set<DExprImpl*> DCount;


	struct CCode
	{
		std::string Body;
		struct ExprLess {
			bool operator()(Expr const &l, Expr const &r) const {
				return l.Uid() < r.Uid();
			}
		};
		std::map<Expr, std::string, ExprLess> Names;
	};
	CCode ToCCode(Expr const &expr);

	// replace a variable with a expresion
	// any sub expressions of `expr` can't have a reference to the any sub expression of this
	Expr ReplaceVariable(Expr const &in_what, Expr const &to_be_replaced, ExprOrDouble const &replaced_by);
	Expr ReplaceVariable(Expr const &in_what, std::pair<Expr, ExprOrDouble> const &);

	Num VE(ExprOrDouble const &expr);

	// a functor wrapper
	struct Func1
	{
		Func1(Expr expr, Expr const &x) : fExpr(expr), fX(CastToVar(x)) {
		}

		Func1(Expr expr, Var const &x) : fExpr(expr), fX(x) {
		}

		double operator()(double x) const {
			fX.SetV(x);
			return fExpr.V();
		}
		Var const fX;
		Expr fExpr;
	};

}
