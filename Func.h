#pragma once
#include "Diff.h"

namespace Diff {

	// a wrapper
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
