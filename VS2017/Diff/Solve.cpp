#include "Diff.h"
#include <exception>

namespace Diff {

	double Solve(Expr const &expr, Expr const &var_,
		double xstart, double xmin, double xmax, double tol,
		double learnRate) {

		Var var = CastToVar(var_);

		if (expr.GetTypeName() != "==") {
			throw std::runtime_error("solve only equation expression");
		}

		SubExpressionVector sev;
		expr.GetSubExpressions(sev);
		Expr l = sev[0];
		Expr r = sev[1];
		Expr v = r - l;

		double x = xstart;

		RebindableExpr d;
		for (int i = 0; i < 100; ++i) {
			var.SetV(x);
			double y = v.V();
			if (fabs(y) < tol) {
				return x;
			}
			if (d.Empty()) d = v.D(var);
			double const dy_dx = d.V();
			if (isnan(dy_dx)) {
				return NAN;
			} else if (dy_dx >= 0 && x == xmin) {
				return NAN;
			} else if (dy_dx <= 0 && x == xmax) {
				return NAN;
			}

			double lastX = x;
			x = x - learnRate*y / dy_dx;
			if (x < xmin) x = xmin;
			if (x > xmax) x = xmax;
			if (x == lastX) {
				return x;
			}
		}


		return NAN;

	}

	
}
