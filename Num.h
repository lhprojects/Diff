#pragma once

#include <math.h>
#include <cmath>

namespace Diff {

	/* Num class is like double, but can trace the calculation error */
	struct Num {

		Num(double v, double e1, double e2) { fE1 = e1; fV = v; fE2 = e2; }
		Num & operator=(Num const &) = default;
		double V() const { return fV; }
		/* up bound of error */
		double E1() const { return fE1; }
		/* squre of error, sqrt of this is more near by the real error */
		double E2() const { return fE2; }
		/* the rms of error */
		double SqrtE2() const { return std::sqrt(fE2); }
	private:
		double fV;
		double fE1;
		double fE2;
	};

	inline double intrie1(double e) {
		return (nextafter(e, INFINITY) - e) / 2;
	}

	inline double intrie2(double e) {
		double e1 = intrie1(e);
		return e1*e1;
	}
	
	inline double squr(double x) { return x*x; }

	inline Num operator+(Num const &l, Num const &r) {
		double v = l.V() + r.V();
		return Num(v, l.E1() + r.E1() + intrie1(v), l.E2() + r.E2() + intrie2(v));
	}
	
	inline Num operator-(Num const &l, Num const &r) {
		double v = l.V() - r.V();
		return Num(v, l.E1() + r.E1() + intrie1(v), l.E2() + r.E2() + intrie2(v));
	}
	
	inline Num operator*(Num const &l, Num const &r) {
		double v = l.V() * r.V();
		return Num(v,
			l.E1() * fabs(r.V()) + r.E1() * fabs(l.V())  + intrie1(v),
			l.E2() * squr(r.V()) + r.E2() * squr(l.V()) + intrie2(v));
	}

	inline Num operator/(Num const &l, Num const &r) {
		double v = l.V() / r.V();
		return Num(v,
			l.E1() / fabs(r.V()) + r.E1() * fabs(l.V() / squr(r.V())) + intrie1(v),
			l.E2() / squr(r.V()) + r.E2() * squr(l.V() / squr(r.V())) + intrie2(v));
	}

	inline Num sin(Num const &x) {
		double v = std::sin(x.V());
		double e1 = x.E1()*fabs(std::cos(x.V())) + intrie1(v);
		double e2 = x.E2()*squr(std::cos(x.V())) + intrie2(v);
		return Num(v, e1, e2);
	}

	inline Num sinh(Num const &x) {
		double v = std::sinh(x.V());
		double e1 = x.E1()*fabs(std::cosh(x.V())) + intrie1(v);
		double e2 = x.E2()*squr(std::cosh(x.V())) + intrie2(v);
		return Num(v, e1, e2);
	}

	inline Num cos(Num const &x) {
		double v = std::cos(x.V());
		double e1 = x.E1()*fabs(std::sin(x.V())) + intrie1(v);
		double e2 = x.E2()*squr(std::sin(x.V())) + intrie2(v);
		return Num(v, e1, e2);
	}

	inline Num tan(Num const &x)
	{
		double v = std::tan(x.V());
		double e1 = x.E1()*fabs(1/std::cos(x.V())) + intrie1(v);
		double e2 = x.E2()*squr(1/std::cos(x.V())) + intrie2(v);
		return Num(v, e1, e2);
	}

	inline Num cosh(Num const &x) {
		double v = std::cosh(x.V());
		double e1 = x.E1()*fabs(std::sinh(x.V())) + intrie1(v);
		double e2 = x.E2()*squr(std::sinh(x.V())) + intrie2(v);
		return Num(v, e1, e2);
	}

	inline Num pow(Num const &x, double n) {
		if (n == 0) return Num(1, 0, 0);
		else if (n == 0.5) {
			double v = std::sqrt(x.V());
			double e1 = x.E1()*fabs(0.5 / v);
			double e2 = x.E2()*squr(0.5 / v);
			return Num(v, e1, e2);
		} else if (n == 1) return x;
		else if (n == 2) {
			double v = x.V() * x.V();
			double e1 = x.E1()*fabs(n*x.V());
			double e2 = x.E2()*squr(n*x.V());
			return Num(v, e1, e2);
		}  else {
			double v = std::pow(x.V(), n);
			double e1 = x.E1()*fabs(n*std::pow(x.V(), n - 1));
			double e2 = x.E2()*squr(n*std::pow(x.V(), n - 1));
			return Num(v, e1, e2);
		}
	}

	inline Num sqrt(Num const &x) {
		return pow(x, 0.5);
	}

	inline Num log(Num const &x) {
		double v = std::log(x.V());
		double e1 = x.E1() / fabs(x.V()) + intrie1(v);
		double e2 = x.E2() / squr(x.V()) + intrie2(v);
		return Num(v, e1, e2);
	}

	inline Num exp(Num const &x) {
		double v = std::exp(x.V());
		double e1 = x.E1() * fabs(v) + intrie1(v);
		double e2 = x.E2() * squr(v) + intrie2(v);
		return Num(v, e1, e2);
	}

}
