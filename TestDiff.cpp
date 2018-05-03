#include "Diff.h"
#include "Func.h"
#include "Num.h"
#include "Quad.h"
#include <chrono>
#include <math.h>
#include <cmath>
#include <random>
#include <algorithm>
using namespace Diff;
//#define TEST_CCODE


#define TPrintf printf
const double PI = 3.1415926535897932384626433;

static int n_total = 0;
static int n_failed = 0;
#define TEST_SAME(x, y) do { n_total++; double x_ = (x); double y_ = (y);\
	if(isnan(x_) || fabs(x_ - y_) > abs(y_)*1E-3)  { n_failed++; TPrintf("FAILED: %s: %3d:"  #x " (%f) == " #y " (%f)\n", __func__, __LINE__, x_, y_); }   } while(0)
#define TEST_TRUE(x) do { n_total++; double x_ = (x); if(!x_)  { n_failed++; TPrintf("FAILED: %s: %3d:"  #x "\n", __func__, __LINE__);}   } while(0)

void test_V() {

	{
		Var three = 3;
		TEST_SAME(three.V(), 3);
		TEST_SAME((three * three).V(), 9);
		TEST_SAME((three + three).V(), 6);
		TEST_SAME((three - three).V(), 0);
		TEST_SAME((three / three).V(), 1);
		TEST_SAME(log(three).V(), log(3));

		TEST_SAME(three.VE().V(), 3);
		TEST_SAME((three * three).VE().V(), 9);
		TEST_SAME((three + three).VE().V(), 6);
		TEST_SAME((three - three).VE().V(), 0);
		TEST_SAME((three / three).VE().V(), 1);
		TEST_SAME(log(three).VE().V(), log(3));

		printf("log(3) %e", log(three).VE().E1());
	}

}
void test_Diff() {



	{ // basics test
		Var x = 4;

		TEST_SAME(D(x, x).V(), 1);
		TEST_SAME(D(x * x, x).V(), 2 * x.V());
		TEST_SAME(D(x + x, x).V(), 2);
		TEST_SAME(D(x - x, x).V(), 0);
		TEST_SAME(D(x / x, x).V(), 0);
		TEST_SAME(D(sqrt(x), x).V(), 0.5 / sqrt(x.V()));
		TEST_SAME(D(sqrt(x), {x, 2}).V(), -0.25 / pow(sqrt(x.V()), 3));
		TEST_SAME(D(sqrt(x), {x, 3}).V(), 3 * 0.125 / pow(sqrt(x.V()), 5));
		TEST_SAME(D(log(x), x).V(), 1 / x.V());
		TEST_SAME(D(pow(x, 2), x).V(), 2 * x.V());
	}
	
	{ // sqrt
	
		Var x = 4;
		TEST_SAME(D(sqrt(x), { x, 1 }).V(), 0.5 / sqrt(x.V()));
		TEST_SAME(D(sqrt(x), { x, 2 }).V(), -0.25 / pow(sqrt(x.V()), 3));
		TEST_SAME(D(sqrt(x), { x, 3 }).V(), 3 * 0.125 / pow(sqrt(x.V()), 5));
	}

	{ // mul
		Var x = 4;
		TEST_SAME(D(x*x*x, {x, 2}).V(), 3 * 2 * x.V());
		TEST_SAME(D(x*x*x*x, { x, 3 }).V(), 4 * 3 * 2 * x.V());
	}

	{ // exp
		Var x("x", 4);
		TEST_SAME(exp(2 * x).V(), exp(2 * x.V()));
		TEST_SAME(D(exp(2 * x), x).V(), 2 * exp(2 * x.V()));
		TEST_SAME(D(exp(2 * x), { x, 2 }).V(), 2 * 2 * exp(2 * x.V()));
		TEST_SAME(D(exp(2 * x), { x, 3 }).V(), 2 * 2 * 2 * exp(2 * x.V()));

		TPrintf("exp(2x)          : %s\n", D(exp(2 * x), x).ToString().c_str());
		TPrintf("exp(2x)'         : %s\n", D(exp(2 * x), x).ToString().c_str());
		TPrintf("exp(2x)''        : %s\n", D(exp(2 * x), x).ToString().c_str());
	}

	{ // sin cos
		Var x("x", 4);
		TEST_SAME(sin(2 * x).V(), sin(2 * x.V()));
		TEST_SAME(sin(2 * x).D(x).V(), 2 * cos(2 * x.V()));
		TEST_SAME(sin(2 * x).D(x).D(x).V(), -4 * sin(2 * x.V()));

		TEST_SAME(sinh(2 * x).V(), sinh(2 * x.V()));
		TEST_SAME(sinh(2 * x).D(x).V(), 2 * cosh(2 * x.V()));
		TEST_SAME(sinh(2 * x).D(x).D(x).V(), 4 * sinh(2 * x.V()));

		TEST_SAME(cos(2 * x).V(), cos(2 * x.V()));
		TEST_SAME(cos(2 * x).D(x).V(), -2 * sin(2 * x.V()));
		TEST_SAME(cos(2 * x).D(x).D(x).V(), -4 * cos(2 * x.V()));

		TEST_SAME(cosh(2 * x).V(), cosh(2 * x.V()));
		TEST_SAME(cosh(2 * x).D(x).V(), 2 * sinh(2 * x.V()));
		TEST_SAME(cosh(2 * x).D(x).D(x).V(), 4 * cosh(2 * x.V()));

		TPrintf("sin(2x)          : %s\n", sin(2 * x).ToString().c_str());
		TPrintf("sin(2x)'         : %s\n", sin(2 * x).D(x).ToString().c_str());
		TPrintf("sin(2x)''        : %s\n", sin(2 * x).D(x).D(x).ToString().c_str());
		TPrintf("sinh(2x)         : %s\n", sinh(2 * x).ToString().c_str());
		TPrintf("sinh(2x)'        : %s\n", sinh(2 * x).D(x).ToString().c_str());
		TPrintf("sinh(2x)''       : %s\n", sinh(2 * x).D(x).D(x).ToString().c_str());
		TPrintf("cos(2x)          : %s\n", cos(2 * x).ToString().c_str());
		TPrintf("cos(2x)'         : %s\n", cos(2 * x).D(x).ToString().c_str());
		TPrintf("cos(2x)''        : %s\n", cos(2 * x).D(x).D(x).ToString().c_str());
		TPrintf("cosh(2x)         : %s\n", cosh(2 * x).ToString().c_str());
		TPrintf("cosh(2x)'        : %s\n", cosh(2 * x).D(x).ToString().c_str());
		TPrintf("cosh(2x)''       : %s\n", cosh(2 * x).D(x).D(x).ToString().c_str());
	}
	{ // sqrt
		Var p("p", 1);
		Const mass = 1;
		Expr energy = sqrt(mass*mass + p*p);
		TEST_SAME(energy.V(), sqrt(2));
		TEST_SAME(energy.D(p).V(), (p / energy).V());
		TPrintf("sqrt(m^2+p^2)        : %s\n", energy.ToString().c_str());
		TPrintf("sqrt(m^2+p^2)'       : %s\n", energy.D(p).ToString().c_str());
	}

	{ // constants fold
		Const mass = 1;
		Const p = 1;
		Expr energy = exp(pow(log(1 + cos(sin(sqrt(mass*mass + p*p)) - 1)), 2));
		TEST_SAME(energy.V(), exp(pow(log(1 + cos(sin(sqrt(mass*mass + p*p)) - 1)).V(), 2)));
		TPrintf("exp(pow(log(1+cos(sin(sqrt(mass*mass + p*p))-1)), 2))        : %s\n", energy.ToString().c_str());
	}

	{ // numerical test
		Const mass = 1;
		Var p = 1;
		Expr energy = exp(pow(log(1 + cos(sin(sqrt(mass*mass + p*p)) - 1)), 2));

		p.SetV(0.9999);
		double e1 = energy.V();
		p.SetV(1.0001);
		double e2 = energy.V();
		p.SetV(1);
		TEST_SAME(energy.D(p).V(), (e2-e1)/0.0002);
	}

	{
		
		Var x("x", 2);
		TPrintf("x+x          : %s\n", (x + x).ToString().c_str());
		TPrintf("x-x          : %s\n", (x - x).ToString().c_str());
		TPrintf("x*x          : %s\n", (x * x).ToString().c_str());
		TPrintf("x/x          : %s\n", (x / x).ToString().c_str());
		TPrintf("x+1+x          : %s\n", (x + 1 + x).ToString().c_str());
		TPrintf("x+1-x          : %s\n", (x + 1 - x).ToString().c_str());
		TPrintf("x+1*x          : %s\n", (x + 1 * x).ToString().c_str());
		TPrintf("x*1/x          : %s\n", (x * 1 / x).ToString().c_str());

		TPrintf("x^2          : %s\n", POW2(x).ToString().c_str());
		TPrintf("x^2'         : %s\n", POW2(x).D(x).ToString().c_str());
		TPrintf("x^2''        : %s\n", POW2(x).D(x).D(x).ToString().c_str());
		TPrintf("x^2'''        : %s\n", POW2(x).D(x).D(x).D(x).ToString().c_str());

		TPrintf("x^2          : %s\n", (x*x).ToString().c_str());
		TPrintf("x^2'         : %s\n", (x*x).D(x).ToString().c_str());
		TPrintf("x^2''        : %s\n", (x*x).D(x).D(x).ToString().c_str());
		TPrintf("x^2'''        : %s\n", (x*x).D(x).D(x).D(x).ToString().c_str());

		TPrintf("sqrt(x)          : %s\n", sqrt(x).ToString().c_str());
		TPrintf("sqrt(x)'         : %s\n", sqrt(x).D(x).ToString().c_str());
		TPrintf("sqrt(x)''        : %s\n", sqrt(x).D(x).D(x).ToString().c_str());
		TPrintf("sqrt(x)'''       : %s\n", sqrt(x).D(x).D(x).D(x).ToString().c_str());
		TPrintf("sqrt(x)''''      : %s\n", sqrt(x).D(x).D(x).D(x).D(x).ToString().c_str());

		TPrintf("x+x          : %s\n", (x + x).ToString().c_str());
		TPrintf("x/sqrt(1+x)          : %s\n", (x / sqrt(1 + x)).ToString().c_str());
		TPrintf("x/sqrt(1+x)'         : %s\n", (x / sqrt(1 + x)).D(x).ToString().c_str());
	}

	for (auto p : DCount) {
		printf("live object: %p\n", (void*)p);
	}
	TEST_TRUE(DCount.size() <= 3);

}

void fix_var() {

	{
		Const mass = 1;
		Var p("p", 1);
		Expr energy = exp(pow(log(1 + cos(sin(sqrt(mass*mass + 1*p*p)) - 1)), 2) / 2);

		auto energy_vars = energy.GetVariablesList();
		for (auto &v : energy_vars) {
			printf("var in energy %s\n", v.GetName().c_str());
		}
		TEST_TRUE(energy_vars.size() == 1);

		Expr e = energy.FixVariable(p);
		auto e_vars = e.GetVariablesList();
		for (auto &v : e_vars) {
			printf("var in energy %s\n", v.GetName().c_str());
		}
		TEST_TRUE(e_vars.size() == 0);

		TEST_SAME(e.V(), energy.V());

		TPrintf("exp(pow(log(1 + cos(sin(sqrt(mass*mass + p*p)) - 1)), 2)/2) fix var: %s\n", e.ToString().c_str());

	}

	for (auto p : DCount) {
		printf("live object: %p\n", (void*)p);
	}
	TEST_TRUE(DCount.size() <= 3);

}
inline double POW2(double x) { return x*x; }
inline double POW4(double x) { return x*x*x*x; }

typedef double Float;
struct Constants
{
	Constants(Float sqrts_,
		Float M_H_,
		Float M_Z_,
		Float w_Z_,
		Float M_W_,
		Float sinThetaW2_,
		Float G_F_,
		Float ZDecayToNeutrinoBranchingRatio) :
		sqrts(sqrts_),
		s(sqrts_*sqrts_),
		m_H(M_H_),
		m_Z(M_Z_),
		m_W(M_W_),
		w_Z(w_Z_),
		a_e(-1),
		sinThetaW2(sinThetaW2_),
		v_e(-1 + 4 * sinThetaW2_),
		G_F(G_F_),
		Br_Z_v(ZDecayToNeutrinoBranchingRatio)
	{
	}
	Float const G_F;
	Float const m_Z;
	Float const w_Z;
	Float const m_W;
	Float const sqrts;
	Float const m_H;

	Float const sinThetaW2;
	Float const a_e;
	Float const v_e;
	Float const s;
	Float const Br_Z_v;
};


struct ConstantsFactory
{
	ConstantsFactory(double sqrts_) {
		G_F = 1.1553787E-5;
		m_Z = 91.1876;
		m_W = 80.385;
		w_Z = 2.4952;
		sqrts = sqrts_;
		m_H = 125;
		sinThetaW2 = 1 - (m_W / m_Z)*(m_W / m_Z);
		Br_Z_v = 0.2;
	}

	Constants Get() {
		return Constants(sqrts, m_H, m_Z, w_Z, m_W, sinThetaW2, G_F, Br_Z_v);
	}

	Float G_F;
	Float m_Z;
	Float w_Z;
	Float m_W;
	Float sqrts;
	Float m_H;

	Float sinThetaW2;
	Float Br_Z_v;

};

struct Formula {

	Var costheta;
	Var eH;
	std::vector<Expr> X_Is;
	std::vector<Expr> X_Ss;
	std::vector<Expr> X_Ws;
	std::vector<Expr> X_Is0;
	std::vector<Expr> X_Ss0;
	std::vector<Expr> X_Ws0;

	void Init(int n)
	{
		for (; (int)X_Is.size() <= n; ) {
			X_Is.push_back(D(X_Is.at(X_Is.size() - 1), costheta));
			costheta.SetV(0);
			X_Is0.push_back(X_Is.back().FixVariable(costheta));
		}
		for (; (int)X_Ws.size() <= n; ) {
			X_Ws.push_back(D(X_Ws.at(X_Ws.size() - 1), costheta));
			costheta.SetV(0);
			X_Ws0.push_back(X_Ws.back().FixVariable(costheta));
		}
		for (; (int)X_Ss.size() <= n; ) {
			X_Ss.push_back(D(X_Ss.at(X_Ss.size() - 1), costheta));
			costheta.SetV(0);
			X_Ss0.push_back(X_Ss.back().FixVariable(costheta));
		}
	}

	Formula(Constants C, double XsectionGlobalFactor) : costheta("costheta", 0), eH("eH", 140)
	{
		Const const m_Z = C.m_Z;
		Const const w_Z = C.w_Z;
		Const const sqrts = C.sqrts;
		Const const mH = C.m_H;
		Const const s = C.s;
		Expr const costheta2 = POW2(costheta);
		Expr const p = sqrt(eH*eH - mH*mH);
		Expr const epsilon_v = sqrts - eH;
		Expr const s_v = epsilon_v*epsilon_v - p*p;
		Expr const s1 = sqrts*(epsilon_v + p*costheta);
		Expr const s2 = sqrts*(epsilon_v - p*costheta);
		Expr const s1s2 = sqrts*sqrts*(POW2(epsilon_v) - POW2(p)*costheta2);
		Expr const c_x = 1 - 2 * s*s_v / (s1s2);
		Expr const s2_x = 1 - c_x*c_x;
		Expr const h1 = 1 + 2 * POW2(C.m_W) / s1;
		Expr const h2 = 1 + 2 * POW2(C.m_W) / s2;
		Expr const h1h2 = 1 + 4 * POW4(C.m_W) / s1s2 + 4 * POW2(C.m_W)*sqrts*epsilon_v / s1s2;

		Expr const r = h1*h1 + h2*h2 + 2 * c_x*h1h2 - s2_x;
		Expr const t1 = h1 + c_x*h2;
		Expr const t2 = h2 + c_x*h1;
		Expr const L = log((h1h2 + c_x + sqrt(r)) / (h1h2 + c_x - sqrt(r)));

		{
			Const const cosThetaW2 = 1 - C.sinThetaW2;
			Expr const f1 = (C.v_e + C.a_e)*POW2(cosThetaW2) / 8;
			Expr const f2 = (s_v - m_Z*m_Z) / ((s - POW2(m_Z)) * (POW2(s_v - m_Z*m_Z) + POW2(m_Z*w_Z)));
			Expr const f3 = 2 - (h1 + 1)*log((h1 + 1) / (h1 - 1)) - (h2 + 1)*log((h2 + 1) / (h2 - 1)) + (h1 + 1)*(h2 + 1)*L / sqrt(r);
			Expr G_I = f1*f2*f3;
			Expr X_I = XsectionGlobalFactor*p*G_I;
			X_Is.push_back(X_I);
			X_Is0.push_back(X_Is.back().FixVariable(costheta));
		}
		{
			Const const f1 = (POW2(C.v_e) + POW2(C.a_e)) / 96;
			Expr const f2 = (s*s_v + s1*s2) / (POW2(s - POW2(m_Z)) * (POW2(s_v - m_Z*m_Z) + POW2(m_Z*w_Z)));
			Expr const Gs = f1*f2;
			Expr X_S = XsectionGlobalFactor*p * 3 * Gs;
			X_Ss.push_back(X_S);
			X_Ss0.push_back(X_Ss.back().FixVariable(costheta));

		}
		{
			Const const cosThetaW2 = 1 - C.sinThetaW2;
			Expr const cosThetaW8 = POW4(cosThetaW2);
			Expr const f1 = cosThetaW8 / (s1s2*r);
			Expr const f2 = (h1 + 1)*(h2 + 1)*(2 / (POW2(h1) - 1) + 2 / (POW2(h2) - 1) - 6 * s2_x / r + (3 * t1*t2 / r - c_x)*L / sqrt(r));
			Expr const f3 = -(2 * t1 / (h2 - 1) + 2 * t2 / (h1 - 1) + (t1 + t2 + s2_x)*L / sqrt(r));
			Expr const G_W = f1*(f2 + f3);

			Expr X_W = XsectionGlobalFactor*p * G_W;
			X_Ws.push_back(X_W);
			X_Ws0.push_back(X_Ws.back().FixVariable(costheta));

		}

		if (1) {
			for (auto x : X_Ws.at(0).GetVariablesList()) {
				printf("X_Ws var: %s\n", x.GetName().c_str());
			}
			for (auto x : X_Is.at(0).GetVariablesList()) {
				printf("X_Is var: %s\n", x.GetName().c_str());
			}
			for (auto x : X_Ss.at(0).GetVariablesList()) {
				printf("X_Ss var: %s\n", x.GetName().c_str());
			}
			for (auto x : X_Ws0.at(0).GetVariablesList()) {
				printf("X_Ws0 var: %s\n", x.GetName().c_str());
			}
			for (auto x : X_Is0.at(0).GetVariablesList()) {
				printf("X_Is0 var: %s\n", x.GetName().c_str());
			}
			for (auto x : X_Ss0.at(0).GetVariablesList()) {
				printf("X_Ss0 var: %s\n", x.GetName().c_str());
			}
		}
	}
};

void test_for() {

	{
		ConstantsFactory cf(250);
		auto c = cf.Get();
		Formula f(c, 1);
		int n = 10;
		f.Init(10);
		f.eH.SetV(140);

		f.costheta.SetV(0.01);
		double y1 = f.X_Ws.at(0).V();
		f.costheta.SetV(0);
		double y0 = f.X_Ws.at(0).V();
		double ypp = (y1 - y0) / (0.5*0.01*0.01);

		// test fixvariable
		TEST_SAME(f.X_Ws0.at(2).V(), ypp);
		TEST_SAME(f.X_Ws.at(2).V(), ypp);
		TEST_SAME(f.X_Ws.at(2).VE().V(), f.X_Ws.at(2).V());

		for (int i = 0; i < 10; ++i) {
			printf("d%d %+e +- %e +- %e\n", i, 
				f.X_Ws.at(i).VE().V(), f.X_Ws.at(i).VE().E1(), f.X_Ws.at(i).VE().SqrtE2());
		}

		printf("nodes %d\n", (int)f.X_Ws.at(10).Nodes());
		printf("nodes %d\n", (int)f.X_Ss.at(10).Nodes());
		printf("nodes %d\n", (int)f.X_Is.at(10).Nodes());
		printf("nodes %d\n", (int)f.X_Ws.at(0).Nodes());
		printf("nodes %d\n", (int)f.X_Ss.at(0).Nodes());
		printf("nodes %d\n", (int)f.X_Is.at(0).Nodes());
		printf("nodes %d\n", (int)f.X_Ws.at(0).FixVariable(f.eH).Nodes());
		printf("nodes %d\n", (int)f.X_Ss.at(0).FixVariable(f.eH).Nodes());
		printf("nodes %d\n", (int)f.X_Is.at(0).FixVariable(f.eH).Nodes());

	}

	for (auto p : DCount) {
		printf("live object: %p\n", (void*)p);
	}
	TEST_TRUE(DCount.size() <= 3);

}

void testfunc() {
	Var x("x", 0);
	Func1 f(1 + x, x);
	TEST_SAME(f(1), 2);
}

void test_int() {

	{ // basic test
		Var x = 0;
		TEST_SAME(Integrate(x*x, { x, Const(0), Const(1) }).V(), 1. / 3);
		TEST_SAME(GaussLegendre64PointsIntegrate(x, Const(0), Const(1), x*x).V(), 1. / 3);
		printf("Integrate(x, 0, 1, x*x)-1/3: %.20f\n", Integrate(x, 0, 1, x*x).V() - 1 / 3.0);
		printf("Integrate(x, 0, 1, x*x)-1/3: %.20f\n", GaussLegendre64PointsIntegrate(x, Const(0), Const(1), x*x).V() - 1 / 3.0);
		printf("Integrate(x, 0, 1, x*x)    : %s\n", GaussLegendre64PointsIntegrate(x, Const(0), Const(1), x*x).ToString().c_str());

		TEST_SAME(Integrate(x, 0, 1, exp(-x)).V(), 1 - exp(-1));
		printf("Integrate(x, 0, 1, exp(-x))-(1 - e^-1): %+ef\n", Integrate(x, 0, 1, exp(-x)).V() - (1 - exp(-1)));
		TEST_SAME(Integrate(x, 0, 1, exp(x)).V(), exp(1) - 1);
		printf("Integrate(x, 0, 1, exp(x))-(e^-1): %+ef\n", Integrate(x, 0, 1, exp(x)).V() - (exp(1) - 1));
		TEST_SAME(Integrate(x, 0, PI, sin(x)).V(), 2);
	}

	{ // test for large range
		Var x = 0;
		printf("Integrate(x, 0,   PI, sin(x))-2: %+e\n", Integrate(x, 0, PI, sin(x)).V()-2);
		printf("Integrate(x, 0,  3PI, sin(x))-2: %+e\n", Integrate(x, 0, 3*PI, sin(x)).V() - 2);
		printf("Integrate(x, 0,  9PI, sin(x))-2: %+e\n", Integrate(x, 0, 9*PI, sin(x)).V() - 2);
		printf("Integrate(x, 0, 17PI, sin(x))-2: %+e\n", Integrate(x, 0, 17*PI, sin(x)).V() - 2);
		printf("Integrate(x, 0, 33PI, sin(x))-2: %+e\n", Integrate(x, 0, 33*PI, sin(x)).V() - 2);
		printf("Integrate(x, 0, 65PI, sin(x))-2: %+e\n", Integrate(x, 0, 65 * PI, sin(x)).V() - 2);
		printf("Integrate(x, 0,129PI, sin(x))-2: %+e\n", Integrate(x, 0, 129* PI, sin(x)).V() - 2);
	}

	{ // sqrt
		Var x = 0;
		TEST_SAME(Integrate(x, 0, 1, sqrt(x)).V(), 2 / 3.0);
		TEST_SAME(Integrate(x, 0, 1, pow(x, 3. / 2)).V(), 2 / 5.0);
		TEST_SAME(Integrate(x, 0, 1, pow(x, 5. / 2)).V(), 2 / 7.0);
		printf("Integrate(x, 0, 1, sqrt(x))-(2/3.)'   : %.20f\n", Integrate(x, 0, 1, sqrt(x)).V() - (2 / 3.0));
		printf("Integrate(x, 0, 1, pow(x,3/2))-(2/5.)': %.20f\n", Integrate(x, 0, 1, pow(x, 3. / 2)).V() - (2 / 5.0));
		printf("Integrate(x, 0, 1, pow(x,5/2))-(2/7.)': %.20f\n", Integrate(x, 0, 1, pow(x, 5. / 2)).V() - (2 / 7.0));
	}

	{ // basic D
		Var t("t", 1);
		Var x("x", 0);
		Expr y = Integrate(x, Const(0), t, x*x*t);

		Expr yopen = GaussLegendre64PointsIntegrate(x, Const(0), t, x*x*t);

		TEST_SAME(y.V(), 1 / 3.0);
		TEST_SAME(D(y, t).V(), 1 / 3.0 + 1);
		TEST_SAME(D(yopen, t).V(), 1 / 3.0 + 1);
		printf("dy/dt       - 4/3 %+e\n", D(y, t).V() - 4 / 3.);
		printf("dy/dt(open) - 4/3 %+e\n", D(yopen, t).V() - 4 / 3.);
	}

	{ // sqrt'
		Var x = 0;
		Var t = 1;
		TEST_SAME(Integrate(x, 0, t, sqrt(t - x)).V(), 2 / 3.0);
		TEST_SAME(GaussLegendre64PointsIntegrate(x, Const(0), t, sqrt(t - x)).V(), 2 / 3.0);
		TEST_TRUE(fabs(D(Integrate(x, 0, t, sqrt(t - x)), t).V() - 1) < 0.01);
		TEST_SAME(D(GaussLegendre64PointsIntegrate(x, Const(0), t, sqrt(t - x)), t).V(), 1);
		printf("Integrate(x, 0, t, sqrt(t - x))'                      - 1: %+.6e\n", D(Integrate(x, 0, t, sqrt(t - x)), t).V() - 1);
		printf("GaussLegendre64PointsIntegrate(x, 0, t, sqrt(t - x))' - 1: %+.6e\n", D(GaussLegendre64PointsIntegrate(x, Const(0), t, sqrt(t - x)), t).V() - 1);
	}

}



#ifdef TEST_CCODE
#include <immintrin.h>

double hsum(__m256d x) {
#ifdef _MSC_VER
	__declspec(align(32)) double a[4] = { };//MSVC
#else
	__attribute__((aligned(32))) double a[4] = { };//GCC, ICC
#endif
       	_mm256_store_pd(a, x);
	return a[0] + a[1] + a[2] + a[3];
}

__m256d _mm256_log_pd(__m256d x)
{
#ifdef _MSC_VER
__declspec(align(32)) double a[4] = { };//MSVC
#else
__attribute__((aligned(32))) double a[4] = { };//GCC, ICC
#endif

	_mm256_store_pd(a, x);
	a[0] = log(a[0]);
	a[1] = log(a[1]);
	a[2] = log(a[2]);
	a[3] = log(a[3]);
	return _mm256_load_pd(a);
}

__m256d _mm256_pow_pd(__m256d x,double n)
{
#ifdef _MSC_VER
__declspec(align(32)) double a[4] = { };//MSVC
#else
__attribute__((aligned(32))) double a[4] = { };//GCC, ICC
#endif

	_mm256_store_pd(a, x);
	a[0] = pow(a[0], n);
	a[1] = pow(a[1], n);
	a[2] = pow(a[2], n);
	a[3] = pow(a[3], n);
	return _mm256_load_pd(a);
}




#include "X_Ws_D0.h"
#include "X_Ws_D1.h"
#include "X_Ws_D2.h"
#include "X_Ws_D3.h"
#include "X_Ws_D4.h"
#include "X_Ws_D5.h"
#include "X_Ws_D6.h"
#include "X_Ws_D0_avx.h"
#include "X_Ws_D1_avx.h"
#include "X_Ws_D2_avx.h"
#include "X_Ws_D3_avx.h"
#include "X_Ws_D4_avx.h"
#include "X_Ws_D5_avx.h"
#include "X_Ws_D6_avx.h"
#endif


void test_code()
{
	{
		Const mass = 1;
		Var p("p", 1);
		Expr energy = exp(pow(log(1 + cos(sin(sqrt(mass*mass + p*p)) - 1)), 2));
		printf("%s\n", energy.ToCCode().Body.c_str());
	}

	{
		int NMax = 6;
		ConstantsFactory cf(250);
		Constants c = cf.Get();
		Formula f(c, 1);
		f.costheta.SetV(0);
		f.costheta.SetV(140);
		f.Init(NMax);

		for (int i = 0; i <= NMax; ++i)
		{
			{
				char b[1024];
				sprintf(b, "X_Ws_D%d.h", i);
				FILE *file = fopen(b, "w");

				CCode ccode = f.X_Ws.at(i).ToCCode();
				fprintf(file, "double X_Ws_D%d(double %s, double %s) {\n",
					i,
					ccode.Names[f.eH].c_str(),
					ccode.Names[f.costheta].c_str());
				fprintf(file, "%s", ccode.Body.c_str());
				fprintf(file, "return %s; }", ccode.Names[f.X_Ws.at(i)].c_str());

				fclose(file);
			}
			{
				char b[1024];
				sprintf(b, "X_Ws_D%d_avx.h", i);
				FILE *file = fopen(b, "w");

				CCode ccode = f.X_Ws.at(i).ToAVXCode();
				fprintf(file, "__m256d X_Ws_D%d_avx(__m256d %s, __m256d %s) {\n",
					i,
					ccode.Names[f.eH].c_str(),
					ccode.Names[f.costheta].c_str());
				fprintf(file, "%s", ccode.Body.c_str());
				fprintf(file, "return %s; }", ccode.Names[f.X_Ws.at(i)].c_str());

				fclose(file);
			}
		}

#ifdef TEST_CCODE

		TEST_SAME(X_Ws_D0(140, 0), f.X_Ws.at(0).V());
		TEST_SAME(X_Ws_D1(140, 0), f.X_Ws.at(1).V());
		TEST_SAME(X_Ws_D2(140, 0), f.X_Ws.at(2).V());
		TEST_SAME(X_Ws_D3(140, 0), f.X_Ws.at(3).V());
		TEST_SAME(X_Ws_D4(140, 0), f.X_Ws.at(4).V());
		TEST_SAME(X_Ws_D5(140, 0), f.X_Ws.at(5).V());
		TEST_SAME(X_Ws_D6(140, 0), f.X_Ws.at(6).V());

		double (*ptr[])(double, double) = { X_Ws_D0, X_Ws_D1, X_Ws_D2,
			X_Ws_D3, X_Ws_D4, X_Ws_D5, X_Ws_D6 };
		__m256d(*ptr_avx[])(__m256d, __m256d) = { X_Ws_D0_avx, X_Ws_D1_avx, X_Ws_D2_avx,
			X_Ws_D3_avx, X_Ws_D4_avx, X_Ws_D5_avx, X_Ws_D6_avx };
#endif
		double delta = 0.0001;
		int Nj = (int)(1 / delta) - 10;
		for (int i = 0; i <= NMax; ++i) {
			int Njj = Nj;
			if (i >= 3) Njj = Nj / 10;

#ifdef TEST_CCODE

			{
				auto t0 = std::chrono::high_resolution_clock::now();
				double sum = 0;
				double theta = 0;
				for (int j = 0; j < Njj; ++j) {
					sum += ptr[i](140, theta);
					theta += delta;
				}
				auto d = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - t0);
				if (i >= 3)
					printf("cpp double order %d, time: 10x%10.1fus, value: %f\n", i, (double)d.count(), sum);
				else
					printf("cpp double order %d, time: %10.1fus, value: %f\n", i, (double)d.count(), sum);
			}
			{
				auto t0 = std::chrono::high_resolution_clock::now();
				double sum = 0;
				double theta = 0;
				for (int j = 0; j < Njj; j += 4) {
					__m256d ev = _mm256_set1_pd(140);
					__m256d tv = _mm256_set_pd(theta, theta + delta, theta + 2* delta, theta + 3* delta);
					__m256d sumv = ptr_avx[i](ev, tv);
					sum += hsum(sumv);
					theta += delta;
					theta += delta;
					theta += delta;
					theta += delta;
				}
				auto d = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - t0);
				if (i >= 3)
					printf("avx double order %d, time: 10x%10.1fus, value: %f\n", i, (double)d.count(), sum);
				else
					printf("avx double order %d, time: %10.1fus, value: %f\n", i, (double)d.count(), sum);
			}
#endif
			{
				auto t0 = std::chrono::high_resolution_clock::now();
				double sum = 0;
				double theta = 0;
				for (int j = 0; j < Njj; ++j) {
					f.costheta.SetV(theta);
					sum += f.X_Ws.at(i).V();
					theta += delta;
				}
				auto d = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - t0);
				if(i >= 3)
					printf("tre double order %d, time: 10x%10.1fus, value: %f\n", i, (double)d.count(), sum);
				else
					printf("tre double order %d, time: %10.1fus, value: %f\n", i, (double)d.count(), sum);

			}
			{
				auto t0 = std::chrono::high_resolution_clock::now();
				double sum = 0;
				double theta = 0;
				
				for (int j = 0; j < Njj; ++j) {
					f.costheta.SetV(theta);
					sum += f.X_Ws.at(i).VE().V();
					theta += delta;
				}
				auto d = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - t0);
				if (i >= 3)
					printf("tre num    order %d, time: 10x%10.1fus, value: %f\n", i, (double)d.count(), sum);
				else
					printf("tre num    order %d, time: %10.1fus, value: %f\n", i, (double)d.count(), sum);
			}
		}
	}

}

void test_num() {

	{
		Num x1(1, 0.001, 1E-4);
		Num x2(1, 0.001, 1E-4);
		Num x = x1 + x2;
		TEST_SAME(x.V(), 2);
		TEST_SAME(x.E1(), 0.002);
		TEST_SAME(x.E2(), 2E-4);
	}

	{
		Num x1(1, 0.001, 1E-4);
		Num x2(1, 0.001, 1E-4);
		Num x = x1 - x2;
		TEST_SAME(x.V(), 0);
		TEST_SAME(x.E1(), 0.002);
		TEST_SAME(x.E2(), 2E-4);
	}

	{
		Num x1(2, 0.001, 1E-4);
		Num x2(2, 0.001, 1E-4);
		Num x = x1 * x2;
		TEST_SAME(x.V(), 4);
		TEST_SAME(x.E1(), 0.004);
		TEST_SAME(x.E2(), 8 * 1E-4);
	}

	{
		Num x1(2, 0.001, 1E-4);
		Num x2(2, 0.001, 1E-4);
		Num x = x1 / x2;
		TEST_SAME(x.V(), 1);
		TEST_SAME(x.E1(), 0.001/2 + 0.001*2/4);
		TEST_SAME(x.E2(), 1E-4 / 4 + 1E-4 / 4);

	}

	{
		std::mt19937 mt;
		std::uniform_real_distribution<double> u;
		double sum = 0;
		double sum1 = 0;
		double sum2 = 0;
		int N = 100000;
		int M = 20;
		std::vector<int> ds(M);
		for (int i = 0;i < N; ++i) {
			Num one(1, 0, 0);
			Num two(2, 0, 0);
			Num a(u(mt), 0, 0);
			Num b(u(mt), 0, 0);
			Num c(u(mt), 0, 0);
			Num d(u(mt), 0, 0);
			Num e(u(mt), 0, 0);
			Num f(u(mt), 0, 0);
			Num m =  (a*a + two*a*b + two*a*c + two*a*d + two*a*e + two*a*f + b*b + two*b*c + two*b*d + two*b*e + two*b*f + c*c + two*c*d + two*c*e + two*c*f + d*d + two*d*e + two*d*f + e*e + two*e*f + f*f);
			Num n = pow(a + b + c + d + e + f, 2);
			Num r = m / n;
			double delta = r.V() - 1;
			sum += delta;
			sum1 += fabs(delta);
			sum2 += delta*delta;
			for (int i = 0; i < M; ++i) {
				if (delta < (i - M/2 + 0.5)*DBL_EPSILON && delta >= (i  - M/2 - 0.5)*DBL_EPSILON) ds[i]++;
			}

		}
		Num one(1, 0, 0);
		Num x(u(mt), 0, 0);
		Num y(u(mt), 0, 0);
		Num z(u(mt), 0, 0);

		Num m = (x*y + y*z + z*x) / (x*y*z) * (x + y + z) / (x*y*z);
		Num n = (one / x + one / y + one / z)*(one / x / y + one / y / z + one / z / x);
		Num r = m / n;

		for (int i = 0; i < M; ++i) {
			printf("[%+e, %+e) %d\n", (i - M / 2 - 0.5) * DBL_EPSILON, (i - M / 2 + 0.5) * DBL_EPSILON, ds[i]);
		}

		printf("mean %e\n", sum / N);
		printf("sigma2 %e\n", sum2 / N);
		printf("num E1 %e\n", r.E1());
		printf("num E2 %e\n", r.E2());

	}
	{

		std::mt19937 mt;
		std::uniform_real_distribution<double> u;
		double sum = 0;
		double sum1 = 0;
		double sum2 = 0;
		int N = 100000;
		int M = 20;
		std::vector<int> ds(M);
		for (int i = 0; i < N; ++i) {

			Num one(1, 0, 0);
			Num x(u(mt) * 100, 0, 0);
			Num r = pow(sin(x), 2) + pow(cos(x), 2);

			double delta = r.V() - 1;
			sum += delta;
			sum1 += fabs(delta);
			sum2 += delta*delta;
			for (int i = 0; i < M; ++i) {
				if (delta < (i - M / 2 + 0.5)*DBL_EPSILON && delta >= (i - M / 2 - 0.5)*DBL_EPSILON) ds[i]++;
			}
		}

		for (int i = 0; i < M; ++i) {
			printf("[%+e, %+e) %d\n", (i - M / 2 - 0.5) * DBL_EPSILON, (i - M / 2 + 0.5) * DBL_EPSILON, ds[i]);
		}

		Num one(1, 0, 0);
		Num x(u(mt) * 100, 0, 0);
		Num r = pow(sin(x), 2) + pow(cos(x), 2);
		printf("mean %e\n", sum / N);
		printf("sigma2 %e\n", sum2 / N);
		printf("num E1 %e\n", r.E1());
		printf("num E2 %e\n", r.E2());

	}
}


void test_quad() {

	TEST_SAME(GaussLegendre16Points([](double x) { return x*x; }, 0, 1), 1. / 3);
	TEST_SAME(GaussLegendre64Points([](double x) { return x*x; }, 0, 1), 1. / 3);
	TEST_SAME(TanhSinh65Points([](double x) { return x*x; }, 0, 1), 1. / 3);
	TEST_SAME(Tanh65Points([](double x) { return x*x; }, 0, 1), 1. / 3);
	TEST_SAME(Erf65Points([](double x) { return x*x; }, 0, 1), 1. / 3);
	printf("GL16     %.20f\n", GaussLegendre16Points([](double x) { return x*x; }, 0, 1));
	printf("GL64     %.20f\n", GaussLegendre64Points([](double x) { return x*x; }, 0, 1));
	printf("TanhSinh %.20f\n", TanhSinh65Points([](double x) { return x*x; }, 0, 1));
	printf("Tanh     %.20f\n", Tanh65Points([](double x) { return x*x; }, 0, 1));
	printf("Erf      %.20f\n", Erf65Points([](double x) { return x*x; }, 0, 1));

	TEST_SAME(GaussLegendre16Points([](double x) { return sqrt(x); }, 0, 1), 2. / 3);
	TEST_SAME(GaussLegendre64Points([](double x) { return sqrt(x); }, 0, 1), 2. / 3);
	TEST_SAME(TanhSinh65Points([](double x) { return sqrt(x); }, 0, 1), 2. / 3);
	TEST_SAME(Tanh65Points([](double x) { return sqrt(x); }, 0, 1), 2. / 3);
	TEST_SAME(Erf65Points([](double x) { return sqrt(x); }, 0, 1), 2. / 3);
	printf("GL16     %.20f\n", GaussLegendre16Points([](double x) { return sqrt(x); }, 0, 1));
	printf("GL64     %.20f\n", GaussLegendre64Points([](double x) { return sqrt(x); }, 0, 1));
	printf("TanhSinh %.20f\n", TanhSinh65Points([](double x) { return sqrt(x); }, 0, 1));
	printf("Tanh     %.20f\n", Tanh65Points([](double x) { return sqrt(x); }, 0, 1));
	printf("Erf      %.20f\n", Erf65Points([](double x) { return sqrt(x); }, 0, 1));

	// these two will fail
	//TEST_SAME(GaussianLegendre16Points([](double x) { return 1 / sqrt(x); }, 0, 1), 2. / 3);
	//TEST_SAME(GaussianLegendre64Points([](double x) { return 1 / sqrt(x); }, 0, 1), 2. / 3);
	TEST_SAME(TanhSinh65Points([](double x) { return 1 / sqrt(x); }, 0, 1), 2);
	TEST_SAME(Tanh65Points([](double x) { return 1 / sqrt(x); }, 0, 1), 2);
	TEST_SAME(Erf65Points([](double x) { return 1 / sqrt(x); }, 0, 1), 2);
	printf("GL16     %.20f\n", GaussLegendre16Points([](double x) { return 1 / sqrt(x); }, 0, 1));
	printf("GL64     %.20f\n", GaussLegendre64Points([](double x) { return 1 / sqrt(x); }, 0, 1));
	printf("TanhSinh %.20f\n", TanhSinh65Points([](double x) { return 1 / sqrt(x); }, 0, 1));
	printf("Tanh     %.20f\n", Tanh65Points([](double x) { return 1 / sqrt(x); }, 0, 1));
	printf("Erf      %.20f\n", Erf65Points([](double x) { return 1 / sqrt(x); }, 0, 1));

	printf("GL16     %+.6e\n", GaussLegendre16Points([](double x) { return 1 / (1 + 20 * x*x); }, -1, 1) - 2 * 0.302049929383142873916842364599);
	printf("GL64     %+.6e\n", GaussLegendre64Points([](double x) { return 1 / (1 + 20 * x*x); }, -1, 1) - 2 * 0.302049929383142873916842364599);
	printf("TanhSinh %+.6e\n", TanhSinh65Points([](double x) { return 1 / (1 + 20 * x*x); }, -1, 1) - 2 * 0.302049929383142873916842364599);
	printf("Tanh     %+.6e\n", Tanh65Points([](double x) { return 1 / (1 + 20 * x*x); }, -1, 1) - 2 * 0.302049929383142873916842364599);
	printf("Erf      %+.6e\n", Erf65Points([](double x) { return 1 / (1 + 20 * x*x); }, -1, 1) - 2 * 0.302049929383142873916842364599);

	printf("GL16     %+.6e\n", GaussLegendre16Points([](double x) { return 1 / (1 + 0.1 * x*x); }, -1, 1) - 2 * 0.96853408234038924938038050);
	printf("GL64     %+.6e\n", GaussLegendre64Points([](double x) { return 1 / (1 + 0.1 * x*x); }, -1, 1) - 2 * 0.96853408234038924938038050);
	printf("TanhSinh %+.6e\n", TanhSinh65Points([](double x) { return 1 / (1 + 0.1 * x*x); }, -1, 1) - 2 * 0.96853408234038924938038050);
	printf("Tanh     %+.6e\n", Tanh65Points([](double x) { return 1 / (1 + 0.1 * x*x); }, -1, 1) - 2 * 0.96853408234038924938038050);
	printf("Erf      %+.6e\n", Erf65Points([](double x) { return 1 / (1 + 0.1 * x*x); }, -1, 1) - 2 * 0.96853408234038924938038050);

	{
		printf("tanh_sinh_65points_h %f\n", tanh_sinh_65points_h);
		printf("Methods  Range    Integrand    Error\n");
		printf("TanhSinh [ 0, 1 ] x*x         %+.6e\n", TanhSinh65Points([](double x) { return x*x; }, 0, 1) - 1/3.);
		printf("TanhSinh [ 0, 1 ] sqrt(x)     %+.6e\n", TanhSinh65Points([](double x) { return sqrt(x); }, 0, 1) - 2/3.);
		printf("TanhSinh [ 0, 1 ] log(x)      %+.6e\n", TanhSinh65Points([](double x) { return log(x); }, 0, 1) - -1);
		printf("TanhSinh [ 0, 1 ] 1/sqrt(x)   %+.6e\n", TanhSinh65Points([](double x) { return 1 / sqrt(x); }, 0, 1) - 2.);
		printf("TanhSinh [ 0, 1 ] pow(x,-2/3) %+.6e\n", TanhSinh65Points([](double x) { return pow(x, -2. / 3); }, 0, 1) - 3.);
		printf("TanhSinh [ 0, 1 ] 1/(1+20xx)  %+.6e\n", TanhSinh65Points([](double x) { return 1 / (1 + 20 * x*x); }, 0, 1) - 0.302049929383142873916842364599);
		printf("TanhSinh [-1, 1 ] 1/(1+20xx)  %+.6e\n", TanhSinh65Points([](double x) { return 1 / (1 + 20 * x*x); }, -1, 1) - 2*0.302049929383142873916842364599);
		printf("TanhSinh [ 0, 1 ] 1/(1+xx)    %+.6e\n", TanhSinh65Points([](double x) { return 1 / (1 + x*x); }, 0, 1) - 0.78539816339744830961566);
		printf("TanhSinh [-1, 1 ] 1/(1+xx)    %+.6e\n", TanhSinh65Points([](double x) { return 1 / (1 + x*x); }, -1, 1) - 2*0.78539816339744830961566);
	}


	{
		TEST_SAME(TanhSinh65Points([](double x1, double x2) { 
			return 1/sqrt(x1*x2); },
			0, 1), 3.1415926535897932384);
		TEST_SAME(TanhSinh65Points([](double x) {
			double one_minus_x = 1 - x;
			if (one_minus_x == 0) one_minus_x = 1 - nextafter(1, -1);
			return 1 / sqrt(x*one_minus_x);
		}, 0, 1), 3.1415926535897932384);
		printf("TanhSinh [ 0, 1 ] 1/sqrt(x(1-x)) %+.6e\n", TanhSinh65Points([](double x, double x2) { return 1 / sqrt(x*x2); }, 0, 1) - Pi);
		printf("TanhSinh [ 0, 1 ] 1/sqrt(x(1-x)) %+.6e\n", TanhSinh65Points([](double x) { 
			double one_minus_x = 1 - x;
			if (one_minus_x == 0) one_minus_x = 1 - nextafter(1, -1);
			return 1 / sqrt(x*one_minus_x);
		}, 0, 1) - Pi);
	}
}


int main() {
	TEST_TRUE(DCount.size() == 0);
	test_num();
	test_quad();
	test_V();
	test_Diff();
	fix_var();
	testfunc();
	test_int();
	test_for();
	test_code();

	printf("%d test(s) failed\n", n_failed);
	return n_failed;
}
