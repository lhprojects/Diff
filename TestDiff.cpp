#include "Diff.h"
#include "Func.h"
using namespace Diff;

#include <math.h>
#define TPrintf printf


static int n_failed = 0;
#define TEST_SAME(x, y) do { double x_ = (x); double y_ = (y);\
	if(isnan(x_) || fabs(x_ - y_) > abs(y_)*1E-3)  { n_failed++; TPrintf("FAILED: " #x " (%f) == " #y " (%f)\n", x_, y_);}   } while(0)
#define TEST_TRUE(x) do { double x_ = (x); if(!x_)  { n_failed++; TPrintf("FAILED: " #x "\n");}   } while(0)

void test_Diff() {


	TEST_TRUE(DCount.size() == 0);
	{
		Var x = 4;
		TEST_SAME(sqrt(x).D(x).V(), 0.5 / sqrt(x.V()));
	}

	{
		Var three = 3;
		TEST_SAME(three.V(), 3);
		TEST_SAME((three * three).V(), 9);
		TEST_SAME((three + three).V(), 6);
		TEST_SAME((three - three).V(), 0);
		TEST_SAME((three / three).V(), 1);
	}

	{
		Var x = 4;
		TEST_SAME(x.D(x).V(), 1);
		TEST_SAME((x * x).D(x).V(), 2 * x.V());
		TEST_SAME((x + x).D(x).V(), 2);
		TEST_SAME((x - x).D(x).V(), 0);
		TEST_SAME((x / x).D(x).V(), 0);
		TEST_SAME(sqrt(x).D(x).V(), 0.5 / sqrt(x.V()));
		TEST_SAME(sqrt(x).D(x).D(x).V(), -0.25 / pow(sqrt(x.V()), 3));
		TEST_SAME(sqrt(x).D(x).D(x).D(x).V(), 3 * 0.125 / pow(sqrt(x.V()), 5));
		TEST_SAME(log(x).D(x).V(), 1 / x.V());
		TEST_SAME(POW2(x).D(x).V(), 2 * x.V());
	}


	{
		Var x = 4;
		TEST_SAME((x*x*x).D(x).D(x).V(), 3 * 2 * x.V());
		TEST_SAME((x*x*x*x).D(x).D(x).D(x).V(), 4 * 3 * 2 * x.V());
	}
	{
		Var x("x", 4);
		TEST_SAME(POW2(x).V(), (x*x).V());
		TEST_SAME(POW2(x).D(x).V(), (x*x).D(x).V());
		TEST_SAME(POW2(x).D(x).D(x).V(), (x*x).D(x).D(x).V());
		TEST_SAME(POW2(x).D(x).D(x).D(x).V(), (x*x).D(x).D(x).D(x).V());

		TEST_SAME(POW4(x).V(), (x*x*x*x).V());
		TEST_SAME(POW4(x).D(x).V(), (x*x*x*x).D(x).V());

	}

	{ // exp
		Var x("x", 4);
		TEST_SAME(exp(2 * x).V(), exp(2 * x.V()));
		TEST_SAME(exp(2 * x).D(x).V(), 2 * exp(2 * x.V()));
		TEST_SAME(exp(2 * x).D(x).D(x).V(), 2 * 2 * exp(2 * x.V()));

		TPrintf("exp(2x)          : %s\n", exp(2 * x).ToString().c_str());
		TPrintf("exp(2x)'         : %s\n", exp(2 * x).D(x).ToString().c_str());
		TPrintf("exp(2x)''        : %s\n", exp(2 * x).D(x).D(x).ToString().c_str());
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
double POW2(double x) { return x*x; }
double POW4(double x) { return x*x*x*x; }

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
			X_Is.push_back(X_Is.at(X_Is.size() - 1).D(costheta));
			costheta.SetV(0);
			X_Is0.push_back(X_Is.back().FixVariable(costheta));
		}
		for (; (int)X_Ws.size() <= n; ) {
			X_Ws.push_back(X_Ws.at(X_Ws.size() - 1).D(costheta));
			costheta.SetV(0);
			X_Ws0.push_back(X_Ws.back().FixVariable(costheta));
		}
		for (; (int)X_Ss.size() <= n; ) {
			X_Ss.push_back(X_Ss.at(X_Ss.size() - 1).D(costheta));
			costheta.SetV(0);
			X_Ss0.push_back(X_Ss.back().FixVariable(costheta));
		}
	}

	Formula(Constants C, double XsectionGlobalFactor) : costheta("costheta", 0), eH("e", 140)
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

		if (0) {
			printf("L %f\n", L.V());
			printf("r %f\n", r.V());
			printf("s1 %f\n", s1.V());
			printf("p %f\n", p.V());
			printf("mH %f\n", mH.V());
			printf("eH %f\n", eH.V());
		}

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

		TEST_SAME(f.X_Ws0.at(2).V(), ypp);
		TEST_SAME(f.X_Ws.at(2).V(), ypp);

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
	Var t("t", 1);
	Var x("x", 0);
	Expr y = Integrate(x, Const(0), t, x*x*t);
	TEST_SAME(y.V(), 1/3.0);
	TEST_SAME(y.D(t).V(), 1 / 3.0 + 1);
	printf("%s\n", y.D(t).ToString().c_str());
}


int main() {

	test_Diff();
	fix_var();
	testfunc();
	test_int();
	test_for();

	printf("%d test(s) failed\n", n_failed);
	return 0;
}
