#include "Quad.h"
#include "Constants.h"
#include <math.h>

namespace Diff {

	double const gl_x_16points[] = {
		0.0950125098376374,
		0.2816035507792589,
		0.4580167776572274,
		0.6178762444026438,
		0.7554044083550030,
		0.8656312023878318,
		0.9445750230732326,
		0.9894009349916499,
	};

	double const gl_w_16points[] = {
		0.1894506104550685,
		0.1826034150449236,
		0.1691565193950025,
		0.1495959888165767,
		0.1246289712555339,
		0.0951585116824928,
		0.0622535239386479,
		0.0271524594117541,
	};


	double GaussLegendre16Points(std::function<double(double x)> const &f, double x0, double x1)
	{
		double delta = (x1 - x0) / 2;
		double median = (x1 + x0) / 2;
		double h = 0;
		for (int i = 0; i < 8; ++i) {
			double x = gl_x_16points[i] * delta + median;
			h += f(x)*gl_w_16points[i];
			x = -gl_x_16points[i] * delta + median;
			h += f(x)*gl_w_16points[i];
		}
		h *= delta;
		return h;
	}

	double const gl_x_64points[] = {
		0.0243502926634244,
		0.0729931217877990,
		0.1214628192961206,
		0.1696444204239928,
		0.2174236437400071,
		0.2646871622087674,
		0.3113228719902110,
		0.3572201583376681,
		0.4022701579639916,
		0.4463660172534641,
		0.4894031457070530,
		0.5312794640198946,
		0.5718956462026340,
		0.6111553551723933,
		0.6489654712546573,
		0.6852363130542333,
		0.7198818501716109,
		0.7528199072605319,
		0.7839723589433414,
		0.8132653151227975,
		0.8406292962525803,
		0.8659993981540928,
		0.8893154459951141,
		0.9105221370785028,
		0.9295691721319396,
		0.9464113748584028,
		0.9610087996520538,
		0.9733268277899110,
		0.9833362538846260,
		0.9910133714767443,
		0.9963401167719553,
		0.9993050417357722,
	};

	double const gl_w_64points[] = {
		0.0486909570091397203833654,
		0.0485754674415034269347991,
		0.0483447622348029571697695,
		0.0479993885964583077281262,
		0.0475401657148303086622822,
		0.0469681828162100173253263,
		0.0462847965813144172959532,
		0.0454916279274181444797710,
		0.0445905581637565630601347,
		0.0435837245293234533768279,
		0.0424735151236535890073398,
		0.0412625632426235286101563,
		0.0399537411327203413866569,
		0.0385501531786156291289625,
		0.0370551285402400460404151,
		0.0354722132568823838106931,
		0.0338051618371416093915655,
		0.0320579283548515535854675,
		0.0302346570724024788679741,
		0.0283396726142594832275113,
		0.0263774697150546586716918,
		0.0243527025687108733381776,
		0.0222701738083832541592983,
		0.0201348231535302093723403,
		0.0179517157756973430850453,
		0.0157260304760247193219660,
		0.0134630478967186425980608,
		0.0111681394601311288185905,
		0.0088467598263639477230309,
		0.0065044579689783628561174,
		0.0041470332605624676352875,
		0.0017832807216964329472961,
	};

	double GaussLegendre64Points(std::function<double(double x)> const &f, double x0, double x1)
	{
		double delta = (x1 - x0) / 2;
		double median = (x1 + x0) / 2;
		double h = 0;

		for (int i = 0; i < 32; ++i) {
			double x = -gl_x_64points[i] * delta + median;
			//printf("%d %f %f %f %f\n", i, x, gl_w_64points[i], f(x), h);
			h += f(x)*gl_w_64points[i];
			x = gl_x_64points[i] * delta + median;
			//printf("%d %f %f %f %f\n", i, x, gl_w_64points[i], f(x), h);
			h += f(x)*gl_w_64points[i];
		}
		h = h * delta;
		return h;
	}


	double const ts_65points_h = 0.13;
	double const tanh_sinh_65points_h = ts_65points_h;

#define SQURE(x) (x*x)

	double const tanh_sinh_w_65points[] = {
		0.5*Pi*ts_65points_h*std::cosh(0 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(0 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(1 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(1 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(2 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(2 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(3 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(3 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(4 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(4 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(5 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(5 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(6 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(6 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(7 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(7 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(8 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(8 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(9 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(9 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(10 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(10 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(11 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(11 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(12 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(12 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(13 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(13 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(14 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(14 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(15 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(15 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(16 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(16 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(17 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(17 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(18 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(18 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(19 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(19 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(20 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(20 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(21 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(21 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(22 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(22 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(23 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(23 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(24 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(24 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(25 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(25 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(26 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(26 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(27 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(27 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(28 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(28 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(29 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(29 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(30 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(30 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(31 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(31 * ts_65points_h))),
		0.5*Pi*ts_65points_h*std::cosh(32 * ts_65points_h) / SQURE(std::cosh(0.5*Pi*std::sinh(32 * ts_65points_h))),
	};


	double const tanh_sinh_1_minus_x_65points[] = {
		1,
		exp(0.5*Pi*std::sinh(-1 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(1 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-2 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(2 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-3 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(3 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-4 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(4 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-5 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(5 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-6 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(6 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-7 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(7 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-8 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(8 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-9 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(9 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-10 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(10 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-11 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(11 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-12 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(12 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-13 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(13 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-14 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(14 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-15 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(15 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-16 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(16 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-17 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(17 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-18 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(18 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-19 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(19 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-20 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(20 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-21 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(21 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-22 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(22 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-23 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(23 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-24 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(24 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-25 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(25 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-26 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(26 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-27 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(27 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-28 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(28 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-29 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(29 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-30 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(30 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-31 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(31 * ts_65points_h)),
		exp(0.5*Pi*std::sinh(-32 * ts_65points_h)) / std::cosh(0.5*Pi*std::sinh(32 * ts_65points_h)),
	};

	double TanhSinh65Points(std::function<double(double x_minus_x0, double x1_minus_x)> const &f, double x0, double x1)
	{

		int n = 32;
		double I = 0;

		for (int k = -32; k <= 32; ++k) {

			double x_minus_x0;
			double x1_minus_x;
			double w;
			if (k >= 0) {
				x1_minus_x = 0.5 * tanh_sinh_1_minus_x_65points[k] * (x1 - x0);
				x_minus_x0 = x1 - x0 - x1_minus_x;
				w = tanh_sinh_w_65points[k];
			} else {
				x_minus_x0 = 0.5 * tanh_sinh_1_minus_x_65points[-k] * (x1 - x0);
				x1_minus_x = x1 - x0 - x_minus_x0;
				w = tanh_sinh_w_65points[-k];
			}

			I += f(x_minus_x0, x1_minus_x) * w;
		}
		I *= (x1 - x0) / 2;
		return I;


	}
	double TanhSinh65Points(std::function<double(double x)> const &f, double x0, double x1)
	{

		int n = 32;
		double I = 0;
		
		for (int k = -32; k <= 32; ++k) {

			double x;
			double w;
			if (k >= 0) x = x1 - 0.5 * tanh_sinh_1_minus_x_65points[k] *(x1 - x0);
			else if (k < 0) x = x0 + 0.5 * tanh_sinh_1_minus_x_65points[-k] * (x1 - x0);

			if (k >= 0) w = tanh_sinh_w_65points[k];
			else if (k < 0) w = tanh_sinh_w_65points[-k];

			I += f(x) * w;
		}
		I *= (x1 - x0) / 2;
		return I;
	}

	double Tanh65Points(std::function<double(double x)> const &f, double x0, double x1)
	{
		double h = Pi / sqrt(65);
		int n = 32;
		double I = 0;
		for (int k = -n; k <= n; ++k) {
			double coshkh = std::cosh(k*h);
			double tanhkh_plus_one = exp(k*h)/std::cosh(k*h);
			double one_minus_tankh = exp(-k*h) / std::cosh(k*h);
			double x;
			if (k < 0) x = 0.5 * tanhkh_plus_one* (x1 - x0) + x0;
			else x = x1 - 0.5 * one_minus_tankh *(x1 - x0);
			I += f(x) / (coshkh*coshkh);
		}
		I *= h;
		I *= (x1 - x0) / 2;
		return I;
	}

	double Erf65Points(std::function<double(double x)> const &f, double x0, double x1)
	{
		double h = 0.5*pow(65, -1. / 3);
		int n = 32;
		double I = 0;
		for (int k = -n; k <= n; ++k) {
			double kh = k*h;
			double erfkh = std::erf(k*h);
			double x = (erfkh + 1) / 2 * (x1 - x0) + x0;
			I += f(x)*exp(-kh*kh);
		}
		I *= 2 / std::sqrt(Pi);
		I *= h;

		I *= (x1 - x0) / 2;
		return I;
	}

}
  