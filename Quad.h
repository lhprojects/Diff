#pragma once

#include <functional>
#include <type_traits>

namespace Diff {
	extern double const gl_x_16points[8];
	extern double const gl_w_16points[8];
	extern double const gl_x_64points[32];
	extern double const gl_w_64points[32];
	double GaussLegendre64Points(std::function<double(double x)> const &f, double x0, double x1);
	double GaussLegendre16Points(std::function<double(double x)> const &f, double x0, double x1);

	extern double const tanh_sinh_65points_h;
	extern double const tanh_sinh_1_minus_x_65points[33];
	extern double const tanh_sinh_w_65points[33];
	double TanhSinh65Points(std::function<double(double x)> const &f, double x0, double x1);
	double TanhSinh65Points(std::function<double(double x_minus_x0, double x1_minus_x)> const &f, double x0, double x1);

	template<class L>
	std::enable_if_t<
		std::is_same< 
			decltype(std::declval<L>()(std::declval<double>())),
			double
		>::value
	, double>
	TanhSinh65Points(L const &l, double x0, double x1) {
		return TanhSinh65Points(std::function<double(double)>(l), x0, x1);
	}

	template<class L, class DoubleArg = int>
	std::enable_if_t<
		std::is_same<
			decltype(std::declval<L>()(std::declval<double>(), std::declval<double>())),
			double
		>::value	
	, double>
	TanhSinh65Points(L const &l, double x0, double x1) {
		return TanhSinh65Points(std::function<double(double, double)>(l), x0, x1);
	}

	double ExpSinh65Points(std::function<double(double x)> const &f, double x0);
	double Tanh65Points(std::function<double(double x)> const &f, double x0, double x1);
	double Erf65Points(std::function<double(double x)> const &f, double x0, double x1);
}
