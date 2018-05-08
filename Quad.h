#pragma once

#include <type_traits>

namespace Diff {
	extern double const gl_x_16points[8];
	extern double const gl_w_16points[8];
	extern double const gl_x_64points[32];
	extern double const gl_w_64points[32];
	double GaussLegendre64Points(double(*stub)(void const* ptr, double x), void const *ptr, double x0, double x1);
	double GaussLegendre16Points(double(*stub)(void const* ptr, double x), void const *ptr, double x0, double x1);

	template<class L>
	std::enable_if_t<
		std::is_same< 
			decltype(std::declval<L>()(std::declval<double>())),
			double
		>::value
	, double>
		GaussLegendre64Points(L const &l, double x0, double x1) {
		double(*stub)(void const*, double) = [](void const *ptr, double x) {
			return reinterpret_cast<L const *>(ptr)->operator()(x);
		};
		return GaussLegendre64Points(stub, reinterpret_cast<void const*>(&l), x0, x1);
	}

	template<class L>
	std::enable_if_t<
		std::is_same< 
			decltype(std::declval<L>()(std::declval<double>())),
			double
		>::value
	, double>
		GaussLegendre16Points(L const &l, double x0, double x1) {
		double(*stub)(void const*, double) = [](void const *ptr, double x) {
			return reinterpret_cast<L const *>(ptr)->operator()(x);
		};
		return GaussLegendre16Points(stub, reinterpret_cast<void const*>(&l), x0, x1);
	}

	extern double const tanh_sinh_65points_h;
	extern double const tanh_sinh_1_minus_x_65points[33];
	extern double const tanh_sinh_w_65points[33];
	double TanhSinh65Points(double(*stub)(void const*, double), void const *, double x0, double x1);
	double TanhSinh65Points(double(*stub)(void const*, double, double), void const *, double x0, double x1);

	template<class L>
	typename std::enable_if<
		std::is_same< 
			decltype(std::declval<L>()(std::declval<double>())),
			double
		>::value
	, double>::type
	TanhSinh65Points(L const &l, double x0, double x1) {
		double(*stub)(void const*, double) = [](void const *ptr, double x) {
			return reinterpret_cast<L const *>(ptr)->operator()(x);
		};
		return TanhSinh65Points(stub, reinterpret_cast<void const*>(&l), x0, x1);
	}

	template<class L, class DoubleArg = int>
	typename std::enable_if<
		std::is_same<
			decltype(std::declval<L>()(std::declval<double>(), std::declval<double>())),
			double
		>::value	
	, double>::type
	TanhSinh65Points(L const &l, double x0, double x1) {
		double(*stub)(void const*, double, double) = [](void const *ptr, double xa, double xb) {
			return reinterpret_cast<L const *>(ptr)->operator()(xa, xb);
		};
		return TanhSinh65Points(stub, reinterpret_cast<void const*>(&l), x0, x1);
	}

	double ExpSinh65Points(double(*stub)(void const* ptr, double x), void const *ptr, double x0);

	template<class L>
	std::enable_if_t<
		std::is_same< 
			decltype(std::declval<L>()(std::declval<double>())),
			double
		>::value
	, double>
		ExpSinh65Points(L const &l, double x0) {
		double(*stub)(void const*, double) = [](void const *ptr, double x) {
			return reinterpret_cast<L const *>(ptr)->operator()(x);
		};
		return ExpSinh65Points(stub, reinterpret_cast<void const*>(&l), x0);
	}


	double Tanh65Points(double(*stub)(void const* ptr, double x), void const *ptr, double x0, double x1);

	template<class L>
	std::enable_if_t<
		std::is_same< 
			decltype(std::declval<L>()(std::declval<double>())),
			double
		>::value
	, double>
		Tanh65Points(L const &l, double x0, double x1) {
		double(*stub)(void const*, double) = [](void const *ptr, double x) {
			return reinterpret_cast<L const *>(ptr)->operator()(x);
		};
		return Tanh65Points(stub, reinterpret_cast<void const*>(&l), x0, x1);
	}

	double Erf65Points(double(*stub)(void const* ptr, double x), void const *ptr, double x0, double x1);

	template<class L>
	std::enable_if_t<
		std::is_same< 
			decltype(std::declval<L>()(std::declval<double>())),
			double
		>::value
	, double>
		Erf65Points(L const &l, double x0, double x1) {
		double(*stub)(void const*, double) = [](void const *ptr, double x) {
			return reinterpret_cast<L const *>(ptr)->operator()(x);
		};
		return Erf65Points(stub, reinterpret_cast<void const*>(&l), x0, x1);
	}

}
