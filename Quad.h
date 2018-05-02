#pragma once

#include <functional>

namespace Diff {
	extern double const gl_x_16points[8];
	extern double const gl_w_16points[8];
	extern double const gl_x_64points[32];
	extern double const gl_w_64points[32];
	double GaussianLegendre64Points(std::function<double(double x)> const &f, double x0, double x1);
	double GaussianLegendre16Points(std::function<double(double x)> const &f, double x0, double x1);
	double TanhSinh65Points(std::function<double(double x)> const &f, double x0, double x1);
	double Tanh65Points(std::function<double(double x)> const &f, double x0, double x1);

	extern double const tanh_sinh_65points_h;
	double Erf65Points(std::function<double(double x)> const &f, double x0, double x1);
}
