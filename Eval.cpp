#include "Diff.h"

namespace Diff {


	struct Evalor {

		std::map<uint64_t, Num> fMem;

		Num Eva(Expr const &x)
		{
			auto it = fMem.find(x.Uid());
			if (it != fMem.end()) {
				return it->second;
			}
			Num num(0,0,0);
			SubExpressionVector subs;
			auto &type_name = x.GetTypeName();
			x.GetSubExpressions(subs);

			if (type_name == "Variable") {
				ParameterVector pars;
				x.GetParameters(pars);
				num = Num(pars.at(0), 0, 0);
			} else if (type_name == "Constant") {
				ParameterVector pars;
				x.GetParameters(pars);
				num = Num(pars.at(0), pars.at(1), 0);
			} else if (type_name == "add") {
				Num t0 = Eva(subs.at(0));
				Num t1 = Eva(subs.at(1));
				num = t0 + t1;
			} else if (type_name == "sub") {
				Num t0 = Eva(subs.at(0));
				Num t1 = Eva(subs.at(1));
				num = t0 - t1;
			} else if (type_name == "mul") {
				Num t0 = Eva(subs.at(0));
				Num t1 = Eva(subs.at(1));
				num = t0 * t1;
			} else if (type_name == "div") {
				Num t0 = Eva(subs.at(0));
				Num t1 = Eva(subs.at(1));
				num = t0 / t1;
			} else if (type_name == "sin") {
				Num t = Eva(subs.at(0));
				num = sin(t);
			} else if (type_name == "cos") {
				Num t = Eva(subs.at(0));
				num = cos(t);
			} else if (type_name == "sinh") {
				Num t = Eva(subs.at(0));
				num = sinh(t);
			} else if (type_name == "cosh") {
				Num t = Eva(subs.at(0));
				num = cosh(t);
			} else if (type_name == "exp") {
				Num t = Eva(subs.at(0));
				num = exp(t);
			} else if (type_name == "log") {
				Num t = Eva(subs.at(0));
				num = log(t);
			} else if (type_name == "tan") {
				Num t = Eva(subs.at(0));
				num = tan(t);
			} else if (type_name == "pow") {
				Num t = Eva(subs.at(0));
				ParameterVector pars;
				x.GetParameters(pars);
				num = pow(t, pars.at(0));
			} else {
				throw std::logic_error("not implemented");
			}
			fMem.insert(std::make_pair(x.Uid(), num));
			return num;
		}

	};

	Num VE(ExprOrDouble const &expr)
	{
		Evalor ev;
		return ev.Eva(expr);
	}


}
