#include "Diff.h"

namespace Diff {


	struct RVWorker {

		Var fVar;
		Expr fBy;
		std::map<uint64_t, RebindableExpr> fMemory;
		RVWorker(Var const &var, Expr const &by) : fVar(var), fBy(by) {
		}

		Expr replace(Expr const &what)
		{
			auto &p = fMemory[what.Uid()];
			if (!p.Empty()) return p;
			std::string const &name = what.GetTypeName();
			SubExpressionVector subs;
			if (name == "Constant") {
				p = what;
			} else if (name == "Variable") {
				if (what.Uid() == fVar.Uid()) {
					p = fBy;
				} else {
					p = what;
				}
			} else if (name == "sin") {
				what.GetSubExpressions(subs);
				auto new_sub0 = replace(subs.at(0));
				if (subs[0].Uid() == new_sub0.Uid()) {
					p = what;
				} else {
					p = sin(new_sub0);
				}
			} else if (name == "cos") {
				what.GetSubExpressions(subs);
				auto new_sub0 = replace(subs.at(0));
				if (subs[0].Uid() == new_sub0.Uid()) {
					p = what;
				} else {
					p = cos(new_sub0);
				}
			} else if (name == "exp") {
				what.GetSubExpressions(subs);
				auto new_sub0 = replace(subs.at(0));
				if (subs[0].Uid() == new_sub0.Uid()) {
					p = what;
				} else {
					p = exp(new_sub0);
				}
			} else if (name == "sinh") {
				what.GetSubExpressions(subs);
				auto new_sub0 = replace(subs.at(0));
				if (subs[0].Uid() == new_sub0.Uid()) {
					p = what;
				} else {
					p = sinh(new_sub0);
				}
			} else if (name == "cosh") {
				what.GetSubExpressions(subs);
				auto new_sub0 = replace(subs.at(0));
				if (subs[0].Uid() == new_sub0.Uid()) {
					p = what;
				} else {
					p = cosh(new_sub0);
				}
			} else if (name == "log") {
				what.GetSubExpressions(subs);
				auto new_sub0 = replace(subs.at(0));
				if (subs[0].Uid() == new_sub0.Uid()) {
					p = what;
				} else {
					p = log(new_sub0);
				}
			} else if (name == "pow") {
				what.GetSubExpressions(subs);
				auto new_sub0 = replace(subs.at(0));
				if (subs[0].Uid() == new_sub0.Uid()) {
					p = what;
				} else {
					ParameterVector pars;
					what.GetParameters(pars);
					p = pow(new_sub0, pars.at(0));
				}
			} else if (name == "add") {
				what.GetSubExpressions(subs);
				auto new_sub0 = replace(subs.at(0));
				auto new_sub1 = replace(subs.at(1));
				if (subs[0].Uid() == new_sub0.Uid() && subs[1].Uid() == new_sub1.Uid()) {
					p = what;
				} else {
					p = new_sub0 + new_sub1;
				}
			} else if (name == "sub") {
				what.GetSubExpressions(subs);
				auto new_sub0 = replace(subs.at(0));
				auto new_sub1 = replace(subs.at(1));
				if (subs[0].Uid() == new_sub0.Uid() && subs[1].Uid() == new_sub1.Uid()) {
					p = what;
				} else {
					p = new_sub0 - new_sub1;
				}
			} else if (name == "mul") {
				what.GetSubExpressions(subs);
				auto new_sub0 = replace(subs.at(0));
				auto new_sub1 = replace(subs.at(1));
				if (subs[0].Uid() == new_sub0.Uid() && subs[1].Uid() == new_sub1.Uid()) {
					p = what;
				} else {
					p = new_sub0 * new_sub1;
				}
			} else if (name == "div") {
				what.GetSubExpressions(subs);
				auto new_sub0 = replace(subs.at(0));
				auto new_sub1 = replace(subs.at(1));
				if (subs[0].Uid() == new_sub0.Uid() && subs[1].Uid() == new_sub1.Uid()) {
					p = what;
				} else {
					p = new_sub0 / new_sub1;
				}
			} else if (name == "Sum") {

				what.GetSubExpressions(subs);
				if (subs.at(1).Uid() == fVar.Uid() && !IsVar(what)) {
					throw std::logic_error("The summation variable must be replaced by exprssion of type `Variable`");
				}
				auto new_sub0 = replace(subs.at(0));
				auto new_sub1 = replace(subs.at(1));
				if (subs[0].Uid() == new_sub0.Uid() && subs[1].Uid() == new_sub1.Uid()) {
					p = what;
				} else {
					ParameterVector pars;
					what.GetParameters(pars);
					Sum(new_sub0, new_sub1, pars.at(0), pars.at(1));
				}
			} else if (name == "Integrate") {

				what.GetSubExpressions(subs);
				if (subs.at(1).Uid() == fVar.Uid() && !IsVar(what)) {
					throw std::logic_error("The integrate variable must be replaced by exprssion of type `Variable`");
				}
				auto new_sub0 = replace(subs.at(0));
				auto new_sub1 = replace(subs.at(1));
				auto new_sub2 = replace(subs.at(2));
				auto new_sub3 = replace(subs.at(3));
				if (subs[0].Uid() == new_sub0.Uid() && subs[1].Uid() == new_sub1.Uid() &&
					subs[2].Uid() == new_sub2.Uid() && subs[3].Uid() == new_sub3.Uid()) {
					p = what;
				} else {
					Integrate(new_sub0, new_sub1, new_sub2, new_sub3);
				}
			} else {
				throw std::logic_error("not implemented!\n");
			}

			return p;
		}

	};


	// replace a variable with a expresion
	// any sub expressions of `expr` can't have a reference to the any sub expression of this
	Expr ReplaceVariable(Expr const &in_what, Expr const &to_be_replaced, ExprOrDouble const &replaced_by) {
		RVWorker worker(CastToVar(to_be_replaced), replaced_by);
		return worker.replace(in_what);
	}

	Expr ReplaceVariable(Expr const &in_what, std::pair<Expr, ExprOrDouble> const & rep) {
		RVWorker worker(CastToVar(rep.first), rep.second);
		return worker.replace(in_what);
	}

}
