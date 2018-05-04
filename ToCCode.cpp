#include "Diff.h"


namespace Diff {

	struct ToCCodeT
	{


		std::string sb;
		std::map<Expr, std::string, CCode::ExprLess> names;

		std::set<std::string> ufs;
		std::map<std::string, std::string> bfs;

		ToCCodeT() {

			static char const *unitary[] = {
				"sin",
				"cos",
				"exp",
				"sinh",
				"cosh",
				"log",
			};
			for (auto x : unitary) {
				ufs.insert(x);
			}
			bfs["add"] = "+";
			bfs["sub"] = "-";
			bfs["mul"] = "*";
			bfs["div"] = "/";
		}

		void add(Expr const &expr)
		{

			if (names.find(expr) != names.end()) {
				return;
			}

			SubExpressionVector subs;
			char n[64];
			char b[256];

			expr.GetSubExpressions(subs);
			for (auto &sub : subs) {
				add(sub);
			}

			auto &type_name = expr.GetTypeName();
			if (type_name == "Constant") {
				sprintf(n, "C_%llu", (unsigned long long)expr.Uid());
				names[expr] = n;

				ParameterVector ps;
				expr.GetParameters(ps);
				sprintf(b, "double const %s = %.20E;\n", n, ps.at(0));

			} else if (type_name == "Variable") {

				Var var = Cast<Var>(expr);
				if (var.GetName().empty()) {
					sprintf(n, "V_%llu", (unsigned long long)expr.Uid());
					names[expr] = n;
				} else {
					names[expr] = var.GetName();
				}
				sprintf(b, "");

			} else if (ufs.find(type_name) != ufs.end()) {
				sprintf(n, "%s_%llu", type_name.c_str(), (unsigned long long)expr.Uid());
				names[expr] = n;

				sprintf(b, "double const %s = %s(%s);\n", n, type_name.c_str(), names[subs.at(0)].c_str());
			} else if (bfs.find(type_name) != bfs.end()) {
				sprintf(n, "%s_%llu", type_name.c_str(), (unsigned long long)expr.Uid());
				names[expr] = n;

				auto &op = bfs.find(type_name)->second;
				sprintf(b, "double const %s = %s %s %s;\n", n,
					names[subs.at(0)].c_str(),
					op.c_str(),
					names[subs.at(1)].c_str());

			} else if (type_name == "pow") {
				sprintf(n, "%s_%llu", type_name.c_str(), (unsigned long long)expr.Uid());
				names[expr] = n;
				ParameterVector ps;
				expr.GetParameters(ps);
				sprintf(b, "double const %s = pow(%s, %.20E);\n", n, names[subs.at(0)].c_str(), ps.at(0));
			} else {
				throw std::logic_error("not implemented");
			}
			sb.append(b);

		}


	};


	CCode ToCCode(Expr const &expr) {
		ToCCodeT t;
		t.add(expr);

		CCode c;
		c.Body = std::move(t.sb);
		c.Names = std::move(t.names);
		return c;
	}
}
