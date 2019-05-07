#include "Diff.h"
#include "Quad.h"
#include "SamllVector.h"
#include <map>
#include <typeinfo>
#include <cmath>

namespace Diff {

	uint64_t guid = 0;
	std::set<DExprImpl*> DCount;
	Const const &Zero();
	Const const &One();
	Const const &Two();

	struct DVariableImpl;

	typedef std::set<DVariableImpl const*> DVariabaleSet;

	struct CCodeState {
		std::string fVarName;
	};

	char const *ToStr(double g) {
		static char b[1000];
		sprintf(b, "%g", g);
		return b;
	}

	/****************************	DExprImpl begin **********************************************/


	enum class ExprType {
		NotSet,
		Const,
		Variable,
	};

	struct DExprImpl
	{
		/* Evaluate the value */
		virtual double EvalV() const = 0;
		virtual Expr EvalD(Var const &s) const = 0;

		virtual std::string const & GetTypeName() const = 0;
		virtual void AddExpressions(std::set<DExprImpl const*> &nodes) const = 0;
		virtual void GetSubExpressions(SubExpressionVector &expr) const = 0;
		virtual void GetParameters(ParameterVector &pars) const { }
		virtual void ToString(std::string &s) const = 0;

		virtual ~DExprImpl();

		DExprImpl();

	private:
		friend struct DConstant;
		friend struct DVariableImpl;
		DExprImpl(ExprType type);
		ExprType const fType;
	public:
		bool IsConst() const { return fType == ExprType::Const; }
		bool IsVariable() const { return fType == ExprType::Variable; }

		uint64_t const fUid;


		// reference counting
	private:
		mutable int fRef;
	public:
		void IncRef() const { ++fRef; }
		void DecRef() const {
			--fRef;
			if (fRef == 0) delete this;
		}

		// differential
	private:
		mutable std::map<uint64_t, RebindableExpr> fDMem;
	public:
		Expr DMem(Var const &s) const;

	private:
		mutable double fVMem;
	public:
		mutable bool fVMemValid;

		/* Get the value, re-evaluate if non-valid */
		double VMem() const
		{
			if (!fVMemValid) {
				fVMem = EvalV();
				fVMemValid = true;
			}
			return fVMem;
		}

		mutable std::vector<DExprImpl const*> fNodesMem;
		// Get the nodes list
		std::vector<DExprImpl const*> const & GetNodesMem() const {
			if (fNodesMem.size() == 0) {
				std::set<DExprImpl const *> nodes;
				AddExpressions(nodes);
				fNodesMem.reserve(nodes.size());
				for (auto &p : nodes) {
					fNodesMem.push_back(p);
				}
			}
			return fNodesMem;
		}


	};

	DExprImpl::DExprImpl(ExprType type) : fType(type), fVMemValid(false), fUid(guid++)
	{
		fRef = 0;
		auto it = DCount.find(this);
		if (it != DCount.end()) {
			printf("FATAL ERROR: you are intializing object <%p> has bee intialized! you can't continue!", (void*)*it);
			exit(1);
		}
		DCount.insert(this);
	}

	DExprImpl::DExprImpl() : DExprImpl(ExprType::NotSet)
	{
	}

	Expr DExprImpl::DMem(Var const & s) const
	{
		auto &p = fDMem[s.Uid()];
		if (p.Empty()) {
			p = EvalD(s);
		}
		return p;
	}


	DExprImpl::~DExprImpl()
	{
		auto it = DCount.find(this);
		if (it == DCount.end()) {
			printf("FATAL ERROR: you are freeing memory <%p>, but should have been freed! You can't continue!", (void*)*it);
			exit(1);
		} else {
			DCount.erase(it);
		}
		
	}


	/****************************	DScalarImpl end **********************************************/


	// put these function as back as possible
	// at lease let thsi function know DExprImpl has virtual functions
	inline std::type_info const &Typeid(Expr const &s) {
		return typeid(*s.fImpl);
	}

	/****************************	DExpr	begin *********************************************/

	Expr::Expr(DExprImpl const &impl) : fImpl(&impl) {
		fImpl->IncRef();
	}

	Expr::~Expr() {
		if (fImpl) {
			fImpl->DecRef();
		}
		fImpl = nullptr;
	}

	Expr::Expr(Expr const &s) {
		fImpl = s.fImpl;
		if (fImpl) fImpl->IncRef();
	}

	Expr::Expr(Expr &&s) {
		fImpl = s.fImpl;
		s.fImpl = nullptr;
	}

	double Expr::V() const {
		// invliad all nodes
		// we just don't know if the variable changed
		auto &nodes = fImpl->GetNodesMem();
		for (auto &p : nodes) {
			p->fVMemValid = false;
		}
		return fImpl->VMem();
	};

	Expr D(Expr const &expr, Expr const &var) {
		return expr.fImpl->DMem(CastToVar(var));
	}

	inline Expr D(Expr const &expr, Var const &var) {
		return expr.fImpl->DMem(var);
	}

	Expr D(Expr const &expr, Expr const &var, int n)
	{
		RebindableExpr result = expr;
		for (int i = 0; i < n; ++i) {
			result = D(result, var);
		}
		return result;
	}

	Expr D(Expr const &expr, std::pair<Expr, int> const &pair)
	{
		return D(expr, pair.first, pair.second);
	}


	size_t Expr::Nodes() const
	{
		return fImpl->GetNodesMem().size();
	}

	uint64_t Expr::Uid() const {
		return fImpl->fUid;
	}

	std::string Expr::ToString() const {
		std::string sb;
		fImpl->ToString(sb);
		return sb;
	}

	Expr Expr::D(Expr const & var) const
	{
		return fImpl->DMem(CastToVar(var));
	}

	/****************************	DExpr	end *********************************************/

	/****************************	RebindableExpr	begin  *********************************************/
	
	RebindableExpr::RebindableExpr() : Expr() {

	}


	RebindableExpr &RebindableExpr::operator=(Expr const &s) {
		if (&s != this) {
			if (fImpl) fImpl->DecRef();
			fImpl = s.fImpl;
			if (fImpl) fImpl->IncRef();
		}
		return *this;
	}

	RebindableExpr &RebindableExpr::operator=(Expr &&s) {
		if (&s != this) {
			if (fImpl) fImpl->DecRef();
			fImpl = s.fImpl;
			s.fImpl = nullptr;
		}
		return *this;
	}
	/****************************	RebindableExpr	end *********************************************/


	/********************	DConstant	begin	*********************************/

	struct DConstant : DExprImpl
	{
		double const fV;
		double const fE;

		DConstant(double v, double e = 0) : DExprImpl(ExprType::Const), fV(v), fE(e) { }

		double EvalV() const override {
			return fV;
		}


		void AddExpressions(std::set<DExprImpl const*> &nodes) const override {
			nodes.insert(this);
		}

		void GetParameters(ParameterVector &v) const override {
			v.push_back(fV);
			v.push_back(fE);
		}

		void GetSubExpressions(SubExpressionVector &) const override {
		}

		static std::string const sTypeName;
		std::string const &GetTypeName() const override {
			return sTypeName;
		}

		Expr EvalD(Var const &s) const override {
			return Zero();
		}

		void ToString(std::string & sb) const override
		{
			char b[100];
			sprintf(b, "%g", fV);
			sb.append(b);
		}

	};

	std::string const DConstant::sTypeName = "Constant";

	/********************	DConstant	end	*********************************/

	/********************	DConst	begin	*********************************/

	DExprImpl const *get_impl(double v)
	{
		DExprImpl const *fImpl;
		if (v == 0) { fImpl = Zero().fImpl; } 
		else if (v == 1) { fImpl = One().fImpl; }
		else if (v == 2) { fImpl = Two().fImpl; }
		else fImpl = new DConstant(v);
		return fImpl;
	}

	// try fold some values
	Const::Const(double v) : Expr(*get_impl(v)){
	}

	// not try fold some values
	Const::Const(DConstant const &p) : Expr(p) {
	}

	Const const &Zero() {
		static Const zero = Const(*new DConstant(0));
		return zero;
	}

	Const const &One() {
		static Const one = Const(*new DConstant(1));
		return one;
	}

	Const const &Two() {
		static Const two = Const(*new DConstant(2));
		return two;
	}

	/********************	DConst	end	*********************************/

	/********************	DVariableImpl	begin	*********************************/
	std::string varName = "Variable";
	struct DVariableImpl : DExprImpl
	{
		double fV;
		double fE;
		std::string const fName;
		DVariableImpl(double v, double e = 0) : DVariableImpl(std::string(), v, e) { }

		DVariableImpl(std::string const &name, double v, double e = 0) : DExprImpl(ExprType::Variable),
			fV(v), fE(e), fName(name) { }

		void SetV(double v, double e = 0) {
			fV = v;
			fE = e;
		}


		double EvalV() const override { return fV; }

		Expr EvalD(Var const &s) const override {
			if (s.fImpl == this) {
				return One();
			} else {
				return Zero();
			}
		}

		void GetSubExpressions(SubExpressionVector &) const override {
		}

		void GetParameters(ParameterVector &pars) const {
			pars.push_back(fV);
		}

		std::string const &GetTypeName() const {
			return varName;
		}

		void AddExpressions(std::set<DExprImpl const*> &nodes) const override {
			nodes.insert(this);
		}

		void ToString(std::string & sb) const override
		{
			if (fName.empty()) {
				sb.append("<").append(std::to_string((size_t)this)).append(">");
			} else {
				sb.append(fName);
			}
		}
	};

	/********************	DVariableImpl	end		*********************************/

	std::vector<Var> Expr::GetVariablesList() const
	{

		std::vector<Var> ret;
		for (auto &k : fImpl->GetNodesMem()) {
			if (k->IsVariable()) {
				auto kk = static_cast<DVariableImpl const *>(k);
				ret.emplace_back(Var(*kk));
			}
		}
		return ret;
	}

	std::string const & Expr::GetTypeName() const
	{
		return fImpl->GetTypeName();
	}

	void Expr::GetSubExpressions(SubExpressionVector& expr) const
	{
		fImpl->GetSubExpressions(expr);
	}

	void Expr::GetParameters(ParameterVector& expr) const
	{
		fImpl->GetParameters(expr);
	}


	/****************************	DVariable	begin *********************************************/


	bool IsVar(Expr const &expr) {
		return expr.fImpl->IsVariable();
	}

	bool IsConst(Expr const &expr) {
		return expr.fImpl->IsConst();
	}

	Var CastToVar(Expr const &expr) {
		// let dynamic_cast throw the exception
		return dynamic_cast<DVariableImpl const&>(*expr.fImpl);
	}

	Const CastToConst(Expr const &expr) {
		return dynamic_cast<DConstant const&>(*expr.fImpl);
	}

	Var::Var(double v) : Expr(*new DVariableImpl(v)) {
	}

	Var::Var(std::string const &name, double v) : Expr(*new DVariableImpl(name, v)) {
	}

	Var::Var(DVariableImpl const &p) : Expr(p) {
	}

	std::string const &Var::GetName() const {
		return static_cast<DVariableImpl const*>(fImpl)->fName;
	}

	void Var::SetV(double v) const {
		static_cast<DVariableImpl*>(const_cast<DExprImpl*>(fImpl))->SetV(v);
	}
	/****************************	DVariable	end *********************************************/

	/****************************	Unitary	Function *********************************************/
	struct DUnitaryFunction : DExprImpl {
		Expr f1;
		DUnitaryFunction(Expr const &s) : f1(s) {
		}
		DUnitaryFunction(Expr &&s) : f1(std::move(s)) {
		}

		void GetSubExpressions(SubExpressionVector &exprs) const override {
			exprs.push_back(f1);
		}

		void AddExpressions(std::set<DExprImpl const*> &nodes) const override {
			if (nodes.insert(this).second) {
				f1.fImpl->AddExpressions(nodes);
			}
		}


	};
	/****************************	Unitary	Function *********************************************/

	std::string const DlogName = "log";
	struct DLog : DUnitaryFunction {

		using DUnitaryFunction::DUnitaryFunction;
		double EvalV() const override {
			return std::log(f1.fImpl->VMem());
		}

		std::string const & GetTypeName() const {
			return DlogName;
		}

		Expr EvalD(Var const &s) const override {
			return D(f1, s) / f1;
		}

		void ToString(std::string & sb) const override
		{
			sb.append("log(");
			f1.fImpl->ToString(sb);
			sb.append(")");
		}

	};

	std::string const DtanName = "tan";
	struct TanImpl : DUnitaryFunction
	{

		using DUnitaryFunction::DUnitaryFunction;
		double EvalV() const override
		{
			return std::tan(f1.fImpl->VMem());
		}

		std::string const & GetTypeName() const
		{
			return DtanName;
		}

		Expr EvalD(Var const &s) const override
		{
			return pow(cos(f1), -2);
		}

		void ToString(std::string & sb) const override
		{
			sb.append("tan(");
			f1.fImpl->ToString(sb);
			sb.append(")");
		}

	};

	std::string const DpowName = "pow";

	struct DPowN : DUnitaryFunction {
		double const fN;

		DPowN(Expr const &s, double n) : DUnitaryFunction(s), fN(n) {
		}
		DPowN(Expr &&s, double n) : DUnitaryFunction(std::move(s)), fN(n) {
		}

		std::string const & GetTypeName() const {
			return DpowName;
		}

		void GetParameters(ParameterVector &v) const override {
			v.push_back(fN);
		}

		double EvalV() const override {
			if (fN == 0) return 1;
			else if (fN == 1) return f1.fImpl->VMem();
			else if (fN == 2) return f1.fImpl->VMem()*f1.fImpl->VMem();
			else if (fN == 0.5) return std::sqrt(f1.fImpl->VMem());
			else if (fN == -0.5) return 1 / std::sqrt(f1.fImpl->VMem());
			return std::pow(f1.fImpl->VMem(), fN);
		}


		Expr EvalD(Var const &s) const override {
			if (fN == 0) return Zero();
			else if (fN == 1) return D(f1, s);
			else if (fN == 2) return Two()*f1*D(f1, s);
			return Const(fN)*pow(f1, fN - 1)*D(f1, s);
		}

		void ToString(std::string & sb) const override
		{
			if (fabs(fN - 1) <= 0.000001) {
				printf("%.18f", fN);
			}
			if (f1.fImpl->IsConst() || f1.fImpl->IsVariable()) {
				f1.fImpl->ToString(sb);
			} else {
				sb.append("(");
				f1.fImpl->ToString(sb);
				sb.append(")");
			}
			sb.append("^");
			sb.append(ToStr(fN));
		}

	};

	std::string const DexpName = "exp";
	struct DExp : DUnitaryFunction
	{
		using DUnitaryFunction::DUnitaryFunction;

		double EvalV() const override {
			return std::exp(f1.fImpl->VMem());
		}

		std::string const & GetTypeName() const {
			return DexpName;
		}

		Expr EvalD(Var const &s) const override
		{
			// make a copy
			// don't refer this self
			return D(f1, s)*exp(f1);
		}

		void ToString(std::string & sb) const override
		{
			sb.append("exp(");
			f1.fImpl->ToString(sb);
			sb.append(")");
		}

	};

	std::string const DsinName = "sin";
	struct DSin : DUnitaryFunction
	{
		using DUnitaryFunction::DUnitaryFunction;


		double EvalV() const override {
			return std::sin(f1.fImpl->VMem());
		}

		std::string const & GetTypeName() const {
			return DsinName;
		}

		Expr EvalD(Var const &s) const override
		{
			return D(f1, s)*cos(f1);
		}

		void ToString(std::string & sb) const override
		{
			sb.append("sin(");
			f1.fImpl->ToString(sb);
			sb.append(")");
		}

	};

	std::string const DcosName = "cos";
	struct DCos : DUnitaryFunction
	{

		using DUnitaryFunction::DUnitaryFunction;

		double EvalV() const override {
			return std::cos(f1.fImpl->VMem());
		}

		std::string const & GetTypeName() const {
			return DcosName;
		}

		Expr EvalD(Var const &s) const override
		{
			return -D(f1, s)*sin(f1);
		}

		void ToString(std::string & sb) const override
		{
			sb.append("cos(");
			f1.fImpl->ToString(sb);
			sb.append(")");
		}

	};

	std::string const DsinhName = "sinh";
	struct DSinh : DUnitaryFunction
	{

		using DUnitaryFunction::DUnitaryFunction;
		
		double EvalV() const override {
			return std::sinh(f1.fImpl->VMem());
		}

		std::string const & GetTypeName() const {
			return DsinhName;
		}

		Expr EvalD(Var const &s) const override
		{
			return D(f1, s)*cosh(f1);
		}

		void ToString(std::string & sb) const override
		{
			sb.append("sinh(");
			f1.fImpl->ToString(sb);
			sb.append(")");
		}

	};

	std::string const DcoshName = "cosh";
	struct DCosh : DUnitaryFunction
	{
		using DUnitaryFunction::DUnitaryFunction;

		double EvalV() const override {
			return std::cosh(f1.fImpl->VMem());
		}

		std::string const & GetTypeName() const {
			return DcoshName;
		}

		Expr EvalD(Var const &s) const override
		{
			return D(f1, s)*sinh(f1);
		}

		void ToString(std::string & sb) const override
		{
			sb.append("cosh(");
			f1.fImpl->ToString(sb);
			sb.append(")");
		}

	};

	/****************************	Binary	Function begin *********************************************/
	struct DBinaryFunction : DExprImpl {
		Expr f1;
		Expr f2;

		DBinaryFunction(Expr const &s1, Expr const &s2) : f1(s1), f2(s2) {
		}
		DBinaryFunction(Expr &&s1, Expr &&s2) : f1(std::move(s1)), f2(std::move(s2)) {
		}
		DBinaryFunction(Expr const &s1, Expr &&s2) : f1(s1), f2(std::move(s2)) {
		}
		DBinaryFunction(Expr &&s1, Expr const &s2) : f1(std::move(s1)), f2(s2) {
		}

		void GetSubExpressions(SubExpressionVector &exprs) const override {
			exprs.push_back(f1);
			exprs.push_back(f2);
		}

		void AddExpressions(std::set<DExprImpl const*> &cleaned) const override {
			if (cleaned.insert(this).second) {
				f1.fImpl->AddExpressions(cleaned);
				f2.fImpl->AddExpressions(cleaned);
			}
		}

	};
	/****************************	Binary	Function end *********************************************/

	std::string const DaddhName = "add";
	struct DAdd : DBinaryFunction {
		using DBinaryFunction::DBinaryFunction;

		std::string const & GetTypeName() const {
			return DaddhName;
		}

		double EvalV() const override {
			return f1.fImpl->VMem() + f2.fImpl->VMem();
		}


		Expr EvalD(Var const &s) const override {
			return D(f1, s) + D(f2, s);
		}

		void ToString(std::string & sb) const override
		{
			f1.fImpl->ToString(sb);
			sb.append("+");
			f2.fImpl->ToString(sb);
		}

	};

	std::string const DsubhName = "sub";
	struct DSub : DBinaryFunction {

		using DBinaryFunction::DBinaryFunction;

		std::string const & GetTypeName() const {
			return DsubhName;
		}

		double EvalV() const override {
			return f1.fImpl->VMem() - f2.fImpl->VMem();
		}

		Expr EvalD(Var const &s) const override {
			return D(f1, s) - D(f2, s);
		}

		void ToString(std::string & sb) const override
		{

			if (1 || typeid(*f1.fImpl) == typeid(DSub) || typeid(*f2.fImpl) == typeid(DAdd)) {
				sb.append("(");
				f1.fImpl->ToString(sb);
				sb.append(")-(");
				f2.fImpl->ToString(sb);
				sb.append(")");
			} else {
				f1.fImpl->ToString(sb);
				sb.append("-");
				f2.fImpl->ToString(sb);
			}
		}

	};

	std::string const DmulhName = "mul";
	struct DMul : DBinaryFunction {
		using DBinaryFunction::DBinaryFunction;

		std::string const & GetTypeName() const {
			return DmulhName;
		}

		double EvalV() const override {
			return f1.fImpl->VMem() * f2.fImpl->VMem();
		}


		Expr EvalD(Var const &s) const override {
			if (f1.fImpl->IsConst()) return f1* D(f2, s);
			else if (f2.fImpl->IsConst()) return f2 * D(f1, s);
			else return D(f1, s) * f2 + f1 * D(f2, s);
		}

		void ToString(std::string & sb) const override
		{
			if (typeid(*f1.fImpl) == typeid(DSub) || typeid(*f2.fImpl) == typeid(DAdd)) {
				sb.append("(");
				f1.fImpl->ToString(sb);
				sb.append(")*(");
				f2.fImpl->ToString(sb);
				sb.append(")");
			} else {
				f1.fImpl->ToString(sb);
				sb.append("*");
				f2.fImpl->ToString(sb);
			}
		}

	};

	std::string const DdivhName = "div";
	struct DDiv : DBinaryFunction {

		using DBinaryFunction::DBinaryFunction;

		std::string const & GetTypeName() const {
			return DdivhName;
		}

		double EvalV() const override {
			return f1.fImpl->VMem() / f2.fImpl->VMem();
		}

		Expr EvalD(Var const &s) const override {
			return D(f1, s) / f2 - f1 * D(f2, s) * pow(f2, -2);
		}

		void ToString(std::string & sb) const override
		{
			if (1 || typeid(*f1.fImpl) == typeid(DSub) || typeid(*f2.fImpl) == typeid(DAdd)) {
				sb.append("(");
				f1.fImpl->ToString(sb);
				sb.append(")/(");
				f2.fImpl->ToString(sb);
				sb.append(")");
			} else {
				f1.fImpl->ToString(sb);
				sb.append("/");
				f2.fImpl->ToString(sb);
			}
		}

	};


	std::string const DinthName = "Integrate";

	struct IntegralImpl : DExprImpl {
		Var fX;
		Expr fX0;
		Expr fX1;
		Expr fY;

		IntegralImpl(Var const &s1, Expr const &s2, Expr const &s3, Expr const &s4) : fX(s1), fX0(s2), fX1(s3), fY(s4) {
		}

		std::string const & GetTypeName() const {
			return DinthName;
		}

		void GetSubExpressions(SubExpressionVector &exprs) const {
			exprs.push_back(fY);
			exprs.push_back(fX);
			exprs.push_back(fX0);
			exprs.push_back(fX1);
		}

		void AddExpressions(std::set<DExprImpl const*> &nodes) const override
		{
			if (nodes.insert(this).second) {
				fY.fImpl->AddExpressions(nodes);
				fX.fImpl->AddExpressions(nodes);
				fX0.fImpl->AddExpressions(nodes);
				fX1.fImpl->AddExpressions(nodes);
			}
		}

		double EvalV() const override
		{
			return GaussLegendre64Points([&](double x) {
				fX.SetV(x);
				return fY.V();
			}, fX0.fImpl->VMem(), fX1.fImpl->VMem());
		}

		Expr EvalD(Var const &s) const override
		{
			RebindableExpr d0;
			Expr dx0 = D(fX0, s);
			if (dx0.fImpl->IsConst() && dx0.fImpl->VMem() == 0) {
				d0 = Const(0);
			} else {
				d0 = -dx0 * ReplaceVariable(fY, fX, fX0);
			}

			RebindableExpr d1;
			Expr dx1 = D(fX1, s);
			if (dx1.fImpl->IsConst() && dx1.fImpl->VMem() == 0) {
				d1 = Const(0);
			} else {
				d1 = dx1 * ReplaceVariable(fY, fX, fX1);
			}

			Expr d2 = Integrate(D(fY, s), { fX, fX0, fX1 });

			return d0 + d1 + d2;
		}

		void ToString(std::string &sb) const override {
			sb.append("int(");
			fX.fImpl->ToString(sb);
			sb.append(",");
			fX0.fImpl->ToString(sb);
			sb.append(",");
			fX1.fImpl->ToString(sb);
			sb.append(",");
			fY.fImpl->ToString(sb);
			sb.append("£©");
		}

	};

	std::string sname = "Sum";
	struct SumImpl : DExprImpl
	{
		Var fI;
		Expr fExpr;
		double fFirst;
		double fLast;
		double fInc;

		SumImpl(Expr const &s1, Var const &i, double first, double last, double inc) : fI(i), fExpr(s1), fFirst(first), fLast(last), fInc(inc) {
		}

		void AddExpressions(std::set<DExprImpl const*> &nodes) const override
		{
			if (nodes.insert(this).second) {
				fI.fImpl->AddExpressions(nodes);
				fExpr.fImpl->AddExpressions(nodes);
			}
		}
		std::string const &GetTypeName() const override {
			return sname;
		}

		void GetSubExpressions(SubExpressionVector &expr) const {
			expr.push_back(fExpr);
			expr.push_back(fI);
		}

		void GetParameters(ParameterVector &pars) const {
			pars.push_back(fFirst);
			pars.push_back(fLast);
			pars.push_back(fInc);
		}

		double EvalV() const override
		{
			double s = 0;
			for (double i = fFirst; i <= fLast; i += fInc) {
				fI.SetV(i);
				s += fExpr.V();
			}
			return s;
		}

		Expr EvalD(Var const &s) const
		{
			return Sum(D(fExpr, s), fI, fFirst, fLast, fInc);
		}


		void ToString(std::string &sb) const override
		{
			sb.append("sum(");
			fExpr.fImpl->ToString(sb);
			sb.append(",");
			fI.fImpl->ToString(sb);
			sb.append(",").append(ToStr(fFirst));
			sb.append(",").append(ToStr(fLast));
			sb.append(",").append(ToStr(fInc));
			sb.append("£©");
		}

	};

	std::string gl64_name = "GL64PointsW";
	struct GL64Weightmpl : DExprImpl {

		Var fIndex;
		GL64Weightmpl(Var const &var) : fIndex(var) {
		}

		void AddExpressions(std::set<DExprImpl const*> &s) const override {
			if (s.insert(this).second) {
				fIndex.fImpl->AddExpressions(s);
			}
		}

		void GetSubExpressions(SubExpressionVector &v) const override {
			v.push_back(fIndex);
		}

		std::string const &GetTypeName() const override {
			return gl64_name;
		}

		double EvalV() const override {
			int idx = (int)fIndex.fImpl->VMem();
			if (idx >= 32) idx = idx - 32;
			return gl_w_64points[idx];
		}

		Expr EvalD(Var const &s) const override {
			if (s.Uid() == fIndex.Uid()) {
				throw std::logic_error("not implemented");
			}
			return Const(0);
		}
		
		void ToString(std::string &sb) const override {
			sb.append("gl64_w[");
			fIndex.fImpl->ToString(sb);
			sb.append("]");
		}
		
	};

	std::string gl64_name_x = "GL64PointsX";
	struct GL64Ximpl : DExprImpl {

		Var fIndex;
		GL64Ximpl(Var const &var) : fIndex(var) {
		}

		void AddExpressions(std::set<DExprImpl const*> &s) const override {
			if (s.insert(this).second) {
				fIndex.fImpl->AddExpressions(s);
			}
		}

		void GetSubExpressions(SubExpressionVector &v) const override {
			v.push_back(fIndex);
		}

		std::string const &GetTypeName() const override {
			return gl64_name_x;
		}

		double EvalV() const override {
			int idx = (int)fIndex.fImpl->VMem();
			if (idx >= 32) {
				idx = idx - 32;
				return -gl_x_64points[idx];
			} else {
				return gl_x_64points[idx];
			}
		}

		Expr EvalD(Var const &s) const override {
			if (s.Uid() == fIndex.Uid()) {
				throw std::logic_error("not implemented");
			}
			return Const(0);
		}
		void ToString(std::string &sb) const override {
			sb.append("gl64_x[");
			fIndex.fImpl->ToString(sb);
			sb.append("]");
		}

	};

	std::string const EqualExprName = "==";
	struct EqualExpr : DExprImpl {
		Expr fLeft;
		Expr fRight;
		EqualExpr(Expr const &a, Expr const &b) : fLeft(a), fRight(b) {
		}

		void AddExpressions(std::set<DExprImpl const*> &s) const override {
			if (s.insert(this).second) {
				fLeft.fImpl->AddExpressions(s);
				fRight.fImpl->AddExpressions(s);
			}
		}

		void GetSubExpressions(SubExpressionVector &v) const override {
			v.push_back(fLeft);
			v.push_back(fRight);
		}

		std::string const &GetTypeName() const override {
			return EqualExprName;
		}

		double EvalV() const override {
			double left = fLeft.fImpl->VMem();
			double right = fLeft.fImpl->VMem();
			return left == right;
		}

		Expr EvalD(Var const &s) const override {
			return Const(0);
		}

		void ToString(std::string &sb) const override {
			fLeft.fImpl->ToString(sb);
			sb.append("==");
			fRight.fImpl->ToString(sb);
		}

	};

	Expr operator==(ExprOrDouble const &a, ExprOrDouble const &b) {
		return *new EqualExpr(a, b);
	}


	Expr operator*(ExprOrDouble const &s1, ExprOrDouble const &s2)
	{
		if (s1.fImpl->IsConst() && s1.fImpl->EvalV() == 1) {
			return s2;
		} else if (s1.fImpl->IsConst() && s1.fImpl->EvalV() == 0) {
			return Zero();
		} else if (s2.fImpl->IsConst() && s2.fImpl->EvalV() == 1) {
			return s1;
		} else if (s2.fImpl->IsConst() && s2.fImpl->EvalV() == 0) {
			return Zero();
		} else if (s1.fImpl->IsConst() && s2.fImpl->IsConst()) {
			return Const(s1.fImpl->EvalV()*s2.fImpl->EvalV());
		} else if (s1.fImpl == s2.fImpl) {
			return pow(s1, 2);
		} else if (s1.fImpl->IsConst()) {
			if (s2.fImpl->GetTypeName() == "mul") {
				SubExpressionVector subs;
				s2.fImpl->GetSubExpressions(subs);
				if (subs.at(0).fImpl->IsConst()) { // merge constant
					double v = subs.at(0).fImpl->EvalV()*s1.fImpl->EvalV();
					return *new DMul(Const(v), subs.at(1));
				}
			}
		} else if (s2.fImpl->IsConst()) {
			if (s1.fImpl->GetTypeName() == "mul") {
				SubExpressionVector subs;
				s1.fImpl->GetSubExpressions(subs);
				if (subs.at(0).fImpl->IsConst()) { // merge constant
					double v = subs.at(0).fImpl->EvalV()*s2.fImpl->EvalV();
					return *new DMul(Const(v), subs.at(1));
				}
			}
		} else if(s2.fImpl->GetTypeName() == "mul"){
			SubExpressionVector subs;
			s2.fImpl->GetSubExpressions(subs);
			if (subs.at(0).fImpl->IsConst()) {
				return subs.at(0) * (s1 * subs.at(1));
			} else if (subs.at(1).fImpl->IsConst()) {
				return subs.at(1) * (s1 * subs.at(0));
			}
		} else if (s1.fImpl->GetTypeName() == "mul") {
			SubExpressionVector subs;
			s1.fImpl->GetSubExpressions(subs);
			if (subs.at(0).fImpl->IsConst()) {
				return subs.at(0) * (subs.at(1) * s2);
			} else if (subs.at(1).fImpl->IsConst()) {
				return subs.at(1) * (subs.at(0) * s2);
			}
		}
		return *new DMul(s1, s2);
	}

	Expr operator+(ExprOrDouble const &s1, ExprOrDouble const &s2)
	{
		if (s1.fImpl->IsConst() && s1.fImpl->EvalV() == 0) {
			return s2;
		} else if (s2.fImpl->IsConst() && s2.fImpl->EvalV() == 0) {
			return s1;
		} else if (s1.fImpl->IsConst() && s2.fImpl->IsConst()) {
			return Const(s1.fImpl->EvalV() + s2.fImpl->EvalV());
		} else if (s1.fImpl == s2.fImpl) {
			return Two()*s1;
		}
		return *new DAdd(s1, s2);
	}

	Expr operator-(ExprOrDouble const &s1, ExprOrDouble const &s2)
	{
		if (s2.fImpl->IsConst() && s2.fImpl->EvalV() == 0) {
			return s1;
		} else if (s1.fImpl->IsConst() && s2.fImpl->IsConst()) {
			return Const(s1.fImpl->EvalV() - s2.fImpl->EvalV());
		} else if (s1.fImpl == s2.fImpl) {
			return Zero();
		}
		return *new DSub(s1, s2);
	}

	Expr operator/(ExprOrDouble const &s1, ExprOrDouble const &s2)
	{
		if (s1.fImpl->IsConst() && s1.fImpl->EvalV() == 0) {
			return Zero();
		} else if (s2.fImpl->IsConst() && s2.fImpl->EvalV() == 1) {
			return s1;
		} else if (s1.fImpl->IsConst() && s2.fImpl->IsConst()) {
			return Const(s1.fImpl->EvalV() / s2.fImpl->EvalV());
		} else if (DPowN const *pw = dynamic_cast<DPowN const*>(s2.fImpl)) {
			return s1*pow(pw->f1, -pw->fN);
		} else if (s1.fImpl == s2.fImpl) {
			return One();
		}
		return *new DDiv(s1, s2);
	}


	Expr sqrt(Expr const &s)
	{
		return pow(s, 0.5);
	}

	Expr log(Expr const &s)
	{
		if (s.fImpl->IsConst()) {
			return Const(std::log(s.fImpl->EvalV()));
		}
		return *new DLog(s);
	}

	Expr exp(Expr const &s)
	{
		if (s.fImpl->IsConst()) {
			return Const(std::exp(s.fImpl->EvalV()));
		}
		return *new DExp(s);
	}

	Expr sin(Expr const &s)
	{
		if (s.fImpl->IsConst()) {
			return Const(std::sin(s.fImpl->EvalV()));
		}
		return *new DSin(s);
	}

	Expr cos(Expr const &s)
	{
		if (s.fImpl->IsConst()) {
			return Const(std::cos(s.fImpl->EvalV()));
		}
		return *new DCos(s);
	}

	Expr tan(Expr const &s)
	{
		if (s.fImpl->IsConst()) {
			return Const(std::tan(s.fImpl->EvalV()));
		}
		return *new TanImpl(s);
	}

	Expr sinh(Expr const &s)
	{
		if (s.fImpl->IsConst()) {
			return Const(std::sinh(s.fImpl->EvalV()));
		}
		return *new DSinh(s);
	}

	Expr cosh(Expr const &s)
	{
		if (s.fImpl->IsConst()) {
			return Const(std::cosh(s.fImpl->EvalV()));
		}
		return *new DCosh(s);
	}

	Expr Integrate(ExprOrDouble const &y, Expr const &x, ExprOrDouble const &from, ExprOrDouble const &to)
	{
		Var var = CastToVar(x);
		return *new IntegralImpl(var, from, to, y);
	}


	Expr Integrate(ExprOrDouble const & y, std::tuple<Expr, ExprOrDouble, ExprOrDouble> const & x)
	{
		Var var = CastToVar(std::get<0>(x));
		return *new IntegralImpl(var, std::get<1>(x), std::get<2>(x), y);
	}

	Expr Sum(Expr const &expr, Expr const &var, double first, double last, double inc)
	{
		return *new SumImpl(expr, CastToVar(var), first, last, inc);
	}

	Expr Sum(Expr const &expr, SumSecondArg const &arg)
	{
		return *new SumImpl(expr, CastToVar(arg.fExpr), arg.fFist, arg.fSecond, arg.fInc);
	}


	Expr GaussLegendre64PointsIntegrate(ExprOrDouble const &y, Expr const &x_, ExprOrDouble const &from, ExprOrDouble const &to)
	{
		Var original_x = CastToVar(x_);
		Var i("i", 0);
		Expr xi = *new GL64Ximpl(i);
		Expr x = 0.5*(1 - xi)*from + 0.5*(1 + xi)*to;
		Expr w = *new GL64Weightmpl(i);
		Expr sd = ReplaceVariable(y, original_x, x) * w;
		Expr h = 0.5 * (to - from) * Sum(sd, i, 0, 63, 1);
		return h;
	}

	Expr GaussLegendre64PointsIntegrate(ExprOrDouble const &y, std::tuple<Expr, ExprOrDouble, ExprOrDouble> const &x)
	{
		return GaussLegendre64PointsIntegrate(y, std::get<0>(x), std::get<1>(x), std::get<2>(x));
	}


	Expr pow(Expr const &s, double n)
	{
		if (s.fImpl->IsConst()) {
			return Const(std::pow(s.fImpl->EvalV(), n));
		} else if (n == 1) {
			return s;
		} else if (DPowN const * pw = dynamic_cast<DPowN const*>(s.fImpl)) {
			return pow(pw->f1, pw->fN*n);
		}
		return *new DPowN(s, n);
	}

	ExprOrDouble::ExprOrDouble(double v) : Expr(*get_impl(v))
	{
	}
}
