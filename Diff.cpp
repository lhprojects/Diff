#include "Diff.h"
#include "Quad.h"
#include "Num.h"
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
		virtual double DoV() const = 0;
		virtual Num DoVE() const = 0;
		virtual Expr DoD(Var const &s) const = 0;
		virtual void ToString(std::string &s) const = 0;
		// expr can't have a reference to any parent of s
		// TODO: remove this in the futrue
		virtual Expr DoReplaceVariable(Var const &s, Expr const &expr) const = 0;

		virtual std::string const & GetTypeName() const = 0;
		virtual void AddExpressions(std::set<DExprImpl const*> &nodes) const = 0;
		virtual void GetSubExpressions(SubExpressionVector &expr) const = 0;
		virtual void GetParameters(SmallVector<double, 2> &pars) const { }

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
		mutable int fRef;
		void IncRef() const { ++fRef; }
		void DecRef() const {
			--fRef;
			if (fRef == 0) delete this;
		}

		Expr ReplaceVariable(Var const &s, Expr const &expr) const;

		mutable double fVMem;
		mutable Num fVEMem;
		mutable bool fVMemValid;
		mutable std::vector<DExprImpl const*> fNodesMem;
		mutable std::map<uint64_t, RebindableExpr> fDoDMem;
		mutable std::map<std::pair<uint64_t, uint64_t>, std::pair<int, RebindableExpr> > fReplaceMem;
		Expr DoDMem(Var const &s) const;
		double VMem() const;

		Num VEMem() const {
			if (!fVMemValid) {
				fVEMem = DoVE();
				fVMemValid = true;
			}
			return fVEMem;
		}

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

	DExprImpl::DExprImpl(ExprType type) : fType(type), fVMemValid(false), fUid(guid++), fVEMem(0, 0, 0)
	{
		fRef = 0;
		DCount.insert(this);
	}

	DExprImpl::DExprImpl() : DExprImpl(ExprType::NotSet)
	{
	}

	Expr DExprImpl::DoDMem(Var const & s) const
	{
#if 1
		auto &p = fDoDMem[s.Uid()];
		if (p.Empty()) {
			p = DoD(s);
		}
		return p;
#else
		return DoD(s);
#endif
	}

	double DExprImpl::VMem() const
	{
		if (!fVMemValid) {
			fVMem = DoV();
			fVMemValid = true;
		}
		return fVMem;
	}

	Expr DExprImpl::ReplaceVariable(Var const & s, Expr const &expr_) const
	{
		if (IsConst() || IsVariable()) {
			return DoReplaceVariable(s, expr_);
		}
		auto &p = fReplaceMem[std::make_pair(s.Uid(), expr_.Uid())];
		if (p.first == 0) { // ok first time to call this
			Expr expr = DoReplaceVariable(s, expr_);
			if (expr.fImpl == this) {
				// just know the result is myself, but don't save a smart pointer
				p.first = 2;
			} else {
				p.first = 1;
				p.second = expr;
			}
		}
		if (p.first == 2) return *this;
		else return p.second;
	}


	DExprImpl::~DExprImpl()
	{
		auto it = DCount.find(this);
		if (it == DCount.end()) {
			printf("FATAL ERROR: you are free memory <%p>, but should have been freed! You can't continue!", (void*)*it);
			exit(1);
		} else {
			DCount.erase(it);
		}
		
		//fReplaceMem.clear();

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
		auto &nodes = fImpl->GetNodesMem();
		for (auto &p : nodes) {
			p->fVMemValid = false;
		}
		return fImpl->VMem();
	};

	Num Expr::VE() const {
		auto &nodes = fImpl->GetNodesMem();
		for (auto &p : nodes) {
			p->fVMemValid = false;
		}
		return fImpl->VEMem();
	};


	Expr D(Expr const &expr, Expr const &var) {
		return expr.fImpl->DoDMem(CastToVar(var));
	}

	inline Expr D(Expr const &expr, Var const &var) {
		return expr.fImpl->DoDMem(var);
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


	Expr Expr::ReplaceVariable(Var const &s, Expr const &expr) const
	{
		return fImpl->ReplaceVariable(s, expr);
	}

	Expr Expr::FixVariable(Var const & s) const
	{
		return fImpl->ReplaceVariable(s, Const(s.V()));
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
		return fImpl->DoDMem(CastToVar(var));
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
		Num const fVE;

		DConstant(double v) : DExprImpl(ExprType::Const), fV(v), fVE(v, 0, 0) { }

		double DoV() const override {
			return fV;
		}

		Num DoVE() const override {
			return fVE;
		}


		void AddExpressions(std::set<DExprImpl const*> &nodes) const override {
			nodes.insert(this);
		}

		void GetParameters(SmallVector<double, 2> &v) const override {
			v.push_back(fV);
		}

		void GetSubExpressions(SmallVector<Expr, 2> &) const override {
		}

		static std::string const sTypeName;
		std::string const &GetTypeName() const override {
			return sTypeName;
		}

		Expr DoD(Var const &s) const override {
			return Zero();
		}

		Expr DoReplaceVariable(Var const &, Expr const &) const override {
			return *this;
		}

		void ToString(std::string & sb) const override
		{
			char b[100];
			sprintf(b, "%f", fV);
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
		Num fVE;
		double fV;
		std::string const fName;
		DVariableImpl(double v) : DVariableImpl(std::string(), v) { }

		DVariableImpl(std::string const &name, double v) : DExprImpl(ExprType::Variable),
			fV(v), fVE(v, 0, 0), fName(name) { }

		void SetV(double v) {
			fV = v;
			fVE = Num(v, 0, 0);
		}


		double DoV() const override { return fV; }
		Num DoVE() const override { return fVE; }

		Expr DoD(Var const &s) const override;

		void GetSubExpressions(SubExpressionVector &) const override {
		}

		std::string const &GetTypeName() const {
			return varName;
		}

		void AddExpressions(std::set<DExprImpl const*> &nodes) const override {
			nodes.insert(this);
		}

		Expr DoReplaceVariable(Var const &s, Expr const &expr) const override {
			if (this == s.fImpl) {
				return expr;
			}
			return *this;
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


	Expr DVariableImpl::DoD(Var const &s) const
	{
		if (s.fImpl == this) {
			return One();
		} else {
			return Zero();
		}
	}

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
		double DoV() const override {
			return std::log(f1.fImpl->VMem());
		}

		Num DoVE() const override {
			return log(f1.fImpl->VEMem());
		}

		Expr DoReplaceVariable(Var const &s, Expr const &expr) const override {
			auto f = f1.ReplaceVariable(s, expr);
			if (f1.fImpl == f.fImpl) {
				return *this;
			}
			return log(f);
		}

		std::string const & GetTypeName() const {
			return DlogName;
		}

		Expr DoD(Var const &s) const override {
			return D(f1, s) / f1;
		}

		void ToString(std::string & sb) const override
		{
			sb.append("log(");
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

		Expr DoReplaceVariable(Var const &s, Expr const &expr) const override {
			auto f = f1.ReplaceVariable(s, expr);
			if (f1.fImpl == f.fImpl) {
				return *this;
			}
			return pow(f, fN);
		}


		double DoV() const override {
			if (fN == 0) return 1;
			else if (fN == 1) return f1.fImpl->VMem();
			else if (fN == 2) return f1.fImpl->VMem()*f1.fImpl->VMem();
			else if (fN == 0.5) return std::sqrt(f1.fImpl->VMem());
			else if (fN == -0.5) return 1 / std::sqrt(f1.fImpl->VMem());
			return std::pow(f1.fImpl->VMem(), fN);
		}

		Num DoVE() const override {
			return pow(f1.fImpl->VEMem(), fN);
		}


		Expr DoD(Var const &s) const override {
			if (fN == 0) return Zero();
			else if (fN == 1) return D(f1, s);
			else if (fN == 2) return Two()*f1*D(f1, s);
			return Const(fN)*pow(f1, fN - 1)*D(f1, s);
		}

		void ToString(std::string & sb) const override
		{
			if (fabs(fN - 1) <= 0.000001) {
				printf("%.20f", fN);
			}
			sb.append("(");
			f1.fImpl->ToString(sb);
			sb.append(")^(").append(ToStr(fN)).append(")");
		}

	};

	std::string const DexpName = "exp";
	struct DExp : DUnitaryFunction
	{
		using DUnitaryFunction::DUnitaryFunction;

		double DoV() const override {
			return std::exp(f1.fImpl->VMem());
		}

		Num DoVE() const override {
			return exp(f1.fImpl->VEMem());
		}

		std::string const & GetTypeName() const {
			return DexpName;
		}

		Expr DoReplaceVariable(Var const &s, Expr const &expr) const override {
			auto f = f1.ReplaceVariable(s, expr);
			if (f1.fImpl == f.fImpl) {
				return *this;
			}
			return exp(f);
		}

		Expr DoD(Var const &s) const override
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


		double DoV() const override {
			return std::sin(f1.fImpl->VMem());
		}

		Num DoVE() const override {
			return sin(f1.fImpl->VEMem());
		}

		std::string const & GetTypeName() const {
			return DsinName;
		}

		Expr DoReplaceVariable(Var const &s, Expr const &expr) const override {
			auto f = f1.ReplaceVariable(s, expr);
			if (f1.fImpl == f.fImpl) {
				return *this;
			}
			return sin(f);
		}

		Expr DoD(Var const &s) const override
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

		double DoV() const override {
			return std::cos(f1.fImpl->VMem());
		}

		Num DoVE() const override {
			return cos(f1.fImpl->VEMem());
		}

		std::string const & GetTypeName() const {
			return DcosName;
		}

		Expr DoReplaceVariable(Var const &s, Expr const &expr) const override {
			auto f = f1.ReplaceVariable(s, expr);
			if (f1.fImpl == f.fImpl) {
				return *this;
			}
			return cos(f);
		}

		Expr DoD(Var const &s) const override
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
		
		double DoV() const override {
			return std::sinh(f1.fImpl->VMem());
		}

		std::string const & GetTypeName() const {
			return DsinhName;
		}

		Expr DoReplaceVariable(Var const &s, Expr const &expr) const override {
			auto f = f1.ReplaceVariable(s, expr);
			if (f1.fImpl == f.fImpl) {
				return *this;
			}
			return sin(f);
		}

		Num DoVE() const override {
			return sinh(f1.fImpl->VEMem());
		}

		Expr DoD(Var const &s) const override
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

		double DoV() const override {
			return std::cosh(f1.fImpl->VMem());
		}

		Num DoVE() const override {
			return cosh(f1.fImpl->VEMem());
		}

		std::string const & GetTypeName() const {
			return DcoshName;
		}

		Expr DoReplaceVariable(Var const &s, Expr const &expr) const override {
			auto f = f1.ReplaceVariable(s, expr);
			if (f1.fImpl == f.fImpl) {
				return *this;
			}
			return cosh(f);
		}

		Expr DoD(Var const &s) const override
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

		double DoV() const override {
			return f1.fImpl->VMem() + f2.fImpl->VMem();
		}

		Num DoVE() const override {
			return f1.fImpl->VEMem() + f2.fImpl->VEMem();
		}

		Expr DoReplaceVariable(Var const &s, Expr const &expr) const override {
			auto f1_ = f1.ReplaceVariable(s, expr);
			auto f2_ = f2.ReplaceVariable(s, expr);
			if (f1_.fImpl == f1.fImpl && f2_.fImpl == f2.fImpl) {
				return *this;
			}
			return f1_ + f2_;
		}


		Expr DoD(Var const &s) const override {
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

		double DoV() const override {
			return f1.fImpl->VMem() - f2.fImpl->VMem();
		}

		Num DoVE() const override {
			return f1.fImpl->VEMem() - f2.fImpl->VEMem();
		}

		Expr DoReplaceVariable(Var const &s, Expr const &expr) const override {
			auto f1_ = f1.ReplaceVariable(s, expr);
			auto f2_ = f2.ReplaceVariable(s, expr);
			if (f1_.fImpl == f1.fImpl && f2_.fImpl == f2.fImpl) {
				return *this;
			}
			return f1_ - f2_;
		}


		Expr DoD(Var const &s) const override {
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

		double DoV() const override {
			return f1.fImpl->VMem() * f2.fImpl->VMem();
		}

		Num DoVE() const override {
			return f1.fImpl->VEMem() * f2.fImpl->VEMem();
		}


		Expr DoD(Var const &s) const override {
			if (f1.fImpl->IsConst()) return f1* D(f2, s);
			else if (f2.fImpl->IsConst()) return f2 * D(f1, s);
			else return D(f1, s) * f2 + f1 * D(f2, s);
		}

		Expr DoReplaceVariable(Var const &s, Expr const &expr) const override {
			auto f1_ = f1.ReplaceVariable(s, expr);
			auto f2_ = f2.ReplaceVariable(s, expr);
			if (f1_.fImpl == f1.fImpl && f2_.fImpl == f2.fImpl) {
				return *this;
			}
			return f1_ * f2_;
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

		double DoV() const override {
			return f1.fImpl->VMem() / f2.fImpl->VMem();
		}


		Num DoVE() const override {
			return f1.fImpl->VEMem() / f2.fImpl->VEMem();
		}


		Expr DoD(Var const &s) const override {
			return D(f1, s) / f2 - f1 * D(f2, s) * pow(f2, -2);
		}

		Expr DoReplaceVariable(Var const &s, Expr const &expr) const override {
			auto f1_ = f1.ReplaceVariable(s, expr);
			auto f2_ = f2.ReplaceVariable(s, expr);
			if (f1_.fImpl == f1.fImpl && f2_.fImpl == f2.fImpl) {
				return *this;
			}
			return f1_ / f2_;
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

		double DoV() const
		{
			return GaussLegendre64Points([&](double x) {
				fX.SetV(x);
				return fY.V();
			}, fX0.fImpl->VMem(), fX1.fImpl->VMem());
		}


		Num DoVE() const override {
			double v = GaussLegendre64Points([&](double x) {
				fX.SetV(x);
				return fY.V();
			}, fX0.fImpl->VMem(), fX1.fImpl->VMem());
			return Num(v, 0, 0);
		}


		Expr DoReplaceVariable(Var const &s, Expr const &expr) const override {
			auto x0_ = fX0.ReplaceVariable(s, expr);
			auto x1_ = fX1.ReplaceVariable(s, expr);
			auto y_ = fY.ReplaceVariable(s, expr);
			if (x0_.fImpl == fX0.fImpl && x1_.fImpl == fX1.fImpl && y_.fImpl == fY.fImpl) {
				return *this;
			}
			return Integrate(y_, { fX, x0_, x1_ });
		}


		Expr DoD(Var const &s) const
		{
			RebindableExpr d0;
			Expr dx0 = D(fX0, s);
			if (dx0.fImpl->IsConst() && dx0.fImpl->VMem() == 0) {
				d0 = Const(0);
			} else {
				d0 = -dx0 * fY.ReplaceVariable(fX, fX0);
			}

			RebindableExpr d1;
			Expr dx1 = D(fX1, s);
			if (dx1.fImpl->IsConst() && dx1.fImpl->VMem() == 0) {
				d1 = Const(0);
			} else {
				d1 = dx1 * fY.ReplaceVariable(fX, fX1);
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

		double DoV() const
		{
			double s = 0;
			for (double i = fFirst; i <= fLast; i += fInc) {
				fI.SetV(i);
				s += fExpr.V();
			}
			return s;
		}


		Num DoVE() const override {
			Num s(0, 0, 0);
			for (double i = fFirst; i <= fLast; i += fInc) {
				fI.SetV(i);
				s = s + fExpr.VE();
			}
			return s;
		}

		Expr DoReplaceVariable(Var const &s, Expr const &expr) const override {
			auto x0_ = fExpr.ReplaceVariable(s, expr);
			if (x0_.fImpl == fExpr.fImpl) {
				return *this;
			}
			return Sum(x0_, fI, fFirst, fLast, fInc);
		}


		Expr DoD(Var const &s) const
		{
			return Sum(D(fExpr, s), fI, fFirst, fLast, fInc);
		}

		void DoToCCode(std::string &sb) const
		{
			throw std::logic_error("not implemented");
		}

		void DoToAVXCode(std::string &sb) const
		{
			throw std::logic_error("not implemented");
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

		double DoV() const override {
			int idx = (int)fIndex.fImpl->VMem();
			if (idx >= 32) idx = idx - 32;
			return gl_w_64points[idx];
		}

		Num DoVE() const { throw std::logic_error("not implemented"); }
		Expr DoD(Var const &s) const { 
			if (s.Uid() == fIndex.Uid()) {
				throw std::logic_error("not implemented");
			}
			return Const(0);
		}
		
		void ToString(std::string &sb) const override {
			sb.append("gl_w[i]");
		}
		
		// expr can't have a reference to any parent of s
		Expr DoReplaceVariable(Var const &s, Expr const &expr) const {
			if (fIndex.Uid() == expr.Uid()) {
				return *this;
			}
			return *new GL64Weightmpl(CastToVar(expr));
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

		double DoV() const override {
			int idx = (int)fIndex.fImpl->VMem();
			if (idx >= 32) {
				idx = idx - 32;
				return -gl_x_64points[idx];
			} else {
				return gl_x_64points[idx];
			}
		}

		Num DoVE() const override { throw std::logic_error("not implemented"); }
		Expr DoD(Var const &s) const {
			if (s.Uid() == fIndex.Uid()) {
				throw std::logic_error("not implemented");
			}
			return Const(0);
		}
		void ToString(std::string &sb) const override {
			sb.append("gl_xi[i]");
		}
		virtual void DoToCCode(std::string &sb) const {
			throw std::logic_error("not implemented");
		}
		virtual void DoToAVXCode(std::string &sb) const {
			throw std::logic_error("not implemented");
		}
		// expr can't have a reference to any parent of s
		virtual Expr DoReplaceVariable(Var const &s, Expr const &expr) const {
			if (fIndex.Uid() == expr.Uid()) {
				return *this;
			}
			return *new GL64Weightmpl(CastToVar(expr));
		}

	};

	Expr operator*(ExprOrDouble const &s1, ExprOrDouble const &s2)
	{
		if (s1.fImpl->IsConst() && s1.fImpl->DoV() == 1) {
			return s2;
		} else if (s1.fImpl->IsConst() && s1.fImpl->DoV() == 0) {
			return Zero();
		} else if (s2.fImpl->IsConst() && s2.fImpl->DoV() == 1) {
			return s1;
		} else if (s2.fImpl->IsConst() && s2.fImpl->DoV() == 0) {
			return Zero();
		} else if (s1.fImpl->IsConst() && s2.fImpl->IsConst()) {
			return Const(s1.fImpl->DoV()*s2.fImpl->DoV());
		} else if (s1.fImpl == s2.fImpl) {
			return pow(s1, 2);
		}
		return *new DMul(s1, s2);
	}

	Expr operator+(ExprOrDouble const &s1, ExprOrDouble const &s2)
	{
		if (s1.fImpl->IsConst() && s1.fImpl->DoV() == 0) {
			return s2;
		} else if (s2.fImpl->IsConst() && s2.fImpl->DoV() == 0) {
			return s1;
		} else if (s1.fImpl->IsConst() && s2.fImpl->IsConst()) {
			return Const(s1.fImpl->DoV() + s2.fImpl->DoV());
		} else if (s1.fImpl == s2.fImpl) {
			return Two()*s1;
		}
		return *new DAdd(s1, s2);
	}

	Expr operator-(ExprOrDouble const &s1, ExprOrDouble const &s2)
	{
		if (s2.fImpl->IsConst() && s2.fImpl->DoV() == 0) {
			return s1;
		} else if (s1.fImpl->IsConst() && s2.fImpl->IsConst()) {
			return Const(s1.fImpl->DoV() - s2.fImpl->DoV());
		} else if (s1.fImpl == s2.fImpl) {
			return Zero();
		}
		return *new DSub(s1, s2);
	}

	Expr operator/(ExprOrDouble const &s1, ExprOrDouble const &s2)
	{
		if (s1.fImpl->IsConst() && s1.fImpl->DoV() == 0) {
			return Zero();
		} else if (s2.fImpl->IsConst() && s2.fImpl->DoV() == 1) {
			return s1;
		} else if (s1.fImpl->IsConst() && s2.fImpl->IsConst()) {
			return Const(s1.fImpl->DoV() / s2.fImpl->DoV());
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
			return Const(std::log(s.fImpl->DoV()));
		}
		return *new DLog(s);
	}

	Expr exp(Expr const &s)
	{
		if (s.fImpl->IsConst()) {
			return Const(std::exp(s.fImpl->DoV()));
		}
		return *new DExp(s);
	}

	Expr sin(Expr const &s)
	{
		if (s.fImpl->IsConst()) {
			return Const(std::sin(s.fImpl->DoV()));
		}
		return *new DSin(s);
	}

	Expr cos(Expr const &s)
	{
		if (s.fImpl->IsConst()) {
			return Const(std::cos(s.fImpl->DoV()));
		}
		return *new DCos(s);
	}

	Expr sinh(Expr const &s)
	{
		if (s.fImpl->IsConst()) {
			return Const(std::sinh(s.fImpl->DoV()));
		}
		return *new DSinh(s);
	}

	Expr cosh(Expr const &s)
	{
		if (s.fImpl->IsConst()) {
			return Const(std::cosh(s.fImpl->DoV()));
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
		Var i = 0;
		Expr xi = *new GL64Ximpl(i);
		Expr x = 0.5*(1 - xi)*from + 0.5*(1 + xi)*to;
		Expr w = *new GL64Weightmpl(i);
		Expr sd = y.ReplaceVariable(original_x, x) * w;
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
			return Const(std::pow(s.fImpl->DoV(), n));
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
