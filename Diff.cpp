#include "Diff.h"
#include <map>
#include <typeinfo>
#include <cmath>

namespace Diff {

	std::set<DExprImpl*> DCount;
	Const const &Zero();
	Const const &One();
	Const const &Two();

	struct DVariableImpl;
	typedef std::set<DVariableImpl const*> DVariabaleSet;



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
		virtual void AddNode(std::set<DExprImpl const*> &nodes) const = 0;
		virtual double DoV() const = 0;
		virtual Expr DoD(Var const &s) const = 0;
		virtual void ToString(std::string &s) const = 0;
		// expr can't have a reference to any parent of s
		Expr ReplaceVariable(Var const &s, Expr const &expr) const;
		virtual Expr DoReplaceVariable(Var const &s, Expr const &expr) const = 0;
		virtual ~DExprImpl();

		DExprImpl();
		DExprImpl(ExprType type);
		ExprType const fType;
		bool IsConst() const { return fType == ExprType::Const; }
		bool IsVariable() const { return fType == ExprType::Variable; }

		mutable int fRef;
		void IncRef() const { ++fRef; }
		void DecRef() const {
			--fRef;
			if (fRef == 0) delete this;
		}

		mutable double fVMem;
		mutable bool fVMemValid;
		mutable std::vector<DExprImpl const*> fNodesMem;
		mutable std::map<DExprImpl const *, Expr*> fDoDMem;
		mutable std::map<std::pair<DVariableImpl const *, DExprImpl const *>, std::pair<int, Expr> > fReplaceMem;
		Expr DoDMem(Var const &s) const;
		double VMem() const;

		std::vector<DExprImpl const*> const & GetNodesMem() const {
			if (fNodesMem.size() == 0) {
				std::set<DExprImpl const *> nodes;
				AddNode(nodes);
				fNodesMem.reserve(nodes.size());
				for (auto &p : nodes) {
					fNodesMem.push_back(p);
				}
			}
			return fNodesMem;
		}


	};

	DExprImpl::DExprImpl(ExprType type) : fType(type), fVMemValid(false)
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
		auto &p = fDoDMem[s.fImpl];
		if (p == nullptr) {
			p = new Expr(DoD(s));
		}
		return *p;
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
		auto &p = fReplaceMem[std::make_pair((DVariableImpl*)s.fImpl, expr_.fImpl)];
		if (p.first == 0) { // ok first time to call this
			Expr expr = DoReplaceVariable(s, expr_);
			if (expr.fImpl == this) {
				// just know the result is myself, but don't have a smart pointer
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

		for (auto &p : fDoDMem) {
			delete p.second;
		}
		fDoDMem.clear();
		
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

	Expr const &Expr::operator=(Expr const &s) {
		if (&s != this) {
			if (fImpl) fImpl->DecRef();
			fImpl = s.fImpl;
			if (fImpl) fImpl->IncRef();
		}
		return *this;
	}

	Expr const &Expr::operator=(Expr &&s) {
		if (&s != this) {
			if (fImpl) fImpl->DecRef();
			fImpl = s.fImpl;
			s.fImpl = nullptr;
		}
		return *this;
	}

	double Expr::V() const {
		auto &nodes = fImpl->GetNodesMem();
		for (auto &p : nodes) {
			p->fVMemValid = false;
		}
		return fImpl->VMem();
	};

	double Expr::W() const {
		return fImpl->VMem();
	};


	Expr Expr::D(Var const &s) const {
		return fImpl->DoDMem(s);
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

	std::string Expr::ToString() const {
		std::string sb;
		fImpl->ToString(sb);
		return sb;
	}

	/****************************	DExpr	end *********************************************/

	/********************	DConstant	begin	*********************************/

	struct DConstant : DExprImpl
	{
		double const fV;
		DConstant(double v) : DExprImpl(ExprType::Const), fV(v) { }

		virtual double DoV() const override {
			return fV;
		}

		void AddNode(std::set<DExprImpl const*> &nodes) const override {
			nodes.insert(this);
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
	struct DVariableImpl : DExprImpl
	{
		double fV;
		std::string const fName;
		DVariableImpl(double v) : DExprImpl(ExprType::Variable), fV(v) { }
		DVariableImpl(std::string const &name, double v) : DExprImpl(ExprType::Variable), fV(v), fName(name) { }
		void SetV(double v) {
			fV = v;
		}

		double DoV() const override { return fV; }
		Expr DoD(Var const &s) const override;

		void AddNode(std::set<DExprImpl const*> &nodes) const override {
			nodes.insert(this);
		}

		Expr DoReplaceVariable(Var const &s, Expr const &expr) const override {
			if (this == s.fImpl) {
				return expr;
			}
			return *this;
		}

		void ToString(std::string &sb) const override;
	};


	Expr DVariableImpl::DoD(Var const &s) const
	{
		if (s.fImpl == this) {
			return One();
		} else {
			return Zero();
		}
	}

	void DVariableImpl::ToString(std::string & sb) const
	{
		if (fName.empty()) {
			sb.append("<").append(std::to_string((size_t)this)).append(">");
		} else {
			sb.append(fName);
		}
	}
	/********************	DVariableImpl	end		*********************************/

	Expr Expr::D(Expr const &s) const {
		if (!IsVar(s)) {
			printf("ERROR differential with %s, %s expected, 0 returned\n", Typeid(s).name(), typeid(DVariableImpl).name());
			return Zero();
		}

		return fImpl->DoDMem(CastToVar(s));
	}

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

		void AddNode(std::set<DExprImpl const*> &nodes) const override {
			if (nodes.insert(this).second) {
				f1.fImpl->AddNode(nodes);
			}
		}

	};
	/****************************	Unitary	Function *********************************************/


	struct DSqrt : DUnitaryFunction {

		using DUnitaryFunction::DUnitaryFunction;
		double DoV() const override {
			return std::sqrt(f1.W());
		}

		Expr DoReplaceVariable(Var const &s, Expr const &expr) const override {
			auto f = f1.ReplaceVariable(s, expr);
			if (f1.fImpl == f.fImpl) {
				return *this;
			}
			return sqrt(f);
		}

		Expr DoD(Var const &s) const override {
			// don't ref this self
			return Const(0.5)*f1.D(s) / sqrt(f1);
		}

		void ToString(std::string & sb) const override
		{
			sb.append("sqrt(");
			f1.fImpl->ToString(sb);
			sb.append(")");
		}
	};

	struct DLog : DUnitaryFunction {

		using DUnitaryFunction::DUnitaryFunction;
		double DoV() const override {
			return std::log(f1.W());
		}

		Expr DoReplaceVariable(Var const &s, Expr const &expr) const override {
			auto f = f1.ReplaceVariable(s, expr);
			if (f1.fImpl == f.fImpl) {
				return *this;
			}
			return log(f);
		}

		Expr DoD(Var const &s) const override {
			return f1.D(s) / f1;
		}

		void ToString(std::string & sb) const override
		{
			sb.append("log(");
			f1.fImpl->ToString(sb);
			sb.append(")");
		}

	};

	struct DPowN : DUnitaryFunction {
		double const fN;

		DPowN(Expr const &s, double n) : DUnitaryFunction(s), fN(n) {
		}
		DPowN(Expr &&s, double n) : DUnitaryFunction(std::move(s)), fN(n) {
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
			else if (fN == 1) return f1.W();
			else if (fN == 2) return f1.W()*f1.W();
			else if (fN == 0.5) return std::sqrt(f1.W());
			else if (fN == -0.5) return 1 / std::sqrt(f1.W());
			return std::pow(f1.W(), fN);
		}

		Expr DoD(Var const &s) const override {
			if (fN == 0) return Zero();
			else if (fN == 1) return f1.D(s);
			else if (fN == 2) return Two()*f1*f1.D(s);
			return Const(fN)*pow(f1, fN - 1)*f1.D(s);
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

	struct DExp : DUnitaryFunction
	{
		using DUnitaryFunction::DUnitaryFunction;

		double DoV() const override {
			return std::exp(f1.W());
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
			return f1.D(s)*exp(f1);
		}

		void ToString(std::string & sb) const override
		{
			sb.append("exp(");
			f1.fImpl->ToString(sb);
			sb.append(")");
		}

	};

	struct DSin : DUnitaryFunction
	{
		using DUnitaryFunction::DUnitaryFunction;


		virtual double DoV() const override {
			return std::sin(f1.W());
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
			return f1.D(s)*cos(f1);
		}

		void ToString(std::string & sb) const override
		{
			sb.append("sin(");
			f1.fImpl->ToString(sb);
			sb.append(")");
		}

	};

	struct DCos : DUnitaryFunction
	{

		using DUnitaryFunction::DUnitaryFunction;

		virtual double DoV() const override {
			return std::cos(f1.W());
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
			return -f1.D(s)*sin(f1);
		}

		void ToString(std::string & sb) const override
		{
			sb.append("cos(");
			f1.fImpl->ToString(sb);
			sb.append(")");
		}

	};

	struct DSinh : DUnitaryFunction
	{

		using DUnitaryFunction::DUnitaryFunction;
		virtual double DoV() const override {
			return std::sinh(f1.W());
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
			return f1.D(s)*cosh(f1);
		}

		void ToString(std::string & sb) const override
		{
			sb.append("sinh(");
			f1.fImpl->ToString(sb);
			sb.append(")");
		}

	};

	struct DCosh : DUnitaryFunction
	{
		using DUnitaryFunction::DUnitaryFunction;

		double DoV() const override {
			return std::cosh(f1.W());
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
			return f1.D(s)*sinh(f1);
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

		void AddNode(std::set<DExprImpl const*> &cleaned) const override {
			if (cleaned.insert(this).second) {
				f1.fImpl->AddNode(cleaned);
				f2.fImpl->AddNode(cleaned);
			}
		}

	};
	/****************************	Binary	Function end *********************************************/

	struct DAdd : DBinaryFunction {
		using DBinaryFunction::DBinaryFunction;

		double DoV() const override {
			return f1.W() + f2.W();
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
			return f1.D(s) + f2.D(s);
		}

		void ToString(std::string & sb) const override
		{
			f1.fImpl->ToString(sb);
			sb.append("+");
			f2.fImpl->ToString(sb);
		}

	};

	struct DSub : DBinaryFunction {

		using DBinaryFunction::DBinaryFunction;

		double DoV() const override {
			return f1.W() - f2.W();
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
			return f1.D(s) - f2.D(s);
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

	struct DMul : DBinaryFunction {
		using DBinaryFunction::DBinaryFunction;

		double DoV() const override {
			return f1.W() * f2.W();
		}

		Expr DoD(Var const &s) const override {
			if (f1.fImpl->IsConst()) return f1* f2.D(s);
			else if (f2.fImpl->IsConst()) return f2*f1.D(s);
			else return f1.D(s) * f2 + f1 * f2.D(s);
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

	struct DDiv : DBinaryFunction {

		using DBinaryFunction::DBinaryFunction;

		double DoV() const override {
			return f1.W() / f2.W();
		}

		Expr DoD(Var const &s) const override {
			return f1.D(s) / f2 - f1*f2.D(s) *pow(f2, -2);
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

	struct IntegralImpl : DExprImpl {
		Var fX;
		Expr fX0;
		Expr fX1;
		Expr fY;

		IntegralImpl(Var const &s1, Expr const &s2, Expr const &s3, Expr const &s4) : fX(s1), fX0(s2), fX1(s3), fY(s4) {
		}

		void AddNode(std::set<DExprImpl const*> &nodes) const override
		{
			nodes.insert(this);
			fX.fImpl->AddNode(nodes);
			fX0.fImpl->AddNode(nodes);
			fX1.fImpl->AddNode(nodes);
			fX.fImpl->AddNode(nodes);
		}

		double DoV() const
		{
			double x0 = fX0.W();
			double y0 = fX1.W();
			double dx = (y0 - x0)/ 1000;
			double h = 0;

			for (int i = 0; i < 1000; ++i) {
				double x = x0 + (i + 0.5)*dx;
				fX.SetV(x);
				double y = fY.V();
				h += y;
			}
			h *= dx;
			return h;
		}


		Expr DoReplaceVariable(Var const &s, Expr const &expr) const override {
			auto x0_ = fX0.ReplaceVariable(s, expr);
			auto x1_ = fX1.ReplaceVariable(s, expr);
			auto y_ = fY.ReplaceVariable(s, expr);
			if (x0_.fImpl == fX0.fImpl && x1_.fImpl == fX1.fImpl && y_.fImpl == fY.fImpl) {
				return *this;
			}
		}


		Expr DoD(Var const &s) const
		{
			Expr d0;
			Expr dx0 = fX0.D(s);
			if (dx0.fImpl->IsConst() && dx0.W() == 0) {
				d0 = Const(0);
			} else {
				d0 = -dx0 * fY.ReplaceVariable(fX, fX0);
			}

			Expr d1;
			Expr dx1 = fX1.D(s);
			if (dx1.fImpl->IsConst() && dx1.W() == 0) {
				d1 = Const(0);
			} else {
				d1 = dx1 * fY.ReplaceVariable(fX, fX1);
			}

			Expr d2 = Integrate(fX, fX0, fX1, fY.D(s));

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


	Expr operator*(Expr const &s1, Expr const &s2) {
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

	Expr operator+(Expr const &s1, Expr const &s2) {
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

	Expr operator-(Expr const &s1, Expr const &s2) {
		if (s2.fImpl->IsConst() && s2.fImpl->DoV() == 0) {
			return s1;
		} else if (s1.fImpl->IsConst() && s2.fImpl->IsConst()) {
			return Const(s1.fImpl->DoV() - s2.fImpl->DoV());
		} else if (s1.fImpl == s2.fImpl) {
			return Zero();
		}
		return *new DSub(s1, s2);
	}

	Expr operator/(Expr const &s1, Expr const &s2) {
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

	Expr Integrate(Expr const &x, Expr const &from, Expr const &to, Expr const &y)
	{
		return *new IntegralImpl(CastToVar(x), from, to, y);
	}

	Expr Integrate(Var const &x, Expr const &from, Expr const &to, Expr const &y)
	{
		return *new IntegralImpl(x, from, to, y);
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
}
