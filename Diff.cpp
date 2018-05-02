#include "Diff.h"
#include "Quad.h"
#include "Num.h"
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
		virtual void AddNode(std::set<DExprImpl const*> &nodes) const = 0;
		virtual double DoV() const = 0;
		virtual Num DoVE() const = 0;
		virtual Expr DoD(Var const &s) const = 0;
		virtual void ToString(std::string &s) const = 0;
		virtual void DoToCCode(std::string &sb) const = 0;
		virtual void DoToAVXCode(std::string &sb) const = 0;
		// expr can't have a reference to any parent of s
		virtual Expr DoReplaceVariable(Var const &s, Expr const &expr) const = 0;
		virtual ~DExprImpl();

		DExprImpl();
		DExprImpl(ExprType type);
		ExprType const fType;
		bool IsConst() const { return fType == ExprType::Const; }
		bool IsVariable() const { return fType == ExprType::Variable; }

		uint64_t const fUid;
		mutable int fRef;
		void IncRef() const { ++fRef; }
		void DecRef() const {
			--fRef;
			if (fRef == 0) delete this;
		}

		mutable std::string fCCodeName;
		mutable int fCCodeID;
		mutable bool fCCodeValid;
		void ToCCode(std::string &sb) const
		{
			if (!fCCodeValid) {
				DoToCCode(sb);
				fCCodeValid = true;
			}
		}
		void ToAVXCode(std::string &sb) const
		{
			if (!fCCodeValid) {
				DoToAVXCode(sb);
				fCCodeValid = true;
			}
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
				AddNode(nodes);
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

	Expr D(Expr const &expr, std::pair<Expr, int> const &pair) {
		RebindableExpr result = expr;
		for (int i = 0; i < pair.second; ++i) {
			result = D(result, pair.first);
		}
		return result;
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

	CCode Expr::ToCCode() const
	{
		auto &nodes = fImpl->GetNodesMem();
		int id = 0;
		for (auto & node : nodes) {
			node->fCCodeID = id;
			node->fCCodeValid = false;
			++id;
		}
		std::string sb;
		fImpl->ToCCode(sb);

		CCode ccode;
		ccode.Body = std::move(sb);
		for (auto &node : nodes) {
			ccode.Names[*node] = std::move(node->fCCodeName);
		}
		return std::move(ccode);
	}

	CCode Expr::ToAVXCode() const
	{
		auto &nodes = fImpl->GetNodesMem();
		int id = 0;
		for (auto & node : nodes) {
			node->fCCodeID = id;
			node->fCCodeValid = false;
			++id;
		}
		std::string sb;
		fImpl->ToAVXCode(sb);

		CCode ccode;
		ccode.Body = std::move(sb);
		for (auto &node : nodes) {
			ccode.Names[*node] = std::move(node->fCCodeName);
		}
		return std::move(ccode);
	}

	/****************************	DExpr	end *********************************************/

	/****************************	RebindableExpr	begin  *********************************************/

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


		void AddNode(std::set<DExprImpl const*> &nodes) const override {
			nodes.insert(this);
		}

		Expr DoD(Var const &s) const override {
			return Zero();
		}

		Expr DoReplaceVariable(Var const &, Expr const &) const override {
			return *this;
		}

		void DoToAVXCode(std::string &sb) const override {
			char n[1024];
			sprintf(n, "C_%d", fCCodeID);
			fCCodeName = n;

			char b[1024];
			sprintf(b, "__m256d const %s = _mm256_set1_pd(%.20E);\n",
				fCCodeName.c_str(), fV);
			sb.append(b);
		}

		void DoToCCode(std::string &sb) const override {
			char n[1024];
			sprintf(n, "C_%d", fCCodeID);
			fCCodeName = n;

			char b[1024];
			sprintf(b, "double const %s = %.20E;\n",
				fCCodeName.c_str(), fV);
			sb.append(b);
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

		void AddNode(std::set<DExprImpl const*> &nodes) const override {
			nodes.insert(this);
		}

		void DoToCCode(std::string &sb) const override {
			if (fName.empty()) {
				char b[1024];
				sprintf(b, "v_%d", fCCodeID);
				fCCodeName = b;
			} else {
				fCCodeName = fName;
			}
		}

		void DoToAVXCode(std::string &sb) const override {
			if (fName.empty()) {
				char b[1024];
				sprintf(b, "v_%d", fCCodeID);
				fCCodeName = b;
			} else {
				fCCodeName = fName;
			}
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
			return std::sqrt(f1.fImpl->VMem());
		}

		Num DoVE() const override {
			return sqrt(f1.fImpl->VEMem());
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
			return Const(0.5)*D(f1, s) / sqrt(f1);
		}
		
		void DoToCCode(std::string &sb) const
		{
			f1.fImpl->ToCCode(sb);

			char n[1024];
			sprintf(n, "sqrt_%d", fCCodeID);
			fCCodeName = n;

			char b[1024];
			sprintf(b, "double const %s = sqrt(%s);\n",
				fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			sb.append(b);
		}

		void DoToAVXCode(std::string &sb) const
		{
			f1.fImpl->ToAVXCode(sb);

			char n[1024];
			sprintf(n, "sqrt_%d", fCCodeID);
			fCCodeName = n;

			char b[1024];
			sprintf(b, "__m256d const %s = _mm256_sqrt_pd(%s);\n",
				fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			sb.append(b);
			
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

		Expr DoD(Var const &s) const override {
			return D(f1, s) / f1;
		}

		void DoToCCode(std::string &sb) const override
		{
			f1.fImpl->ToCCode(sb);

			char n[1024];
			sprintf(n, "log_%d", fCCodeID);
			fCCodeName = n;

			char b[1024];
			sprintf(b, "double const %s = log(%s);\n",
				fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			sb.append(b);
		}

		void DoToAVXCode(std::string &sb) const override
		{
			f1.fImpl->ToAVXCode(sb);

			char n[1024];
			sprintf(n, "log_%d", fCCodeID);
			fCCodeName = n;

			char b[1024];
			sprintf(b, "__m256d const %s = _mm256_log_pd(%s);\n",
				fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			
			sb.append(b);
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

		void DoToCCode(std::string &sb) const override
		{
			f1.fImpl->ToCCode(sb);

			char n[1024];
			sprintf(n, "pow_%d", fCCodeID);
			fCCodeName = n;


			char b[1024];
			if (fN == 0) {
				sprintf(b, "double const %s = 1;\n",
					fCCodeName.c_str());
			} else if (fN == 1) {
				sprintf(b, "double const %s = %s;\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == 2) {
				sprintf(b, "double const %s = %s * %s;\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == 3) {
				sprintf(b, "double const %s = %s * %s * %s;\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == 4) {
				sprintf(b, "double const %s = %s * %s * %s * %s;\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == 5) {
				sprintf(b, "double const %s = %s * %s * %s * %s * %s;\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == 0.5) {
				sprintf(b, "double const %s = sqrt(%s);\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == 1.5) {
				sprintf(b, "double const %s = sqrt(%s) * %s;\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == 2.5) {
				sprintf(b, "double const %s = sqrt(%s) * %s * %s;\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == -0.5) {
				sprintf(b, "double const %s = 1 / sqrt(%s);\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == -1.5) {
				sprintf(b, "double const %s = 1 / (sqrt(%s) * %s);\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == -2) {
				sprintf(b, "double const %s = 1 / (%s * %s);\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == -2.5) {
				sprintf(b, "double const %s = 1 / (sqrt(%s) * %s * %s);\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == -3) {
				sprintf(b, "double const %s = 1 / (%s * %s * %s);\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == -3.5) {
				sprintf(b, "double const %s = 1 / (sqrt(%s) * %s * %s * %s);\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == -4) {
				sprintf(b, "double const %s = 1 / (%s * %s * %s * %s);\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == -4.5) {
				sprintf(b, "double const %s = 1 / (sqrt(%s) * %s * %s * %s * %s);\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == -5) {
				sprintf(b, "double const %s = 1 / (%s * %s * %s * %s * %s);\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == -6) {
				sprintf(b, "double const %s = 1 / (%s * %s * %s * %s * %s * %s);\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(),
					f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == -7) {
				sprintf(b, "double const %s = 1 / (%s * %s * %s * %s * %s * %s * %s);\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(),
					f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else {
				sprintf(b, "double const %s = pow(%s, %.20E);\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), fN);
			}

			sb.append(b);
		}

		void DoToAVXCode(std::string &sb) const override
		{
			f1.fImpl->ToAVXCode(sb);

			char n[1024];
			sprintf(n, "pow_%d", fCCodeID);
			fCCodeName = n;


			char b[1024];
			if (fN == 0) {
				sprintf(b, "__m256d const %s = _mm256_set1_pd(1.);\n",
					fCCodeName.c_str());
			} else if (fN == 1) {
				sprintf(b, "__m256d const %s = %s;\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == 2) {
				sprintf(b, "__m256d const %s = _mm256_mul_pd(%s, %s);\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == 3) {
				sprintf(b, "__m256d const %s = _mm256_mul_pd(_mm256_mul_pd(%s, %s), %s);\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == 4) {
				sprintf(b, "__m256d const %s = _mm256_mul_pd(_mm256_mul_pd(%s, %s), _mm256_mul_pd(%s, %s));\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == 5) {
				sprintf(b, "__m256d const %s = _mm256_mul_pd(%s, _mm256_mul_pd(_mm256_mul_pd(%s, %s), _mm256_mul_pd(%s, %s)));\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == 0.5) {
				sprintf(b, "__m256d const %s = _mm256_sqrt_pd(%s);\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == 1.5) {
				sprintf(b, "__m256d const %s = _mm256_mul_pd(_mm256_sqrt_pd(%s), %s);\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == 2.5) {
				sprintf(b, "__m256d const %s = _mm256_mul_pd(_mm256_sqrt_pd(%s), _mm256_mul_pd(%s, %s));\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == -0.5) {
				sprintf(b, "__m256d const %s = _mm256_div_pd(_mm256_set1_pd(1.), _mm256_sqrt_pd(%s));\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == -1.5) {
				sprintf(b, "__m256d const %s = _mm256_div_pd(_mm256_set1_pd(1.), _mm256_mul_pd(_mm256_sqrt_pd(%s), %s));\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == -2) {
				sprintf(b, "__m256d const %s = _mm256_div_pd(_mm256_set1_pd(1.), _mm256_mul_pd(%s, %s));\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == -2.5) {
				sprintf(b, "__m256d const %s = _mm256_div_pd(_mm256_set1_pd(1.), _mm256_mul_pd(_mm256_mul_pd(%s, %s), _mm256_sqrt_pd(%s)));\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == -3) {
				sprintf(b, "__m256d const %s = _mm256_div_pd(_mm256_set1_pd(1.), _mm256_mul_pd(_mm256_mul_pd(%s, %s), %s));\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == -3.5) {
				sprintf(b, "__m256d const %s = _mm256_div_pd(_mm256_set1_pd(1.), _mm256_mul_pd(_mm256_mul_pd(%s, %s), _mm256_mul_pd(%s, _mm256_sqrt_pd(%s))));\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == -4) {
				sprintf(b, "__m256d const %s = _mm256_div_pd(_mm256_set1_pd(1.), _mm256_mul_pd(_mm256_mul_pd(%s, %s), _mm256_mul_pd(%s, %s)));\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == -4.5) {
				sprintf(b, "__m256d const %s = _mm256_div_pd(_mm256_set1_pd(1.), _mm256_mul_pd(_mm256_sqrt_pd(%s), _mm256_mul_pd(_mm256_mul_pd(%s, %s), _mm256_mul_pd(%s, %s))));\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == -5) {
				sprintf(b, "__m256d const %s = _mm256_div_pd(_mm256_set1_pd(1.), _mm256_mul_pd(%s, _mm256_mul_pd(_mm256_mul_pd(%s, %s), _mm256_mul_pd(%s, %s))));\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == -6) {
				sprintf(b, "__m256d const %s = _mm256_div_pd(_mm256_set1_pd(1.), "
					"_mm256_mul_pd(       _mm256_mul_pd(%s, %s),        _mm256_mul_pd(_mm256_mul_pd(%s, %s), _mm256_mul_pd(%s, %s))        )"
						");\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(),
					f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else if (fN == -7) {
				sprintf(b, "__m256d const %s = _mm256_div_pd(_mm256_set1_pd(1.), "
					"_mm256_mul_pd(            _mm256_mul_pd(_mm256_mul_pd(%s, %s), %s),                  _mm256_mul_pd(_mm256_mul_pd(%s, %s), _mm256_mul_pd(%s, %s))      )"
					");\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(),
					f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			} else {
				sprintf(b, "__m256d const %s = _mm256_pow_pd(%s, %.20E);\n",
					fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str(), fN);
			}

			sb.append(b);
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
			return std::exp(f1.fImpl->VMem());
		}

		Num DoVE() const override {
			return exp(f1.fImpl->VEMem());
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

		void DoToCCode(std::string &sb) const override
		{
			f1.fImpl->ToCCode(sb);

			char n[1024];
			sprintf(n, "exp_%d", fCCodeID);
			fCCodeName = n;

			char b[1024];
			sprintf(b, "double const %s = exp(%s);\n",
				fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			sb.append(b);
		}

		void DoToAVXCode(std::string &sb) const override
		{
			f1.fImpl->ToAVXCode(sb);

			char n[1024];
			sprintf(n, "exp_%d", fCCodeID);
			fCCodeName = n;

			char b[1024];
			sprintf(b, "__m256d const %s = _mm256_exp_pd(%s);\n",
				fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());

			sb.append(b);
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


		double DoV() const override {
			return std::sin(f1.fImpl->VMem());
		}

		Num DoVE() const override {
			return sin(f1.fImpl->VEMem());
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

		void DoToAVXCode(std::string &sb) const override
		{
			f1.fImpl->ToAVXCode(sb);

			char n[1024];
			sprintf(n, "sin_%d", fCCodeID);
			fCCodeName = n;

			char b[1024];
			sprintf(b, "__m256d const %s = _mm256_sin_pd(%s);\n",
				fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());

			sb.append(b);
		}

		void DoToCCode(std::string &sb) const override
		{
			f1.fImpl->ToCCode(sb);

			char n[1024];
			sprintf(n, "sin_%d", fCCodeID);
			fCCodeName = n;

			char b[1024];
			sprintf(b, "double const %s = sin(%s);\n",
				fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			sb.append(b);
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

		double DoV() const override {
			return std::cos(f1.fImpl->VMem());
		}

		Num DoVE() const override {
			return cos(f1.fImpl->VEMem());
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

		void DoToAVXCode(std::string &sb) const override
		{
			f1.fImpl->ToAVXCode(sb);

			char n[1024];
			sprintf(n, "cos_%d", fCCodeID);
			fCCodeName = n;

			char b[1024];
			sprintf(b, "__m256d const %s = _mm256_cos_pd(%s);\n",
				fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());

			sb.append(b);
		}


		void DoToCCode(std::string &sb) const 
		{
			f1.fImpl->ToCCode(sb);

			char n[1024];
			sprintf(n, "cos_%d", fCCodeID);
			fCCodeName = n;

			char b[1024];
			sprintf(b, "double const %s = cos(%s);\n",
				fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			sb.append(b);
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
		
		double DoV() const override {
			return std::sinh(f1.fImpl->VMem());
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

		void DoToCCode(std::string &sb) const override
		{
			f1.fImpl->ToCCode(sb);

			char n[1024];
			sprintf(n, "sinh_%d", fCCodeID);
			fCCodeName = n;

			char b[1024];
			sprintf(b, "double const %s = sinh(%s);\n",
				fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			sb.append(b);
		}

		void DoToAVXCode(std::string &sb) const override
		{
			f1.fImpl->ToAVXCode(sb);

			char n[1024];
			sprintf(n, "sinh_%d", fCCodeID);
			fCCodeName = n;

			char b[1024];
			sprintf(b, "__m256d const %s = _mm256_sinh_pd(%s);\n",
				fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());

			sb.append(b);
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
			return std::cosh(f1.fImpl->VMem());
		}

		Num DoVE() const override {
			return cosh(f1.fImpl->VEMem());
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

		void DoToAVXCode(std::string &sb) const
		{
			f1.fImpl->ToAVXCode(sb);

			char n[1024];
			sprintf(n, "cosh_%d", fCCodeID);
			fCCodeName = n;

			char b[1024];
			sprintf(b, "__m256d const %s = _mm256_cosh_pd(%s);\n",
				fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());

			sb.append(b);
		}

		void DoToCCode(std::string &sb) const
		{
			f1.fImpl->ToCCode(sb);

			char n[1024];
			sprintf(n, "cosh_%d", fCCodeID);
			fCCodeName = n;

			char b[1024];
			sprintf(b, "double const %s = cosh(%s);\n",
				fCCodeName.c_str(), f1.fImpl->fCCodeName.c_str());
			sb.append(b);
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

		void DoToCCode(std::string &sb) const
		{
			f1.fImpl->ToCCode(sb);
			f2.fImpl->ToCCode(sb);

			char n[1024];
			sprintf(n, "add_%d", fCCodeID);
			fCCodeName = n;

			char b[1024];
			sprintf(b, "double const %s = %s + %s;\n",
				fCCodeName.c_str(),
				f1.fImpl->fCCodeName.c_str(),
				f2.fImpl->fCCodeName.c_str()
				);
			sb.append(b);
		}

		void DoToAVXCode(std::string &sb) const
		{
			f1.fImpl->ToAVXCode(sb);
			f2.fImpl->ToAVXCode(sb);

			char n[1024];
			sprintf(n, "add_%d", fCCodeID);
			fCCodeName = n;

			char b[1024];
			sprintf(b, "__m256d const %s = _mm256_add_pd(%s, %s);\n",
				fCCodeName.c_str(),
				f1.fImpl->fCCodeName.c_str(),
				f2.fImpl->fCCodeName.c_str()
				);
			sb.append(b);
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

		void DoToCCode(std::string &sb) const
		{
			f1.fImpl->ToCCode(sb);
			f2.fImpl->ToCCode(sb);

			char n[1024];
			sprintf(n, "sub_%d", fCCodeID);
			fCCodeName = n;

			char b[1024];
			sprintf(b, "double const %s = %s - %s;\n",
				fCCodeName.c_str(),
				f1.fImpl->fCCodeName.c_str(),
				f2.fImpl->fCCodeName.c_str()
				);
			sb.append(b);
		}

		void DoToAVXCode(std::string &sb) const  override
		{
			f1.fImpl->ToAVXCode(sb);
			f2.fImpl->ToAVXCode(sb);

			char n[1024];
			sprintf(n, "sub_%d", fCCodeID);
			fCCodeName = n;

			char b[1024];
			sprintf(b, "__m256d const %s = _mm256_sub_pd(%s, %s);\n",
				fCCodeName.c_str(),
				f1.fImpl->fCCodeName.c_str(),
				f2.fImpl->fCCodeName.c_str()
				);
			sb.append(b);
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

		void DoToAVXCode(std::string &sb) const override
		{
			f1.fImpl->ToAVXCode(sb);
			f2.fImpl->ToAVXCode(sb);

			char n[1024];
			sprintf(n, "mul_%d", fCCodeID);
			fCCodeName = n;

			char b[1024];
			sprintf(b, "__m256d const %s = _mm256_mul_pd(%s, %s);\n",
				fCCodeName.c_str(),
				f1.fImpl->fCCodeName.c_str(),
				f2.fImpl->fCCodeName.c_str()
				);
			sb.append(b);
		}


		void DoToCCode(std::string &sb) const
		{
			f1.fImpl->ToCCode(sb);
			f2.fImpl->ToCCode(sb);

			char n[1024];
			sprintf(n, "mul_%d", fCCodeID);
			fCCodeName = n;

			char b[1024];
			sprintf(b, "double const %s = %s * %s;\n",
				fCCodeName.c_str(),
				f1.fImpl->fCCodeName.c_str(),
				f2.fImpl->fCCodeName.c_str()
				);
			sb.append(b);
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

		void DoToCCode(std::string &sb) const override
		{
			f1.fImpl->ToCCode(sb);
			f2.fImpl->ToCCode(sb);

			char n[1024];
			sprintf(n, "div_%d", fCCodeID);
			fCCodeName = n;

			char b[1024];
			sprintf(b, "double const %s = %s / %s;\n",
				fCCodeName.c_str(),
				f1.fImpl->fCCodeName.c_str(),
				f2.fImpl->fCCodeName.c_str()
				);
			sb.append(b);
		}

		void DoToAVXCode(std::string &sb) const override
		{
			f1.fImpl->ToAVXCode(sb);
			f2.fImpl->ToAVXCode(sb);

			char n[1024];
			sprintf(n, "div_%d", fCCodeID);
			fCCodeName = n;

			char b[1024];
			sprintf(b, "__m256d const %s = _mm256_div_pd(%s, %s);\n",
				fCCodeName.c_str(),
				f1.fImpl->fCCodeName.c_str(),
				f2.fImpl->fCCodeName.c_str()
				);
			sb.append(b);
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

		IntegralImpl(Var const &s1, Expr const &s2, Expr const &s3, Expr const &s4) : fX(0.0), fX0(s2), fX1(s3), fY(s4.ReplaceVariable(s1, fX)) {
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
			return Integrate(fX, x0_, x1_, y_);
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

			Expr d2 = Integrate(fX, fX0, fX1, D(fY, s));

			return d0 + d1 + d2;
		}

		void DoToCCode(std::string &sb) const
		{
			throw std::logic_error("not implemented");
		}
		void DoToAVXCode(std::string &sb) const
		{
			throw std::logic_error("not implemented");
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

	Expr Integrate(Var const &x, Expr const &from, Expr const &to, Expr const &y)
	{
		return *new IntegralImpl(x, from, to, y);
	}

	Expr GaussLegendre64PointsIntegrate(Expr const &x_, Expr const &from, Expr const &to, Expr const &y)
	{
		Var x = CastToVar(x_);
		RebindableExpr h = Const(0);
		for (int i = 0; i < 32; ++i) {
			{
				Expr xi = 0.5*(1 + gl_x_64points[i])*from + 0.5*(1 - gl_x_64points[i])*to;
				h = h + y.fImpl->ReplaceVariable(x, xi) * gl_w_64points[i];
				//printf("%d %f %f %f %f %f\n", i, xi.V(), gl_w_64points[i], xi.V()*xi.V(), y.fImpl->ReplaceVariable(x, xi).V(), h.V());
			}
			{
				Expr xi = 0.5*(1 - gl_x_64points[i])*from + 0.5*(1 + gl_x_64points[i])*to;
				h = h + y.fImpl->ReplaceVariable(x, xi) * gl_w_64points[i];
				//printf("%d %f %f %f %f %f\n", i, xi.V(), gl_w_64points[i], xi.V()*xi.V(), y.fImpl->ReplaceVariable(x, xi).V(), h.V());
			}
		}
		h = h * (to - from) / 2;

		return h;


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
