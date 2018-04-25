#include "Diff.h"
#include <map>
#include <typeinfo>
#include <math.h>

std::set<DExprImpl*> DCount;

struct DVariableImpl;
typedef std::set<DVariableImpl const*> DVariabaleSet;

DConst const &DZero() {
	static DConst zero = DConst(0);
	return zero;
}

DConst const &DOne() {
	static DConst one = DConst(1);
	return one;
}

DConst const &DTwo() {
	static DConst zero = DConst(2);
	return zero;
}


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
	virtual DExpr DoD(DVariable const &s) const = 0;
	virtual void ToString(std::string &s) const = 0;
	virtual DExpr DoFixVariable(DVariable const &s) const = 0;
	virtual void CountVariables(DVariabaleSet &vs) const = 0;
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
	mutable std::set<DExprImpl const*> fNodesMem;
	mutable std::map<DExprImpl const *, DExpr*> fDoDMem;
	mutable std::map<DExprImpl const *, DExpr*> fFixMem;
	DExpr DoDMem(DVariable const &s) const;
	DExpr FixMem(DVariable const &s) const;
	double VMem() const;

	std::set<DExprImpl const*> const & GetNodesMem() const {
		if (fNodesMem.size() == 0) {
			AddNode(fNodesMem);
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

DExpr DExprImpl::DoDMem(DVariable const & s) const
{
#if 1
	auto &p = fDoDMem[s.fImpl];
	if (p == nullptr) {
		p = new DExpr(DoD(s));
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

DExpr DExprImpl::FixMem(DVariable const & s) const
{
#if 1
	// these two function will return this
	// do not memory it
	if (!this->IsConst() && !this->IsVariable()) {
		auto &p = fFixMem[s.fImpl];
		if (p == nullptr) {
			p = new DExpr(DoFixVariable(s));
		}
		return *p;
	}
	return DoFixVariable(s);
#else
	return FixVariable(s);
#endif
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
	for (auto &p : fFixMem) {
		delete p.second;
	}
	fFixMem.clear();

}


/****************************	DScalarImpl end **********************************************/


// put these function as back as possible
// at lease let thsi function know DExprImpl has virtual functions
inline std::type_info const &Typeid(DExpr const &s) {
	return typeid(*s.fImpl);
}

/****************************	DExpr	begin *********************************************/

DExpr::DExpr(DExprImpl const*impl) : fImpl(impl) {
	impl->IncRef();
}

DExpr::~DExpr() {
	if (fImpl) {
		fImpl->DecRef();
	}
	fImpl = nullptr;
}

DExpr::DExpr(DExpr const &s) {
	fImpl = s.fImpl;
	if (fImpl) fImpl->IncRef();
}

DExpr::DExpr(DExpr &&s) {
	fImpl = s.fImpl;
	s.fImpl = nullptr;
}

DExpr const &DExpr::operator=(DExpr const &s) {
	if (&s != this) {
		if (fImpl) fImpl->DecRef();
		fImpl = s.fImpl;
		if (fImpl) fImpl->IncRef();
	}
	return *this;
}

DExpr const &DExpr::operator=(DExpr &&s) {
	if (&s != this) {
		if (fImpl) fImpl->DecRef();
		fImpl = s.fImpl;
		s.fImpl = nullptr;
	}
	return *this;
}

double DExpr::V() const {
	auto &nodes  = fImpl->GetNodesMem();
	for (auto &p : nodes) {
		p->fVMemValid = false;
	}
	return fImpl->VMem();
};

double DExpr::W() const {
	return fImpl->VMem();
};


DExpr DExpr::D(DVariable const &s) const {
	return fImpl->DoDMem(s);
}

size_t DExpr::Nodes() const
{
	return fImpl->GetNodesMem().size();
}

std::vector<DVariable> DExpr::GetVariablesList() const
{
	std::vector<DVariable> ret;
	DVariabaleSet vs;
	fImpl->CountVariables(vs);
	for (auto &k : vs) {
		ret.emplace_back(k);
	}
	return ret;
}

DExpr DExpr::FixVariable(DVariable const & s) const
{
	return fImpl->FixMem(s);
}

std::string DExpr::ToString() const {
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

	DExpr DoD(DVariable const &s) const override {
		return DZero();
	}

	void CountVariables(DVariabaleSet &) const override {
	}

	DExpr DoFixVariable(DVariable const &) const override {
		return this;
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

inline DConst::DConst(double v) : DExpr(new DConstant(v)) {
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
	DExpr DoD(DVariable const &s) const override;

	void AddNode(std::set<DExprImpl const*> &nodes) const override {
		nodes.insert(this);
	}

	void CountVariables(DVariabaleSet &vs) const override {
		vs.insert(this);
	}

	DExpr DoFixVariable(DVariable const &sb) const override {
		if (this == sb.fImpl) {
			return DConst(fV);
		} else {
			return this;
		}
	}
	void ToString(std::string &sb) const override;
};


DExpr DVariableImpl::DoD(DVariable const &s) const
{
	if (s.fImpl == this) {
		return DOne();
	} else {
		return DZero();
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

DExpr DExpr::D(DExpr const &s) const {
	if (!s.fImpl->IsVariable()) {
		printf("ERROR differential with %s, %s expected, 0 returned\n", Typeid(s).name(), typeid(DVariableImpl).name());
		return DZero();
	}

	return fImpl->DoDMem((DVariable const &)s);
}

/****************************	DVariable	begin *********************************************/

DVariable::DVariable(double v) : DExpr(new DVariableImpl(v)) {
}

DVariable::DVariable(std::string const &name, double v) : DExpr(new DVariableImpl(name, v)) {
}

DVariable::DVariable(DVariableImpl const *p) : DExpr(p) { }

std::string const &DVariable::GetName() const {
	return static_cast<DVariableImpl const*>(fImpl)->fName;
}

void DVariable::SetV(double v) {
	static_cast<DVariableImpl*>(const_cast<DExprImpl*>(fImpl))->SetV(v);
}
/****************************	DVariable	end *********************************************/

/****************************	Unitary	Function *********************************************/
struct DUnitaryFunction : DExprImpl {
	DExpr f1;
	DUnitaryFunction(DExpr const &s) : f1(s) {
	}
	DUnitaryFunction(DExpr &&s) : f1(std::move(s)) {
	}

	void AddNode(std::set<DExprImpl const*> &nodes) const override {
		if (nodes.insert(this).second) {
			f1.fImpl->AddNode(nodes);
		}
	}

	void CountVariables(DVariabaleSet &vs) const override {
		f1.fImpl->CountVariables(vs);
	}

};
/****************************	Unitary	Function *********************************************/


struct DSqrt : DUnitaryFunction {

	using DUnitaryFunction::DUnitaryFunction;
	double DoV() const override {
		return sqrt(f1.W());
	}

	DExpr DoFixVariable(DVariable const &s) const override {
		return sqrt(f1.FixVariable(s));
	}

	DExpr DoD(DVariable const &s) const override {
		// don't ref this self
		return DConst(0.5)*f1.D(s)/sqrt(f1);
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
		return log(f1.W());
	}

	DExpr DoFixVariable(DVariable const &s) const override {
		return log(f1.FixVariable(s));
	}

	DExpr DoD(DVariable const &s) const override {
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

	DPowN(DExpr const &s, double n) : DUnitaryFunction(s), fN(n) {
	}
	DPowN(DExpr &&s, double n) : DUnitaryFunction(std::move(s)), fN(n) {
	}

	DExpr DoFixVariable(DVariable const &s) const override {
		return pow(f1.FixVariable(s), fN);
	}

	double DoV() const override {
		if (fN == 0) return 1;
		else if (fN == 1) return f1.W();
		else if (fN == 2) return f1.W()*f1.W();
		else if (fN == 0.5) return sqrt(f1.W());
		else if (fN == -0.5) return 1/sqrt(f1.W());
		return pow(f1.W(), fN);
	}

	DExpr DoD(DVariable const &s) const override {
		if (fN == 0) return DZero();
		else if (fN == 1) return f1.D(s);
		else if (fN == 2) return DTwo()*f1*f1.D(s);
		return DConst(fN)*pow(f1, fN - 1)*f1.D(s);
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
		return exp(f1.W());
	}

	DExpr DoFixVariable(DVariable const &s) const override {
		return exp(f1.FixVariable(s));
	}

	DExpr DoD(DVariable const &s) const override
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
		return sin(f1.W());
	}

	DExpr DoFixVariable(DVariable const &s) const override {
		return sin(f1.FixVariable(s));
	}

	DExpr DoD(DVariable const &s) const override
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
		return cos(f1.W());
	}

	DExpr DoFixVariable(DVariable const &s) const override {
		return cos(f1.FixVariable(s));
	}

	DExpr DoD(DVariable const &s) const override
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
		return sinh(f1.W());
	}

	DExpr DoFixVariable(DVariable const &s) const override {
		return sinh(f1.FixVariable(s));
	}

	DExpr DoD(DVariable const &s) const override
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
		return cosh(f1.W());
	}

	DExpr DoFixVariable(DVariable const &s) const override {
		return cosh(f1.FixVariable(s));
	}

	DExpr DoD(DVariable const &s) const override
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

/****************************	Binary	Function *********************************************/
struct DBinaryFunction : DExprImpl {
	DExpr f1;
	DExpr f2;

	DBinaryFunction(DExpr const &s1, DExpr const &s2) : f1(s1), f2(s2) {
	}
	DBinaryFunction(DExpr &&s1, DExpr &&s2) : f1(std::move(s1)), f2(std::move(s2)) {
	}
	DBinaryFunction(DExpr const &s1, DExpr &&s2) : f1(s1), f2(std::move(s2)) {
	}
	DBinaryFunction(DExpr &&s1, DExpr const &s2) : f1(std::move(s1)), f2(s2) {
	}

	void AddNode(std::set<DExprImpl const*> &cleaned) const override {
		if (cleaned.insert(this).second) {
			f1.fImpl->AddNode(cleaned);
			f2.fImpl->AddNode(cleaned);
		}
	}

	void CountVariables(DVariabaleSet &vs) const override {
		f1.fImpl->CountVariables(vs);
		f2.fImpl->CountVariables(vs);
	}

};
/****************************	Binary	Function *********************************************/

struct DAdd : DBinaryFunction {
	using DBinaryFunction::DBinaryFunction;

	double DoV() const override {
		return f1.W() + f2.W();
	}


	DExpr DoFixVariable(DVariable const &s) const override {
		return f1.FixVariable(s) + f2.FixVariable(s);
	}

	DExpr DoD(DVariable const &s) const override {
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

	DExpr DoFixVariable(DVariable const &s) const override {
		return f1.FixVariable(s) - f2.FixVariable(s);
	}


	DExpr DoD(DVariable const &s) const override {
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

	DExpr DoD(DVariable const &s) const override {
		if (f1.fImpl->IsConst()) return f1* f2.D(s);
		else if (f2.fImpl->IsConst()) return f2*f1.D(s);
		else return f1.D(s) * f2 + f1 * f2.D(s);
	}

	DExpr DoFixVariable(DVariable const &s) const override {
		return f1.FixVariable(s) * f2.FixVariable(s);
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

	DExpr DoD(DVariable const &s) const override {
		return f1.D(s) / f2 - f1*f2.D(s) *pow(f2, -2);
	}

	DExpr DoFixVariable(DVariable const &s) const override {
		return f1.FixVariable(s) / f2.FixVariable(s);
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


DExpr operator*(DExpr const &s1, DExpr const &s2) {
	if (s1.fImpl->IsConst() && s1.fImpl->DoV() == 1) {
		return s2;
	} else if (s1.fImpl->IsConst() && s1.fImpl->DoV() == 0) {
		return DZero();
	} else if (s2.fImpl->IsConst() && s2.fImpl->DoV() == 1) {
		return s1;
	} else if (s2.fImpl->IsConst() && s2.fImpl->DoV() == 0) {
		return DZero();
	} else if (s1.fImpl->IsConst() && s2.fImpl->IsConst()) {
		return DConst(s1.fImpl->DoV()*s2.fImpl->DoV());
	} else if (s1.fImpl == s2.fImpl) {
		return pow(s1, 2);
	}
	return new DMul(s1, s2);
}

DExpr operator+(DExpr const &s1, DExpr const &s2) {
	if (s1.fImpl->IsConst() && s1.fImpl->DoV() == 0) {
		return s2;
	} else if (s2.fImpl->IsConst() && s2.fImpl->DoV() == 0) {
		return s1;
	} else if (s1.fImpl->IsConst()  && s2.fImpl->IsConst()) {
		return DConst(s1.fImpl->DoV() + s2.fImpl->DoV());
	} else if (s1.fImpl == s2.fImpl) {
		return DTwo()*s1;
	}
	return new DAdd(s1, s2);
}

DExpr operator-(DExpr const &s1, DExpr const &s2) {
	if (s2.fImpl->IsConst() && s2.fImpl->DoV() == 0) {
		return s1;
	} else if (s1.fImpl->IsConst() && s2.fImpl->IsConst()) {
		return DConst(s1.fImpl->DoV() - s2.fImpl->DoV());
	} else if (s1.fImpl == s2.fImpl) {
		return DZero();
	}
	return new DSub(s1, s2);
}

DExpr operator/(DExpr const &s1, DExpr const &s2) {
	if (s1.fImpl->IsConst() && s1.fImpl->DoV() == 0) {
		return DZero();
	} else if (s2.fImpl->IsConst() && s2.fImpl->DoV() == 1) {
		return s1;
	} else if (s1.fImpl->IsConst() && s2.fImpl->IsConst()) {
		return DConst(s1.fImpl->DoV() / s2.fImpl->DoV());
	} else if (DPowN const *pw = dynamic_cast<DPowN const*>(s2.fImpl)) {
		return s1*pow(pw->f1, -pw->fN);
	} else if (s1.fImpl == s2.fImpl) {
		return DOne();
	}
	return new DDiv(s1, s2);
}


DExpr sqrt(DExpr const &s)
{
	return pow(s, 0.5);
}

DExpr log(DExpr const &s)
{
	if (s.fImpl->IsConst()) {
		return DConst(log(s.fImpl->DoV()));
	}
	return new DLog(s);
}

DExpr exp(DExpr const &s)
{
	if (s.fImpl->IsConst()) {
		return DConst(exp(s.fImpl->DoV()));
	}
	return new DExp(s);
}

DExpr sin(DExpr const &s)
{
	if (s.fImpl->IsConst()) {
		return DConst(sin(s.fImpl->DoV()));
	}
	return new DSin(s);
}

DExpr cos(DExpr const &s)
{
	if (s.fImpl->IsConst()) {
		return DConst(cos(s.fImpl->DoV()));
	}
	return new DCos(s);
}

DExpr sinh(DExpr const &s)
{
	if (s.fImpl->IsConst()) {
		return DConst(sinh(s.fImpl->DoV()));
	}
	return new DSinh(s);
}

DExpr cosh(DExpr const &s)
{
	if (s.fImpl->IsConst()) {
		return DConst(cosh(s.fImpl->DoV()));
	}
	return new DCosh(s);
}

DExpr pow(DExpr const &s, double n)
{
	if (s.fImpl->IsConst()) {
		return DConst(pow(s.fImpl->DoV(), n));
	} else if (n == 1) {
		return s;
	} else if (DPowN const * pw = dynamic_cast<DPowN const*>(s.fImpl)) {
		return pow(pw->f1, pw->fN*n);
	}
	return new DPowN(s, n);
}
