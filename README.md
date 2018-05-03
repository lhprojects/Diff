# Diff

The diff is a c++ library that can calculate the differential analytically and integral numerically.

## First Impression
Let's show you a sample of codes uses Diff
```
Var x = 3;
printf("%f", x.V()); // -> 3

Const k = 2;
Expr dkx2_dx = D(k*x*x, x); 

printf("%f", dkx2_dx.V()); // -> 12
v.SetV(4)
printf("%f", dkx2_dx.V()); // -> 16
```

## Basic Concept

The main classes in Diff are
* `Expr`，a handler of a `General` expression.
* `Var`，a handler of `Variable` expression. `Variable` expression `is a` `General` expression, and `Var` `is a` `Expr`.
* `Const`, a handler of `Constant` object. `Constant` expression `is a` `General` expression, and `Const` `is a` `Expr`.
* `ExprOrDoube`, a handler of `General` expression. `ExprOrDoube` is an `Expr`. `ExprOrDoube` has additional constructors.

For examples:
```
using namespace Diff;

// Create a `Variable`, and `v` is the handler of this expression
Var v = 0;
// `v2` is a handler of the expression that `v` handles.
Var v2 = v;
// Error. The value of `Variable` must be intialized.
Var v3;
// Ok. You can change the value of a `Variable`.
v.SetV(1)

// create a `Constant` expression, and `c` is the handler of this expression.
Const c = 0;
// `c2` is a handler of the expression that `c` handles.
Const c2 = 0;
// Error. You can't change the value of a `Constant`.
// c.SetV(1)

// `Expr` is a handler of `General` expression.
Expr ev = v;
Expr ec = c;

// ExprOrDouble is has a constructor that takes a `double`.
// This constructor will create a `Const` epxression,
// and this handler will handle the `Const` expression
void foo(ExprOrDouble x);
Var w = 0;
// OK.
foo(1);
// OK.
foo(w);
```

## Get More Types of Expressions

You can make more types of expressions using operators or functions. For examples: 
```
using namespace Diff;

// `x` is a handler of `Variable`.
Var x = 0;

// `x_plus_x` is a `Add` expression, `And` expression `is a` `General` expression.
Expr x_plus_x = x + x;

// `exp_x` is a `Exp` expression, `Exp` expression `is a` `General` expression.
Expr exp_x = exp(x);
```

Supported operators or functions are
```
operator+
operator-(binary)
operator-(unitary)
operator*
operator/
sin
cos
sinh
cosh
exp
log
pow(expr, double)
```
## Evaluation
Any expression is tree of expressions. For example:
```
Var x = 1;
Var y = 1;
Var z = 1;
Expr w = x * y + z;
//            w
//           / \
//         add  z
//        /   \
//       x     y
w.V()
```
The evaluation is performed at the time of invoking the `.V()`.

## Functional Programming

* Very `Expr` instance is a handler of immutable expression, or in c++ terms, a pointer to a const object. The structure of the tree can't be changed. 
* `Expr`, `Const`, `Var`, `ExprOrDouble` are immutable handlers. It means they can't rebind the other expressions after they are initialized.


## `RebindableExpr`
`Expr`, `Const`, and `Var` are  `const` handlers. It means they can't rebind the other expressions, after they are intialized. It enforces the functional programming. For example:
```
Var x1 = 1;
Var x2 = 2;
Expr y = x1;
// Error.
// y = x2;

// Error.
Expr z;
{
     Expr w = x*x;
     z = w*w;
}
// Ok with lambda.
Expr z = [&](){
  Expr w = x*x;
  return w*w;
}();
```
But functional programming is not very convenient compared to structured programming. So the `RebindableExpr` is provided. `RebindableExpr ` `is a` `Expr`. Relation between of `RebindableExpr` and `Expr` is just like that of `Expression const *` and `Expression const * const`. For exmaple:
```
Var x = ...;
RebinableExpr z;
for(int i = 0; i < n; ++i) {
     z = z*x;
}
```

## `D`
The function `D` returns a expression represents the differential. Chain rule is used to create the new expressions tree. Usage example:
```
Var x = 0;
Expr d = D(x*x, x);

x.SetV(1);
d.V() // -> 2
```

## `Integrate`
The function `Integrate` returns a expression represents the definite integral. The Integral is evaluated numerically by 65 points Gauss–Legendre algorithm. Usage exmaple:
```
Var x = 0;
Var t = 0;
Expr I = Integrate(x*x*t, {x, 0, 1});
t.SetV(1);
I.V(); // -> 0.5
```

And the differential tree is created by the chain rule:
```
d Integrate(y, {x, x_low, x_up}) / dt = y(x_up) d x_up / dt - y(x_low) d x_low / dt + Integrate(dy / dt, { x, x_low, d_up})
````

The function `GaussLegendre64PointsIntegrate` returns a expression which is sum of 65 points of integrand. i.e
```
Sum(y(x_low*(1-gs_x[i])+x_up*(1-gs[i])) * gs_w[i], {i, 0, 63})
```
The the expression tree of differntial is crated by above formula. This is useful in the case of the integration itself is evaluated precsisely, but the derivation can't be.


