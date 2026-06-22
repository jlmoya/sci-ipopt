// Copyright (C) 2010 - DIGITEO - Michael Baudin


////////////////////////////////////////////////////////////////////////
// Definition of the problem
// This is Hock - Schittkowski Problem #71
// http://www.math.uni-bayreuth.de/~kschittkowski/tp_coll1.htm
// This is used as an example for Ipopt / C++ :
// http://www.coin-or.org/Ipopt/documentation/node34.html

function f = objfun ( x,x_is_new )
  f = prod(x([1 4])) * sum(x(1:3)) + x(3)
endfunction
function c = confun(x,x_is_new)
  c = [
    prod(x) - 25
    sum(x.^2)-40
  ]
endfunction

////////////////////////////////////////////////////////////////////////
// Define gradient and Hessian matrix

function gf = gradf ( x,x_is_new )
  s13 = sum(x(1:3))
  p14 = prod(x([1 4]))
  gf(1) = p14 + x(4) * s13
  gf(2) = p14
  gf(3) = p14 + 1
  gf(4) = x(1) * s13
endfunction

function H = hessf ( x)
    H = zeros(4,4)
    // Fill lower triangular part
    H(1,1) = 2*x(4)
    H(2,1) = x(4)
    H(3,1) = x(4)
    H(4,1) = 2*x(1) + x(2) + x(3)
    H(4,2) = x(1)
    H(4,3) = x(1)
    // Set upper triangular part
    H(1,2) = H(2,1)
    H(1,3) = H(3,1)
    H(1,4) = H(4,1)
    H(2,4) = H(4,2)
    H(3,4) = H(4,3)
endfunction

function y = dg1(x)
    y(1) = prod(x([2 3 4]))
    y(2) = prod(x([1 3 4]))
    y(3) = prod(x([1 2 4]))
    y(4) = prod(x([1 2 3]))
endfunction

function y = dg2(x)
    y = 2*x
endfunction

function H = Hg1(x)
    H = zeros(4,4)
    // Fill lower triangular part
    H(2,1) = x(3) * x(4)
    H(3,1) = x(2) * x(4)
    H(3,2) = x(1) * x(4)
    H(4,1) = x(2) * x(3)
    H(4,2) = x(1) * x(3)
    H(4,3) = x(1) * x(2)
    // Set upper triangular part
    H(1,2) = H(2,1)
    H(1,3) = H(3,1)
    H(2,3) = H(3,2)
    H(1,4) = H(4,1)
    H(2,4) = H(4,2)
    H(3,4) = H(4,3)
endfunction

function H = Hg2(x)
  H = zeros(4,4)
  H(1,1) = 2
  H(2,2) = 2
  H(3,3) = 2
  H(4,4) = 2
endfunction

// The Lagrangian
function y = dh(x,x_new,obj_weight,lambda,lambda_new)
  y = obj_weight * hessf ( x ) + lambda(1) * Hg1(x) + lambda(2) * Hg2(x)
endfunction

// The constraints
function y=dg(x,x_is_new)
  nbvar = size(x,"*")
  //
  i = 1
  k = (i-1)*nbvar
  y(k+1:k+nbvar) = dg1(x)
  //
  i = 2
  k = (i-1)*nbvar
  y(k+1:k+nbvar) = dg2(x)
endfunction

////////////////////////////////////////////////////////////////////////
// The starting point is feasible w.r.t. the first constraint
x0 = [1,5,5,1]'
fx0 = 16

// The sparsity structure of the constraints
// this particular Jacobian is dense
sparse_dg = [
  1 1
  1 2
  1 3
  1 4
  2 1
  2 2
  2 3
  2 4
];


// The sparsity structure of the Lagrangian
// the Hessian for this problem is actually dense
sparse_dh = [
  1 1
  1 2
  1 3
  1 4
  2 1
  2 2
  2 3
  2 4
  3 1
  3 2
  3 3
  3 4
  4 1
  4 2
  4 3
  4 4
];

lower = [1;1;1;1];
upper = [5;5;5;5];
nb_constr = 2;
var_lin_type(1) = 1; // Non-Linear
var_lin_type(2) = 1; // Non-Linear (this variable appears nonlinearly in the objective function or at least in one constraint)
var_lin_type(3) = 1; // Non-Linear
var_lin_type(4) = 1; // Non-Linear (this variable appears nonlinearly in the objective function or at least in one constraint)
constr_lin_type (1) = 1; // Non-Linear
constr_lin_type (2) = 1; // Non-Linear
constr_lhs(1) = 25;
constr_rhs(1) = %inf;
constr_lhs(2) = 40;
constr_rhs(2) = 40;

params = struct()
params.hessian_approximation="exact";
//params.print_level=12;
//    params = add_param(params,"journal_level",12);

[x_sol, f_sol, extra] = ipopt(x0, objfun, gradf, confun, dg, sparse_dg, dh, sparse_dh, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, [],params);


