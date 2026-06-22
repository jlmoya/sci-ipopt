// Copyright (C) 2010 - DIGITEO - Michael Baudin

// Test IPOPT In Scilab on a Rosenbrock constrained test case

function f = rosenbrock_f ( x )
    f = 100.0 *(x(2)-x(1)^2)^2 + (1-x(1))^2
endfunction

function df = rosenbrock_df ( x )
    df(1) = - 400. * ( x(2) - x(1)**2 ) * x(1) -2. * ( 1. - x(1) )
    df(2) = 200. * ( x(2) - x(1)**2 )
endfunction

// Nonlinear inequality constraints
function g = rosenbrock_g ( x )
  g = x(1)^2 + x(2)^2 - 1.5
endfunction

function dg = rosenbrock_dg ( x )
  dg(1) = 2 * x(1)
  dg(2) = 2 * x(2)
endfunction

function Hf = rosenbrock_Hf ( x )
  Hf = zeros(2,2)
  Hf(1,1) = diag(-400*x(2) + 1200*x(1).^2 + 2)
  Hf(2,2) = 200
  Hf = Hf - diag(400*x(1),1) - diag(400*x(1),-1)
endfunction

function Hg = rosenbrock_Hg ( x )
  Hg = [
    2 0
    0 2
    ];
endfunction

// The Hessian of the Lagrangian
function y = rosenbrock_hessian ( x , lambda , obj_weight )
  Hf = rosenbrock_Hf ( x )
  Hg = rosenbrock_Hg ( x )
  y = obj_weight * Hf + lambda(1) * Hg
endfunction


// The sparsity structure of the constraints
sparse_dg = [
  1 1
  1 2
];

// The sparsity structure of the Lagrangian
sparse_dh = [
  1 1
  1 2
  2 1
  2 2
];

upper = []
lower = []
// Not Feasible starting point
x0  = [-1.9 2.0]' 
nb_constr = 1;
var_lin_type(1) = 1; // Non-Linear
var_lin_type(2) = 1; // Non-Linear
constr_lin_type (1) = 1; // Non-Linear
constr_rhs(1) = 0;
constr_lhs(1) = -%inf;

////////////////////////////////////////////////////////////////////////

params = init_param();
params = add_param(params,"hessian_approximation","exact");
[x_sol, f_sol, extra] = ipopt(x0, rosenbrock_f, rosenbrock_df, ..
  rosenbrock_g, rosenbrock_dg, sparse_dg, rosenbrock_hessian, sparse_dh, var_lin_type, ..
  constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);

// Expected solution : 0.9072,0.8228


