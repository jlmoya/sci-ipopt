// =============================================================================
// Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
//
// Copyright (C) 2010 - DIGITEO - Michael Baudin
// Copyright (C) 2020 - UTC - Stéphane MOTTELET
//
// This file is hereby licensed under the terms of the GNU GPL v2.0,
// For more information, see the COPYING file which you should have received
//

//<-- CLI SHELL MODE -->
//<-- NO CHECK REF -->

// Check that we can use numerical derivatives to compute the
// derivative of the functions

////////////////////////////////////////////////////////////////////////
// Inputs
function f = objfun ( x, new_x )
//pause
  f = exp(x(1))*(4*x(1)^2 + 2*x(2)^2 + 4*x(1)*x(2) + 2*x(2) + 1);
endfunction
function [c, ceq] = confun(x)
  // Nonlinear inequality constraints
  c = [
    1.5 + x(1)*x(2) - x(1) - x(2)
    -x(1)*x(2) - 10
  ];
  // Nonlinear equality constraints
  ceq = [];
endfunction
////////////////////////////////////////////////////////////////////////
// Transform the problem into a ipopt input


function y = ipopt_df ( x, new_x )
//pause
  y = numderivative ( objfun , x )
endfunction

// The constraints
function y = ipopt_g ( x, new_x )
//pause
  [c, ceq] = confun(x)
  y = c
endfunction

function y = ipopt_dg ( x, new_x )
//pause
  y = numderivative ( ipopt_g , x )
  y = matrix ( y' , 4 , 1 )
endfunction

function y = ipopt_gi ( x , i )
//pause
  [c, ceq] = confun(x)
  y = c(i)
endfunction

// The Lagrangian
function y = ipopt_dh ( x , new_x, obj_weight, lambda, new_lambda )
  [Jf , Hf] = numderivative(objfun,x)
  [Jg1 , Hg1] = numderivative(list(ipopt_gi,1),x)
  [Jg2 , Hg2] = numderivative(list(ipopt_gi,2),x)
  y = obj_weight * Hf + lambda(1) * Hg1 + lambda(2) * Hg2
  y = y([1 2 4]); // lower triangle
endfunction

// Make a starting guess
x0 = [-1 1]';

// The sparsity structure of the constraints
sparse_dg = [
  1 1
  1 2
  2 1
  2 2
];

// The sparsity structure of the Lagrangian
sparse_dh = [
  1 1
  2 1
  2 2
];

upper = [%inf %inf]';
lower = [-%inf -%inf]';
nb_constr = 2;
var_lin_type = [1 1]'; // Non-Linear
constr_lin_type = [1 1]'; // Non-Linear
constr_rhs = [0 0]';
constr_lhs = [-%inf -%inf]';

////////////////////////////////////////////////////////////////////////

[x_sol1, f_sol, extra] = ipopt(x0, objfun, ipopt_df, ipopt_g, ..
  ipopt_dg, sparse_dg, ipopt_dh, sparse_dh, var_lin_type, ..
  constr_lin_type, constr_rhs, constr_lhs, lower, upper,[]);

////////////////////////////////////////////////////////////////////////

params = struct();
params.hessian_approximation = "limited-memory";
params.jacobian_approximation = "finite-difference-values";

[x_sol2, f_sol, extra] = ipopt(x0, objfun, ipopt_df, ipopt_g, ..
  ipopt_dg, sparse_dg, ipopt_dh, sparse_dh, var_lin_type, ..
  constr_lin_type, constr_rhs, constr_lhs, lower, upper,[],params);

assert_checkalmostequal(x_sol1,x_sol2,1e-7,1e-7)
