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

// Definition of the optimization problem
// The objective function
function y=f(x,x_new)
  y=4*x(1) - x(2)^2 - 12
endfunction

function y=df(x,x_new)
  y(1) = 4;
  y(2) = -2*x(2);
endfunction

// The constraints
function y=g(x,x_new)
  y(1) = - 10*x(1) + x(1)^2 - 10*x(2) + x(2)^2 + 34;
  y(2) = 20 - x(1)^2 - x(2)^2
endfunction
function y=dg(x,x_new)
  y(1) = -10 + 2*x(1);
  y(2) = -10 + 2*x(2);
  y(3) = -2*x(1);
  y(4) = -2*x(2);
endfunction

// The sparsity structure of the constraints
sparse_dg = [1 1; ...
             1 2; ...
             2 1; ...
             2 2];

// The Lagrangian
function y = dh(x,x_new,obj_weight,lambda,lambda_new)
  y(1) = lambda(1)*2 - lambda(2)*2;
  y(2) = -obj_weight*2 + lambda(1)*2 - lambda(2)*2;
endfunction

// The sparsity structure of the Lagrangian
sparse_dh = [1 1; ...
             2 2];

upper = [15;15];
lower = [-15;-15];
x0    = [-12;-12]; // Feasible starting point

var_lin_type(1) = 1; // Non-Linear
var_lin_type(2) = 1; // Non-Linear (this variable appears nonlinearly in the objective function or at least in one constraint)

constr_lin_type (1) = 1; // Non-Linear
constr_lin_type (2) = 1; // Non-Linear

constr_rhs(1) = 0;
constr_rhs(2) = 0;
constr_lhs(1) = -10000;
constr_lhs(2) = 0;

////////////////////////////////////////////////////////////////////////

params.hessian_approximation = "limited-memory";

[x_sol, f_sol, extra] = ipopt(x0, f, df, g, dg, sparse_dg, dh, sparse_dh, ..
  var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, [], params);

assert_checkalmostequal ( x_sol , [1.0537922  4.3462078]' , 1.e-7 );
assert_checkalmostequal ( f_sol , -26.67435344 , 1.e-7 );
assert_checkalmostequal ( extra.lambda , [0 ; 0] , 1.4);
assert_checkequal(extra.status , "Solve_Succeeded" );
assert_checkequal(extra.it_count , 17 );
assert_checktrue(extra.cpu_time > 0);
assert_checkequal(extra.fobj_eval, 25);
assert_checkequal(extra.fobj_grad_eval , 18 );
assert_checkequal(extra.constr_eval , 25 );
assert_checkequal(extra.constr_jac_eval , 19 );
assert_checkequal(extra.hess_eval , 0 );
assert_checkalmostequal (extra.dual_inf , 0,0, 1e-10);
assert_checkalmostequal(extra.constr_viol , 0, 0, 1e-10);
assert_checkalmostequal(extra.kkt_error , 0, 0, 1e-8);
assert_checkalmostequal(extra.complementarity , 0,0, 1e-8);

////////////////////////////////////////////////////////////////////////

params = struct();
params.hessian_approximation = "exact";

[x_sol, f_sol, extra] = ipopt(x0, f, df, g, dg, sparse_dg, dh, sparse_dh, ..
  var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, [], params);

assert_checkalmostequal ( x_sol , [1.0537922  4.3462078]' , 1.e-7 );
assert_checkalmostequal ( f_sol , -26.67435344 , 1.e-7 );
assert_checkalmostequal ( extra.lambda , [0 ; 0] , 1.4);
assert_checkequal(extra.status , "Solve_Succeeded" );
assert_checkequal(extra.it_count , 26 );
assert_checktrue(extra.cpu_time > 0);
assert_checkequal(extra.fobj_eval, 35);
assert_checkequal(extra.fobj_grad_eval , 27 );
assert_checkequal(extra.constr_eval , 35 );
assert_checkequal(extra.constr_jac_eval , 28 );
assert_checkequal(extra.hess_eval , 26 );
assert_checkalmostequal (extra.dual_inf , 0,0, 1e-10);
assert_checkalmostequal(extra.constr_viol , 0, 0, 1e-10);
assert_checkalmostequal(extra.kkt_error , 0, 0, 1e-8);
assert_checkalmostequal(extra.complementarity , 0,0, 1e-8);



