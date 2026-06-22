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

function y = rosenbrock_f ( x , x_new )
    y = 100.0 *(x(2)-x(1)^2)^2 + (1-x(1))^2;
endfunction

function y = rosenbrock_df ( x , x_new )
    y(1) = -400*(x(2)-x(1)^2)*x(1) - 2*(1-x(1));
    y(2) = 200*(x(2)-x(1)^2);
endfunction

// Nonlinear inequality constraints
function y = rosenbrock_g ( x , x_new )
  y = x(1)^2 + x(2)^2 - 1.5;
endfunction

function y = rosenbrock_dg ( x , x_new )
  y(1) = 2 * x(1);
  y(2) = 2 * x(2);
endfunction

function y = rosenbrock_Hf ( x , x_new )
  y = zeros(2,2);
  y(1,1) = diag(-400*x(2) + 1200*x(1).^2 + 2);
  y(2,2) = 200;
  y = y - diag(400*x(1),1) - diag(400*x(1),-1);
endfunction

function y = rosenbrock_Hg ( x , x_new )
  y = [2 0; ...
       0 2];
endfunction

// Check derivatives of F
x = [-1.9 2.0]';
[df1 , Hf1] = numderivative(list(rosenbrock_f,%t),x,[],[],H_form="blockmat");
df2 = rosenbrock_df ( x )';
Hf2 = rosenbrock_Hf ( x );
printf('Checking the derivatives of F:\n');
printf('norm(df1-df2) = %f norm(Hf1-Hf2) = %f\n', norm(df1-df2), norm(Hf1-Hf2));



// Check derivatives of G
x = [-1.9 2.0]';
[dg1 , Hg1] = numderivative(rosenbrock_g,x,[],[],H_form="blockmat");
dg2 = rosenbrock_dg ( x )';
Hg2 = rosenbrock_Hg ( x );

printf('Checking the derivatives of G:\n');
printf('norm(dg1-dg2) = %f norm(Hg1-Hg2) = %f\n', norm(dg1-dg2), norm(Hg1-Hg2));

// The sparsity structure of the constraints
sparse_dg = [1 1; ..
              1 2];

// The Hessian of the Lagrangian
function y = rosenbrock_hessian ( x , new_x, obj_weight, lambda, new_lambda )
  Hf = rosenbrock_Hf(x);
  Hg = rosenbrock_Hg(x);
  y  = obj_weight * Hf + lambda(1) * Hg;
  y=y([1 2 4]);
endfunction

// The sparsity structure of the Lagrangian
sparse_dh = [1 1; ...
             2 1; ...
             2 2];

upper = [];
lower = [];

// Not Feasible starting point
x0                  = [-1.9 2.0]';
nb_constr           = 1;
var_lin_type        = [1 1]; // Non-Linear
constr_lin_type     = 1; // Non-Linear
constr_rhs          = 0;
constr_lhs          = -%inf;

////////////////////////////////////////////////////////////////////////

params = struct();
params.hessian_approximation = "exact"
// We use a limited-bfgs approximation for the Hessian.
//params = add_param(params,"hessian_approximation","limited-memory");

[x_sol, f_sol, extra] = ipopt(x0, rosenbrock_f, rosenbrock_df,  ...
                              rosenbrock_g, rosenbrock_dg, sparse_dg, ...
                              rosenbrock_hessian, sparse_dh, ...
                              var_lin_type, constr_lin_type, ...
                              constr_rhs, constr_lhs, ...
                              lower, upper, [],params);
                              
assert_checkalmostequal(x_sol,[0.90723379674169202;0.82275515858492032])
