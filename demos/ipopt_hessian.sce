function y=f(x,x_new)
    y=4*x(1) - x(2)^2 - 12;
end
  
function y=df(x,x_new)
    y(1) = 4;
    y(2) = -2*x(2);
end

// The constraints
function y=g(x,x_new)
    y(1) = - 10*x(1) + x(1)^2 - 10*x(2) + x(2)^2 + 34;
    y(2) = 20 - x(1)^2 - x(2)^2;
end

function y=dg(x,x_new)
    y(1) = -10 + 2*x(1)
    y(2) = -10 + 2*x(2);
    y(3) = -2*x(1);
    y(4) = -2*x(2);
end

// The sparsity structure of the constraints
sparse_dg = [1 1; 1 2; 2 1; 2 2];
//The Hessian of the Lagrangian (lower triangle non-zero terms)
//function y = dh(x,x_new,obj_weight,lambda,lambda_new)
//    y(1) = lambda(1)*2 - lambda(2)*2;
//    y(2) = -obj_weight*2 + lambda(1)*2 - lambda(2)*2;
//end
function H = dh(x,x_new,obj_weight,lambda,lambda_new)
    y(1) = lambda(1)*2 - lambda(2)*2;
    y(2) = -obj_weight*2 + lambda(1)*2 - lambda(2)*2;
    H = diag(y);
end

// The sparsity structure of the lower triangle of the Hessian of the Lagrangian
//sparse_dh = [1 1;2 2];
sparse_dh=[]; // full hessian

upper = [15;15];
lower = [-15;-15];
x0    = [-12;-12]; // Feasible starting point

var_lin_type    = [1 1]; // Non-Linear
constr_lin_type = [1 1]; // Non-Linear
constr_rhs      = [0 0];
constr_lhs      = [-10000 0];
params = struct();
params.hessian_approximation = "exact";
params.derivative_test = "second-order";
params.derivative_test_print_all = "yes";

[x_sol, f_sol, extra] = ipopt(x0, f, df, g, dg, sparse_dg, dh, sparse_dh, var_lin_type, constr_lin_type, ...
                              constr_rhs, constr_lhs, lower, upper, [], params);
disp(x_sol)
disp(extra)
