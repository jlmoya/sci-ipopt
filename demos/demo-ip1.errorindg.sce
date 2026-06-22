
function f = objfun ( x )
  f = exp(x(1))*(4*x(1)^2 + 2*x(2)^2 + 4*x(1)*x(2) + 2*x(2) + 1)
endfunction
function c = confun(x)
  c = [
    1.5 + x(1)*x(2) - x(1) - x(2)
    -x(1)*x(2) - 10
  ]
endfunction

function gf = gradf ( x )
  gf(1) = exp(x(1))*(4*x(1)^2 + 2*x(2)^2 + 4*x(1)*x(2) + 6*x(2) + 1+ 8*x(1) );
  gf(2) = exp(x(1))*(4*x(2) + 4*x(1) + 2);
endfunction

function H = hessf ( x )
  H(1,1) = exp(x(1))*(4*x(1)^2 + 2*x(2)^2 + 4*x(1)*x(2) + 10*x(2) + 9 + 16*x(1) );
  H(1,2) = exp(x(1))*(4*x(2) + 4*x(1) + 6);
  H(2,1) = H(1,2);
  H(2,2) = exp(x(1))*4;
endfunction

function y = dg1(x)
  y = [x(2)-1;x(1)-1]
endfunction

function y = dg2(x)
  y = [-x(2);-x(1)]
endfunction

function H = Hg1(x)
  H = [
    0 1
    1 0
    ]
endfunction

function H = Hg2(x)
  H = [
    0 -1
    -1 0
    ]
endfunction


function y = dh(x,lambda,obj_weight)
  y = obj_weight * hessf ( x ) + lambda(1) * Hg1(x) + lambda(2) * Hg2(x)
endfunction

function y=dg(x)
// What a wrong line !!!
y = x.dg
endfunction

x0 = [-1,1];

sparse_dg = [
  1 1
  1 2
  2 1
  2 2
];


sparse_dh = [
  1 1
  2 2
];

upper = [15;15];
lower = [-15;-15];
x0 = [-1;1]; 
nb_constr = 2;
var_lin_type(1) = 1; 
var_lin_type(2) = 1;
constr_lin_type (1) = 1; 
constr_lin_type (2) = 1; 
constr_rhs(1) = 0;
constr_rhs(2) = 0;
constr_lhs(1) = -10000;
constr_lhs(2) = -10000;

params = init_param();
params = add_param(params,"hessian_approximation","exact");
[x_sol, f_sol, extra] = ipopt(x0, objfun, gradf, confun, dg, sparse_dg, dh, sparse_dh, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);

// The value of %pi is not known anymore !
sin(%pi)

