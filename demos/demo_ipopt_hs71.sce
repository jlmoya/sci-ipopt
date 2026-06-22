// Copyright (C) 2010 - DIGITEO - Michael Baudin


////////////////////////////////////////////////////////////////////////
// Definition of the problem
// This is Hock - Schittkowski Problem #71
// http://www.math.uni-bayreuth.de/~kschittkowski/tp_coll1.htm
// This is used as an example for Ipopt / C++ :
// http://www.coin-or.org/Ipopt/documentation/node34.html

function f = objfun ( x, x_new )
	f = prod(x([1 4])) * sum(x(1:3)) + x(3)
endfunction
function c = confun(x, x_new)
	c = [
	prod(x) - 25
	sum(x.^2)-40
	]
endfunction

////////////////////////////////////////////////////////////////////////
// Define gradient and Hessian matrix

function gf = gradf ( x, x_new )
	s13 = sum(x(1:3))
	p14 = prod(x([1 4]))
	gf(1) = p14 + x(4) * s13
	gf(2) = p14
	gf(3) = p14 + 1
	gf(4) = x(1) * s13
endfunction

function H = hessf ( x, x_new )
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
function y = dh(x,lambda,obj_weight)
	//pause
	y = obj_weight * hessf ( x ) + lambda(1) * Hg1(x) + lambda(2) * Hg2(x)
endfunction

// The constraints
function y=dg(x, x_new)
	//pause
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
x0 = [1,5,5,1]';
fx0 = 16;

////////////////////////////////////////////////////////////////////////
// The expected optimum, according to HS
xopt =  [1,4.7429994,3.8211503,1.3794082]'
fopt = 17.0140173
//objfun(xopt)
//confun(xopt)

////////////////////////////////////////////////////////////////////////
// Check gradient/Hessian

if ( %f ) then

	// Check the gradient/Hessian of f
	df_e = gradf ( x0 )
	Hf_e = hessf ( x0 )
	[J_c,H_c] = derivative ( objfun , x0 , order = 4 , H_form = "blockmat" )
	norm ( df_e - J_c' )
	norm ( Hf_e - H_c )

	function y = g1 ( x )
		c = confun(x)
		y = c(1)
	endfunction

	function y = g2 ( x )
		c = confun(x)
		y = c(2)
	endfunction

	// Check the gradient/Hessian of g1
	dg1_e = dg1 ( x0 )
	Hg1_e = Hg1 ( x0 )
	[J_c,H_c] = derivative ( g1 , x0 , order = 4 , H_form = "blockmat" )
	norm ( dg1_e - J_c' )
	norm ( Hg1_e - H_c )

	// Check the gradient/Hessian of g2
	dg2_e = dg2 ( x0 )
	Hg2_e = Hg2 ( x0 )
	[J_c,H_c] = derivative ( g2 , x0 , order = 4 , H_form = "blockmat" )
	norm ( dg2_e - J_c' )
	norm ( Hg2_e - H_c )

end

////////////////////////////////////////////////////////////////////////

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

// the variables have lower bounds of 1
lower = [1;1;1;1];
// the variables have upper bounds of 5
upper = [5;5;5;5];
var_lin_type(1) = 1; // Non-Linear
var_lin_type(2) = 1; // Non-Linear (this variable appears nonlinearly in the objective function or at least in one constraint)
var_lin_type(3) = 1; // Non-Linear (this variable appears nonlinearly in the objective function or at least in one constraint)
var_lin_type(4) = 1; // Non-Linear (this variable appears nonlinearly in the objective function or at least in one constraint)
constr_lin_type (1) = 1; // Non-Linear
constr_lin_type (2) = 1; // Non-Linear
// the first constraint g1 has a lower bound of 25
constr_lhs(1) = 25;
// the first constraint g1 has NO upper bound, here we set it to %inf.
constr_rhs(1) = %inf;
// the second constraint g2 is an equality constraint, so we set the
// upper and lower bound to the same value
constr_lhs(2) = 40;
constr_rhs(2) = 40;

params = init_param();
// We use the given Hessian
//params = add_param(params,"hessian_approximation","exact");
// We use a limited-bfgs approximation for the Hessian.
params = add_param(params,"hessian_approximation","limited-memory");

params = add_param(params,"print_level",12);
params = add_param(params,"journal_level",12);


//[x_sol, f_sol, extra] = ipopt(x0, objfun, gradf, confun, dg, sparse_dg, dh, sparse_dh, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);
[x_sol, f_sol, extra] = ipopt(x0, objfun, gradf, confun, dg, sparse_dg, [], [], var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);

x_sol
f_sol

////////////////////////////////////////////////////////////////////////
// Draw ipopt while progressing

function f = objfunD ( x )
	f = objfun ( x )
	//plot(x(1),x(2),"r+")
	xs = strcat(string(x)," ")
	mprintf("x=[%s]\n",xs)
endfunction
[x_sol, f_sol, extra] = ipopt(x0, objfunD, gradf, confun, dg, sparse_dg, dh, sparse_dh, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);


