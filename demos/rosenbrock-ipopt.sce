// Copyright (C) 2010 - DIGITEO - Michael Baudin
// Copyright (C) 2020-2021 - UTC - Stéphane MOTTELET
//
// This file is hereby licensed under the terms of the GNU GPL v3.0,
// For more information, see the COPYING file which you should have received

// Test IPOPT In Scilab on Rosenbrock "banana" function

function rosenbrock_ipopt()

    function f = rosenbrock_f ( x , x_new )
        f = 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
    endfunction

    function df = rosenbrock_df ( x , x_new )
        df = [ -400*(x(2)-x(1)**2)*x(1)-2*(1-x(1))
                200*(x(2)-x(1)**2)];
    endfunction

    // Nonlinear inequality constraints
    function g = rosenbrock_g ( x , x_new )
      g = x(1)^2 + x(2)^2 - 1;
    endfunction

    function dg = rosenbrock_dg ( x , x_new )
      dg = [2*x(1)
            2*x(2)];
    endfunction

    // Hessian of f
    function Hf = rosenbrock_Hf ( x , x_new )
      Hf = zeros(2,2);
      Hf(1,1) = diag(-400*x(2) + 1200*x(1).^2 + 2);
      Hf(2,2) = 200;
      Hf = Hf - diag(400*x(1),1) - diag(400*x(1),-1);
    endfunction

    // Hessian of g
    function Hg = rosenbrock_Hg ( x , x_new )
      Hg = [2 0
            0 2];
    endfunction

    function f = rosenbrock_fC ( x1 , x2 )
        x = [x1 x2];
        f = rosenbrock_f(x,%t);
    endfunction

    function g = rosenbrock_gC ( x1 , x2 )
      x = [x1 x2];
      g = rosenbrock_g(x,%t);
    endfunction

    // The Hessian of the Lagrangian
    function y = rosenbrock_hessian ( x , lambda , obj_weight , x_new , lambda_new )
      Hf = rosenbrock_Hf ( x );
      Hg = rosenbrock_Hg ( x );
      y = obj_weight * Hf + lambda(1) * Hg;
      y = y([1 2 4]);
    endfunction

    function out = callback(params)
        global ipopt_hdl
        if isfield(params,"x")
            if isempty(ipopt_hdl) || ~is_handle_valid(ipopt_hdl)
                ipopt_hdl = plot(params.x(1),params.x(2),"-or")
            else
                ipopt_hdl.data = [ipopt_hdl.data;params.x'];
            end
        end
        out = %t;
    endfunction

    ////////// Main program //////////////

    drawlater
    clf
    gcf().color_map = jet(10)
    gca().background=color("gray95")

    x = linspace ( -2 , 2 , 100 );
    y = linspace ( -1 , 3 , 100 );

    // Plot contours of G
    [xc,yc]=contour2di( x , y , feval(x,y,rosenbrock_gC) , [0 0] );
    xc = xc(2:1+yc(1));
    yc = yc(2:1+yc(1));
    xfpoly(xc,yc,color("white"))
    gce().clip_state = "on"

    // Plot contours of F
    contour ( x , y , rosenbrock_fC , [1 10 100 500 1000] )
    h = gce().children.children;
    delete(h(h.type=="Text"))
    isoview on
    // The sparsity structure of the constraints
    sparse_dg = [
      1 1
      1 2
    ];

    // The sparsity structure of the Lagrangian
    sparse_dh = [
      1 1
      2 1
      2 2
    ];

    upper = [%inf %inf]';
    lower = [-%inf -%inf]';
    // Not Feasible starting point
    x0  = [grand("unf",-1.5,1.5) grand("unf",-0.5,2.5)];
    nb_constr = 1;
    var_lin_type = [1 1]'; // Non-Linear
    constr_lin_type = 1; // Non-Linear
    constr_rhs = 0;
    constr_lhs = -%inf;

    ////////////////////////////////////////////////////////////////////////

    params = struct();
    // We use the given Hessian
    params.hessian_approximation = "exact";
    // We use a limited-bfgs approximation for the Hessian.
    //params = add_param(params,"hessian_approximation","limited-memory");

    xstring(x0(1),x0(2),"$\,\LARGE \mathbf{x_0}$")

    drawnow

    [x_sol, f_sol, extra] = ipopt(x0, rosenbrock_f, rosenbrock_df, ..
      rosenbrock_g, rosenbrock_dg, sparse_dg, rosenbrock_hessian, sparse_dh, var_lin_type, ..
      constr_lin_type, constr_rhs, constr_lhs, lower, upper, callback, params);

    xstring(x_sol(1),x_sol(2),"$\,\LARGE \mathbf{\hat x}$")

    //plot(1,1,"ob")
    //xstring(1,1,["Unconstrained";"solution"])

    // Expected solution :  0.7864151   0.6176982

    disp(x_sol)
    disp(extra)

    demo_viewCode("rosenbrock-ipopt.sce")
end

rosenbrock_ipopt()
clear rosenbrock_ipopt 
clearglobal ipopt_hdl
