// =============================================================================
// Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
//
// Copyright (C) 2020-2023 - UTC - Stéphane MOTTELET
//
// This file is hereby licensed under the terms of the GNU GPL v3.0,
// For more information, see the COPYING file which you should have received

function ipopt_quad()

    function [f,df] = costf(x,new_x,a,b)
        ax = a*x;
        f = .5*x'*ax-b'*x;
        // gradient is computed at almost zero cost
        df = ax-b;
    endfunction

    function h = hessian(x,new_x,obj,lambda,new_lambda,a)
        // hessian function is called only once
        [ij,values] = spget(a);
        h = obj*values;
    endfunction

    function out = callback(param,n)
        global ipopt_hdl
        out = %t;
        if isfield(param,"x")
            x = matrix(param.x,-1,2);
            if isempty(ipopt_hdl) || ~is_handle_valid(ipopt_hdl)
                ipopt_hdl = plot(x(:,1),x(:,2),'o');
                ipopt_hdl.mark_size = 1;
                gca().data_bounds=[-1.1,1.1,-1.1,1.1];
             else
                ipopt_hdl.data=x;
             end
             title([msprintf("Large sparse quadratic program with %2.0e unknowns",n)
             "Plot of x as pairs of coordinates"
             msprintf("iter = %d, cost = %f",param.iter, param.obj_value)])
        end
    endfunction

    // a is sparse, tridiagonal and positive definite
    n = 46340; // limitation is due to https://codereview.scilab.org/#/c/21675/
    ij = [1:n 1:n-1
          1:n 2:n]';
    v = grand(2*n-1,1,"nor",0,1);
    a = sparse(ij,v);
    a = a'*a+speye(n,n);
    b = grand(n,1,"nor",0,1);

    // ipopt parameters


    // optimization problem
    problem = struct();
    problem.x0 = zeros(n,1);
    problem.x_lower = -ones(n,1);
    problem.x_upper = ones(n,1);
    problem.f = list(costf,a,b);
    problem.sparse_dh = spget(a);
    problem.dh = list(hessian,a);
    problem.int_cb = list(callback,n);

    problem.print_level = 3;
    problem.hessian_approximation = "exact"
    problem.hessian_constant = "yes"
    clf

    x = ipopt(problem);

    demo_viewCode("ipopt_quad.sce");
end

ipopt_quad()
clear ipopt_quad
clearglobal ipopt_hdl
