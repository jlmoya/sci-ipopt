// =============================================================================
// Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
//
// Copyright (C) 2020-2021 - UTC - Stéphane MOTTELET
//
// This file is hereby licensed under the terms of the GNU GPL v3.0,
// For more information, see the COPYING file which you should have received

function ipopt_simonescu()

    function [f,df] = costf(x, new_x)
        f = 0.1*x(1)*x(2);
        if argn(1) == 2 then
            df = 0.1*[x(2);x(1)]
        end
    endfunction

    function out =  g(x,y)
        out = constr([x;y], %t)
    endfunction

    function out =  f(x,y)
        out = costf([x;y], %t)
    endfunction

    function g = constr(x,new_x)
        rt = 1;
        rs = 0.2;
        n = 8;
        g = x(1)*x(1)+x(2)*x(2) - (rt+rs*cos(n*atan(x(1),x(2))))^2;
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

    constrLhs = -%inf;
    constrRhs = 0;
    constrJacPattern = [ 
       1.   1.
       1.   2.];

    binf = [-1.5;-1.5]
    bsup = [ 1.5; 1.5]

    opt = struct();
    //opt.start_with_resto = "yes";
    opt.jacobian_approximation = "finite-difference-values";
    opt.hessian_approximation = "limited-memory";
    opt.constr_viol_tol=1e-8;
    opt.print_level = 5;
    opt.tol = 1e-13;

    drawlater
    clf
    gcf().color_map = jet(30);
    gca().background=color("gray95")

    x = linspace ( -1.5 , 1.5 , 100 );
    y = linspace ( -1.5 , 1.5 , 100 );

    // Plot contours of G

    [xc,yc]=contour2di( x , y , feval(x,y,g) , [0 0] );
    xc = xc(2:1+yc(1));
    yc = yc(2:1+yc(1));
    xfpoly(xc,yc,color("white"))
    gce().clip_state = "on"


    //// Plot contours of F
    //contour ( x , y , rosenbrock_fC , [1 10 100 500 1000] )
    contour ( x , y , f,[-0.2:0.02:.2 ])
    h = gce().children.children;
    delete(h(h.type=="Text"))

    h = gce().children.children;
    delete(h(h.type=="Text"))
    isoview on
    //1.4739762
    //   1.1443051
    x0 = grand(1,2,"unf",-1.45,1.45)
    xstring(x0(1),x0(2),"$\,\LARGE \mathbf{x_0}$")
    drawnow

    [x_sol,f_sol,extra] = ipopt(x0, costf, [], constr, [], [], [], [], [], [], ...
      constrRhs, constrLhs, binf, bsup, callback, opt);

    xstring(x_sol(1),x_sol(2),"$\,\LARGE \mathbf{\hat x}$")

    demo_viewCode("ipopt_simonescu.sce");
end

ipopt_simonescu()
clear ipopt_simonescu
clearglobal ipopt_hdl
