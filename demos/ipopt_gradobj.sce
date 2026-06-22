// =============================================================================
// Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
//
// Copyright (C) 2020-2024 - UTC - Stéphane MOTTELET
//
// This file is hereby licensed under the terms of the GNU GPL v3.0,
// For more information, see the COPYING file which you should have received

function ipopt_gradobj()

    function f=residual(x,_t,_y)
        z=x(1)*cos(x(2)*_t+x(3)).*exp(x(4)*_t.^2);
        r=z-_y;
        f=sum(r.*r);    
    end

    function [f,df] = costf(x, new_x, _t, _y)
        f = residual(x);
        if argn(1) == 2 then
            n=length(x)    
            df=zeros(n,1);
            xi=x;
            //gradient with complex step
            dx = 1e-100;
            for k=1:n
                xi(k)=complex(x(k),dx);
                df(k) = imag(residual(xi,_t,_y))/dx;
                xi(k)=x(k);
            end
        end
    end

    function [g,dg] = constr(x,new_x)
        g = norm(x)^2-6;
        if argn(1) == 2 then
            dg =  2*x;
        end    
    end

    function out = intcb(param,_t,_y)
        global ipopt_hdl
        if isfield(param,"x")
            x = param.x;
            z=x(1)*cos(x(2)*_t+x(3)).*exp(x(4)*_t.^2);
            if isempty(ipopt_hdl) || ~is_handle_valid(ipopt_hdl)
                drawlater
                ipopt_hdl = plot(_t,z,_t,_y,'o')
                gca().data_bounds(3:4)=[-1.5,1.5];
                legend model data
                drawnow
             else
                ipopt_hdl(1).data(:,2) = z;
             end     
        end
        out = %t;
    end

    t=linspace(-10,10,100);
    y=cos(2*t+1).*exp(-0.1*t.^2);
    y=y+rand(y,'normal')/10;

    problem = struct();
    problem.x0 = rand(4,1);
    problem.f = list(costf,t,y);
    problem.g = constr;
    problem.constr_rhs = 0;
    problem.constr_lhs = -%inf;
    problem.x_lower = [0 0 0 -%inf];
    problem.x_upper = [%inf %inf %inf 0];
    problem.int_cb =  list(intcb,t,y);

    clf

    x=ipopt(problem);

    demo_viewCode("ipopt_gradobj.sce")

end

ipopt_gradobj()
clear ipopt_gradobj
clearglobal ipopt_hdl






