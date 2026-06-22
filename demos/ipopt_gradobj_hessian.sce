// =============================================================================
// Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
//
// Copyright (C) 2020-2021 - UTC - Stéphane MOTTELET
//
// This file is hereby licensed under the terms of the GNU GPL v3.0,
// For more information, see the COPYING file which you should have received

function ipopt_gradobj_hessian()
    
function [r,dr] = residual(x,_t,_y)
    _exp = exp(x(4)*_t.^2);
    _cosexp = cos(x(2)*_t+x(3)).*_exp;
    z=x(1)*_cosexp;
    r=z-_y;
    if argn(1) == 2 then
        _sinexp = -x(1)*sin(x(2)*_t+x(3)).*_exp;
        dr = [_cosexp,  _t.*_sinexp,  _sinexp,  x(1)*_t.^2.*_cosexp];
    end
endfunction


function [f,df] = costf(x, new_x, _t, _y)
    if argn(1) == 1 then
        r = residual(x,_t,_y);
        f = r.'*r; // use dot transpose for complex step
    else
        [r,dr] = residual(x,_t,_y);
        f = r.'*r;
        df = 2*dr.'*r; // use dot transpose for complex step
    end
endfunction

function [g,dg] = constr(x,new_x)
    g = x.'*x-6; // use dot transpose for complex step
    if argn(1) == 2 then
        dg = 2*x;
    end    
endfunction

function h = hessian(x,new_x,obj,lambda,new_lambda,_t,_y)
    n=length(x)
    nc=length(lambda)
    Hx=zeros(n,n);
    Hl=zeros(n,n);
    xi=x;
    //gradient with complex step
    dx = 1e-100;
    for k=1:n
        xi(k)=complex(x(k),dx);
        [f,df] = costf(xi, new_x, _t, _y)
        [g,dg] = constr(xi, new_x)
        Hx(:,k) = imag(df)/dx;
        Hl(:,k) = imag(dg)/dx;
        xi(k)=x(k);
    end
    h = obj*Hx+lambda*Hl;
    h = (h + h')/2;
endfunction

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
                ipopt_hdl(2).data(:,2) = z;
             end     
        end
        out = %t;
    end

t=linspace(-10,10,100)';
y=cos(2*t+1).*exp(-0.1*t^2);
y=y+rand(y,'normal')/10;

params=struct();
params.derivative_test = "second-order"
params.derivative_test_print_all = "yes"
params.hessian_approximation = "exact"

problem = struct();
problem.x0 = [1;1;1;1];
problem.f = list(costf,t,y);
problem.g = constr;
problem.dh = list(hessian,t,y);

problem.constr_rhs = 0;
problem.constr_lhs = -%inf;
problem.x_lower = [0 0 0 -%inf];
problem.x_upper = [%inf %inf %inf 0];
problem.int_cb =  list(intcb,t,y);
problem.params = params;
                
clf
       
[x,f,extra]=ipopt(problem);

end

ipopt_gradobj_hessian()
clear ipopt_gradobj_hessian
clearglobal ipopt_hdl







