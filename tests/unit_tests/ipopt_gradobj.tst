// =============================================================================
// Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
//
// Copyright (C) 2020 - UTC - Stéphane MOTTELET
//
// This file is hereby licensed under the terms of the GNU GPL v2.0,
// For more information, see the COPYING file which you should have received
//


//<-- NO CHECK REF -->

function f=residu(x)
    z=x(1)*cos(x(2)*t+x(3)).*exp(x(4)*t.^2);
    r=z-yb;
    f=sum(r.*r);    
endfunction

function [f,df] = costf(x, new_x)
    f = residu(x);
    if argn(1) == 2 then
        n=length(x)    
        df=zeros(n,1);
        xi=x;
        //gradient with complex step
        dx = 1e-100;
        for k=1:n
            xi(k)=complex(x(k),dx);
            df(k) = imag(residu(xi))/dx;
            xi(k)=x(k);
        end
    end
endfunction

function [g,dg] = constr(x,new_x)
    g = norm(x)^2-6;
    if argn(1) == 2 then
        dg =  2*x;
    end    
endfunction

function out = intcb(param)
//    disp(param)
//    out = param.iter < 10;
    out = %t;
endfunction

rand("seed",19);

t=linspace(-10,10,50);
a=1;
b=2;
c=1;
d=-.1;
y=a*cos(b*t+c).*exp(d*t^2);
yb=y+rand(y,'normal')/10;

x0=rand(4,1);
x0(4)=-x0(4)

//load data.sod

constrLhs = -%inf;
constrRhs = 0;
constrLinearity = 0;
constrJacPattern = [ 
   1.   1.
   1.   2.
   1.   3.
   1.   4];

binf = zeros(x0);
bsup = %inf*ones(x0);

binf(4)=-%inf;
bsup(4)=0;

opt = struct();
opt.hessian_approximation = "limited-memory";
opt.print_level = 5;

[x,f,extra] = ipopt(x0, costf, [], constr, [], [], [], [], [], [], ...
  constrRhs, constrLhs, binf, bsup, intcb, opt);

clf
tt=linspace(-10,10,200);
z=x(1)*cos(x(2)*tt+x(3)).*exp(x(4)*tt.^2);
drawlater
plot(tt,z,t,yb,'o')
gca().data_bounds(3:4)=[-1.5,1.5];
drawnow

assert_checkalmostequal(x,[1.010524;1.9913363;1.0023487;-0.0933706],1e-6);
assert_checkalmostequal(extra.lambda,0.2394931,1e-6);
