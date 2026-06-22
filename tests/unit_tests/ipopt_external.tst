// =============================================================================
// Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
//
// Copyright (C) 2010 - Yann Collette
// Copyright (C) 2020 - UTC - Stéphane MOTTELET
//
// This file is hereby licensed under the terms of the GNU GPL v2.0,
// For more information, see the COPYING file which you should have received
//


// <-- CLI SHELL MODE -->
// <-- NO CHECK REF -->

f  = [];
df = [];
g  = [];
dg = [];
dh = [];
sparse_dg = [];
sparse_dh = [];

current_dir = pwd();

function out = f(x,x_new)
    out = 5*x(1)^2 - 3*x(2)^2;
end

function out = df(x,x_new)
    out = [10*x(1)
           -6*x(2)];
end

function out = g(x,x_new)
    out = -x;
end

function out = dg(x,x_new)
    out = [-1;-1];
end


f_C = ['#include <math.h>'
        'int f_C(double * x, double * f, int n_size_x, double x_new)'
        '{'
        '  f[0] = 5*pow(x[0],2) - 3*pow(x[1],2);'
        '  return 1;'
        '}'];

df_C = ['#include <math.h>'
        'int df_C(double * x, double * f, int n_size_x, double x_new)'
        '{'
        '  f[0] = 10*x[0];'
        '  f[1] = -6*x[1];'
        '  return 1;'
        '}'];

g_C = ['#include <math.h>'
       'int g_C(double * x, int n_size_x, double * g, int n_size_g, double x_new)'
       '{'
       '  g[0] = -x[0];'
       '  g[1] = -x[1];'
       '  return 1;'
       '}'];

dg_C = ['#include <math.h>'
        '#include <stdio.h>'
        'int dg_C(double * x, int n_size_x, double new_x, int n_size_g, double * values)'
        '{'
        '  values[0] = -1;'
        '  values[1] = -1;'
        '  return 1;'
        '}'];

cd TMPDIR;
mputl(f_C, TMPDIR+'/ipopt_demo_f_C.c');
mputl(df_C,TMPDIR+'/ipopt_demo_df_C.c');
mputl(g_C, TMPDIR+'/ipopt_demo_g_C.c');
mputl(dg_C,TMPDIR+'/ipopt_demo_dg_C.c');

// compile the C code
printf('Compilation of the f_C function\n');
f_C_link_handle  = ilib_for_link('f_C', 'ipopt_demo_f_C.c', [],'c');
printf('Compilation of the df_C function\n');
df_C_link_handle = ilib_for_link('df_C','ipopt_demo_df_C.c',[],'c');
printf('Compilation of the g_C function\n');
g_C_link_handle  = ilib_for_link('g_C', 'ipopt_demo_g_C.c', [],'c');
printf('Compilation of the dg_C function\n');
dg_C_link_handle = ilib_for_link('dg_C','ipopt_demo_dg_C.c',[],'c');

// incremental linking
link(f_C_link_handle, 'f_C', 'c');
link(df_C_link_handle,'df_C','c');
link(g_C_link_handle, 'g_C', 'c');
link(dg_C_link_handle,'dg_C','c');

cd(current_dir);

upper = [4;4];
lower = [-4;-4];
x0    = [1;1];
  
var_lin_type(1) = 1; // Non-Linear
var_lin_type(2) = 1; // Non-Linear
constr_lin_type (1) = 0; // Linear
constr_lin_type (2) = 0; // Linear
constr_rhs(1) = 0;
constr_rhs(2) = 0;
constr_lhs(1) = -%inf;
constr_lhs(2) = -%inf;

sparse_dg = [1 1
             2 2];
             
////////////////////////////////////////////////////////////////////////

params = struct();
params.hessian_approximation = "limited-memory";

[x_sol, f_sol, extra] = ipopt(x0, 'f_C', 'df_C', 'g_C', 'dg_C', sparse_dg, dh, sparse_dh, ...
                              var_lin_type, [], constr_rhs, constr_lhs, lower, upper, [], params);

// unlink the optimization problem

//ulink([f_C_link_handle, df_C_link_handle, g_C_link_handle, dg_C_link_handle]);

assert_checkalmostequal(x_sol,[0.0000204;4],0,1e-6);
assert_checkalmostequal(extra.lambda,[0.000204;0],0,1e-6);
