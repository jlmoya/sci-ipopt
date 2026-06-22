// ====================================================================
// Copyright - UTC - Stéphane MOTTELET 2020
//
// This file is released into the public domain
// ====================================================================
demopath = get_absolute_file_path("sci_ipopt.dem.gateway.sce");

subdemolist = ["Fitting data", "ipopt_gradobj.sce"
               "Constrained Rosenbrock", "rosenbrock-ipopt.sce"
               "Non-convex", "ipopt_simonescu.sce"
               "Large scale quadratic program", "ipopt_quad.sce"];

subdemolist(:,2) = demopath + subdemolist(:,2);
// ====================================================================
