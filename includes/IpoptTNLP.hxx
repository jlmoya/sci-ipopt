//
// Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
// Copyright (C) 2020 - UTC - Stéphane MOTTELET
//
// This file is hereby licensed under the terms of the GNU GPL v2.0,
// For more information, see the COPYING file which you should have received
//
//--------------------------------------------------------------------------

#include <IpIpoptApplication.hpp>
#include <IpSolveStatistics.hpp>
#include <IpTNLP.hpp>
#include <IpIpoptCalculatedQuantities.hpp>
#include <IpSmartPtr.hpp>
#include "IpIpoptData.hpp"
#include "IpTNLPAdapter.hpp"
#include "IpOrigIpoptNLP.hpp"

#include "OptimizationManager.hxx"

using namespace Ipopt;

class IpoptTNLP : public TNLP
{
public:

    IpoptTNLP() : m_OptimizationManager(NULL) {}
    virtual ~IpoptTNLP() {}

    virtual bool get_variables_linearity(Index n, TNLP::LinearityType* var_types);
    virtual bool get_constraints_linearity(Index m, TNLP::LinearityType* const_types);
    virtual bool get_nlp_info(Index& n, Index&m, Index& nnz_jac_g,Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style);
    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,Index m, Number* g_l, Number* g_u);
    virtual bool get_starting_point(Index n, bool init_x, Number* x,bool init_z, Number* z_L, Number* z_U,
    Index m, bool init_lambda,Number* lambda);
    virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);
    virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);
    virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);
    virtual bool eval_jac_g(Index n, const Number* x, bool new_x,Index m, Index nele_jac, Index* iRow, Index *jCol,Number* values);
    virtual bool eval_h(Index n, const Number* x, bool new_x,Number obj_factor, Index m, const Number* lambda,bool new_lambda,
    Index nele_hess, Index* iRow, Index* jCol, Number* values);
    virtual Index get_number_of_nonlinear_variables();
    virtual bool  get_list_of_nonlinear_variables(Index num_nonlin_vars, Index* pos_nonlin_vars);
    virtual void finalize_solution(SolverReturn status,Index n, const Number *x, const Number *z_L, const Number *z_U, Index m,
    const Number *g, const Number *lambda, Number obj_value, const IpoptData *ip_data,
    IpoptCalculatedQuantities *ip_cq);
    virtual bool intermediate_callback(AlgorithmMode mode, Index iter, Number obj_value, Number inf_pr, Number inf_du, Number mu, 
    Number d_norm, Number regularization_size, Number alpha_du, Number alpha_pr, Index ls_trials, 
    const IpoptData *ip_data, IpoptCalculatedQuantities *ip_cq);

    virtual OptimizationManager* getOptimizationManager();
    virtual void setOptimizationManager(OptimizationManager *);

    typedef int (*call_f_objective)(double *, double *, int, double);
    typedef int (*call_f_objective_grad)(double *, double *, int, double);
    typedef int (*call_f_constraints)(double *, int, double *, int, double);
    typedef int (*call_f_constraints_jac)(double *, int, double, int, double *);
    typedef int (*call_f_hessian)(double *, int, double, double, int, double *, double, double *);

private:

    OptimizationManager* m_OptimizationManager;
};