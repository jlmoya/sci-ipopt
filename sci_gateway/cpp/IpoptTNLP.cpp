//
// Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
// Copyright (C) 2020 - UTC - Stéphane MOTTELET
//
// This file is hereby licensed under the terms of the GNU GPL v2.0,
// For more information, see the COPYING file which you should have received
//
//--------------------------------------------------------------------------

#include "IpoptTNLP.hxx"

// IpoptTNLP methods

void IpoptTNLP::setOptimizationManager(OptimizationManager *manager)
{
    m_OptimizationManager = manager;
}

OptimizationManager* IpoptTNLP::getOptimizationManager()
{
    return m_OptimizationManager;
}

bool IpoptTNLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
    OptimizationManager* manager = getOptimizationManager();
    manager->setIsComputed(new_x);    
    manager->increaseFunctionCount(OptimizationManager::FOBJ);

    dynlibFunPtr pFunc = manager->getEntryPointFunction(OptimizationManager::FOBJ);
    if (pFunc != NULL)
    {
        // Computation is made by user C/C++/Fortran dynamically lindked entrypoint
        return (bool) ((call_f_objective)pFunc)((double *) x, (double *) &obj_value, (int) n, (double) new_x);
    }    
    // Computation is made by Scilab user function
    if (!manager->computeFunctionsAndDerivatives((double *)x, new_x, OptimizationManager::FOBJ))
    {
        return false;
    }
    return manager->getObjective((double *)&obj_value);
}

bool IpoptTNLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
    OptimizationManager* manager = getOptimizationManager();
    manager->setIsComputed(new_x);
    manager->increaseFunctionCount(OptimizationManager::DFOBJ);

    dynlibFunPtr pFunc = manager->getEntryPointFunction(OptimizationManager::DFOBJ);
    if (pFunc != NULL)
    {
        // Computation is made by user C/C++/Fortran dynamically lindked entrypoint
        return (bool) ((call_f_objective_grad)pFunc)((double *) x, (double *) grad_f, (int) n, (double) new_x);
    }
    // Computation is made by Scilab user function
    if (!manager->computeFunctionsAndDerivatives((double *)x, new_x, OptimizationManager::DFOBJ))
    {
        return false;
    }
    return manager->getComputed((double *)grad_f, OptimizationManager::DFOBJ);
}

bool IpoptTNLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
    OptimizationManager* manager = getOptimizationManager();
    manager->setIsComputed(new_x);
    manager->increaseFunctionCount(OptimizationManager::GOBJ);

    dynlibFunPtr pFunc = manager->getEntryPointFunction(OptimizationManager::GOBJ);
    if (pFunc != NULL)
    {
        // Computation is made by user C/C++/Fortran dynamically lindked entrypoint
        return (bool) ((call_f_constraints)pFunc)((double *) x, (int)n, (double *) g, (int) m, (double) new_x);
    }
    // Computation is made by Scilab user function
    if (!manager->computeFunctionsAndDerivatives((double *)x, new_x, OptimizationManager::GOBJ))
    {
        return false;
    }
    return manager->getComputed((double *)g, OptimizationManager::GOBJ);
}

bool IpoptTNLP::eval_jac_g(Index n, const Number* x, bool new_x,Index m, Index nele_jac, Index* iRow, Index *jCol,Number* values)
{
    OptimizationManager* manager = getOptimizationManager();
    
    if (x == NULL) // get jacobian pattern
    {
        return manager->getConstraintsJacobianPattern((int *)iRow, (int *)jCol);
    }
    manager->setIsComputed(new_x);
    manager->increaseFunctionCount(OptimizationManager::DGOBJ);
    
    dynlibFunPtr pFunc = manager->getEntryPointFunction(OptimizationManager::DGOBJ);
    if (pFunc != NULL)
    {
        // Computation is made by user C/C++/Fortran dynamically lindked entrypoint
        return (bool) ((call_f_constraints_jac)pFunc)((double *) x, (int) n, (double) new_x, (int) m, (double *) values );
    }
    // Computation is made by Scilab user function
    if (!manager->computeFunctionsAndDerivatives((double *)x, new_x, OptimizationManager::DGOBJ))
    {
        return false;
    }
    return manager->getComputed((double *)values, OptimizationManager::DGOBJ);
}

bool IpoptTNLP::eval_h(Index n, const Number* x, bool new_x,Number obj_factor, Index m, const Number* lambda,bool new_lambda,
Index nele_hess, Index* iRow, Index* jCol, Number* values)
{
    OptimizationManager* manager = getOptimizationManager();

    if (x == NULL) // get Hessian *lower triangular* part pattern
    {
        return manager->getHessianLowerTriPattern((int *)iRow, (int *)jCol);
    }
    manager->setIsComputed(new_x);
    manager->increaseFunctionCount(OptimizationManager::DHOBJ);

    dynlibFunPtr pFunc = manager->getEntryPointFunction(OptimizationManager::DHOBJ);
    if (pFunc != NULL)
    {
        // Computation is made by user C/C++/Fortran dynamically lindked entrypoint
        return (bool) ((call_f_hessian)pFunc)( (double *)x, (int)n, (double)new_x, (double)obj_factor, (int)m, (double *)lambda, (double)new_lambda, (double *)values);
    }
    // Computation is made by Scilab user function
    if (!manager->computeHessian((double *)x, new_x, (double)obj_factor, (double *)lambda, new_lambda, OptimizationManager::DHOBJ))
    {
        return false;
    }
    return manager->getComputed((double *)values, OptimizationManager::DHOBJ);
}

bool IpoptTNLP::get_variables_linearity(Index n, TNLP::LinearityType* var_types)
{
    OptimizationManager* manager = getOptimizationManager();
    bool ret = manager->getVariablesLinearity(var_types);
    return ret;
}

bool IpoptTNLP::get_constraints_linearity(Index m, TNLP::LinearityType* const_lin_types)
{
    OptimizationManager* manager = getOptimizationManager();
    bool ret = manager->getConstraintsLinearity(const_lin_types);
    return ret;
}

Index IpoptTNLP::get_number_of_nonlinear_variables()
{
    OptimizationManager* manager = getOptimizationManager();
    return manager->getNumberOfNonLinearVariables();
}

bool IpoptTNLP::get_list_of_nonlinear_variables(Index num_nonlin_vars, Index* pos_nonlin_vars)
{
    OptimizationManager* manager = getOptimizationManager();
    return manager->getNonLinearVariables(pos_nonlin_vars);
}


bool IpoptTNLP::get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u)
{
    OptimizationManager* manager = getOptimizationManager();
    bool ret = manager->getBoundsInfo((double *)x_l, (double *)x_u, (double *)g_l, (double *)g_u);
    return ret;
}

bool IpoptTNLP::get_nlp_info(Index& n, Index&m, Index& nnz_jac_g, Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style)
{
    OptimizationManager* manager = getOptimizationManager();
    n =          manager->getNumberOfVariables();
    m =          manager->getNumberOfConstraints();
    nnz_jac_g =  manager->getNumberOfConstraintsJacobianNonZeroTerms();
    nnz_h_lag =  manager->getNumberOfHessianNonZeroTerms();
    index_style = TNLP::C_STYLE;
    return true;
}

bool IpoptTNLP::get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda)
{
    OptimizationManager* manager = getOptimizationManager();
    return  manager->getX0((double *)x);
}

void IpoptTNLP::finalize_solution(SolverReturn status,
Index n, const Number* x, const Number* z_L,
const Number* z_U, Index m, const Number* g,
const Number* lambda, Number obj_value,
const IpoptData* ip_data, IpoptCalculatedQuantities* ip_cq)                                  
{
    OptimizationManager* manager = getOptimizationManager();
    SmartPtr<IpoptApplication> ipopt_app = manager->getIpoptApp();
    
    int int_fobj_eval, 
        int_constr_eval, 
        int_fobj_grad_eval,
        int_constr_jac_eval, 
        int_hess_eval;
        
    double dual_inf, constr_viol, complementarity, kkt_error;

    manager->setSolution((double *)x, (double *)lambda, (double *) z_L, (double *) z_U, (double)obj_value);
    manager->setSolverReturn(status);
    manager->setdblExtra(OptimizationManager::status, (double) status);

    if (Ipopt::IsValid(ipopt_app->Statistics()))
    {
        manager->setdblExtra(OptimizationManager::it_count, (double)ipopt_app->Statistics()->IterationCount());
        manager->setdblExtra(OptimizationManager::cpu_time, (double)ipopt_app->Statistics()->TotalCPUTime());
        ipopt_app->Statistics()->NumberOfEvaluations(int_fobj_eval,
        int_constr_eval,
        int_fobj_grad_eval,
        int_constr_jac_eval,
        int_hess_eval);

        manager->setdblExtra(OptimizationManager::fobj_eval, (double)int_fobj_eval);
        manager->setdblExtra(OptimizationManager::constr_eval, (double)int_constr_eval);
        manager->setdblExtra(OptimizationManager::fobj_grad_eval, (double)int_fobj_grad_eval);
        manager->setdblExtra(OptimizationManager::constr_jac_eval, (double)int_constr_jac_eval);
        manager->setdblExtra(OptimizationManager::hess_eval, (double)int_hess_eval);

        ipopt_app->Statistics()->ScaledInfeasibilities(dual_inf, constr_viol, complementarity, kkt_error);

        manager->setdblExtra(OptimizationManager::dual_inf, dual_inf);
        manager->setdblExtra(OptimizationManager::constr_viol, constr_viol);
        manager->setdblExtra(OptimizationManager::complementarity, complementarity);
        manager->setdblExtra(OptimizationManager::kkt_error, kkt_error);
    }
    else
    {
        //failsafe mode, in case of optimization stopped because systematic Nan or Inf returned by user function (or other problem)
        manager->setdblExtra(OptimizationManager::it_count, (double)ip_data->iter_count());
        manager->setdblExtra(OptimizationManager::kkt_error, ip_cq->curr_nlp_error());
        manager->setdblExtra(OptimizationManager::constr_viol, ip_cq->curr_constraint_violation());          
    }
}

bool IpoptTNLP::intermediate_callback(AlgorithmMode mode, Index iter, Number obj_value, Number inf_pr, Number inf_du, Number mu, 
    Number d_norm, Number regularization_size, Number alpha_du, Number alpha_pr, Index ls_trials, 
    const IpoptData *ip_data, IpoptCalculatedQuantities *ip_cq)
{
    OptimizationManager* manager = getOptimizationManager();

    std::vector<std::pair<std::wstring, double>> scalarQuantities = { 
    {L"algorithm_mode", (double)mode},
    {L"iter", (double)iter},
    {L"obj_count", (double)manager->getFunctionCount(OptimizationManager::FOBJ)},
    {L"obj_value", obj_value},
    {L"prev_obj_value", (double)manager->getPrevFobjValue()},
    {L"inf_pr", inf_pr},
    {L"inf_du", inf_du},
    {L"mu", mu},
    {L"d_norm",d_norm},
    {L"regularization_size",regularization_size},
    {L"alpha_du", alpha_du},
    {L"alpha_pr", alpha_pr},
    {L"ls_trials", (double)ls_trials},
    {L"curr_nlp_error", ip_cq->curr_nlp_error()}};

    manager->setPrevFobjValue(obj_value);
    bool bRet = manager->intermediateCallback(scalarQuantities, ip_data, ip_cq);
    
    return bRet;
}