//
// Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
// Copyright (C) 2020-2023 - UTC - Stéphane MOTTELET
//
// This file is hereby licensed under the terms of the GNU GPL v2.0,
// For more information, see the COPYING file which you should have received
//
//--------------------------------------------------------------------------

#include <IpIpoptApplication.hpp>
#include <IpSmartPtr.hpp>

#include "function.hxx"
#include "double.hxx"
#include "string.hxx"
#include "struct.hxx"
#include "configvariable.hxx"
#include "commentexp.hxx"

#include "parameters.hxx"
#include "argumentsPositions.hxx"
#include "OptimizationManager.hxx"
#include "IpoptTNLP.hxx"
#include "manage_params.hpp"
#include "scilabjournal.hpp"

extern "C"
{
#include "localization.h"
#include "Scierror.h"
#include "sciprint.h"
}

/*------------------------------------ SCILAB GATEWAY -----------------------------------------------------------------------*/

types::Function::ReturnValue sci_ipopt(types::typed_list &in, int _iRetCount, types::typed_list &out)
{
    // [x_sol, f_sol, extra] = ipopt(x0, f, df, g, dg, sparse_dg, dh, sparse_dh, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);
    types::Function::ReturnValue ret = types::Function::Error;
    OptimizationManager* manager = NULL;
    bool bRet;
    int iPrintLevel = 0;
    ApplicationReturnStatus ipopt_status;

    std::map<Ipopt::SolverReturn,std::wstring> IpoptStatus {
    {SUCCESS, L"Optimal solution found"},                  
    {MAXITER_EXCEEDED , L"Maximum number of iterations exceeded"},       
    {CPUTIME_EXCEEDED , L"Maximum CPU time exceeded"},       
    {STOP_AT_TINY_STEP , L"Search direction is becoming Too Small"},
    {STOP_AT_ACCEPTABLE_POINT , L"Solved To acceptable level"},                
    {LOCAL_INFEASIBILITY , L"Converged to a point of local infeasibility. Problem may be infeasible"},               
    {USER_REQUESTED_STOP , L" Stopping optimization at current point as requested by user"},              
    {FEASIBLE_POINT_FOUND , L"Feasible point for square problem found"},       
    {DIVERGING_ITERATES , L"Iterates diverging; problem might be unbounded"},                
    {RESTORATION_FAILURE , L" Restoration Failed"},         
    {ERROR_IN_STEP_COMPUTATION , L"Error in step computation (regularization becomes too large?)"},     
    {INVALID_NUMBER_DETECTED , L"Invalid number in NLP function or derivative detected"},        
    {TOO_FEW_DEGREES_OF_FREEDOM , L"Problem has too few degrees of freedom"},                    
    {INVALID_OPTION , L"Invalid option encountered"}, 
    {OUT_OF_MEMORY , L" Not enough memory"},
    {INTERNAL_ERROR , L"Internal error"}};
    
    std::vector<std::wstring> wstrExtraFields = {
    L"status",      
    L"it_count",        
    L"cpu_time",        
    L"fobj_eval",       
    L"fobj_grad_eval",  
    L"constr_eval",     
    L"constr_jac_eval", 
    L"hess_eval",       
    L"dual_inf",        
    L"constr_viol",     
    L"complementarity", 
    L"kkt_error"}; 
    
    wchar_t pwstrArgName[16][16] = {L"x0",L"f",L"df",L"g",L"dg",L"sparse_dg",L"dh",L"sparse_dh",
       L"var_lin_type",L"constr_lin_type",L"constr_rhs",L"constr_lhs",L"x_lower",
        L"x_upper",L"int_cb"};
    
    
    if (_iRetCount > 3)
    {
      Scierror(999, _("%s: Wrong number of output argument(s): at most %d expected.\n"), "ipopt", 3);
      return types::Function::Error;
    }

    if (in.size() == 1 && in[0]->isStruct())
    {
        types::Struct *pStruct = in[0]->getAs<types::Struct>()->clone();

        if (pStruct->getSize() != 1)
        {
            Scierror(999, _("%s: Wrong size for input argument #%d : A single struct expected.\n"), "ipopt",1);
            return types::Function::Error;            
        }

        types::typed_list in2;
        for (int i=0; i<PARAMS_POS; i++)
        {
            if (pStruct->exists(pwstrArgName[i]))
            {
                in2.push_back(pStruct->extractField(pwstrArgName[i]));
                pStruct->removeField(pwstrArgName[i]);
            }
            else
            {
                in2.push_back(new types::Double(0,0));
            }
        }

        types::String *pStrRemaining=pStruct->getFieldNames();
        if (pStrRemaining->getSize() > 0)
        {
            in2.push_back(pStruct);
        }
        pStrRemaining->killMe();

        types::Function::ReturnValue ret = sci_ipopt(in2, _iRetCount, out);

        pStruct->killMe();
        return ret;
    }

    if (in.size() != 15 && in.size() != 16)
    {
        Scierror(999, _("%s: Wrong number of input argument(s): %d to %d expected.\n"), "ipopt", 15,16);
        return types::Function::Error;
    }

    if ( in.size()==PARAMS_POS+1 && (in[PARAMS_POS]->isStruct() == false 
        && (in[PARAMS_POS]->isDouble()==false || in[PARAMS_POS]->getAs<types::Double>()->isEmpty() == false )) )
    {
        Scierror(999, _("%s: Wrong type for input argument #%d: A structure expected.\n"), "ipopt", 16);
        return types::Function::Error;    
    }

    // Create an instance of the OptimizationManager
    manager = new OptimizationManager(L"ipopt");

    // Create an instance of the IpoptApplication
    SmartPtr<IpoptApplication> ipopt_app = new IpoptApplication(false); // true because we want messages to be printed in the bash console
    manager->setIpoptApp(ipopt_app);

    // Manage options (last input argument in ipopt() call)

    if (in.size() == PARAMS_POS+1)
    {
        try
        {
            types::Struct *pStruct = in[PARAMS_POS]->getAs<types::Struct>()->clone();
            manage_ipopt_params(pStruct, ipopt_app->Options(), true);
            types::String *pStrRemaining=pStruct->getFieldNames();
            pStruct->killMe();
            if (pStrRemaining->getSize() > 0)
            {
                Scierror(999, _("%s: Unknown field name \"%ls\" in problem or parameter structure.\n"), "ipopt", pStrRemaining->get(0));
                pStrRemaining->killMe();
                return types::Function::Error;
            }
            pStrRemaining->killMe();
        }
        catch (ast::InternalError &ie)
        {
            throw(ie);
        }
    }

    // Suppress printing IPOpt startup message
    ipopt_app->Options()->SetStringValue("sb","yes");
    // Default should be "limited memory"
    ipopt_app->Options()->SetStringValueIfUnset("hessian_approximation","limited-memory");

    std::string strJacOpt, strHessOpt, strDerTest;
    ipopt_app->Options()->GetStringValue("jacobian_approximation", strJacOpt, "");
    ipopt_app->Options()->GetStringValue("hessian_approximation", strHessOpt, "");
    ipopt_app->Options()->GetStringValue("derivative_test", strDerTest, "");

    if (strJacOpt == "finite-difference-values" && (strDerTest== "first-order" || strDerTest== "second-order"))
    {
        delete manager;
        Scierror(999, _("%s: jacobian_approximation must be ""exact"" to check first-order derivatives.\n"), "ipopt", 16);
        return types::Function::Error;    
        return types::Function::Error;            
    }

    if (strHessOpt == "limited-memory" && strDerTest== "second-order")
    {
        delete manager;
        Scierror(999, _("%s: hessian_approximation must be ""exact"" to check second-order derivatives.\n"), "ipopt", 16);
        return types::Function::Error;    
        return types::Function::Error;            
        
    }
    
    bRet = manager->parseMatrices(in); IPOPT_ERROR;
    
    // parse the functions arguments (string, list, or macro)
    // please don't change the order of first 4 calls below:
    bRet = manager->parseFunction(in, FOBJ_POS,  OptimizationManager::FOBJ); IPOPT_ERROR;
    bRet = manager->parseFunction(in, DFOBJ_POS, OptimizationManager::DFOBJ); IPOPT_ERROR;
    bRet = manager->parseFunction(in, GOBJ_POS,  OptimizationManager::GOBJ); IPOPT_ERROR;
    bRet = manager->parseFunction(in, DGOBJ_POS, OptimizationManager::DGOBJ); IPOPT_ERROR;
    bRet = manager->parseFunction(in, DHOBJ_POS, OptimizationManager::DHOBJ); IPOPT_ERROR;
    bRet = manager->parseFunction(in, INTCB_POS, OptimizationManager::INTCB); IPOPT_ERROR;
            
    // allocate some members (i.e. vectors of double) of the OptimizationManager object
    manager->setAllStuff();

    // Add a new journal to print all the ipopt informations in the scilab console
    ipopt_app->Options()->GetIntegerValue("print_level", iPrintLevel, "");
    
    SmartPtr<Journal> console = new ScilabJournal((EJournalLevel)iPrintLevel);
    ipopt_app->Jnlst()->AddJournal(console);
    ipopt_app->RethrowNonIpoptException(true);    
    
    // Initialize

    try
    {
        ipopt_status = ipopt_app->Initialize();
    }
    catch(IpoptException &E)
    {
        //There has been a failure to solve a problem with Ipopt.
        sciprint("ipopt: Exception - %s\n",E.Message().c_str());
    }
    catch(...)
    {
        sciprint("ipopt: An unmanaged exception has been thrown during initialization\n");
    }

    if (ipopt_status != Solve_Succeeded)
    {
        Scierror(999,"%s: Error during initialization!\n", "ipopt");
        delete manager;
        return types::Function::Error;    
    }
    
    // create the nonlinear program
    SmartPtr<IpoptTNLP> mynlp = new IpoptTNLP();
    mynlp->setOptimizationManager(manager);

    try
    {
        // launch the optimization
        ipopt_status = ipopt_app->OptimizeTNLP(mynlp);
    }
    catch(IpoptException &E)
    {
        //There has been a failure to solve a problem with Ipopt.
        sciprint("ipopt: Exception - %s\n",E.Message().c_str());
    }
    catch(std::bad_alloc &E)
    {
        sciprint("ipopt: allocation problem - %s\n", E.what());
    }
    catch(const ast::InternalError& ie)
    {
        delete manager;
        char *message = wide_string_to_UTF8(ie.GetErrorMessage().c_str());        
        Scierror(999, message);
        FREE(message);
        return types::Function::Error;            
    }

    // create Scilab variables for output arguments
    types::Double* pDblX = new types::Double(manager->getNumberOfVariables(), 1);
    types::Double* pDblF = new types::Double(1, 1);
    // these ones as fields of "extra" structure
    types::Double* pDblLambda = new types::Double(manager->getNumberOfConstraints(), 1);
    types::Double* pDblZL = new types::Double(manager->getNumberOfVariables(), 1);
    types::Double* pDblZU = new types::Double(manager->getNumberOfVariables(), 1);
    
    // get primal and dual variables
    manager->getSolution(pDblX->get(), pDblLambda->get(), pDblZL->get(), pDblZU->get(), pDblF->get());

    types::Struct *pStructExtra = new types::Struct(1,1);
    // add "message" field
    pStructExtra->addField(L"message");
    // as a human readable string
    types::String *pStr = new types::String(IpoptStatus[manager->getSolverReturn()].c_str());
    pStructExtra->get(0)->set(L"message", pStr);           

    // add Lagrange multipliers
    pStructExtra->addField(L"lambda");
    pStructExtra->get(0)->set(L"lambda", pDblLambda);           
    pStructExtra->addField(L"zL");
    pStructExtra->get(0)->set(L"zL", pDblZL);           
    pStructExtra->addField(L"zU");
    pStructExtra->get(0)->set(L"zU", pDblZU);           

    // scalar fields
    for (int i=OptimizationManager::status; i <= OptimizationManager::kkt_error; i++)
    {
        types::Double *pDbl = new types::Double(manager->getdblExtra(i));
        pStructExtra->addField(wstrExtraFields[i]);
        pStructExtra->get(0)->set(wstrExtraFields[i], pDbl);           
    }
    
    out.push_back(pDblX);

    if (_iRetCount >= 2)
    {
      out.push_back(pDblF);
    }
    else
    {
      delete pDblF;
    }
    if (_iRetCount == 3)
    {
      out.push_back(pStructExtra);
    }
    else
    {
      delete pStructExtra;
    }

    delete manager;

    return types::Function::OK;    
}
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
