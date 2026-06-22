//
// Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
// Copyright (C) 2020 - UTC - Stéphane MOTTELET
//
// This file is hereby licensed under the terms of the GNU GPL v2.0,
// For more information, see the COPYING file which you should have received
//
//--------------------------------------------------------------------------

#ifndef _OPTIMIZATIONMANAGER_HXX_
#define _OPTIMIZATIONMANAGER_HXX_

#include <array>
#include <IpIpoptApplication.hpp>
#include <IpSolveStatistics.hpp>
#include <IpTNLP.hpp>
#include <IpIpoptCalculatedQuantities.hpp>
#include <IpSmartPtr.hpp>
#include "IpIpoptData.hpp"
#include "IpTNLPAdapter.hpp"
#include "IpOrigIpoptNLP.hpp"

#include "function.hxx"
#include "double.hxx"
#include "sparse.hxx"
#include <Eigen/Sparse>
#include "string.hxx"
#include "struct.hxx"
#include "configvariable.hxx"
#include "commentexp.hxx"

extern "C"
{
#include "localization.h"
#include "Scierror.h"
#include "sciprint.h"
}

typedef void(*dynlibFunPtr)();

class OptimizationManager
{
    public :

    OptimizationManager(const std::wstring& callerName);
    virtual ~OptimizationManager();
    enum functionKind {FOBJ=0,DFOBJ,GOBJ,DGOBJ,DHOBJ,INTCB};

    char pstrArgName[16][16] = {"x0","f","df","g","dg","sparse_dg","dh","sparse_dh",
        "var_lin_type","constr_lin_type","constr_rhs","constr_lhs","x_lower",
         "x_upper","int_cb","params"};
        
    virtual bool computeFunctionsAndDerivatives(double *, bool, functionKind); // for all functions but Hessian
    virtual bool computeHessian(double *, bool, double, double *, bool, functionKind);
    virtual void checkHessianOutput(types::InternalType *, functionKind);
    virtual bool intermediateCallback(std::vector<std::pair<std::wstring, double>>, 
         const Ipopt::IpoptData* , Ipopt::IpoptCalculatedQuantities*);
    virtual bool parseFunction(types::typed_list &, int, functionKind);
    virtual bool parseMatrices(types::typed_list &);
    template <typename T>  bool parseHessianPattern(T*,int,bool);


    virtual bool getSparseValuesFromIncompletePattern(double *, types::Sparse *, types::Sparse *);
    virtual bool getComputed(double *, functionKind);
    virtual bool getConstraintsJacobianPattern(int *, int *);
    virtual bool getHessianLowerTriPattern(int *, int *);
    virtual bool getConstraintsLinearity(Ipopt::TNLP::LinearityType*);
    virtual bool getBoundsInfo(double *, double *, double *, double *);
    virtual bool getX0(double *);
    virtual bool getVariablesLinearity(Ipopt::TNLP::LinearityType*);
    virtual bool getNonLinearVariables(Ipopt::Index*);
    virtual bool setAllStuff();
    virtual bool setSolution(double *, double *, double *, double *, double);
    virtual bool getSolution(double *, double *, double *, double *, double *);
    
    enum extraFields {status=0,          
    it_count,        
    cpu_time,        
    fobj_eval,       
    fobj_grad_eval,  
    constr_eval,     
    constr_jac_eval, 
    hess_eval,       
    dual_inf,        
    constr_viol,     
    complementarity, 
    kkt_error};

    void setPrevFobjValue(double dblvalue)
    {
        m_dblPrevFobjValue = dblvalue;
    }
    double getPrevFobjValue()
    {
        return m_dblPrevFobjValue;
    }
    void increaseFunctionCount(functionKind what)
    {
        m_iNbFunctionCount[what]++;
    }
    int getFunctionCount(functionKind what)
    {
        return m_iNbFunctionCount[what];
    }
    void setdblExtra (int field, double value) 
    {
        m_dblExtra[field] = value;
    }
    double getdblExtra (int field) 
    {
        return m_dblExtra[field];
    }
    char *getCallerStr()
    {
        return m_pstrCaller;
    }
    int getNumberOfVariables()
    {
        return m_iNbVariables;
    }
    int getNumberOfNonLinearVariables()
    {
        return m_iNbNonLinearVariables;
    }
    int getNumberOfConstraints()
    {
        return m_iNbConstraints;
    }
    int getNumberOfConstraintsJacobianNonZeroTerms()
    {
        return m_iNbConstraintsJacobianNonZeroTerms;
    }
    int getNumberOfHessianNonZeroTerms()
    {
        return m_iNbLowerTriHessianNonZeroTerms;
    }
    bool getObjective(double *pdblValue)
    {
        *pdblValue = *m_pdblIpoptInput[FOBJ];
        return true;
    }
    void setIsComputed(bool new_x)
    {
        if (new_x)
        {
            m_bIsComputed.fill(false);
        }
    }
    void setIpoptApp(Ipopt::SmartPtr<Ipopt::IpoptApplication> ipopt_app)
    {
        m_ipopt_app = ipopt_app;
    }
    Ipopt::SmartPtr<Ipopt::IpoptApplication> getIpoptApp()
    {
        return m_ipopt_app;
    }
    Ipopt::SolverReturn getSolverReturn()
    {
        return m_iSolverReturn;
    }
    void setSolverReturn(Ipopt::SolverReturn status)
    {
        m_iSolverReturn = status;
    }
    
    dynlibFunPtr getEntryPointFunction(functionKind  what)
    {
        return m_pEntryPointFunction[what];
    }
    
    private :
    std::wstring m_wstrCaller;
    char *m_pstrCaller;
    int m_iNbVariables;
    int m_iNbConstraints;
    int m_iNbIneqConstraints;
    int m_iNbEqConstraints;
    int m_iNbConstraintsJacobianNonZeroTerms;
    int m_iNbNonLinearVariables;

    int m_iNbUserHessianNonZeroTerms;
    int m_iNbLowerTriHessianNonZeroTerms;
    int *m_pIndexOfLowerTriHessianTerms;
    int *m_piHessianRow;
    int *m_piHessianCol;
        
    types::Double *m_pDblX0;
    types::Double *m_pDblVariablesLinearity;
    types::Double *m_pDblXLower;
    types::Double *m_pDblXUpper;
    types::Double *m_pDblConstraintsLhs;
    types::Double *m_pDblConstraintsRhs;
    types::Double *m_pDblConstraintsJacPattern;
    types::Sparse *m_pSpConstraintsJacPattern;
    types::Double *m_pDblHessianPattern;
    types::Sparse *m_pSpHessianPattern;
    types::Double *m_pDblConstraintsLinearity;
    
    std::array<functionKind,5> m_actualWhat = {FOBJ,DFOBJ,GOBJ,DGOBJ,DHOBJ};
    std::array<bool,5> m_bIsComputed;    
    std::array<int,5>  m_iNbFunctionCount;
    std::array<types::Callable *,6> m_pCallFunction;

    std::array<dynlibFunPtr,5> m_pEntryPointFunction;

    std::array<double *,5> m_pdblIpoptInput;
    std::array<int,5> m_iSizeOfScilabOutput;
    std::array<char *,6> m_pCallFunctionName;
    std::array<types::typed_list,6> m_pParameters;
    std::array<double,12> m_dblExtra;

    Ipopt::SmartPtr<Ipopt::IpoptApplication> m_ipopt_app; // a pointer to the current solver to retrieve the final statistics
    Ipopt::SolverReturn m_iSolverReturn;

    int    *m_piNonLinearVariables;
    double *m_pdblSolution;
    double *m_pdblLambda;
    double *m_pdblZL;
    double *m_pdblZU;
    double m_dblPrevFobjValue;
    double m_dblFobjSolution;
};

#endif
