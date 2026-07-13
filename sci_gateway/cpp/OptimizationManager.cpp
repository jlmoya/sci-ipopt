//
// Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
// Copyright (C) 2020 - UTC - Stéphane MOTTELET
//
// This file is hereby licensed under the terms of the GNU GPL v2.0,
// For more information, see the COPYING file which you should have received
//
//--------------------------------------------------------------------------

#include "OptimizationManager.hxx"
#include "argumentsPositions.hxx"

extern "C"
{
// #include "Sciprint.h"
    void C2F(dscal)(int*, double*, double*, int*);
}

#define OPTMNGR_DELETE(p) if (p != NULL) delete[] p

// OptimizationManager methods

OptimizationManager::OptimizationManager(const std::wstring& callerName)
{
    m_wstrCaller = callerName;
    m_pstrCaller = wide_string_to_UTF8(callerName.c_str());
    m_iNbVariables = 0;
    m_iNbConstraints = 0;
    m_iNbIneqConstraints = 0;
    m_iNbEqConstraints = 0;
    m_iNbConstraintsJacobianNonZeroTerms = 0;
    m_iNbUserHessianNonZeroTerms = 0;
    
    m_pdblSolution = NULL;
    m_pdblLambda = NULL;
    m_pdblZL = NULL;
    m_pdblZU = NULL;
    m_piNonLinearVariables = NULL;
    m_pIndexOfLowerTriHessianTerms = NULL;
    m_piHessianRow = NULL;
    m_piHessianCol = NULL;
      
    // Scilab internal types
    m_pDblX0 = NULL;
    m_pDblVariablesLinearity = NULL;
    m_pDblXLower = NULL;
    m_pDblXUpper = NULL;
    m_pDblConstraintsLhs = NULL;
    m_pDblConstraintsRhs = NULL;
    m_pDblConstraintsJacPattern = NULL;
    m_pSpConstraintsJacPattern = NULL;
    m_pDblHessianPattern = NULL;
    m_pSpHessianPattern = NULL;
    m_pDblConstraintsLinearity = NULL;

    // arrays
    m_pCallFunction.fill(NULL);
    m_pEntryPointFunction.fill(NULL);
    m_pdblIpoptInput.fill(NULL);
    m_iSizeOfScilabOutput.fill(0);
    m_pCallFunctionName.fill(NULL);
    m_iNbFunctionCount.fill(0);
    m_dblExtra.fill(std::numeric_limits<double>::quiet_NaN());    
    
    //doubles
    m_dblFobjSolution = std::numeric_limits<double>::quiet_NaN();
    m_dblPrevFobjValue = -std::numeric_limits<double>::infinity();
}

OptimizationManager::~OptimizationManager()
{
    for (int i=0; i<5; i++)
    {
        OPTMNGR_DELETE(m_pdblIpoptInput[i]);
        OPTMNGR_DELETE(m_pCallFunctionName[i]);
    }
    
    OPTMNGR_DELETE(m_pIndexOfLowerTriHessianTerms);
    OPTMNGR_DELETE(m_piHessianRow);
    OPTMNGR_DELETE(m_piHessianCol);
    OPTMNGR_DELETE(m_piNonLinearVariables);
    OPTMNGR_DELETE(m_pdblSolution);
    OPTMNGR_DELETE(m_pdblLambda);
    OPTMNGR_DELETE(m_pdblZU);
    OPTMNGR_DELETE(m_pdblZL);
}

bool OptimizationManager::parseMatrices(types::typed_list &in)
{
    std::list<int> iFullMatrices = {X_POS, VAR_LIN_TYPE_POS,CONSTR_LIN_TYPE_POS, CONSTR_RHS_POS,CONSTR_LHS_POS,X_LOWER_POS,X_UPPER_POS};
  
    std::list<int> iFullOrSparseMatrices = {SPARSE_DGOBJ_POS,SPARSE_DHOBJ_POS};

    // basic tests
    for (int i : iFullMatrices)
    {
        if (in[i]->isDouble() == false || in[i]->getAs<types::Double>()->isComplex())
        {
            Scierror(999, _("%s: Wrong type for input argument #%d (%s): A real matrix expected.\n"), getCallerStr(), i+1, pstrArgName[i]);
            return false;    
        }
    }

    for (int i : iFullOrSparseMatrices)
    {
        if ((in[i]->isDouble() == false || in[i]->getAs<types::Double>()->isComplex()) && in[i]->isSparse() == false)
        {
            Scierror(999, _("%s: Wrong type for input argument #%d (%s): A real matrix expected.\n"), getCallerStr(), i+1, pstrArgName[i]);
            return false;    
        }
    }

    if (in[X_POS]->getAs<types::Double>()->isEmpty())
    {
        Scierror(999, _("%s: Wrong size for input argument #%d (%s): Non-empty matrix expected.\n"), getCallerStr(), X_POS, pstrArgName[X_POS]);
        return false;    
    }
    else
    {
        m_pDblX0 = in[X_POS]->getAs<types::Double>();
        m_iNbVariables = m_pDblX0->getSize();
    }
    
    for (int i : {VAR_LIN_TYPE_POS, X_LOWER_POS, X_UPPER_POS})
    {
        if (in[i]->getAs<types::Double>()->isEmpty() == false && in[i]->getAs<types::Double>()->getSize() != m_iNbVariables)
        {
            Scierror(999,_("%s: Wrong size for input argument #%d (%s): A vector of size %d expected.\n"), getCallerStr(), i+1, pstrArgName[i], m_iNbVariables);
            return false;    
        }  
    }
    
    // Verify vectors and matrices
    
    m_pDblXLower                 = in[X_LOWER_POS]->getAs<types::Double>();
    m_pDblXUpper                 = in[X_UPPER_POS]->getAs<types::Double>();
    m_pDblVariablesLinearity     = in[VAR_LIN_TYPE_POS]->getAs<types::Double>();
    m_pDblConstraintsLhs         = in[CONSTR_LHS_POS]->getAs<types::Double>();
    m_pDblConstraintsRhs         = in[CONSTR_RHS_POS]->getAs<types::Double>();
    m_pDblConstraintsLinearity   = in[CONSTR_LIN_TYPE_POS]->getAs<types::Double>();
        
    int iCLhs = m_pDblConstraintsLhs->getSize();
    int iCRhs = m_pDblConstraintsRhs->getSize();
        
    m_iNbConstraints = std::max(iCLhs, iCRhs);
    if (m_iNbConstraints > 0)
    {
        if (std::min(iCLhs, iCRhs) > 0 && iCLhs != iCRhs)
        {
            Scierror(999,_("%s: Arguments #%d (%s) and #%d (%s) must have the same sizes.\n"), getCallerStr(), 
            pstrArgName[CONSTR_LHS_POS], CONSTR_LHS_POS+1, pstrArgName[CONSTR_RHS_POS], CONSTR_RHS_POS+1);
            return false;    
        }
        if (in[GOBJ_POS]->isDouble()) // i.e. is an empty matrix
        {
            Scierror(999,_("%s: parameter %d (%s) should be a function or a list.\n"), getCallerStr(), GOBJ_POS+1, pstrArgName[GOBJ_POS]);
            return false;    
        }
        m_iNbIneqConstraints = m_iNbConstraints;
        if (iCLhs != 0 && iCRhs != 0) // constraints may not be all inequalities
        {
            for (int i = 0; i < m_iNbConstraints; i++)
            {
                m_iNbEqConstraints += (int) (m_pDblConstraintsLhs->get(i) == m_pDblConstraintsRhs->get(i));
            }
             m_iNbIneqConstraints = m_iNbConstraints - m_iNbEqConstraints;
        }
    }
    else if (in[GOBJ_POS]->isDouble() == false) // i.e. is not a empty matrix hence is a function or a list
    {
        Scierror(999,_("%s: Arguments #%d (%s) and #%d (%s) cannot be both empty.\n"), getCallerStr(), CONSTR_LHS_POS+1, 
        pstrArgName[CONSTR_LHS_POS], CONSTR_RHS_POS+1, pstrArgName[CONSTR_RHS_POS]);
        return false;    
    }

    // Verify constraints jacobian
    if (m_iNbConstraints > 0)
    {
        if (in[SPARSE_DGOBJ_POS]->isDouble())
        {
            m_pDblConstraintsJacPattern  = in[SPARSE_DGOBJ_POS]->getAs<types::Double>();
                
            int iNbMaxTerms = (int) m_iNbVariables*m_iNbConstraints;

            if (m_pDblConstraintsJacPattern->isEmpty())
            {
                m_iNbConstraintsJacobianNonZeroTerms = iNbMaxTerms; // full pattern will be set in OptimizationManager::getConstraintsJacobianPattern
            }
            else if (m_pDblConstraintsJacPattern->getCols() != 2)
            {
                Scierror(999,_("%s: Wrong size for input argument #%d (%s): A matrix with %d columns expected.\n"), getCallerStr(), 
                SPARSE_DGOBJ_POS+1, pstrArgName[SPARSE_DGOBJ_POS], 2);            
                return false;    
            }
            else if (m_pDblConstraintsJacPattern->getRows() > iNbMaxTerms)
            {
                Scierror(999,_("%s: Wrong size for input argument #%d (%s): A matrix with at most %d rows expected.\n"), getCallerStr(), 
                SPARSE_DGOBJ_POS+1, pstrArgName[SPARSE_DGOBJ_POS], iNbMaxTerms);            
                return false;    
            }
            else
            {
                for (int i=0; i<m_pDblConstraintsJacPattern->getRows(); i++)
                {
                    int iCstr = m_pDblConstraintsJacPattern->get(i,0);
                    int iVar = m_pDblConstraintsJacPattern->get(i,1);
                    if (iCstr < 1 || iCstr > m_iNbConstraints || iVar < 1 || iVar > m_iNbVariables)
                    {
                        Scierror(999,_("%s: invalid value in argument #%d (%s) at row %d.\n"), getCallerStr(), SPARSE_DGOBJ_POS+1, pstrArgName[SPARSE_DGOBJ_POS], i+1);    
                        return false;                            
                    }
                }
                m_iNbConstraintsJacobianNonZeroTerms = m_pDblConstraintsJacPattern->getRows();
            }
        
            if (m_pDblConstraintsLinearity->isEmpty() == false && m_pDblConstraintsLinearity->getSize() != m_iNbConstraints)
            {
                Scierror(999,_("%s: Wrong size for input argument #%d (%s): A vector of size %d expected.\n"), getCallerStr(), 
                CONSTR_LIN_TYPE_POS+1, pstrArgName[CONSTR_LIN_TYPE_POS], m_iNbConstraints);
                return false;     
            }        
        }
        else // Scilab sparse will give the pattern
        {
            m_pSpConstraintsJacPattern  = in[SPARSE_DGOBJ_POS]->getAs<types::Sparse>();
            if (m_pSpConstraintsJacPattern->getRows() != m_iNbConstraints || m_pSpConstraintsJacPattern->getCols() != m_iNbVariables)
            {
                Scierror(999,_("%s: Wrong size for input argument #%d (%s): A matrix with %d rows and %d columns expected.\n"), getCallerStr(), 
                SPARSE_DGOBJ_POS+1, pstrArgName[SPARSE_DGOBJ_POS], m_iNbConstraints, m_iNbVariables);            
                return false;     
            }
            m_iNbConstraintsJacobianNonZeroTerms = m_pSpConstraintsJacPattern->nonZeros(); // empty matrix
        }
    }
    else // m_iNbConstraints == 0
    {
        if (in[SPARSE_DGOBJ_POS]->isDouble() == false || in[SPARSE_DGOBJ_POS]->getAs<types::Double>()->isEmpty() == false)
        {
            Scierror(999,_("%s: Wrong size for input argument #%d (%s): Should be an empty matrix.\n"), getCallerStr(), 
            SPARSE_DGOBJ_POS+1, pstrArgName[SPARSE_DGOBJ_POS]);
            return false;             
        }
        if (m_pDblConstraintsLinearity->isEmpty() == false)
        {
            Scierror(999,_("%s: Wrong size for input argument #%d (%s): Should be an empty matrix.\n"), getCallerStr(), 
            CONSTR_LIN_TYPE_POS+1, pstrArgName[CONSTR_LIN_TYPE_POS]);
            return false;             
        }
        m_pDblConstraintsJacPattern  = in[SPARSE_DGOBJ_POS]->getAs<types::Double>();        
    }
    
    // Verify Hessian lower triangle sparsity pattern.
    // Unlike the constraints Jacobian pattern, we fully verify the individual terms because
    // we need the actual number of non-zero terms *in the lower triangular part* for Ipopt app initialization,
    // which is done before the request for the pattern itself.
    
    std::string strHessOpt;
    getIpoptApp()->Options()->GetStringValue("hessian_approximation", strHessOpt, "");
    
    
//    if (!in[DHOBJ_POS]->isDouble()) // i.e. is a function or a list
    if (strHessOpt == "exact")
    {
        int iNbMaxTerms = (int) m_iNbVariables*(m_iNbVariables+1.0)/2.0;
        
        if (in[SPARSE_DHOBJ_POS]->isDouble())
        {
            m_pDblHessianPattern = in[SPARSE_DHOBJ_POS]->getAs<types::Double>();
            if (m_pDblHessianPattern->isEmpty())
            {
                m_iNbUserHessianNonZeroTerms = iNbMaxTerms; // full pattern will be set in OptimizationManager::getHessianLowerTriPattern
                m_iNbLowerTriHessianNonZeroTerms = iNbMaxTerms;
            }
            else if (m_pDblHessianPattern->getCols() != 2)
            {
                Scierror(999,_("%s: Wrong size for input argument #%d (%s): A matrix with %d columns expected.\n"), getCallerStr(),
                SPARSE_DHOBJ_POS+1, pstrArgName[SPARSE_DHOBJ_POS],2);            
                return false;    
            }
            else
            {
                // check pattern and build index of lower triangle terms in vector of non-zero terms
                return parseHessianPattern(m_pDblHessianPattern->get(), m_pDblHessianPattern->getRows(), false);
            }            
        }
        else // pattern given by a sparse matrix
        {
            m_pSpHessianPattern = in[SPARSE_DHOBJ_POS]->getAs<types::Sparse>();
            if (m_pSpHessianPattern->getRows() != m_iNbVariables || m_pSpHessianPattern->getCols() != m_iNbVariables)
            {
                Scierror(999,_("%s: Wrong size for input argument #%d (%s): a sparse matrix with %d rows and %d columns expected.\n"), getCallerStr(), 
                SPARSE_DGOBJ_POS+1, pstrArgName[SPARSE_DHOBJ_POS], m_iNbVariables, m_iNbVariables);            
                return false;     
            }
            int iNonZeros = m_pSpHessianPattern->nonZeros();
            int *piRowCols = new int[iNonZeros*2]; // deleted in OptimizationManager::parseHessianPattern

            m_pSpHessianPattern->outputRowCol(piRowCols);
            // check pattern, build index of lower triangle terms in vector of non-zero terms
            // and ask for deletion of piRowCols before exit
            return parseHessianPattern(piRowCols, iNonZeros, true);            
        }        
    }
        
    return true;
}


bool OptimizationManager::parseFunction(types::typed_list &in, int iPos, functionKind  what)
{
    types::InternalType *pIn = in[iPos];

    std::string strJacOpt, strHessOpt, strDerTest, strJacCConstOpt, strJacDConstOpt;
    getIpoptApp()->Options()->GetStringValue("jacobian_approximation", strJacOpt, "");
    getIpoptApp()->Options()->GetStringValue("jac_d_constant", strJacDConstOpt, "");
    getIpoptApp()->Options()->GetStringValue("jac_c_constant", strJacCConstOpt, "");
    getIpoptApp()->Options()->GetStringValue("jacobian_approximation", strJacOpt, "");
    getIpoptApp()->Options()->GetStringValue("hessian_approximation", strHessOpt, "");
    getIpoptApp()->Options()->GetStringValue("derivative_test", strDerTest, "");
    
    if (pIn->isDouble() && pIn->getAs<types::Double>()->isEmpty())
    {
        m_pCallFunction[what] = NULL;

        if (what == FOBJ || 
            (what == DGOBJ && strJacOpt == "exact" && in[GOBJ_POS]->isString()) ||
            (what == DHOBJ && strHessOpt == "exact" &&  ((strJacDConstOpt == "no" && m_iNbIneqConstraints > 0) || (strJacCConstOpt == "no" && m_iNbEqConstraints > 0))) ||
            (what == DFOBJ && in[GOBJ_POS]->isString()) ) 
        {
            Scierror(999, _("%s: parameter %d (%s) should be a string, a function or a list.\n"), getCallerStr(), iPos+1, pstrArgName[iPos]);
            return false;
        }
        if (what == DFOBJ) 
        {
            // user function f proto is [f,grad_f] = fun(x)
            m_actualWhat[DFOBJ] = FOBJ;
        }
        else if (what == DHOBJ)
        {
            // user function f prototype will be supposed  [f,grad_f,hess_f] = fun(x)
            // if user further defines hessian_approximation == "exact"
            // Check of the actual proto will be done during call
            if (strHessOpt == "exact" && in[DFOBJ_POS]->isDouble() == false)
            {
                Scierror(999, _("%s: parameter %d (%s) must be empty.\n"), getCallerStr(), iPos+1, pstrArgName[iPos]);
                return false;                
            }
            m_actualWhat[DHOBJ] = FOBJ;
        }
        else if (what == DGOBJ && in[GOBJ_POS]->isString() == false)
        {
            // user function g proto is [f,jac_g] = fun(x)
            m_actualWhat[DGOBJ] = GOBJ;
        }
    }
    else if (pIn->isCallable())
    {
        m_pCallFunction[what] = pIn->getAs<types::Callable>();
    }
    else if (pIn->isList())
    {
        types::List *pList = pIn->getAs<types::List>();
        types::InternalType *pIFirst =  pList->get(0);
        if (pIFirst->isCallable())
        {
            m_pCallFunction[what] = pIFirst->getAs<types::Callable>();            
        }
        else
        {
            Scierror(999,_("%s: first element of parameter %d (%s) should be a function.\n"), getCallerStr(), iPos+1, pstrArgName[iPos]);
            return false;    
        }
        for (int i=1; i<pList->getSize(); i++)
        {
            m_pParameters[what].push_back(pList->get(i));
        }        
    }
    else if (pIn->isString())
    {
        ConfigVariable::EntryPointStr* func = ConfigVariable::getEntryPoint(pIn->getAs<types::String>()->get(0));
        if (func == NULL)
        {
            Scierror(999,_("%s: unable to find entry point %ls.\n"), getCallerStr(), pIn->getAs<types::String>()->get(0));
            return false;
        }
        else
        {
            m_pEntryPointFunction[what] = func->functionPtr;
        }
    }
    else
    {
        Scierror(999, _("%s: parameter %d (%s) should be a string, a function or a list.\n"), getCallerStr(), iPos+1, pstrArgName[iPos]);
        return false;    
    }

    if (m_pCallFunction[what] != NULL)
    {
        m_pCallFunctionName[what] = wide_string_to_UTF8(m_pCallFunction[what]->getName().c_str());        
    }
        
    return true;
}

bool OptimizationManager::getX0(double *pdblX0)
{
    if (m_pDblX0->isEmpty() == false)
    {
        memcpy(pdblX0, m_pDblX0->get(), sizeof(double)*m_iNbVariables);
        return true;
    }
    return false;
}

bool OptimizationManager::getVariablesLinearity(Ipopt::TNLP::LinearityType* piVarLin)
{
    if (m_pDblVariablesLinearity->isEmpty() == false)
    {
        for (int k=0; k < m_iNbVariables; k++)
        {
            piVarLin[k] = m_pDblVariablesLinearity->get(k) == 0. ? Ipopt::TNLP::LINEAR : Ipopt::TNLP::NON_LINEAR;
        }
        return true;
    }
    return false;
}

bool OptimizationManager::getNonLinearVariables(Ipopt::Index *pos_nonlin_vars)
{
    if (m_piNonLinearVariables != NULL)
    {
        for (int i = 0; i < m_iNbNonLinearVariables; i++)
        {
            pos_nonlin_vars[i] = m_piNonLinearVariables[i];
        }
        return true;        
    }
    return false;
}

bool OptimizationManager::getSolution(double *pdblX, double *pdblLambda, double *pdblZL, double *pdblZU, double *pdblF)
{
    if (m_pdblSolution != NULL)
    {
        memcpy(pdblX, m_pdblSolution, sizeof(double)*m_iNbVariables);
        memcpy(pdblLambda, m_pdblLambda, sizeof(double)*m_iNbConstraints);        
        memcpy(pdblZL, m_pdblZL, sizeof(double)*m_iNbVariables);        
        memcpy(pdblZU, m_pdblZU, sizeof(double)*m_iNbVariables);        
        *pdblF = m_dblFobjSolution;
    }
        
    return true;
}

bool OptimizationManager::getBoundsInfo(double *pdblXLower, double *pdblXUpper, double *pdblConstraintsLhs, double *pdblConstraintsRhs)
{
    int iOne = 1;
    if (m_pDblXLower->isEmpty() == false)
    {
        memcpy(pdblXLower, m_pDblXLower->get(), sizeof(double)*m_iNbVariables);
    }
    else
    {
        std::fill(pdblXLower, pdblXLower+m_iNbVariables, -std::numeric_limits<double>::infinity());
    }
    
    if (m_pDblXUpper->isEmpty() == false)
    {
        memcpy(pdblXUpper, m_pDblXUpper->get(), sizeof(double)*m_iNbVariables);
    }
    else
    {
        std::fill(pdblXUpper, pdblXUpper+m_iNbVariables, std::numeric_limits<double>::infinity());
    }
    
    if (m_pDblConstraintsLhs->isEmpty() == false)
    {
        memcpy(pdblConstraintsLhs, m_pDblConstraintsLhs->get(), sizeof(double)*m_iNbConstraints);
    }
    else
    {
        std::fill(pdblConstraintsLhs, pdblConstraintsLhs+m_iNbConstraints, -std::numeric_limits<double>::infinity());
    }
    
    if (m_pDblConstraintsRhs->isEmpty() == false)
    {
        memcpy(pdblConstraintsRhs, m_pDblConstraintsRhs->get(), sizeof(double)*m_iNbConstraints);
    }
    else
    {
        std::fill(pdblConstraintsRhs, pdblConstraintsRhs+m_iNbConstraints, std::numeric_limits<double>::infinity());
    }
    return true;
}

bool OptimizationManager::getComputed(double *pdblOut, functionKind  what)
{    
    if (m_pdblIpoptInput[what] != NULL)
    {
        // WARNING:
        // for the Hessian m_iNbLowerTriHessianNonZeroTerms is less than m_iSizeOfScilabOutput[DHOBJ]
        // if user gives redondant (i,j) amd (j,i) terms
        int iSize = (what == DHOBJ ? m_iNbLowerTriHessianNonZeroTerms : m_iSizeOfScilabOutput[what]);
        // check NaNs and Infs
        double *pdbl = m_pdblIpoptInput[what];
        for (int k = 0; k < iSize; k++, pdbl++)
        {
            if (!std::isfinite(*pdbl))
            {
                return false;
            }
        }       
        memcpy(pdblOut, m_pdblIpoptInput[what], sizeof(double)*iSize);
        return true;            
    }
    return false;
}

bool OptimizationManager::getConstraintsJacobianPattern(int *piRow, int *piCol)
{
    if (m_pSpConstraintsJacPattern != NULL)
    {
        // Pattern given by Scilab Sparse
        int iNonZeros = m_pSpConstraintsJacPattern->nonZeros();
        int *piRowCols = new int[iNonZeros*2];
        m_pSpConstraintsJacPattern->outputRowCol(piRowCols);
        for (int k=0; k < iNonZeros; k++)
        {
            //pattern is given with Scilab indices (starting from 1)
            piRow[k] = (int)piRowCols[k] - 1;
            piCol[k] = (int)piRowCols[k+iNonZeros] - 1;;                    
        }
        delete[] piRowCols;
    }
    else if (m_pDblConstraintsJacPattern->isEmpty() == false)
    {
        // Pattern given by (nnz,2) sized full matrix
        for (int k=0; k < m_iNbConstraintsJacobianNonZeroTerms; k++)
        {
            //pattern is given with Scilab indices (starting from 1)
            piRow[k] = (int) m_pDblConstraintsJacPattern->get(k,0) - 1;
            piCol[k] = (int) m_pDblConstraintsJacPattern->get(k,1) - 1;
        }
    }
    else // full jacobian (column-major order)
    {
        int k = 0;
        for (int j = 0; j < m_iNbVariables; j++)
        {
            for (int i = 0; i < m_iNbConstraints; i++)
            {
                piRow[k] = i; 
                piCol[k] = j;
                ++k;
            }                
        }
    }
    return true;
}

template<typename T> 
bool OptimizationManager::parseHessianPattern(T *piRowCols, int iSize, bool bDelete)
{
    m_iNbUserHessianNonZeroTerms = iSize;
    m_iNbLowerTriHessianNonZeroTerms = 0;

    m_piHessianRow = new int[iSize]; // OptimizationManager class member, deleted at object delete
    m_piHessianCol = new int[iSize]; // OptimizationManager class member, deleted at object delete
    m_pIndexOfLowerTriHessianTerms = new int[iSize];  // OptimizationManager class member, deleted at object delete

    for (int k=0; k < iSize; k++)
    {
        //pattern is given with Scilab indices (starting from 1)
        int iRow = (int)piRowCols[k];
        int iCol = (int)piRowCols[k+iSize];;                    
        if (iRow < 1 || iRow > m_iNbVariables || iCol < 1 || iCol > m_iNbVariables)
        {
            Scierror(999,_("%s: invalid position (%d,%d) in argument #%d (%s).\n"), getCallerStr(), 
            iRow, iCol, SPARSE_DHOBJ_POS+1, pstrArgName[SPARSE_DHOBJ_POS]);
            if (bDelete) // should not be necessary
            {
                delete[] piRowCols;
            }
            return false;                            
        }
        if (iRow >= iCol) // keep only lower triangular part of pattern
        {
            m_piHessianRow[m_iNbLowerTriHessianNonZeroTerms] = iRow-1;
            m_piHessianCol[m_iNbLowerTriHessianNonZeroTerms] = iCol-1;
            m_pIndexOfLowerTriHessianTerms[m_iNbLowerTriHessianNonZeroTerms] = k; // remember index of non-zero terms to be copied in OptimizationManager::computeHessian
            m_iNbLowerTriHessianNonZeroTerms++;
        }
    }
    if (bDelete)
    {
        delete[] piRowCols;
    }
    if (m_iNbLowerTriHessianNonZeroTerms == 0)
    {
        Scierror(999,_("%s: invalid Hessian pattern in argument #%d (%s).\n"), getCallerStr(), 
        SPARSE_DHOBJ_POS+1, pstrArgName[SPARSE_DHOBJ_POS]);        
        return false;
    }
    return true;
}

bool OptimizationManager::getHessianLowerTriPattern(int *piRow, int *piCol)
{
    if (m_pSpHessianPattern != NULL || m_pDblHessianPattern->isEmpty() == false)
    {
        memcpy(piRow, m_piHessianRow, m_iNbLowerTriHessianNonZeroTerms*sizeof(int));
        memcpy(piCol, m_piHessianCol, m_iNbLowerTriHessianNonZeroTerms*sizeof(int));
    }
    else // full hessian 
    {
        int k = 0;
        for (int j = 0; j < m_iNbVariables; j++)
        {
            for (int i = j; i < m_iNbVariables; i++)
            {
                piRow[k] = i; 
                piCol[k] = j;
                ++k;
            }                
        }        
    }
        
    return true;
}

bool OptimizationManager::getConstraintsLinearity(Ipopt::TNLP::LinearityType* piConstrLin)
{
    if (m_pDblConstraintsLinearity->isEmpty() == false)
    {
        for (int k=0; k < m_iNbConstraints; k++)
        {
            piConstrLin[k] = m_pDblConstraintsLinearity->get(k) == 1. ? Ipopt::TNLP::LINEAR : Ipopt::TNLP::NON_LINEAR;
        }
        return true;
    }
    return false;
}

bool OptimizationManager::setAllStuff()
{
    m_pdblIpoptInput[FOBJ] = new double[1];
    m_iSizeOfScilabOutput[FOBJ] = 1;
    
    m_pdblIpoptInput[DFOBJ] = new double[m_iNbVariables];
    m_iSizeOfScilabOutput[DFOBJ] = m_iNbVariables;

    // set variables linearity vector + number of non-linear variables
    if (m_pDblVariablesLinearity->isEmpty())
    {
        // all are supposed non-linear
        m_iNbNonLinearVariables = -1; 
    }
    else
    {
        m_piNonLinearVariables = new int[m_iNbVariables];
        m_iNbNonLinearVariables = 0;
        for (int i=0; i<m_iNbVariables; i++)
        {
            if (m_pDblVariablesLinearity->get(i) != 0) // non-linear
            {
                m_piNonLinearVariables[m_iNbNonLinearVariables] = i;
                m_iNbNonLinearVariables++;
            }
        }
    }

    if (m_iNbConstraints > 0)
    {
        m_pdblIpoptInput[GOBJ] = new double[m_iNbConstraints];
        m_iSizeOfScilabOutput[GOBJ] = m_iNbConstraints;
        m_pdblIpoptInput[DGOBJ] = new double[m_iNbConstraintsJacobianNonZeroTerms];
        m_iSizeOfScilabOutput[DGOBJ] = m_iNbConstraintsJacobianNonZeroTerms;
    }
    
    std::string strHessOpt;
    getIpoptApp()->Options()->GetStringValue("hessian_approximation", strHessOpt, "");
    
//    if (m_pCallFunction[DHOBJ] != NULL)
    if (strHessOpt == "exact")
    {
        m_pdblIpoptInput[DHOBJ] = new double[m_iNbLowerTriHessianNonZeroTerms];        
        m_iSizeOfScilabOutput[DHOBJ] = m_iNbUserHessianNonZeroTerms;
    }

    m_pdblSolution  = new double[m_iNbVariables];
    m_pdblLambda  = new double[m_iNbConstraints];
    m_pdblZL = new double[m_iNbVariables];
    m_pdblZU = new double[m_iNbVariables];
    
    return true;
}

bool OptimizationManager::computeFunctionsAndDerivatives(double *x, bool new_x, functionKind origWhat)
{
    int iOne = 1;
    int iRetCount = 1;
    char errorMsg[256] = "";

    types::typed_list in;
    types::typed_list out;
    types::optional_list opt;
    types::Double* pDblX = NULL;
    types::Bool *pBNewX = NULL;
    
    functionKind  whatArray[3];
    functionKind  what = m_actualWhat[origWhat];
    types::Callable *pCall = m_pCallFunction[what];
    
    if (pCall == NULL)
    {
        return false;
    }
    
    if (m_bIsComputed[origWhat]) // if already computed (can occur with gradObj or gradConstr)
    {
        return true;
    }

    pDblX = new types::Double(m_iNbVariables, 1);
    memcpy(pDblX->get(), x, sizeof(double)*m_iNbVariables);
    in.push_back(pDblX);

    pBNewX = new types::Bool(new_x);
    in.push_back(pBNewX);
    
    // optional user input parameters 
    for (auto pIn : m_pParameters[what])
    {
        in.push_back(pIn);
    }

    for (auto pIn : in)
    {
        pIn->IncreaseRef();
    }

    // when grad_f (resp jac_g) is required but not provided (empty matrix in call) then user function has to return also f (resp g)
    // when hessian is required then user function has to return f,grad_f,h
    
    iRetCount = 0;
    for (int i = FOBJ; i <= origWhat; i++) // enum is {FOBJ,DFOBJ,GOBJ,DGOBJ,DHOBJ}
    {
        if (m_actualWhat[i] == what)
        {
             m_bIsComputed[i] = true;
             whatArray[iRetCount++] = (functionKind) i;
        }
    }
    
    try
    {
        // new std::wstring(L"") is deleted in destructor of ast::CommentExp
        ConfigVariable::clearLastError();
        pCall->invoke(in, opt, iRetCount, out, ast::CommentExp(Location(), new std::wstring(L"")));
    }
    catch(ast::InternalError& ie)
    {
        for (auto pIn : m_pParameters[what])
        {
            pIn->DecreaseRef();
        }
        if (ConfigVariable::getLastErrorFunction() == L"")
        {
            // necessary for errors such as "Wrong number of input arguments"
            std::wostringstream ostr;
            ostr << m_pCallFunctionName[what] << L": " <<  ie.GetErrorMessage();
            ie.SetErrorMessage(ostr.str());
        }
        throw(ie);
    }

    for (auto pOut : out) // safety IncreaseRef
    {
        pOut->IncreaseRef();
    }

    for (auto pIn : in)
    {
        pIn->DecreaseRef();
        if (pIn->isDeletable())
        {
            delete pIn;
        }
    }

    if (out.size() != iRetCount)
    {
        snprintf(errorMsg, sizeof(errorMsg), _("%s: Wrong number of output argument(s): %d expected.\n"), m_pCallFunctionName[what], iRetCount);
        throw ast::InternalError(errorMsg);
    }
        
    // check type of output and copy to destination if type and dimensions are ok.
    for (int i = 0; i < iRetCount; i++)
    {
        if (whatArray[i] == DGOBJ && m_pSpConstraintsJacPattern != NULL) 
        {
            // Sparse matrix expected
            if (out[i]->isSparse() == false || out[i]->getAs<types::Sparse>()->isComplex())
            {
                snprintf(errorMsg, sizeof(errorMsg), _("%s: Wrong type for output argument #%d: sparse real matrix expected.\n"), m_pCallFunctionName[what], i+1);
                throw ast::InternalError(errorMsg);
            }
            types::Sparse *pSp = out[i]->getAs<types::Sparse>();
            if (pSp->getRows() != m_iNbConstraints || pSp->getCols() != m_iNbVariables)
            {
                snprintf(errorMsg, sizeof(errorMsg), _("%s: Wrong size for output argument #%d: a matrix with %d rows and %d columns expected.\n"), m_pCallFunctionName[what], i+1, m_iNbConstraints, m_iNbVariables);
                throw ast::InternalError(errorMsg);
            }
            if (pSp->nonZeros() > m_iNbConstraintsJacobianNonZeroTerms)
            {
                snprintf(errorMsg, sizeof(errorMsg), _("%s: Wrong sparsity pattern for output argument #%d: at most %d non-zero elements expected.\n"), m_pCallFunctionName[what], i+1, m_iNbConstraintsJacobianNonZeroTerms);
                throw ast::InternalError(errorMsg);                    
            }
            if (pSp->nonZeros() == m_iNbConstraintsJacobianNonZeroTerms)
            {
                // if the number of non-zeros match, we don't check the sparsity and cross fingers...
                memcpy(m_pdblIpoptInput[whatArray[i]], pSp->matrixReal->valuePtr(), sizeof(double)*m_iNbConstraintsJacobianNonZeroTerms);                                       
            }
            else
            {
                // monotonicity of pattern is not achieved, let's try to fix this
                if (getSparseValuesFromIncompletePattern(m_pdblIpoptInput[whatArray[i]], m_pSpConstraintsJacPattern, pSp) == false)
                {
                    snprintf(errorMsg, sizeof(errorMsg), _("%s: Wrong sparsity pattern for output argument #%d.\n"), m_pCallFunctionName[what], i+1);
                    throw ast::InternalError(errorMsg);                                            
                }                 
            }            
        }
        else if (whatArray[i] == DHOBJ) // in case of [f,g,h] = fun(x) syntax, i.e. options.HessianFcn = "objective"
        {
            checkHessianOutput(out[i], DHOBJ);
        }
        else // all other cases, typically full matrix output (objective, gradient, constraints, non-zero terms of constraints Jacobian) 
        {
            if (out[i]->isDouble() == false || out[i]->getAs<types::Double>()->isComplex() || out[i]->getAs<types::Double>()->isEmpty())
            {
                snprintf(errorMsg, sizeof(errorMsg), _("%s: Wrong type for output argument #%d: Real matrix expected.\n"), m_pCallFunctionName[what], i+1);
                throw ast::InternalError(errorMsg);            
            }
            if (out[i]->getAs<types::Double>()->getSize() != m_iSizeOfScilabOutput[whatArray[i]])
            {
                snprintf(errorMsg, sizeof(errorMsg), _("%s: Wrong size for output argument #%d: A matrix of size %d expected.\n"), m_pCallFunctionName[what], i+1,m_iSizeOfScilabOutput[whatArray[i]]);   
                throw ast::InternalError(errorMsg);
            }
            memcpy(m_pdblIpoptInput[whatArray[i]], out[i]->getAs<types::Double>()->get(), sizeof(double)*m_iSizeOfScilabOutput[whatArray[i]]);                                
        }
    }

    for (auto pOut : out)
    {
        pOut->DecreaseRef();
        if (pOut->isDeletable())
        {
            delete pOut;
        }
    }

    return true;
}

bool OptimizationManager::computeHessian(double *x, bool new_x, double obj_factor, double *lambda, bool new_lambda, functionKind what)
{
    int iOne = 1;
    int iRetCount = 1;
    char errorMsg[256] = "";

    types::typed_list in;
    types::typed_list out;
    types::optional_list opt;
    types::Double* pDblX = NULL;
    types::Bool *pBNewX = NULL;
    types::Double* pDblObjFactor = NULL;
    types::Double* pDblLambda = NULL;
    types::Bool* pBNewLambda = NULL;
    // what is always DHOBJ
    types::Callable *pCall = m_pCallFunction[what];

    if (m_actualWhat[what] != what) // can be the case if Hessian is computed in objective function (implies only linear constraints, i.e. Hessian of Lagrangian is Hessian of objective)
    {
        if (OptimizationManager::computeFunctionsAndDerivatives(x, new_x, what))
        {
            // don't forget to post-multiply Hessian of objective by obj_factor
            if (obj_factor != 1.0)
            {
                C2F(dscal)(&m_iNbLowerTriHessianNonZeroTerms, &obj_factor, m_pdblIpoptInput[what], &iOne);
            }
            return true;
        }
        return false;
    }

    if (pCall == NULL)
    {
        return false;
    }
    
    pDblX = new types::Double(m_iNbVariables, 1);
    memcpy(pDblX->get(), x, sizeof(double)*m_iNbVariables);
    in.push_back(pDblX);

    pBNewX = new types::Bool(new_x);
    in.push_back(pBNewX);

    pDblObjFactor = new types::Double(obj_factor);
    in.push_back(pDblObjFactor);
    pDblLambda = new types::Double(m_iNbConstraints, 1);
    memcpy(pDblLambda->get(), lambda, sizeof(double)*m_iNbConstraints);
    in.push_back(pDblLambda);
    pBNewLambda = new types::Bool(new_lambda);
    in.push_back(pBNewLambda);
    
    // optional user input parameters 
    for (auto pIn : m_pParameters[what])
    {
        in.push_back(pIn);
    }

    for (auto pIn : in)
    {
        pIn->IncreaseRef();
    }

    try
    {
        // new std::wstring(L"") is deleted in destructor of ast::CommentExp
        pCall->invoke(in, opt, 1, out, ast::CommentExp(Location(), new std::wstring(L"")));
    }
    catch(ast::InternalError& ie)
    {
        for (auto pIn : m_pParameters[what])
        {
            pIn->DecreaseRef();
        }
        if (ConfigVariable::getLastErrorFunction() == L"")
        {
            // necessary for errors such as "Wrong number of input arguments"
            std::wostringstream ostr;
            ostr << m_pCallFunctionName[what] << L": " <<  ie.GetErrorMessage();
            ie.SetErrorMessage(ostr.str());
        }
        throw(ie);
    }

    for (auto pOut : out) // safety IncreaseRef
    {
        pOut->IncreaseRef();
    }

    for (auto pIn : in)
    {
        pIn->DecreaseRef();
        if (pIn->isDeletable())
        {
            delete pIn;
        }
    }

    if (out.size() != iRetCount)
    {
        snprintf(errorMsg, sizeof(errorMsg), _("%s: Wrong number of output argument(s): %d expected.\n"), m_pCallFunctionName[what], iRetCount);
        throw ast::InternalError(errorMsg);
    }
    
    checkHessianOutput(out[0], what);

    for (auto pOut : out)
    {
        pOut->DecreaseRef();
        if (pOut->isDeletable())
        {
            delete pOut;
        }
    }

    return true;
}

void OptimizationManager::checkHessianOutput(types::InternalType *pIn, functionKind what)
{
    char errorMsg[256] = "";

    if (m_pSpHessianPattern != NULL)
    {
        // Hessian is supposed to be given as a Scilab sparse matrix
        if (pIn->isSparse() == false || pIn->getAs<types::Sparse>()->isComplex())
        {
            snprintf(errorMsg, sizeof(errorMsg), _("%s: Wrong type for output argument #%d: sparse real matrix expected.\n"), m_pCallFunctionName[what], 1);
            throw ast::InternalError(errorMsg);
        }
        types::Sparse *pSp = pIn->getAs<types::Sparse>();
        if (pSp->getRows() != m_iNbVariables || pSp->getCols() != m_iNbVariables)
        {
            snprintf(errorMsg, sizeof(errorMsg), _("%s: Wrong size for output argument #%d: a matrix with %d rows and %d columns expected.\n"), m_pCallFunctionName[what], 1, m_iNbVariables, m_iNbVariables);
            throw ast::InternalError(errorMsg);
        }
        if (pSp->nonZeros() > m_iNbUserHessianNonZeroTerms)
        {
            snprintf(errorMsg, sizeof(errorMsg), _("%s: Wrong sparsity pattern for output argument #%d: at most %d non-zero elements expected.\n"), m_pCallFunctionName[what], 1, m_iNbConstraintsJacobianNonZeroTerms);
            throw ast::InternalError(errorMsg);                    
        }
        if (m_iNbLowerTriHessianNonZeroTerms == m_iNbUserHessianNonZeroTerms)
        {
            // User specified lower triangular part of Hessian, hence we (almost) copy data as is,
            // i.e. we don't use the m_pIndexOfLowerTriHessianTerms array.            
            if (pSp->nonZeros() == m_iNbLowerTriHessianNonZeroTerms)
            {
                // if the number of non-zeros match, we don't check the sparsity pattern and cross fingers...
                memcpy(m_pdblIpoptInput[what], pSp->matrixReal->valuePtr(), sizeof(double)*m_iNbLowerTriHessianNonZeroTerms);                                       
            }
            else
            {
                // stability of pattern was not achieved, let's try to fix this
                if (getSparseValuesFromIncompletePattern(m_pdblIpoptInput[what], m_pSpHessianPattern, pSp) == false)
                {
                    snprintf(errorMsg, sizeof(errorMsg), _("%s: Wrong sparsity pattern for output argument #%d.\n"), m_pCallFunctionName[what], 1);
                    throw ast::InternalError(errorMsg);                                            
                }          
            }      
        }
        else
        {
            // User specified lower AND upper triangular parts and of Hessian    
            double *pdblSrc;
            if (pSp->nonZeros() == m_iNbUserHessianNonZeroTerms)
            {
                // if the number of non-zeros match, we don't check the sparsity pattern and cross fingers,
                pdblSrc = pSp->matrixReal->valuePtr();
            }
            else
            {
                // stability of pattern was not achieved, let's try to fix this
                pdblSrc = new double[m_iNbUserHessianNonZeroTerms];
                if (getSparseValuesFromIncompletePattern(pdblSrc, m_pSpHessianPattern, pSp) == false)
                {
                    delete[] pdblSrc;
                    snprintf(errorMsg, sizeof(errorMsg), _("%s: Wrong sparsity pattern for output argument #%d.\n"), m_pCallFunctionName[what], 1);
                    throw ast::InternalError(errorMsg);                                            
                }                 
            }
            // Then we just have to copy the lower triangle terms by using m_pIndexOfLowerTriHessianTerms
            double *pdblDest = m_pdblIpoptInput[what];
            for (int i = 0; i < m_iNbLowerTriHessianNonZeroTerms; i++)
            {
                pdblDest[i] = pdblSrc[m_pIndexOfLowerTriHessianTerms[i]];
            }            
            if (pdblSrc != pSp->matrixReal->valuePtr())
            {
                delete[] pdblSrc;
            }
        }
    }
    else
    {
        types::Double *pDbl =  pIn->getAs<types::Double>();
        if (pIn->isDouble() == false || pIn->getAs<types::Double>()->isComplex() || pIn->getAs<types::Double>()->isEmpty())
        {
            snprintf(errorMsg, sizeof(errorMsg), _("%s: Wrong type for output argument #%d: Real matrix expected.\n"), m_pCallFunctionName[what], 1);
            throw ast::InternalError(errorMsg);            
        }

        int iDims = pDbl->getDims();

        if (m_pDblHessianPattern->isEmpty())
        {
            // Hessian is full, expected returned as such from user function
            if (iDims != 2 || pDbl->getSize() != (m_iNbVariables*m_iNbVariables) || pDbl->getRows() != pDbl->getCols())
            {
                snprintf(errorMsg, sizeof(errorMsg), _("%s: Wrong size for output argument #%d: A matrix of size %d x %d expected.\n"), m_pCallFunctionName[what], 1,m_iNbVariables,m_iNbVariables);           
                throw ast::InternalError(errorMsg);            
            }
            // fast copy of lower triangular part of matrix, column compressed, matching pattern set in OptimizationManager::parseMatrices
            int iOne = 1;
            int iSize = m_iNbVariables;
            double *pdblSrc = pDbl->get();
            double *pdblDest = m_pdblIpoptInput[what];
            
            for (int j=0; j < m_iNbVariables; j++)
            {
                memcpy(pdblDest, pdblSrc, sizeof(double)*iSize);
                pdblSrc += iSize + j+1;
                pdblDest += iSize;
                iSize--;
            }            
        }
        else if (pDbl->getSize() != m_iSizeOfScilabOutput[what])
        {
            snprintf(errorMsg, sizeof(errorMsg), _("%s: Wrong size for output argument #%d: A matrix of size %d expected.\n"), m_pCallFunctionName[what], 1,m_iSizeOfScilabOutput[what]);   
            throw ast::InternalError(errorMsg);
        }
        else // user gave Hessian pattern
        {
            if (m_iNbLowerTriHessianNonZeroTerms == m_iNbUserHessianNonZeroTerms) // User specified lower triangular part of Hessian, hence we copy data as is.
            {
                memcpy(m_pdblIpoptInput[what], pDbl->get(), sizeof(double)*pDbl->getSize());        
            }
            else // user specified both parts of Hessian in Ipopt call, hence we copy only previously identified lower tri terms.
            {
                double *pdblDest = m_pdblIpoptInput[what];
                double *pdblSrc = pDbl->get();
                for (int i = 0; i < m_iNbLowerTriHessianNonZeroTerms; i++)
                {
                    pdblDest[i] = pdblSrc[m_pIndexOfLowerTriHessianTerms[i]];
                }
            }
        }
    }
}

bool OptimizationManager::setSolution(double *x, double *lambda, double *z_L, double *z_U, double fobj)
{
    if (m_pdblSolution != NULL)
    {
        memcpy(m_pdblSolution, x, sizeof(double)*m_iNbVariables);
        memcpy(m_pdblLambda, lambda, sizeof(double)*m_iNbConstraints);
        memcpy(m_pdblZL, z_L, sizeof(double)*m_iNbVariables);
        memcpy(m_pdblZU, z_U, sizeof(double)*m_iNbVariables);
        m_dblFobjSolution = fobj;        
    }
    return true;
}

bool OptimizationManager::intermediateCallback(std::vector<std::pair<std::wstring, double>> scalarQuantities, 
const Ipopt::IpoptData *ip_data, Ipopt::IpoptCalculatedQuantities *ip_cq)
{
    types::typed_list in;
    types::typed_list out;
    types::optional_list opt;
    int iRetCount = 1;
    char errorMsg[256] = "";   
    bool bRet;

    if (m_pCallFunction[INTCB] == NULL)
    {
        return true;
    }
    
    types::Struct *pStructParam = new types::Struct(1,1);

    // get primal and dual variables (if available)
    
    Ipopt::TNLPAdapter* tnlp_adapter = NULL;
    if( ip_cq != NULL )
    {
        Ipopt::OrigIpoptNLP* orignlp;
        orignlp = dynamic_cast<Ipopt::OrigIpoptNLP*>(GetRawPtr(ip_cq->GetIpoptNLP()));
        if( orignlp != NULL)
        {
            tnlp_adapter = dynamic_cast<Ipopt::TNLPAdapter*>(GetRawPtr(orignlp->nlp()));
          
            types::Double *pDblPrimal = new types::Double(m_iNbVariables, 1); 
            types::Double *pDbDualEqs = new types::Double(m_iNbConstraints, 1);
            types::Double *pDbDualLbs = new types::Double(m_iNbVariables, 1);
            types::Double *pDbDualUbs = new types::Double(m_iNbVariables, 1);
          
            tnlp_adapter->ResortX(*ip_data->curr()->x(), pDblPrimal->get());
            tnlp_adapter->ResortG(*ip_data->curr()->y_c(), *ip_data->curr()->y_d(), pDbDualEqs->get());
            tnlp_adapter->ResortBnds(*ip_data->curr()->z_L(), pDbDualLbs->get(),
            *ip_data->curr()->z_U(), pDbDualUbs->get());

            pStructParam->addField(L"x");
            pStructParam->get(0)->set(L"x", pDblPrimal);                               
            pStructParam->addField(L"lambda");
            pStructParam->get(0)->set(L"lambda", pDbDualEqs);                               
            pStructParam->addField(L"zL");
            pStructParam->get(0)->set(L"zL", pDbDualLbs);                               
            pStructParam->addField(L"zU");
            pStructParam->get(0)->set(L"zU", pDbDualUbs);                               
        }
    }
    
    // standard quantities
    
    for (int i=0; i<scalarQuantities.size(); i++)
    {
        pStructParam->addField(scalarQuantities[i].first.c_str());        
        types::Double *pDbl = new types::Double(scalarQuantities[i].second);
        pStructParam->get(0)->set(scalarQuantities[i].first.c_str(), pDbl);           
    }
    in.push_back(pStructParam);
    
    // optional user input parameters 
    for (auto pIn : m_pParameters[INTCB])
    {
        in.push_back(pIn);
    }

    for (auto pIn : in)
    {
        pIn->IncreaseRef();
    }    
    
    try
    {
        // new std::wstring(L"") is deleted in destructor of ast::CommentExp
        m_pCallFunction[INTCB]->invoke(in, opt, iRetCount, out, ast::CommentExp(Location(), new std::wstring(L"")));
    }
    catch(ast::InternalError& ie)
    {
        for (auto pIn : in)
        {
            pIn->DecreaseRef();
            if (pIn->isDeletable())
            {
                delete pIn;
            }
        }    
        throw(ie);
    }
    
    out[0]->IncreaseRef(); // safety IncreaseRef
    
    for (auto pIn : in)
    {
        pIn->DecreaseRef();
        if (pIn->isDeletable())
        {
            delete pIn;
        }
    }    

    if (out.size() != 1)
    {
        snprintf(errorMsg, sizeof(errorMsg), _("%s: Wrong number of output argument(s): %d expected.\n"), m_pCallFunctionName[INTCB], 1);
        throw ast::InternalError(errorMsg);
    }

    if (out[0]->isBool() == false)
    {
        snprintf(errorMsg, sizeof(errorMsg), _("%s: Wrong type for output argument #%d: boolean expected.\n"), m_pCallFunctionName[INTCB], 1);
        throw ast::InternalError(errorMsg);            
    }
    
    bRet = out[0]->getAs<types::Bool>()->get(0);
        
    out[0]->DecreaseRef();
    if (out[0]->isDeletable())
    {
        delete out[0];
    }
    
    return bRet;        
}

bool OptimizationManager::getSparseValuesFromIncompletePattern(double *pdblValues, types::Sparse *pSp0, types::Sparse *pSp1)
{
    // Extract values of pSp1 following pattern of pSp0.
    // Sparsity pattern of pSp1 is supposed to be a subset of the one of pSp0.
    // This function is called if and only if pSp0->nonZeros() != pSp1->nonZeros()
    //
    // Adapted from IPOPT MATLAB-interface provided by
    // Enrico Bertolazzi (enrico.bertolazzi@unitn.it) and Peter Carbonetto
    
    types::Sparse::RealSparse_t* sp0 = pSp0->matrixReal;
    types::Sparse::RealSparse_t* sp1 = pSp1->matrixReal;
    
    memset(pdblValues, 0.0, sizeof(double) * sp0->nonZeros());
    
    int *iCol0 =  sp0->innerIndexPtr(); // column index of values
    int *iRowBegin0 =  sp0->outerIndexPtr(); // index of each row begining in values vector

    int *iCol1 =  sp1->innerIndexPtr(); // column index of values (matrix to be inserted)
    int *iRowBegin1 =  sp1->outerIndexPtr(); // index of each row begining in values vector (matrix to be inserted)

    double *pdbl1 = sp1->valuePtr();
    
    for (int iRow = 0; iRow < sp0->rows(); iRow++)
    {
        for (int iSrc = iRowBegin1[iRow], iDest = iRowBegin0[iRow]; iSrc < iRowBegin1[iRow+1]; iSrc++, iDest++ ) 
        {
            while (iCol0[iDest] < iCol1[iSrc] && iDest < iRowBegin0[iRow+1])
            {
                iDest++;
            }
            if (iCol0[iDest] != iCol1[iSrc])
            {
                // pattern of pSp0 is not a subset of the one of pSp1
                return false;
            }
            pdblValues[iDest] = pdbl1[iSrc];
        }
    }
    return true;
}





























