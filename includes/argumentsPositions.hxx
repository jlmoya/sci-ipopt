//
// Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
// Copyright (C) 2020 - UTC - Stephane MOTTELET
//
// This file is hereby licensed under the terms of the GNU GPL v2.0,
// For more information, see the COPYING file which you should have received
//
//--------------------------------------------------------------------------

#ifndef __ARGUMENTSPOSITIONS_HXX__
#define __ARGUMENTSPOSITIONS_HXX__

#define X_POS                   0
#define FOBJ_POS                1
#define DFOBJ_POS               2
#define GOBJ_POS                3
#define DGOBJ_POS               4
#define SPARSE_DGOBJ_POS        5
#define DHOBJ_POS               6
#define SPARSE_DHOBJ_POS        7
#define VAR_LIN_TYPE_POS        8
#define CONSTR_LIN_TYPE_POS     9
#define CONSTR_RHS_POS          10
#define CONSTR_LHS_POS          11
#define X_LOWER_POS             12
#define X_UPPER_POS             13
#define INTCB_POS               14
#define PARAMS_POS              15

#define IPOPT_ERROR \
    if (bRet == false)\
    {\
        delete manager;\
        return types::Function::Error;\
    }
     
#endif