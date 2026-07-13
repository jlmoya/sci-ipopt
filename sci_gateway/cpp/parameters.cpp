/*
 *  Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
 *  Copyright (C) 2020-2023 - UTC - Stéphane MOTTELET
 *
 * This file is hereby licensed under the terms of the GNU GPL v2.0,
 * For more information, see the COPYING file which you should have received
 * along with this program.
 *
/*--------------------------------------------------------------------------*/

#include "function.hxx"
#include "double.hxx"
#include "string.hxx"
#include "struct.hxx"
#include "parameters.hxx"

extern "C"
{
#include "localization.h"
#include "Scierror.h"
#include "sciprint.h"
}

bool getDoubleInPList(types::Struct *pStruct, const wchar_t * _pwstLabel, double * _pdblValue,
                        double _dblDefaultValue, bool _bLog, enum type_check _eCheck, ...)
{
    char message[256];
    bool bRet = false;
    
    if (pStruct->exists(_pwstLabel) == true)
    {
        types::InternalType *pI;
        bRet = true;
        pStruct->extract(_pwstLabel, pI);
        if (pI->isDouble() && pI->getAs<types::Double>()->isEmpty()==false)
        {
            *_pdblValue = pI->getAs<types::Double>()->get(0);
            pStruct->removeField(_pwstLabel);
        }
        else
        {
            snprintf(message, sizeof(message),_("%s: wrong value type for parameter \"%ls\": %s expected.\n"), "ipopt",_pwstLabel, "double");
            throw ast::InternalError(message);
        }
    }
    else
    {
        *_pdblValue = _dblDefaultValue;        
    }

/* Now check parameters */

    if (_eCheck != CHECK_NONE)
    {
        va_list vl;
        int nb_value_to_check = 0;
        double value_to_check = 0;
        int check_res = 0;

        va_start(vl, _eCheck);

        switch (_eCheck)
        {
            case CHECK_MIN:
                    value_to_check = va_arg(vl, double);
                    va_end(vl);
                    if (value_to_check > *_pdblValue)
                    {
                        snprintf(message, sizeof(message),_("%s: wrong value %f for parameter %ls: lower bound is %f.\n"), "ipopt", *_pdblValue, _pwstLabel, value_to_check);
                        throw ast::InternalError(message);
                    }
                    break;
            case CHECK_MAX:
                    value_to_check = va_arg(vl, double);
                    va_end(vl);
                    if (value_to_check < *_pdblValue)
                    {
                        snprintf(message, sizeof(message),_("%s: wrong value %f for parameter %ls: upper bound is %f.\n"), "ipopt", *_pdblValue, _pwstLabel, value_to_check);
                        throw ast::InternalError(message);
                    }
                    break;
            case CHECK_BOTH:
                    /* First value is the min bound */
                    value_to_check = va_arg(vl, double);
                    if (value_to_check > *_pdblValue)
                    {
                        snprintf(message, sizeof(message),_("%s: wrong value %f for parameter %ls: lower bound is %f.\n"), "ipopt", *_pdblValue, _pwstLabel, value_to_check);
                        throw ast::InternalError(message);
                    }
                    /* Second value is the max bound */
                    value_to_check = va_arg(vl, double);
                    va_end(vl);
                    if (value_to_check < *_pdblValue)
                    {
                        snprintf(message, sizeof(message),_("%s: wrong value %f for parameter %ls: upper bound is %f.\n"), "ipopt", *_pdblValue, _pwstLabel, value_to_check);
                        throw ast::InternalError(message);
                    }
                    break;
            // case CHECK_VALUES:
            //         /* First parameters is int and contains the number of values to check */
            //         nb_value_to_check = va_arg(vl, int);
            //         check_res = 0;
            //         for (int i = 0; i < nb_value_to_check; i++)
            //         {
            //             value_to_check = va_arg(vl, double);
            //             check_res = check_res || (value_to_check == *_pdblValue);
            //         }
            //
            //         if (!check_res)
            //         {
            //             if ( (_bLog))
            //             {
            //                 sciprint(_("%s: wrong value for parameter %ls: value %f\n"), "getDoubleInPList", _pwstLabel, *_pdblValue);
            //                 sciprint(_("%s: awaited parameters: "), "getDoubleInPList");
            //                 va_start(vl, _eCheck);
            //                 nb_value_to_check = va_arg(vl, int);
            //                 for (int i = 0; i < nb_value_to_check; i++)
            //                 {
            //                     value_to_check = va_arg(vl, double);
            //                     sciprint(" %f", value_to_check);
            //                 }
            //                 sciprint("\n");
            //             }
            //
            //             *_pdblValue = _dblDefaultValue;
            //
            //             va_end(vl);
            //              Scierror(999, _("%s: wrong value for parameter %ls: value %f\n"), "getDoubleInPList", _pwstLabel, *_pdblValue);
            //             return types::Function::Error;
            //         }
            //         va_end(vl);
            //         break;
                    default:
                va_end(vl);
        }
    }
    return bRet;    
}

bool getIntInPList(types::Struct *pStruct, const wchar_t * _pwstLabel, int * _piValue,
                        int _iDefaultValue, bool _bLog, enum type_check _eCheck, ...)
{
    bool bRet = false;
    char message[256];
    
    if (pStruct->exists(_pwstLabel) == true)
    {
        types::InternalType *pI;
        bRet = true;
        pStruct->extract(_pwstLabel, pI);
        if (pI->isDouble())
        {
            *_piValue = (int)pI->getAs<types::Double>()->get(0);
            pStruct->removeField(_pwstLabel);
        }
        else
        {
            snprintf(message, sizeof(message),_("%s: wrong value type for parameter \"%ls\": %s expected.\n"), "ipopt",_pwstLabel, "double");
            throw ast::InternalError(message);
        }
    }
    else
    {
        *_piValue = _iDefaultValue;        
    }

/* Now check parameters */

    if (_eCheck != CHECK_NONE)
    {
        va_list vl;
        int nb_value_to_check = 0;
        int value_to_check = 0;
        int check_res = 0;

        va_start(vl, _eCheck);

        switch (_eCheck)
        {
            case CHECK_MIN:
                value_to_check = va_arg(vl, int);
                va_end(vl);
                if (value_to_check > *_piValue)
                {
                    snprintf(message, sizeof(message),_("%s: wrong value %d for parameter %ls: lower bound is %d.\n"), "ipopt", *_piValue, _pwstLabel, value_to_check);
                    throw ast::InternalError(message);
                }
                break;
            case CHECK_MAX:
                value_to_check = va_arg(vl, int);
                va_end(vl);
                if (value_to_check < *_piValue)
                {
                    snprintf(message, sizeof(message),_("%s: wrong value %d for parameter %ls: upper bound is %d.\n"), "ipopt", *_piValue, _pwstLabel, value_to_check);
                    throw ast::InternalError(message);
                }
                break;
            case CHECK_BOTH:
                /* First value is the min bound */
                value_to_check = va_arg(vl, int);
                if (value_to_check > *_piValue)
                {
                    snprintf(message, sizeof(message),_("%s: wrong value %d for parameter %ls: lower bound is %d.\n"), "ipopt", *_piValue, _pwstLabel, value_to_check);
                    throw ast::InternalError(message);
                }
                /* Second value is the max bound */
                value_to_check = va_arg(vl, int);
                va_end(vl);
                if (value_to_check < *_piValue)
                {
                    snprintf(message, sizeof(message),_("%s: wrong value %d for parameter %ls: upper bound is %d.\n"), "ipopt", *_piValue, _pwstLabel, value_to_check);
                    throw ast::InternalError(message);
                }
                break;
            // case CHECK_VALUES:
            //     /* First parameters is int and contains the number of values to check */
            //     nb_value_to_check = va_arg(vl, int);
            //     check_res = 0;
            //     for (int i = 0; i < nb_value_to_check; i++)
            //     {
            //         value_to_check = va_arg(vl, int);
            //         check_res = check_res || (value_to_check == *_piValue);
            //     }
            //
            //     if (!check_res)
            //     {
            //         if ( (_bLog))
            //         {
            //             sciprint(_("%s: wrong value for parameter %ls: value %d\n"), "getIntInPList", _pwstLabel, *_piValue);
            //             sciprint(_("%s: awaited parameters: "), "getIntInPList");
            //             va_start(vl, _eCheck);
            //             nb_value_to_check = va_arg(vl, int);
            //             for (int i = 0; i < nb_value_to_check; i++)
            //             {
            //                 value_to_check = va_arg(vl, int);
            //                 sciprint(" %d", value_to_check);
            //             }
            //             sciprint("\n");
            //         }
            //
            //         *_piValue = _iDefaultValue;
            //
            //         va_end(vl);
            //          Scierror(999, _("%s: wrong value for parameter %ls: value %d\n"), "getIntInPList", _pwstLabel, *_piValue);
            //         return types::Function::Error;
            //     }
            //     va_end(vl);
            //     break;
            default:
                va_end(vl);
        }
    }

    return bRet;    
}

bool getStringInPList(types::Struct *pStruct, const wchar_t * _pwstLabel, char ** _pstValue,
                        const char * _pstDefaultValue, bool _bLog, enum type_check _eCheck, ...)
{
    bool bRet = false;
    char message[256];

    if (pStruct->exists(_pwstLabel) == true)
    {
        types::InternalType *pI;
        bRet = true;
        pStruct->extract(_pwstLabel, pI);
        if (pI->isString())
        {
            *_pstValue = wide_string_to_UTF8(pI->getAs<types::String>()->get(0));            
            pStruct->removeField(_pwstLabel);
        }
        else
        {
            snprintf(message, sizeof(message),_("%s: wrong value type for parameter \"%ls\": %s expected.\n"), "ipopt",_pwstLabel, "string");
            throw ast::InternalError(message);
        }
    }
    else
    {
        *_pstValue = os_strdup(_pstDefaultValue);
    }

    /* Now check parameters */

    if (_eCheck != CHECK_NONE)
    {
        va_list vl;
        int nb_value_to_check = 0;
        char * value_to_check = 0;
        int check_res = 0;

        va_start(vl, _eCheck);

        switch (_eCheck)
        {
            case CHECK_VALUES:
                /* First parameters is int and contains the number of values to check */
                nb_value_to_check = va_arg(vl, int);
                check_res = 0;
                for (int i = 0; i < nb_value_to_check; i++)
                {
                    value_to_check = va_arg(vl, char *);
                    check_res = check_res || (strcmp(value_to_check, *_pstValue) == 0);
                }

                if (!check_res)
                {
                    std::ostringstream os;
                    snprintf(message, sizeof(message),_("%s: wrong value \"%s\" for parameter %ls. Expected values are "), "ipopt", *_pstValue, _pwstLabel );
                    os << message;
                    va_end(vl);
                    va_start(vl, _eCheck);
                    nb_value_to_check = va_arg(vl, int);
                    for (int i = 0; i < nb_value_to_check; i++)
                    {
                        value_to_check = va_arg(vl, char *);
                        os << "\"" << value_to_check << "\"";
                        if (i < nb_value_to_check-1)
                        {
                            os << ", ";                            
                        }
                    }
                    va_end(vl);

                    if (*_pstValue)
                    {
                        FREE(*_pstValue);
                        *_pstValue = NULL;
                    }
                    throw ast::InternalError(os.str());
                }

                va_end(vl);
                break;
            default:
                va_end(vl);
        }
    }

    return bRet; 
}