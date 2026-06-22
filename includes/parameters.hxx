//
// Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
// Copyright (C) 2020 - UTC - Stéphane MOTTELET
//
// This file is hereby licensed under the terms of the GNU GPL v2.0,
// For more information, see the COPYING file which you should have received
//
//--------------------------------------------------------------------------

enum type_check
{
    CHECK_NONE,
    CHECK_SIZE,
    CHECK_MIN,
    CHECK_MAX,
    CHECK_BOTH,
    CHECK_VALUES
};

bool getDoubleInPList(types::Struct *, const wchar_t * _pwstLabel, double * _pdblValue,
                        double _dblDefaultValue, bool _bLog, enum type_check _eCheck, ...);


bool getIntInPList(types::Struct *, const wchar_t * _pwstLabel, int * _piValue,
                        int _iDefaultValue, bool _bLog, enum type_check _eCheck, ...);


bool getStringInPList(types::Struct *, const wchar_t * _pwstLabel, char ** _pstValue,
                        const char * _pstDefaultValue, bool _bLog, enum type_check _eCheck, ...);
