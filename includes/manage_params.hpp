//
// Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
// Copyright (C) 2020 - UTC - Stéphane MOTTELET
//
// This file is hereby licensed under the terms of the GNU GPL v2.0,
// For more information, see the COPYING file which you should have received
//
//--------------------------------------------------------------------------

#ifndef MANAGE_PARAMS_HPP
#define MANAGE_PARAMS_HPP
#include <IpIpoptApplication.hpp>
#include "struct.hxx"

int manage_ipopt_params(types::Struct *pStructOptions, Ipopt::SmartPtr<Ipopt::OptionsList> options, int Log);
#endif
