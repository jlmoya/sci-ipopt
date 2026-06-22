//  Copyright (C) 2008-2010 Yann Collette
//  Copyright (C) 2020 - UTC - Stéphane MOTTELET
//
//  SCIIPOPT is free software; you can redistribute it and/or modify it
//  under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 2, or (at your option) any
//  later version.
//
//  This part of code is distributed with the FURTHER condition that it 
//  can be compiled and linked with the Scilab libraries and it can be 
//  used within the Scilab environment.
//
//  SCIIPOPT is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.

mode(-1);
lines(0);
try
    version = getversion("scilab");
    assert_checktrue(version(1)+version(2)/10 >= 6.1)
catch
    error(gettext("Scilab 6.1 or more is required."));  
end;

// ====================================================================
if ~with_module("development_tools") then
  error(msprintf(gettext("%s module not installed."),"development_tools"));
end
// ====================================================================
TOOLBOX_NAME = "sci_ipopt";
TOOLBOX_TITLE = "SciIPOpt";
// ====================================================================
toolbox_dir = get_absolute_file_path("builder.sce");
// ====================================================================

tbx_builder_macros(toolbox_dir);
tbx_builder_gateway(toolbox_dir);
tbx_builder_help(toolbox_dir);
tbx_build_loader(toolbox_dir);
tbx_build_cleaner(toolbox_dir);

clear toolbox_dir TOOLBOX_NAME TOOLBOX_TITLE;
