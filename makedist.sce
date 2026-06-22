// =============================================================================
// Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
//
// Copyright (C) 2024 - UTC - Stéphane MOTTELET
//
// This file is hereby licensed under the terms of the GNU GPL v3.0,
// For more information, see the COPYING file which you should have received

path_makedist = get_absolute_file_path("makedist.sce");
cd(path_makedist);

modulename = "sci_ipopt"
version = mgetl("VERSION")
http_get(msprintf("%s/%s/%s/DESCRIPTION","https://atoms.scilab.org/toolboxes",modulename,version),"DESCRIPTION");

[sci,v]=getversion()
if getos() <> "Windows" then
    filename = modulename+"_"+mgetl("VERSION")+".bin."+v(2)+"."+getos()+".tar.gz"
    cd ..
    compress(filename,"sci-ipopt",format="tar",compression="gzip");
else
    filename =  modulename+"_"+mgetl("VERSION")+".bin."+v(2)+"."+getos()+".zip"
    cd ..
    compress(filename,"sci-ipopt",format="zip");
end
cd(path_makedist);
