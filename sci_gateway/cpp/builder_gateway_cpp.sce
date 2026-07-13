//
// Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
// Copyright (C) 2020-2023 - Stephane MOTTELET - Université de Technologie de Compiègne
//
// This file is hereby licensed under the terms of the GNU GPL v2.0,
// For more information, see the COPYING file which you should have received
//
//--------------------------------------------------------------------------

path_builder = get_absolute_file_path('builder_gateway_cpp.sce');
path_root = fullpath(fullfile(path_builder,'..','..'));
cd(path_builder);

files_in_src = ['parameters.cpp'
'OptimizationManager.cpp'
'IpoptTNLP.cpp'
'manage_ipopt_params.cpp'
'scilabjournal.cpp'
'scilabexception.cpp'];

for i=1:size(files_in_src,"*")
    copyfile(fullfile("..","..","src","cpp",files_in_src(i)),files_in_src(i));
end

files_to_compile = ['sci_ipopt.cpp'
'parameters.cpp'
'OptimizationManager.cpp'
'IpoptTNLP.cpp'
'manage_ipopt_params.cpp'
'scilabjournal.cpp'
'scilabexception.cpp'];

if getos() == "Windows"
    tools_path  = fullfile(path_root,'thirdparty',getos());
elseif getos() == "Darwin"
    // macOS arm64 port: use the Homebrew-installed IPOPT instead of the bundled thirdparty/
    tools_path  = "/opt/homebrew/opt/ipopt";
else
    [v,opt]=getversion();
    arch = opt(2);
    if arch == "x64"
        arch = "x86_64"
    end
    tools_path  = fullfile(path_root,'thirdparty',getos(),arch);
end
include_ipopt = fullfile(tools_path,'include','coin-or');  
cflags = ' -I' + fullfile(path_root,'includes');
cflags = cflags + ' -DUSE_MUMPS'
cflags = cflags + ' -I' + include_ipopt;
if getos() == "Linux" then
    cflags = cflags + ' -I' + fullfile(SCI,'lib','Eigen','include');
elseif getos() == "Darwin"
    // macOS arm64 port: Eigen from Homebrew (SCI/lib/Eigen/includes is absent in a source build)
    cflags = cflags + ' -I/opt/homebrew/opt/eigen/include/eigen3';
    // macOS has no libintl.h on the default include path; Scilab's localization.h needs it.
    cflags = cflags + ' -I/opt/homebrew/opt/gettext/include';
    conda = getenv("CONDA_PREFIX","UNDEFINED");
    if conda <> "UNDEFINED"
        cflags = cflags + " -I" + fullfile(conda,"include");
        cflags = cflags + " -I" + fullfile(conda,"include","eigen3");
    end
elseif getos() == "Windows"
    cflags = cflags + ' -I' + fullfile(SCI,'libs','Eigen','includes');
    ipoptlib = fullfile(tools_path,"lib","ipopt.lib");
    coinmumpslib = fullfile(tools_path,"lib","coinmumps.lib");
    if ~isfile(ipoptlib)
        movefile(fullfile(tools_path,"lib","ipopt.dll.lib"),ipoptlib);
    end
    if ~isfile(coinmumpslib)
        movefile(fullfile(tools_path,"lib","coinmumps.dll.lib"),coinmumpslib);
    end
end

libs = [];

if getos()=="Windows" then
    lib_base_dir = fullfile("..\..\thirdparty\Windows","lib");
    libs = [fullfile(lib_base_dir,"coinmumps")
            fullfile(lib_base_dir,"ipopt")];
    ldflags = ""
elseif getos() == "Darwin"
    // Homebrew IPOPT links MUMPS (libdmumps) transitively — no separate -lcoinmumps
    ldflags  = '-L' + fullfile(tools_path,'lib') + ' -lipopt';
    // Homebrew IPOPT pulls the gcc runtime (-lemutls_w, -lgfortran) transitively; the
    // compiler's built-in gcc lib search path goes stale after a gcc bump, so add the
    // current gcc lib dir explicitly (version-independent 'current' symlink).
    gcc_lib = stripblanks(unix_g("dirname $(ls /opt/homebrew/lib/gcc/current/gcc/*/*/libemutls_w.a 2>/dev/null | head -1)"));
    if gcc_lib <> "" then ldflags = ldflags + " -L" + gcc_lib; end
    if isdir("/opt/homebrew/lib/gcc/current") then ldflags = ldflags + " -L/opt/homebrew/lib/gcc/current"; end
else
    ldflags  = '-L' + fullfile(tools_path,'lib') + ' -lipopt -lcoinmumps';
end
ilib_verbose(2)

tbx_build_gateway('ipopt',['ipopt', 'sci_ipopt','cppsci'],files_to_compile, path_builder, libs, ldflags, cflags);

if getos() == "Windows" then
    loader = mgetl("loader.sce");
    loader = strsubst(loader,libs(1),libs(1)+"-3")
    loader = strsubst(loader,libs(2),libs(2)+"-3")
    loader = strsubst(loader,lib_base_dir,"..\..\thirdparty\Windows\bin")
    mputl(loader,"loader.sce")
end

