#include "context.hxx"
#include "cpp_gateway_prototype.hxx"
#include "libipopt.hxx"
extern "C"
{
#include "libipopt.h"
}

#define MODULE_NAME L"libipopt"

int libipopt(wchar_t* _pwstFuncName)
{
    if(wcscmp(_pwstFuncName, L"ipopt") == 0){ symbol::Context::getInstance()->addFunction(types::Function::createFunction(L"ipopt", &sci_ipopt, MODULE_NAME)); }

    return 1;
}
