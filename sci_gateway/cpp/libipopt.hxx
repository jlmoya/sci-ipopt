#ifndef __LIBIPOPT_GW_HXX__
#define __LIBIPOPT_GW_HXX__

#ifdef _MSC_VER
#ifdef LIBIPOPT_GW_EXPORTS
#define LIBIPOPT_GW_IMPEXP __declspec(dllexport)
#else
#define LIBIPOPT_GW_IMPEXP __declspec(dllimport)
#endif
#else
#define LIBIPOPT_GW_IMPEXP
#endif

extern "C" LIBIPOPT_GW_IMPEXP int libipopt(wchar_t* _pwstFuncName);

CPP_GATEWAY_PROTOTYPE(sci_ipopt);

#endif /* __LIBIPOPT_GW_HXX__ */
