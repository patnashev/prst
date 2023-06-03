#pragma once

//#define BOINC
//#define NETPRST

#define PRST_VERSION "9.0"
#define VERSION_BUILD "704"

inline void print_banner()
{
    printf("PRST version " PRST_VERSION "." VERSION_BUILD ", GWnum library version " GWNUM_VERSION);
#ifdef GMP
    arithmetic::GMPArithmetic* gmp = dynamic_cast<arithmetic::GMPArithmetic*>(&arithmetic::GiantsArithmetic::default_arithmetic());
    if (gmp != nullptr)
        printf(", GMP library version %s", gmp->version().data());
#endif
    printf("\n");
}
