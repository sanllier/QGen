#ifndef QBIT_H
#define QBIT_H

#include "defs.h"
#include <complex>

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

typedef std::complex< BASETYPE > QComplex;

//------------------------------------------------------------

struct QBit
{
    QComplex a;
    QComplex b;

    QBit() {}
    QBit( QComplex a, QComplex b )
        : a(a)
        , b(b)
    {}
};

//------------------------------------------------------------
}

//------------------------------------------------------------
#endif
