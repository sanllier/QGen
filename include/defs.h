#ifndef DEFS_H
#define DEFS_H
//------------------------------------------------------------

#define BASETYPE float
#define MPI_BASETYPE MPI_FLOAT

//------------------------------------------------------------

typedef std::complex< BASETYPE > QComplex;

//-----------------------------------------------------------

enum
{
    ROOT_ID = 0
};

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
#endif
