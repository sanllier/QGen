#ifdef GPU

#include "qindivid_gpu.h"
#include "qgen.h"
#include "qrotoperator.h"

#include "sharedmtrand.h"
#include "mpicheck.h"

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

QGPUIndivid::QGPUIndivid( long long size, MPI_Comm generalComm, MPI_Comm rowComm, int coords[2] )
    : QBaseIndivid( size, generalComm,  rowComm, coords )
{

}

//------------------------------------------------------------

QGPUIndivid::~QGPUIndivid()
{

}

//------------------------------------------------------------

void QGPUIndivid::resize( long long newSize )
{

}

//------------------------------------------------------------

void QGPUIndivid::setInitial()
{

}

//------------------------------------------------------------

void QGPUIndivid::evolve( const QBaseIndivid& bestInd )
{

}

//------------------------------------------------------------

void QGPUIndivid::bcast( int root )
{

}

//-----------------------------------------------------------

QBaseIndivid& QGPUIndivid::operator=( const QBaseIndivid& rInd )
{
    return *this;
}

//-----------------------------------------------------------

BASETYPE QGPUIndivid::getThetaForQBit( const QBaseIndivid& bestInd, long long qbitIndex ) const
{
    return 0.0;
}

//------------------------------------------------------------
}

#endif
