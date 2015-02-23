#ifndef INTERFACES_H
#define INTERFACES_H

#include "qindivid_base.h"
#include "qobservestate.h"

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------
   
class IProcessScreen
{
public:
    virtual ~IProcessScreen() {}
    virtual void operator()( long long cycle, 
                             const int coords[2], 
                             const QBaseIndivid& totalBest, 
                             const QBaseIndivid& iterBest ) = 0;
};

//------------------------------------------------------------

class IFitness
{
public:
    virtual ~IFitness() {}
    virtual BASETYPE operator()( MPI_Comm indComm, const QObserveState& observeState, long long startQBit, int idx ) = 0;
};

//------------------------------------------------------------

class IRepair
{
public:
    virtual ~IRepair() {}
    virtual void operator()( MPI_Comm indComm, QObserveState& observeState, long long startQBit, int idx ) = 0;
};

//-----------------------------------------------------------
}

//------------------------------------------------------------
#endif
