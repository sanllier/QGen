#include "qobservstate.h"
#include "qindivid.h"
#include "sharedmtrand.h"

#include <cstdlib>

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

QObservState::QObservState( const QIndivid& ind ): m_state(0), m_stateSize(0)
{
    process( ind );
}

void QObservState::process( const QIndivid& ind )
{
    requestMemory( ind.qsize() );

    for ( size_t i = 0; i < ind.qsize(); ++i )
    {
        //const BASETYPE randVal = (BASETYPE)(rand()) / RAND_MAX;
        const BASETYPE randVal = ( BASETYPE )( SharedMTRand::getClosedInstance()() );
        const BASETYPE mod = std::abs( ind.at(i).a );
        m_state[i] = randVal >= mod * mod;  // if randVal < |a|^2 then '0', otherwise '1'
    }
}

void QObservState::requestMemory( size_t memSize )
{
    if ( m_stateSize != memSize )
    {
        clear();
        m_state = new bool[ memSize ];
        m_stateSize = memSize;
    }
}

//------------------------------------------------------------
}
