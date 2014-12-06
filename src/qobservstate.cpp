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
    if ( m_stateSize != ind.qsize() )
    {
        clear();
        m_state = new bool[ ind.qsize() ];
        m_stateSize = ind.qsize();
    }

    for ( size_t i = 0; i < ind.qsize(); ++i )
    {
        const BASETYPE randVal = (BASETYPE)(rand()) / RAND_MAX;//( BASETYPE )( SharedMTRand::getClosedInstance()() );
        const BASETYPE mod = std::abs( ind.at(i).a );
        m_state[i] = randVal >= mod * mod;  // if randVal < |a|^2 then '0', otherwise '1'
    }
}

//------------------------------------------------------------
}
