#include "qobservstate.h"

#include <cstdlib>

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

QObservState::QObservState( const QIndivid& ind )
{
    process( ind );
}

void QObservState::process( const QIndivid& ind )
{
    if ( m_state.size() < ind.qsize() )
        m_state.resize( ind.qsize(), false );

    float randVal;
    BASETYPE mod;
    for ( size_t i = 0; i < ind.qsize(); ++i )
    {
        randVal = ( ( float )std::rand() / RAND_MAX );
        mod = std::abs( ind.at(i).a );
        m_state[i] = randVal >= mod * mod;  // if randVal < |a|^2 then '0', otherwise '1'
    }
}

//------------------------------------------------------------
}
