#ifndef QOBSERVSTATE_H
#define QOBSERVSTATE_H

#include <vector>
#include <string>

#include "qindivid.h"

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

class QObservState
{
public:
    QObservState() {}
    QObservState( const QIndivid& ind );
    ~QObservState() {}

    void process( const QIndivid& ind );
    inline bool at( size_t pos ) const { if ( pos < 0 || pos >= m_state.size() ) \
        throw std::string( "QObservState out of bounds" ).append( __FUNCTION__ ); return m_state[ pos ]; }
    
    inline size_t size() const { return m_state.size(); }
    inline void clear() { m_state.clear(); }

private:
    std::vector< bool > m_state;
};

//------------------------------------------------------------
}

//------------------------------------------------------------
#endif
