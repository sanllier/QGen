#ifndef QOBSERVSTATE_H
#define QOBSERVSTATE_H

#include <vector>
#include <string>
#include <cstring>
#include <iostream>

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

class QIndivid;

class QObservState
{
public:
    QObservState(): m_state(0), m_stateSize(0) {}
    QObservState( const QIndivid& ind );
    ~QObservState() 
    {
        delete[] m_state;
    }

    void process( const QIndivid& ind );  // generate bit-vector by input QIndivid
    inline bool at( size_t pos ) const 
    { 
        if ( pos < 0 || pos >= m_stateSize )
            throw std::string( "QObservState out of bounds" ).append( __FUNCTION__ ); 
        return m_state[ pos ]; 
    }

    inline void setBit( size_t pos, bool val )
    {
        if ( pos < 0 || pos >= m_stateSize )
            throw std::string( "QObservState out of bounds" ).append( __FUNCTION__ ); 
        m_state[ pos ] = val;
    }
    
    inline long long size() const { return m_stateSize; }
    inline void clear()
    { 
        delete[] m_state;
        m_state = 0;
    }
    inline bool* data() { return m_state; }
    inline const bool* data() const { return m_state; }

    void requestMemory( long long memSize );

    QObservState& operator=( const QObservState& rState )
    {
        if ( m_stateSize != rState.m_stateSize )
        {
            clear();
            m_state = new bool[ size_t( rState.m_stateSize ) ];
            m_stateSize = rState.m_stateSize;
        }

        std::memcpy( m_state, rState.m_state, size_t( m_stateSize * sizeof( bool ) ) );
        return *this;
    }

    void print( std::ostream& oStr ) const
    {
        for ( int i = 0; i < (int)m_stateSize; ++i )
            oStr << ( m_state[i] ? 1 : 0 );
    }

private:
    bool* m_state;
    long long m_stateSize;
};

//------------------------------------------------------------
}

//------------------------------------------------------------
#endif
