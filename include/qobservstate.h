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
    
    inline size_t size() const { return m_stateSize; }
    inline void clear()
    { 
        delete[] m_state;
        m_state = 0;
    }
    inline bool* data() { return m_state; }

    void requestMemory( size_t memSize );

    QObservState& operator=( const QObservState& rState )
    {
        if ( m_stateSize != rState.m_stateSize )
        {
            clear();
            m_state = new bool[ rState.m_stateSize ];
            m_stateSize = rState.m_stateSize;
        }

        std::memcpy( m_state, rState.m_state, m_stateSize * sizeof( bool ) );
        return *this;
    }

    void print()
    {
        for ( int i = 0; i < (int)m_stateSize; ++i )
            std::cout << (m_state[i] ? 1 : 0) << " ";
        std::cout << "\r\n";
    }

private:
    bool* m_state;
    size_t m_stateSize;
};

//------------------------------------------------------------
}

//------------------------------------------------------------
#endif
