#ifndef QOBSERVSTATE_H
#define QOBSERVSTATE_H

#include <iostream>
#include "mpi.h"

#if defined( GPU ) && defined( CURAND )
    #include "random_curand.h"
#elseif defined( STDRAND )
    #include "random_def.h"
#else 
    #include "random_mtrand.h"
#endif

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

class QBaseIndivid;

//------------------------------------------------------------

class QObserveState
{
public:
    QObserveState( unsigned seed );
    ~QObserveState(); 

    void observe( const QBaseIndivid& ind );  // generate observe-vector by input individ

    void clear();

    inline long long size() const { return m_stateSize; }
    bool at( long long pos ) const;
    void set( long long pos, bool val );

#ifdef GPU
    const bool* gpuBuffer() const { return m_gpuBuf; }
#endif
  
    void bcast( int root, MPI_Comm comm );

    QObserveState& operator=( const QObserveState& rState );
    
private:
    void resize( long long size );

private:
    bool* m_state;
    long long m_stateSize;

    IRandom* m_rand;

#ifdef GPU
    bool* m_gpuBuf;
#endif
};

//------------------------------------------------------------
}

//------------------------------------------------------------
#endif
