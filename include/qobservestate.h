#ifndef QOBSERVSTATE_H
#define QOBSERVSTATE_H

#include <iostream>

#include "mpi.h"

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

class QBaseIndivid;

//------------------------------------------------------------

class QObserveState
{
public:
    QObserveState();
    ~QObserveState(); 

    void observe( const QBaseIndivid& ind );  // generate observe-vector by input individ

    void clear();

    inline long long size() const { return m_stateSize; }
    bool at( long long pos ) const;
    void set( long long pos, bool val );
  
    void bcast( int root, MPI_Comm comm );

    QObserveState& operator=( const QObserveState& rState );

    const bool* gpuBuffer() const { return m_gpuBuf; }

private:
    void resize( long long size );

private:
    bool* m_state;
    long long m_stateSize;

#ifdef GPU
    bool* m_gpuBuf;
#endif
};

//------------------------------------------------------------
}

//------------------------------------------------------------
#endif
