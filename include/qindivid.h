#ifndef QINDIVID_H
#define QINDIVID_H

#include <complex>

#include "mpi.h"

#include "qrotoperator.h"
#include "qobservstate.h"
#include "defs.h"

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

class QFitnessClass
{
public:
    virtual ~QFitnessClass() {}
    virtual BASETYPE operator()( MPI_Comm indComm, const QObservState& observState, int startQBit, int idx ) = 0;
};

class QRepairClass
{
public:
    virtual ~QRepairClass() {}
    virtual void operator()( MPI_Comm indComm, const QObservState& observState, int startQBit, int idx ) = 0;
};

//------------------------------------------------------------

class QIndivid
{
public:
    QIndivid( long long size, MPI_Comm generalComm, MPI_Comm rowComm, int coords[2] );
    virtual ~QIndivid();

    void resize( long long newSize );
    void setInitial();
    
    inline const QBit& localAt( size_t pos ) const 
    { 
        if ( pos < 0 || pos >= m_localLogicSize ) 
            throw std::string( "QIndivid out of bounds" ).append( __FUNCTION__ ); 
        return m_data[ pos ]; 
    }
    inline QBit& localAt( size_t pos ) 
    { 
        if ( pos < 0 || pos >= m_localLogicSize ) 
            throw std::string( "QIndivid out of bounds" ).append( __FUNCTION__ ); 

        m_needRecalcFitness = true;
        return m_data[ pos ]; 
    }

    inline long long qSize() const { return m_globalLogicSize; }
    inline long long locaQSize() const { return m_localLogicSize; }
    inline long long startQBit() const { return m_startQbit; }

    inline const QBit* raw() const { return m_data; }
    inline QBit* raw() { return m_data; }

    BASETYPE getFitness( QFitnessClass* fClass = 0 );
    inline BASETYPE getFitnessUNSAFE() const { return m_fitness; } // unsafe, m_fitness may be dirty

    inline QObservState& getObsState() 
    {
        m_needRecalcFitness = true;
        return m_obsState;
    }
    inline const QObservState& getObsState() const { return m_obsState; }

    QIndivid& operator=( const QIndivid& ind );

    void tick( const QIndivid& bestInd );
    void repair( QRepairClass* repClass );
    inline void updateObsState() 
    { 
        m_obsState.process( *this ); 
        m_needRecalcFitness = true;
    }

    void bcast( int root );

private:
    QBit* m_data;
    long long m_localLogicSize;
    long long m_globalLogicSize;
    long long m_startQbit;

    QObservState m_obsState;

    BASETYPE m_fitness;
    bool m_needRecalcFitness;

    //--------------------------

    MPI_Comm m_indComm;
    MPI_Comm m_rowComm;
    int m_coords[2];
};

//------------------------------------------------------------
}

//------------------------------------------------------------
#endif
