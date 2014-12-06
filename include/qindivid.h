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
    virtual BASETYPE operator()( const QObservState& observState ) = 0;
};

class QRepairClass
{
public:
    virtual ~QRepairClass() {}
    virtual void operator()( QObservState& observState ) = 0;
};

//------------------------------------------------------------

class QIndivid
{
public:
    QIndivid( size_t size = 0 );
    QIndivid( const QIndivid& ind );
    virtual ~QIndivid();

    void resize( size_t size );
    void setInitial();
    
    inline const QBit& at( size_t pos ) const 
    { 
        if ( pos < 0 || pos >= m_dataLogicSize ) 
            throw std::string( "QIndivid out of bounds" ).append( __FUNCTION__ ); 
        return m_data[ pos ]; 
    }
    inline QBit& at( size_t pos ) 
    { 
        if ( pos < 0 || pos >= m_dataLogicSize ) 
            throw std::string( "QIndivid out of bounds" ).append( __FUNCTION__ ); 

        m_needRecalcFitness = true;
        return m_data[ pos ]; 
    }

    inline size_t qsize() const { return m_dataLogicSize; }

    inline const QBit* raw() const { return m_data; }
    inline QBit* raw() { return m_data; }

    inline BASETYPE getFitness( QFitnessClass* fClass = 0 )
    {
        if ( m_needRecalcFitness && !fClass )
            throw std::string( "QIndivid is trying to recalc fitness with (NULL) func" ).append( __FUNCTION__ );  

        if ( m_needRecalcFitness )
        {
            m_fitness = (*fClass)( m_obsState );
            m_needRecalcFitness = false;
        }

        return m_fitness;
    }
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

    void bcast( int root, MPI_Comm comm );

private:
    QBit* m_data;
    size_t m_dataLogicSize;

    QObservState m_obsState;

    BASETYPE m_fitness;
    bool m_needRecalcFitness;
};

//------------------------------------------------------------
}

//------------------------------------------------------------
#endif
