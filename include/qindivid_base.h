#ifndef QINDIVID_BASE_H
#define QINDIVID_BASE_H

#include <complex>

#include "mpi.h"
#include "qobservestate.h"
#include "defs.h"

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

class QFitnessClass
{
public:
    virtual ~QFitnessClass() {}
    virtual BASETYPE operator()( MPI_Comm indComm, const QObserveState& observeState, long long startQBit, int idx ) = 0;
};

//------------------------------------------------------------

class QRepairClass
{
public:
    virtual ~QRepairClass() {}
    virtual void operator()( MPI_Comm indComm, QObserveState& observeState, long long startQBit, int idx ) = 0;
};

//------------------------------------------------------------

enum EIndividType
{
    INDIVID_TYPE_CPU,
    INDIVID_TYPE_GPU
};

//------------------------------------------------------------

class QBaseIndivid
{
public:
    QBaseIndivid( long long size, MPI_Comm generalComm, MPI_Comm rowComm, int coords[2] );
    virtual ~QBaseIndivid();

    virtual EIndividType getType() const = 0;

    virtual void resize( long long newSize );
    virtual void setInitial() = 0;
    
    inline long long qSize() const { return m_globalLogicSize; }
    inline long long locaQSize() const { return m_localLogicSize; }
    inline long long firstQBit() const { return m_firstQbit; }

    BASETYPE calculateFitness( QFitnessClass* fClass );
    BASETYPE getFitness() const { return m_fitness; } // unsafe, m_fitness may be dirty

    void calculateObservState();
    QObserveState& getObservState();
    const QObserveState& getObservState() const;

    virtual void evolve( const QBaseIndivid& bestInd ) = 0;
    void repair( QRepairClass* repClass );

    virtual void bcast( int root );

    virtual QBaseIndivid& operator=( const QBaseIndivid& rInd ) = 0;

protected:
    virtual BASETYPE getThetaForQBit( const QBaseIndivid& bestInd, long long qbitIndex ) const = 0;

protected:
    static MPI_Datatype MPI_QBIT;

    long long m_localLogicSize;
    long long m_globalLogicSize;
    long long m_firstQbit;

    QObserveState m_observeState;

    BASETYPE m_fitness;
    bool m_needRecalcFitness;

    struct SBaseIndividContext
    {
        MPI_Comm indComm;
        MPI_Comm rowComm;
        int coords[2];

        SBaseIndividContext()
            : indComm( MPI_COMM_NULL )
            , rowComm( MPI_COMM_NULL )
        {}
    };
    SBaseIndividContext m_context;
};

//------------------------------------------------------------
}

//------------------------------------------------------------
#endif
