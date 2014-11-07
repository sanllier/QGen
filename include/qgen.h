#ifndef QGEN_H
#define QGEN_H

#include <complex>
#include <vector>

#include "mpi.h"

#include "defs.h"
#include "qobservstate.h"
#include "qindivid.h"

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

class QFitnessClass
{
public:
    virtual ~QFitnessClass() {}
    virtual BASETYPE operator()( const QObservState& observState ) = 0;
};

//------------------------------------------------------------

struct QGenProcessSettings
{
    long long cycThreshold;
    double timeThreshold;
    int individsNum;
    int indSize;
    long long catastropheThreshold;
    float targetFitness;
    float accuracy;
    QFitnessClass* fClass;

    QGenProcessSettings()
        : cycThreshold(0)
        , timeThreshold( 0.0 )
        , individsNum(0)
        , indSize(0)
        , catastropheThreshold(0)
        , targetFitness(0.0f)
        , accuracy(0.0f)
        , fClass(0) {}
};

class QGenProcess
{
public:
    QGenProcess( const QGenProcessSettings& settings, MPI_Comm comm = MPI_COMM_WORLD );
    ~QGenProcess();

    inline bool active() const { return m_comm != MPI_COMM_NULL; }
    void process();

    const QIndivid& getBestIndivid() const { return m_processBest.individ; }

    inline bool isMaster() const { return m_myID == ROOT_ID; }

private:
    bool findBest();

private:
    static int m_instancesCount;
    MPI_Comm m_comm;
    int m_myID;
    int m_commSize;

    std::vector< QIndivid > m_individs;
    QGenProcessSettings m_settings;

    QObservState m_obsStatePreCached;

    struct BestSolution
    {
        int rank;
        int loc;
        float fitness;
        QIndivid individ;

        BestSolution()
            : rank(-1)
            , loc(-1)
            , fitness( 0.0f )
        {}
    };
    BestSolution m_processBest;

    BestSolution m_iterationBest;
    bool m_iterationBestDirty;

};

//-----------------------------------------------------------
}

//------------------------------------------------------------
#endif
