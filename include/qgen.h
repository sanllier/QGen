#ifndef QGEN_H
#define QGEN_H

#include <complex>
#include <vector>

#include "mpi.h"
#include "mpicheck.h"
#include "qindivid_base.h"
#include "defs.h"

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------
   
class QProcessScreen
{
public:
    virtual void operator()( long long cycle, 
                             const int coords[2], 
                             const QBaseIndivid& totalBest, 
                             const QBaseIndivid& iterBest ) = 0;
};

//------------------------------------------------------------

struct SQGenProcessSettings
{
#ifdef GPU
    bool gpu;
#endif
    long long cycThreshold;
    int individsNum;
    int indSize;
    int topoRows;
    int topoCols;
    BASETYPE targetFitness;
    BASETYPE accuracy;
    QFitnessClass*  fClass;
    QRepairClass*   repClass;
    QProcessScreen* screenClass;

    SQGenProcessSettings()
        : cycThreshold(0)
        , individsNum(0)
        , indSize(0)
        , topoRows(1)
        , topoCols(1)
        , targetFitness( BASETYPE(0) )
        , accuracy( BASETYPE(0) )
        , fClass(0) 
        , repClass(0)
        , screenClass(0)
    #ifdef GPU
        , gpu(false)
    #endif
    {}
};

//------------------------------------------------------------

class QGenProcess
{
public:
    QGenProcess( const SQGenProcessSettings& settings, MPI_Comm comm = MPI_COMM_WORLD );
    ~QGenProcess();

    double process();

    const QBaseIndivid* getBestIndivid() const;

    bool isMaster() const;
    bool isMasterInd() const;

private:
    BASETYPE findIterationBestInd();

    bool active() const;

private:
    static int m_instancesCount;

    SQGenProcessSettings m_settings;

    struct SQGenProcessContext;
    SQGenProcessContext* m_ctx;

    std::vector< QBaseIndivid* > m_individs;
    
    struct SBestSolution;
    SBestSolution* m_totalBest;
    SBestSolution* m_iterBest;
};

//-----------------------------------------------------------
}

//------------------------------------------------------------
#endif
