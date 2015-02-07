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

struct SQGenParams
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
    std::string outFile;

    SQGenParams( const char* file = 0, QFitnessClass* fC = 0, QRepairClass* rC = 0, QProcessScreen* sC = 0 );

    void init( const char* file, QFitnessClass* fC, QRepairClass* rC = 0, QProcessScreen* sC = 0 );
};

//------------------------------------------------------------

class QGenProcess
{
public:
    QGenProcess( const SQGenParams& params, MPI_Comm comm = MPI_COMM_WORLD );
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

    SQGenParams m_params;

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
