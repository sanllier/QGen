#ifndef QGEN_H
#define QGEN_H

#include <vector>

#include "mpi.h"
#include "mpicheck.h"
#include "qindivid_base.h"
#include "interfaces.h"
#include "randomizer.h"
#include "sparams.h"
#include "defs.h"

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

class QGenProcess
{
public:
    QGenProcess( const SParams& params, MPI_Comm comm = MPI_COMM_WORLD );
    ~QGenProcess();

    double process();

    QBaseIndivid* getBestIndivid( int coords[2] );

    bool isMaster() const;
    bool isMasterInd() const;

private:
    BASETYPE findIterationBestInd();

    bool active() const;

private:
    SParams m_params;

    struct SQGenProcessContext;
    SQGenProcessContext* m_ctx;

    std::vector< QBaseIndivid* > m_individs;

    QBaseIndivid* m_totalBest;
    QBaseIndivid* m_iterBest;

    Randomizer* m_randomizer;
};

//-----------------------------------------------------------
}

//------------------------------------------------------------
#endif
