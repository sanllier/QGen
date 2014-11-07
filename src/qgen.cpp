#include "qgen.h"
#include "qrotoperator.h"

#include <iostream>
#include <string>

#define CHECK( __ERR_CODE__ ) checkMPIRes( __ERR_CODE__, __FUNCTION__ )
#define ACTIVE_ONLY if ( !active() ) return
#define ACTIVE_ONLY_R( _x_ ) if ( !active() ) return _x_

#define MASTER_PRINT( _x_ ) if ( m_myID == ROOT_ID ) { std::cout << _x_ ; std::cout.flush(); }

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

int QGenProcess::m_instancesCount = 0;

//-----------------------------------------------------------

void checkMPIRes( int errCode, const char* location )
{
    if ( errCode != MPI_SUCCESS )
    {
        char errText[ MPI_MAX_ERROR_STRING ];
        int len;
        MPI_Error_string( errCode, errText, &len );
        throw std::string( errText ).append( location );
    }
}

//-----------------------------------------------------------

QGenProcess::QGenProcess( const QGenProcessSettings& settings, MPI_Comm comm/* = MPI_COMM_WORLD*/ )
    : m_settings( settings ) 
{
    if ( !m_settings.fClass )
        throw std::string( "Invalid fitness class: NULL pointer." ).append( __FUNCTION__ );

    int isMPIInitialized = 0;
    CHECK( MPI_Initialized( &isMPIInitialized ) );
    if ( !isMPIInitialized )
        CHECK( MPI_Init( 0, 0 ) );

    CHECK( MPI_Comm_rank( comm, &m_myID ) );
    int initialCommSize = 0;
    CHECK( MPI_Comm_size( comm, &initialCommSize ) );

    if ( initialCommSize <= 0 )
        throw std::string( "Invalid communicator: comm size <= 0. " ).append( __FUNCTION__ );
         
    if ( m_myID == ROOT_ID )
        ++m_instancesCount;        
    CHECK( MPI_Bcast( &m_instancesCount, 1, MPI_INT, ROOT_ID, comm ) );

    const int activesCount = m_settings.individsNum / initialCommSize <= 0 ? m_settings.individsNum : initialCommSize;
    int *ranks = new int[ activesCount ];
    for ( int i = 0; i < activesCount; ++i )
        ranks[i] = i;

    MPI_Group initialGroup;
    MPI_Group workingGroup;
    CHECK( MPI_Comm_group( comm, &initialGroup ) );
    CHECK( MPI_Group_incl( initialGroup, activesCount, ranks, &workingGroup ) );
    delete[] ranks;
    ranks = 0;

    CHECK( MPI_Comm_create( comm, workingGroup, &m_comm ) ); 
    CHECK( MPI_Group_free( &initialGroup ) );
    CHECK( MPI_Group_free( &workingGroup ) );
    if ( m_comm != MPI_COMM_NULL )
    {
        CHECK( MPI_Comm_rank( m_comm, &m_myID ) );
        CHECK( MPI_Comm_size( m_comm, &m_commSize ) );
    }
    else
    {
        m_myID = -1;
        m_commSize = 0;
        return;
    }

    const int myIndsNum = m_settings.individsNum / m_commSize + ( m_myID < m_settings.individsNum % m_commSize ? 1 : 0 );
    m_individs.resize( myIndsNum < 0 ? 0 : myIndsNum );
    for ( QIndivid& ind : m_individs )
    {
        ind.resize( settings.indSize );
        ind.setInitial();
    }

    m_processBest.individ.resize( settings.indSize );
    m_iterationBest.individ.resize( settings.indSize );
}

QGenProcess::~QGenProcess()
{
    --m_instancesCount;

    if ( active() )
        CHECK( MPI_Comm_free( &m_comm ) );

    if ( m_instancesCount <= 0 )
        CHECK( MPI_Finalize() );
}

//-----------------------------------------------------------

struct FLoatInt
{
    float fitness;
    int rank;
};

bool QGenProcess::findBest()
{
    ACTIVE_ONLY_R( false );

    BASETYPE max = BASETYPE(0);
    BASETYPE cur = BASETYPE(0);
    int localBestIndividIdx = 0;
    for ( int i = 0; i < (int)m_individs.size(); ++i )
    {
        m_obsStatePreCached.process( m_individs[i] );
        cur = ( *m_settings.fClass )( m_obsStatePreCached );
        if ( cur > max || i == 0 )
        {
            max = cur;
            localBestIndividIdx = i;
        }
    }
    
    FLoatInt localBestIndividRef = { max, m_myID };
    FLoatInt globalBestIndividRef;
    CHECK( MPI_Allreduce( &localBestIndividRef, &globalBestIndividRef, 1, MPI_FLOAT_INT, MPI_MAXLOC, m_comm ) );
    CHECK( MPI_Bcast( &localBestIndividIdx, 1, MPI_INT, globalBestIndividRef.rank, m_comm ) );

    const bool dontReplaceFlag = ( m_processBest.rank == globalBestIndividRef.rank && m_processBest.loc == localBestIndividIdx ) \
        || m_processBest.fitness > globalBestIndividRef.fitness ;

    BestSolution& targetSln = dontReplaceFlag ? m_iterationBest : m_processBest;
    QIndivid& targetInd = globalBestIndividRef.rank == m_myID ? m_individs[ localBestIndividIdx ] : targetSln.individ;
    m_iterationBestDirty = !dontReplaceFlag;

    CHECK( MPI_Bcast( targetInd.raw(), targetInd.qsize() * 2, MPI_FLOAT, globalBestIndividRef.rank, m_comm ) );
    if ( globalBestIndividRef.rank == m_myID )
        targetSln.individ = m_individs[ localBestIndividIdx ];

    targetSln.rank    = globalBestIndividRef.rank;
    targetSln.loc     = localBestIndividIdx;
    targetSln.fitness = globalBestIndividRef.fitness;
    return m_iterationBestDirty;
}

//-----------------------------------------------------------

void QGenProcess::process()
{
    ACTIVE_ONLY;

    double elapsedTime = 0;
    long long catastropheGen = 1;
    long long bestSolutionGen = 0;
    float factor = 0.0f;

    for ( long long cycle = 1; cycle <= m_settings.cycThreshold; ++cycle )
    {
        factor = 0.5f * ( ( m_settings.cycThreshold - cycle ) / catastropheGen ) / m_settings.cycThreshold;
        if ( !findBest() )
            ++bestSolutionGen;
        else
            bestSolutionGen = 0;

        MASTER_PRINT( m_processBest.fitness << "\n" );

        if ( m_settings.targetFitness > 0.0f )
        {
            if ( m_settings.accuracy > 0.0f )
            {
                if ( fabs( m_processBest.fitness - m_settings.targetFitness ) <= m_settings.accuracy )
                    return;
            }
            else if ( m_processBest.fitness >= m_settings.targetFitness )
                return;
        }

        if ( m_settings.targetFitness > 0.0f && m_processBest.fitness >= m_settings.targetFitness )
            return;

        for ( QIndivid& ind : m_individs )
            ind.tick( m_processBest.individ, factor );

        if ( m_settings.catastropheThreshold > 0 && m_settings.catastropheThreshold <= bestSolutionGen )
        {
            m_processBest = m_iterationBest;
            catastropheGen = cycle;
        }
    }
}

//------------------------------------------------------------
}
