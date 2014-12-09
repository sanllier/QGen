#include "qgen.h"
#include "qrotoperator.h"
#include "qobservstate.h"
#include "sharedmtrand.h"
#include "mpicheck.h"

#include <iostream>
#include <string>

#define ACTIVE_ONLY if ( !active() ) return
#define ACTIVE_ONLY_R( _x_ ) if ( !active() ) return _x_

#define MASTER_PRINT( _x_ ) if ( m_myID == ROOT_ID ) { std::cout << _x_ ; std::cout.flush(); }

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

int QGenProcess::m_instancesCount = 0;
MPI_Datatype QGenProcess::MPI_QBIT = MPI_DATATYPE_NULL;

//-----------------------------------------------------------

QGenProcess::QGenProcess( const QGenProcessSettings& settings, MPI_Comm comm/* = MPI_COMM_WORLD*/ )
    : m_settings( settings ) 
    , m_comm( MPI_COMM_NULL )
{
    if ( !m_settings.fClass )
        throw std::string( "Invalid fitness class: NULL pointer." ).append( __FUNCTION__ );

    int isMPIInitialized = 0;
    CHECK( MPI_Initialized( &isMPIInitialized ) );
    if ( !isMPIInitialized )
    {
        CHECK( MPI_Init( 0, 0 ) );
        CHECK( MPI_Type_vector( 1, 4, sizeof( BASETYPE ), MPI_BASETYPE, &MPI_QBIT ) );
        CHECK( MPI_Type_commit( &MPI_QBIT ) );
    }

    int initialCommRank = 0;
    int initialCommSize = 0;
    CHECK( MPI_Comm_rank( comm, &initialCommRank ) );
    CHECK( MPI_Comm_size( comm, &initialCommSize ) );

    if ( initialCommSize <= 0 )
        throw std::string( "Invalid communicator: comm size <= 0. " ).append( __FUNCTION__ );
         
    if ( initialCommRank == ROOT_ID )
        ++m_instancesCount;        
    CHECK( MPI_Bcast( &m_instancesCount, 1, MPI_INT, ROOT_ID, comm ) );

    const int activesCount = std::min( m_settings.individsNum, initialCommSize );
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
    for ( int i = 0; i < (int)m_individs.size(); ++i )
    {
        m_individs[i].resize( settings.indSize );
        m_individs[i].setInitial();
    }

    m_totalBest.ind.resize( settings.indSize );
    m_iterBest.ind.resize( settings.indSize );
}

QGenProcess::~QGenProcess()
{
    --m_instancesCount;

    if ( active() )
        CHECK( MPI_Comm_free( &m_comm ) );

    if ( m_instancesCount <= 0 )
    {
        CHECK( MPI_Type_free( &MPI_QBIT ) );
        CHECK( MPI_Finalize() );        
    }
}

//-----------------------------------------------------------

struct NetIndividRef
{
    float fitness;
    int procRank;
};

BASETYPE QGenProcess::findIterationBestInd()
{
    ACTIVE_ONLY_R( BASETYPE(0) );

    BASETYPE maxFit = BASETYPE(0);
    BASETYPE curFit = BASETYPE(0);
    int localBestIndIdx = 0;
    for ( int i = 0; i < (int)m_individs.size(); ++i )
    {
        if ( m_settings.repClass )
            m_individs[i].repair( m_settings.repClass );

        m_individs[i].updateObsState();
        //m_individs[i].getObsState().print();

        curFit = m_individs[i].getFitness( m_settings.fClass );
        if ( curFit > maxFit || i == 0 )
        {
            maxFit = curFit;
            localBestIndIdx = i;
        }
    }
    
    NetIndividRef localBestIndRef = { maxFit, m_myID };
    NetIndividRef globalBestIndRef;
    CHECK( MPI_Allreduce( &localBestIndRef, &globalBestIndRef, 1, MPI_FLOAT_INT, MPI_MAXLOC, m_comm ) );

    int remoteLocalBestIndIdx = localBestIndIdx;
    CHECK( MPI_Bcast( &remoteLocalBestIndIdx, 1, MPI_INT, globalBestIndRef.procRank, m_comm ) );

    if ( globalBestIndRef.procRank == m_myID )
        m_iterBest.ind = m_individs[ localBestIndIdx ];
    m_iterBest.ind.bcast( globalBestIndRef.procRank, m_comm );

    m_iterBest.procRank = globalBestIndRef.procRank;
    m_iterBest.localIdx = remoteLocalBestIndIdx;

    return globalBestIndRef.fitness;
}

//-----------------------------------------------------------

bool QGenProcess::immigration()
{
    ACTIVE_ONLY_R( false );

    //if ( m_settings.immigrationThreshold <= 0 || m_settings.immigrationSize < 0 )
    //    return false;

    //QIndivid& individ = m_iterationBest.individ;
    //const long long indSize = individ.qsize();
    //const long long immSize = m_settings.immigrationSize >= indSize ? indSize : m_settings.immigrationSize;

    //for ( long long i = 0; i < immSize; ++i )
    //{
    //    const size_t pos = ( size_t )( SharedMTRand::get32UnsignedInstance()() ) % indSize;
    //    QBit& temp = individ.at( pos );
    //    temp.a = ( float )( SharedMTRand::getClosedInstance()() );
    //    temp.b = std::sqrt( 1 - temp.a );
    //}
    
    MASTER_PRINT( "IMMIGRATION\n" );

    return true;
}

//-----------------------------------------------------------

void QGenProcess::process()
{
    ACTIVE_ONLY;

    for ( long long cycle = 1; cycle <= m_settings.cycThreshold; ++cycle )
    {
        BASETYPE bestFitness = findIterationBestInd();
        MASTER_PRINT( bestFitness << "\n" );

        if ( bestFitness > m_totalBest.ind.getFitnessUNSAFE() || cycle == 1 )
            m_totalBest = m_iterBest;
        
        if ( m_settings.targetFitness > 0.0f )
        {
            if ( m_settings.accuracy > 0.0f )
            {
                if ( fabs( m_totalBest.ind.getFitnessUNSAFE() - m_settings.targetFitness ) <= m_settings.accuracy )
                    return;
            }
            else if ( m_totalBest.ind.getFitnessUNSAFE() >= m_settings.targetFitness )
                return;
        }

        for ( int i = 0; i < (int)m_individs.size(); ++i )
            m_individs[i].tick( m_totalBest.ind );
    }
}

//------------------------------------------------------------
}
