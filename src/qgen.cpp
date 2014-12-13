#include "qgen.h"
#include "qrotoperator.h"
#include "qobservstate.h"
#include "sharedmtrand.h"
#include "mpicheck.h"

#include <iostream>
#include <string>

#define ACTIVE_ONLY if ( !active() ) return
#define ACTIVE_ONLY_R( _x_ ) if ( !active() ) return _x_

#define MASTER_PRINT( _x_ ) if ( m_ctx.generalRank == ROOT_ID ) { std::cout << _x_ ; std::cout.flush(); }

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

int QGenProcess::m_instancesCount = 0;
MPI_Datatype QGenProcess::MPI_QBIT = MPI_DATATYPE_NULL;

//-----------------------------------------------------------

QGenProcess::QGenProcess( const QGenProcessSettings& settings, MPI_Comm comm/* = MPI_COMM_WORLD*/ )
    : m_settings( settings ) 
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
    const int requestedCommSize = m_settings.topoCols * m_settings.topoRows;

    if ( initialCommSize <= 0 )
        throw std::string( "Invalid communicator: comm size <= 0. " ).append( __FUNCTION__ );
    if ( initialCommSize < requestedCommSize )
        throw std::string( "Initial communicator size is smaller than requsted. " ).append( __FUNCTION__ );

    if ( initialCommRank == ROOT_ID )
        ++m_instancesCount;        
    CHECK( MPI_Bcast( &m_instancesCount, 1, MPI_INT, ROOT_ID, comm ) );

    const int activesCount = std::min( requestedCommSize, initialCommSize );
    int* ranks = new int[ activesCount ];
    for ( int i = 0; i < activesCount; ++i )
        ranks[i] = i;

    MPI_Group initialGroup;
    MPI_Group workingGroup;
    CHECK( MPI_Comm_group( comm, &initialGroup ) );
    CHECK( MPI_Group_incl( initialGroup, activesCount, ranks, &workingGroup ) );

    delete[] ranks;
    ranks = 0;
    CHECK( MPI_Comm_create( comm, workingGroup, &m_ctx.generalComm ) ); 
    CHECK( MPI_Group_free( &initialGroup ) );
    CHECK( MPI_Group_free( &workingGroup ) );

    if ( m_ctx.generalComm != MPI_COMM_NULL )
    {
        CHECK( MPI_Comm_rank( m_ctx.generalComm, &m_ctx.generalRank ) );
        CHECK( MPI_Comm_size( m_ctx.generalComm, &m_ctx.generalSize ) );

        int dims[2] = { m_settings.topoRows, m_settings.topoCols };
        int periods[2] = { 0, 0 };
        MPI_Comm cartComm = MPI_COMM_NULL;
        CHECK( MPI_Cart_create( m_ctx.generalComm, 2, dims, periods, false, &cartComm ) );
        if ( cartComm == MPI_COMM_NULL )
            throw std::string( "Some problems with cart communicator. " ).append( __FUNCTION__ );

        if ( m_settings.topoCols > 1 )
        {
            int remainRowCommDims[2] = { 0, 1 };
            CHECK( MPI_Cart_sub( cartComm, remainRowCommDims, &m_ctx.rowComm ) );
            CHECK( MPI_Comm_rank( m_ctx.rowComm, &m_ctx.rowRank ) );
            CHECK( MPI_Comm_size( m_ctx.rowComm, &m_ctx.rowSize ) );
        }

        CHECK( MPI_Cart_coords( cartComm, m_ctx.generalRank, 2, m_ctx.coords ) );

        const int localIndsNum = m_settings.individsNum / m_settings.topoCols +
            ( m_ctx.coords[1] < m_settings.individsNum % m_settings.topoCols ? 1 : 0 );
        const int maxIters = m_settings.individsNum / m_settings.topoCols +
            ( m_settings.individsNum % m_settings.topoCols ? 1 : 0 );

        m_individs.resize( localIndsNum );
        for ( int i = 0; i < maxIters; ++i )
            if ( i < localIndsNum )
                m_individs[i] = new QIndivid( m_settings.indSize, cartComm, m_ctx.rowComm, m_ctx.coords );
            else
                QIndivid( m_settings.indSize, cartComm, m_ctx.rowComm, m_ctx.coords );

        m_totalBest.ind = new QIndivid( m_settings.indSize, cartComm, MPI_COMM_NULL, m_ctx.coords );
        m_iterBest.ind  = new QIndivid( m_settings.indSize, cartComm, MPI_COMM_NULL, m_ctx.coords );

        CHECK( MPI_Comm_free( &cartComm ) );
    }
    else
    {
        m_ctx.generalRank = -1;
        m_ctx.generalSize = 0;
    }
}

QGenProcess::~QGenProcess()
{
    --m_instancesCount;

    if ( active() )
        CHECK( MPI_Comm_free( &m_ctx.generalComm ) );

    for ( int i = 0; i < int( m_individs.size() ); ++i )
        delete m_individs[i];

    //delete m_totalBest.ind;
    //delete m_iterBest.ind;

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
        m_individs[i]->updateObsState();

        if ( m_settings.repClass )
            m_individs[i]->repair( m_settings.repClass );

        //m_individs[i].getObsState().print();

        curFit = m_individs[i]->getFitness( m_settings.fClass );
        if ( curFit > maxFit || i == 0 )
        {
            maxFit = curFit;
            localBestIndIdx = i;
        }
    }
    
    NetIndividRef localBestIndRef = { maxFit, m_ctx.rowRank };
    NetIndividRef globalBestIndRef;
    if ( m_ctx.rowComm != MPI_COMM_NULL )
        CHECK( MPI_Allreduce( &localBestIndRef, &globalBestIndRef, 1, MPI_FLOAT_INT, MPI_MAXLOC, m_ctx.rowComm ) );
    else
        globalBestIndRef = localBestIndRef;
    
    int remoteLocalBestIndIdx = localBestIndIdx;
    if ( globalBestIndRef.procRank >= 0 )
        CHECK( MPI_Bcast( &remoteLocalBestIndIdx, 1, MPI_INT, globalBestIndRef.procRank, m_ctx.rowComm ) );

    if ( globalBestIndRef.procRank == m_ctx.rowRank )
        *m_iterBest.ind = *m_individs[ localBestIndIdx ];

    m_iterBest.ind->bcast( globalBestIndRef.procRank );   

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

double QGenProcess::process()
{
    ACTIVE_ONLY_R( 0.0 );

    double startTime = MPI_Wtime();

    for ( long long cycle = 1; cycle <= m_settings.cycThreshold; ++cycle )
    {
        BASETYPE bestFitness = findIterationBestInd();
        MASTER_PRINT( bestFitness << "\n" );

        if ( bestFitness > m_totalBest.ind->getFitnessUNSAFE() || cycle == 1 )
            m_totalBest = m_iterBest;
        
        if ( m_settings.targetFitness > 0.0f )
        {
            if ( m_settings.accuracy > 0.0f )
            {
                if ( fabs( m_totalBest.ind->getFitnessUNSAFE() - m_settings.targetFitness ) <= m_settings.accuracy )
                    return MPI_Wtime() - startTime;
            }
            else if ( m_totalBest.ind->getFitnessUNSAFE() >= m_settings.targetFitness )
                return MPI_Wtime() - startTime;
        }

        for ( int i = 0; i < (int)m_individs.size(); ++i )
            m_individs[i]->tick( *m_totalBest.ind );
    }

    return MPI_Wtime() - startTime;
}

//------------------------------------------------------------
}
