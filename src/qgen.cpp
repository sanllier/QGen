#include "qgen.h"
#include "qindivid_cpu.h"
#include "sharedmtrand.h"
#include "mpicheck.h"
#include "pugixml.hpp"
#ifdef GPU
    #include "qindivid_gpu.h"
#endif

#include <string>
#include <string.h>

//------------------------------------------------------------

#define ACTIVE_ONLY if ( !active() ) return
#define ACTIVE_ONLY_R( _x_ ) if ( !active() ) return _x_

#define MASTER_PRINT( _x_ ) if ( m_ctx.generalRank == ROOT_ID ) { std::cout << _x_ ; std::cout.flush(); }

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

int QGenProcess::m_instancesCount = 0;

//------------------------------------------------------------

struct QGenProcess::SQGenProcessContext
{
    MPI_Comm generalComm;
    MPI_Comm rowComm;

    int generalRank;
    int rowRank;

    int generalSize;
    int rowSize;

    int coords[2];

    SQGenProcessContext()
        : generalComm( MPI_COMM_NULL )
        , rowComm( MPI_COMM_NULL )
        , generalRank( -1 )
        , rowRank( -1 )
        , generalSize(0)
        , rowSize(0)
    {}
};

//-----------------------------------------------------------

struct QGenProcess::SBestSolution
{
    int procRank;
    int localIdx;
    MPI_Comm rowComm;
    QBaseIndivid* ind;

    SBestSolution()
        : procRank(-1)
        , localIdx(-1)
        , rowComm( MPI_COMM_NULL )
        , ind(0)
    {}
    ~SBestSolution()
    {
        delete ind;
        // CRAP
        //CHECK( MPI_Comm_free( &rowComm ) );
    }

    SBestSolution& operator=( const SBestSolution& rSol )
    {
        procRank = rSol.procRank;
        localIdx = rSol.localIdx;
        *ind     = *rSol.ind;

        return *this;
    }
};

//-----------------------------------------------------------

SQGenParams::SQGenParams( const char* file/* = 0*/, QFitnessClass* fC/* = 0*/, QRepairClass* rC/* = 0*/, QProcessScreen* sC/* = 0*/)
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
{
    if ( file )
        init( file, fC, rC, sC );
}

//-----------------------------------------------------------

void SQGenParams::init( const char* file, QFitnessClass* fC, QRepairClass* rC/* = 0*/, QProcessScreen* sC/* = 0*/)
{
    if ( !file || !file[0] )
        throw std::string( "Some problems with params file. " ).append( __FUNCTION__ );

    pugi::xml_document doc;
    doc.load_file( file );
    
    pugi::xml_node qgenNode = doc.child( "qgen" );    
    if ( !qgenNode )
        throw std::string( "Some problems with params file. " ).append( __FUNCTION__ );

    for ( pugi::xml_node node = qgenNode.child( "parameter" ); node; node = node.next_sibling() )
    {
        const char* name = node.attribute( "name" ).as_string();
        if ( !name )
            continue;

        if ( 0 == strcmp( "cycle-threshold", name ) )
        {
            cycThreshold = (long long)node.attribute( "value" ).as_uint(0);
        }
        else if ( 0 == strcmp( "individs-num", name ) )
        {
            individsNum = (int)node.attribute( "value" ).as_int(0);
        }
        else if ( 0 == strcmp( "individ-size", name ) )
        {
            indSize = (int)node.attribute( "value" ).as_int(0);
        }
        else if ( 0 == strcmp( "topology-rows", name ) )
        {
            topoRows = (int)node.attribute( "value" ).as_int(0);
        }
        else if ( 0 == strcmp( "topology-cols", name ) )
        {
            topoCols = (int)node.attribute( "value" ).as_int(0);
        }
        else if ( 0 == strcmp( "target-fitness", name ) )
        {
            targetFitness = (BASETYPE)node.attribute( "value" ).as_double(0.0);
        }
        else if ( 0 == strcmp( "target-accuracy", name ) )
        {
            accuracy = (BASETYPE)node.attribute( "value" ).as_double(0.0);
        }
        else if ( 0 == strcmp( "out-file", name ) )
        {
            outFile = node.attribute( "value" ).as_string("");
        }
    #ifdef GPU
        else if ( 0 == strcmp( "use-gpu", name ) )
        {
            gpu = node.attribute( "value" ).as_bool( false );
        }
    #endif
    }

    fClass      = fC;
    repClass    = rC;
    screenClass = sC;
}

//-----------------------------------------------------------

QGenProcess::QGenProcess( const SQGenParams& params, MPI_Comm comm/* = MPI_COMM_WORLD*/ )
    : m_params( params ) 
    , m_ctx(0)
    , m_totalBest(0)
    , m_iterBest(0)
{
    if ( !m_params.fClass )
        throw std::string( "Invalid fitness class: NULL pointer." ).append( __FUNCTION__ );

    int isMPIInitialized = 0;
    CHECK( MPI_Initialized( &isMPIInitialized ) );
    if ( !isMPIInitialized )
        CHECK( MPI_Init( 0, 0 ) );

    int initialCommRank = 0;
    int initialCommSize = 0;
    CHECK( MPI_Comm_rank( comm, &initialCommRank ) );
    CHECK( MPI_Comm_size( comm, &initialCommSize ) );

    SharedMTRand::getClosedInstance( initialCommRank + 1 ); //setting seed

    const int requestedCommSize = m_params.topoCols * m_params.topoRows;

    if ( initialCommSize <= 0 )
        throw std::string( "Invalid communicator: comm size <= 0. " ).append( __FUNCTION__ );
    if ( initialCommSize < requestedCommSize )
        throw std::string( "Initial communicator size is smaller than requsted. " ).append( __FUNCTION__ );

    if ( initialCommRank == ROOT_ID )
        ++m_instancesCount;        
    CHECK( MPI_Bcast( &m_instancesCount, 1, MPI_INT, ROOT_ID, comm ) );

    m_ctx = new SQGenProcessContext();

    const int activesCount = std::min( requestedCommSize, initialCommSize );
    int* ranks = new int[ activesCount ];
    for ( int i = 0; i < activesCount; ++i )
        ranks[i] = i;

    MPI_Group initialGroup;
    MPI_Group workingGroup;
    CHECK( MPI_Comm_group( comm, &initialGroup ) );
    CHECK( MPI_Group_incl( initialGroup, activesCount, ranks, &workingGroup ) );
    CHECK( MPI_Comm_create( comm, workingGroup, &(m_ctx->generalComm) ) ); 

    delete[] ranks;
    CHECK( MPI_Group_free( &initialGroup ) );
    CHECK( MPI_Group_free( &workingGroup ) );

    if ( m_ctx->generalComm != MPI_COMM_NULL )
    {
        CHECK( MPI_Comm_rank( m_ctx->generalComm, &(m_ctx->generalRank) ) );
        CHECK( MPI_Comm_size( m_ctx->generalComm, &(m_ctx->generalSize) ) );

        if ( m_params.topoCols <= 0 || m_params.topoRows <= 0 )
            throw std::string( "Invalid topology settings. " ).append( __FUNCTION__ );

        int dims[2] = { m_params.topoRows, m_params.topoCols };
        int periods[2] = { 0, 0 };
        MPI_Comm cartComm = MPI_COMM_NULL;
        CHECK( MPI_Cart_create( m_ctx->generalComm, 2, dims, periods, false, &cartComm ) );
        if ( cartComm == MPI_COMM_NULL )
            throw std::string( "Some problems with cart communicator. " ).append( __FUNCTION__ );

        if ( m_params.topoCols > 1 )
        {
            int remainRowCommDims[2] = { 0, 1 };
            CHECK( MPI_Cart_sub( cartComm, remainRowCommDims, &(m_ctx->rowComm) ) );
            CHECK( MPI_Comm_rank( m_ctx->rowComm, &(m_ctx->rowRank) ) );
            CHECK( MPI_Comm_size( m_ctx->rowComm, &(m_ctx->rowSize) ) );
        }

        CHECK( MPI_Cart_coords( cartComm, m_ctx->generalRank, 2, m_ctx->coords ) );

        const int localIndsNum = m_params.individsNum / m_params.topoCols +
            ( m_ctx->coords[1] < m_params.individsNum % m_params.topoCols ? 1 : 0 );
        const int maxIters = m_params.individsNum / m_params.topoCols +
            ( m_params.individsNum % m_params.topoCols ? 1 : 0 );

        m_individs.resize( localIndsNum );
        for ( int i = 0; i < maxIters; ++i )
        {
            if ( i < localIndsNum )
            {
                m_individs[i] =
                #ifdef GPU
                    m_params.gpu ? (QBaseIndivid*)new QGPUIndivid( m_params.indSize, cartComm, m_ctx->rowComm, m_ctx->coords ):
                #endif
                    (QBaseIndivid*)new QCPUIndivid( m_params.indSize, cartComm, m_ctx->rowComm, m_ctx->coords );
            }
            else
            {
                #ifdef GPU
                    if ( m_params.gpu )
                        QGPUIndivid( m_params.indSize, cartComm, m_ctx->rowComm, m_ctx->coords );
                    else
                #endif
                        QCPUIndivid( m_params.indSize, cartComm, m_ctx->rowComm, m_ctx->coords );
            }
        }

        m_totalBest = new SBestSolution;
        m_iterBest  = new SBestSolution;

        if ( m_params.topoCols > 1 )
        {
            int remainRowCommDims[2] = { 0, 1 };
            CHECK( MPI_Cart_sub( cartComm, remainRowCommDims, &m_iterBest->rowComm) );
            m_iterBest->ind =
            #ifdef GPU
                m_params.gpu ? (QBaseIndivid*)new QGPUIndivid( m_params.indSize, cartComm, m_iterBest->rowComm, m_ctx->coords ):
            #endif
                (QBaseIndivid*)new QCPUIndivid( m_params.indSize, cartComm, m_iterBest->rowComm, m_ctx->coords );
        }
        else
        {
            m_iterBest->ind =
            #ifdef GPU
                m_params.gpu ? (QBaseIndivid*)new QGPUIndivid( m_params.indSize, cartComm, MPI_COMM_NULL, m_ctx->coords ):
            #endif
                (QBaseIndivid*)new QCPUIndivid( m_params.indSize, cartComm, MPI_COMM_NULL, m_ctx->coords );
        }

        m_totalBest->ind =
        #ifdef GPU
            m_params.gpu ? (QBaseIndivid*)new QGPUIndivid( m_params.indSize, cartComm, MPI_COMM_NULL, m_ctx->coords ):
        #endif
            (QBaseIndivid*)new QCPUIndivid( m_params.indSize, cartComm, MPI_COMM_NULL, m_ctx->coords );
        
        CHECK( MPI_Comm_free( &cartComm ) );
    }
    else
    {
        m_ctx->generalRank = -1;
        m_ctx->generalSize = 0;
    }
}

//-----------------------------------------------------------

QGenProcess::~QGenProcess()
{
    for ( int i = 0; i < int( m_individs.size() ); ++i )
        delete m_individs[i];

    if ( active() )
    {
        CHECK( MPI_Comm_free( &(m_ctx->generalComm) ) );
        delete m_totalBest;
        delete m_iterBest;
        // CRAP
        //CHECK( MPI_Comm_free( &m_ctx.rowComm ) );  
    }
    
    delete m_ctx;

    --m_instancesCount;
    if ( m_instancesCount <= 0 )
        CHECK( MPI_Finalize() );        
}

//-----------------------------------------------------------

const QBaseIndivid* QGenProcess::getBestIndivid() const
{
    return m_totalBest->ind;
}

//-----------------------------------------------------------

bool QGenProcess::isMaster() const
{ 
    return m_ctx->generalRank == ROOT_ID; 
}

//-----------------------------------------------------------

bool QGenProcess::isMasterInd() const
{
    return m_ctx->coords[1] == 0; 
}

//-----------------------------------------------------------

bool QGenProcess::active() const 
{ 
    return m_ctx->generalComm != MPI_COMM_NULL;
}

//-----------------------------------------------------------

double QGenProcess::process()
{
    ACTIVE_ONLY_R( 0.0 );

    if ( m_params.individsNum < 1 )
        throw std::string( "QGenProcess trying to process with too small number of individs(must be >=1). " ).append( __FUNCTION__ ); 

    double startTime = MPI_Wtime();

    for ( long long cycle = 1; cycle <= m_params.cycThreshold; ++cycle )
    {
        BASETYPE bestFitness = findIterationBestInd();

        if ( bestFitness > m_totalBest->ind->getFitness() || cycle == 1 )
            *m_totalBest = *m_iterBest;

        if ( m_params.screenClass )
            ( *m_params.screenClass )( cycle, m_ctx->coords, *( m_totalBest->ind ), *( m_iterBest->ind ) );

        if ( m_params.targetFitness > 0.0f )
        {
            if ( m_params.accuracy > 0.0f )
            {
                if ( fabs( m_totalBest->ind->getFitness() - m_params.targetFitness ) <= m_params.accuracy )
                    return MPI_Wtime() - startTime;
            }
            else if ( m_totalBest->ind->getFitness() >= m_params.targetFitness )
                return MPI_Wtime() - startTime;
        }

        for ( int i = 0; i < (int)m_individs.size(); ++i )
            m_individs[i]->evolve( *( m_iterBest->ind ) );
    }

    return MPI_Wtime() - startTime;
}

//-----------------------------------------------------------

struct SNetIndividRef
{
    float fitness;
    int   procRank;
};

BASETYPE QGenProcess::findIterationBestInd()
{
    ACTIVE_ONLY_R( BASETYPE(0) );

    BASETYPE maxFit = BASETYPE(0);
    BASETYPE curFit = BASETYPE(0);
    int localBestIndIdx = 0;
    for ( int i = 0; i < (int)m_individs.size(); ++i )
    {
        m_individs[i]->calculateObservState();

        if ( m_params.repClass )
            m_individs[i]->repair( m_params.repClass );

        curFit = m_individs[i]->calculateFitness( m_params.fClass );

        if ( curFit > maxFit || i == 0 )
        {
            maxFit = curFit;
            localBestIndIdx = i;
        }
    }
    
    SNetIndividRef localBestIndRef = { maxFit, m_ctx->rowRank };
    SNetIndividRef globalBestIndRef;
    int remoteLocalBestIndIdx = localBestIndIdx;

    if ( m_ctx->rowComm != MPI_COMM_NULL )
    {
        CHECK( MPI_Allreduce( &localBestIndRef, &globalBestIndRef, 1, MPI_FLOAT_INT, MPI_MAXLOC, m_ctx->rowComm ) );
        CHECK( MPI_Bcast( &remoteLocalBestIndIdx, 1, MPI_INT, globalBestIndRef.procRank, m_ctx->rowComm ) );
    }
    else
    {
        globalBestIndRef = localBestIndRef;
    }       

    if ( globalBestIndRef.procRank == m_ctx->rowRank )
        *( m_iterBest->ind ) = *m_individs[ localBestIndIdx ];

    m_iterBest->ind->bcast( globalBestIndRef.procRank );

    m_iterBest->procRank = globalBestIndRef.procRank;
    m_iterBest->localIdx = remoteLocalBestIndIdx;

    return globalBestIndRef.fitness;
}

//------------------------------------------------------------
}
