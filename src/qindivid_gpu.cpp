#ifdef GPU

#include "qindivid_gpu.h"
#include "qgen.h"
#include "defs_gpu.h"

#include "curand.h"
#include "cuda_runtime.h"
#include "cuda_error_handler.h"

#include "mpicheck.h"

#include <time.h>

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

int QGPUIndivid::deviceId = -1;

//------------------------------------------------------------

extern "C" void launchSetInitialKernel( GPUQbit* data, const BASETYPE* rands, long long size );
extern "C" void launchObserveKernel( const GPUQbit* data, const BASETYPE* rands, bool* buf, long long size );
extern "C" void launchThetaBufferKernel( const GPUQbit* data, const GPUQbit* best, const bool* obs,
                                        const bool* bestObs, const BASETYPE* theta, const bool* isBetter,
                                        BASETYPE* out, long long size );
extern "C" void launchEvolveKernel( GPUQbit* data, const BASETYPE* theta, long long size );

//------------------------------------------------------------

BASETYPE* QGPUIndivid::m_gpuThetaFiled = 0;

//------------------------------------------------------------

QGPUIndivid::QGPUIndivid( long long size, MPI_Comm generalComm, MPI_Comm rowComm, int coords[2], BASETYPE thetaFrac/* = BASETYPE(1)*/ )
    : QBaseIndivid( size, generalComm,  rowComm, coords, thetaFrac )
    , m_cpuBuf(0)
    , m_thetaBuf(0)
    , m_isBetterGPUBuf(0)
    , m_isBetterBuf(0)
    , m_gen(0)
    , m_randomBufGPU(0)
{
    if ( deviceId == -1 )
    {
        deviceId = selectDevice();
        SAFE_CALL( cudaSetDevice( deviceId ) );
    }

    if ( !m_gpuThetaFiled )
    {
        SAFE_CALL( cudaMalloc( &m_gpuThetaFiled, 8 * sizeof( BASETYPE ) ) );
        SAFE_CALL( cudaMemcpy( m_gpuThetaFiled, m_thetaField, 8 * sizeof( BASETYPE ), cudaMemcpyHostToDevice ) );
    }

    m_gen = new curandGenerator_t;
    SAFE_CURAND_CALL( curandCreateGenerator( ( curandGenerator_t * )m_gen, CURAND_RNG_PSEUDO_DEFAULT ) );

    resize( size );
    setInitial(); 
}

//------------------------------------------------------------

QGPUIndivid::~QGPUIndivid()
{
    if ( m_data )
        SAFE_CALL( cudaFree( m_data ) );

    if ( m_thetaBuf )
        SAFE_CALL( cudaFree( m_thetaBuf ) );

    if ( m_isBetterGPUBuf )
        SAFE_CALL( cudaFree( m_isBetterGPUBuf ) );

    SAFE_CURAND_CALL( curandDestroyGenerator( *( curandGenerator_t * )m_gen ) );
    SAFE_CALL( cudaFree( m_randomBufGPU ) );

    if ( m_isBetterBuf )
        delete[] m_isBetterBuf;

    if ( m_cpuBuf )
        delete[] m_cpuBuf;
}

//------------------------------------------------------------

bool QGPUIndivid::resize( long long newSize )
{
    if ( !QBaseIndivid::resize( newSize ) )
        return false;

    if ( m_data )
        SAFE_CALL( cudaFree( m_data ) );
    SAFE_CALL( cudaMalloc( &m_data, size_t( m_localLogicSize * sizeof( QBit ) ) ) );

    if ( m_thetaBuf )
        SAFE_CALL( cudaFree( m_thetaBuf ) );
    SAFE_CALL( cudaMalloc( &m_thetaBuf, size_t( m_localLogicSize * sizeof( BASETYPE ) ) ) );

    if ( m_isBetterGPUBuf )
        SAFE_CALL( cudaFree( m_isBetterGPUBuf ) );
    SAFE_CALL( cudaMalloc( &m_isBetterGPUBuf, size_t( m_localLogicSize ) ) );

    if ( m_isBetterBuf )
        delete[] m_isBetterBuf;
    m_isBetterBuf = new bool[ size_t( m_localLogicSize ) ];

    if ( m_cpuBuf )
        delete[] m_cpuBuf;
    m_cpuBuf = new QBit[ size_t( m_localLogicSize ) ];

    if ( m_randomBufGPU )
        SAFE_CALL( cudaFree( m_randomBufGPU ) );
    SAFE_CALL( cudaMalloc( &m_randomBufGPU, size_t( newSize * sizeof( BASETYPE ) ) ) );

    return true;
}

//------------------------------------------------------------

void QGPUIndivid::setInitial()
{
    SAFE_CURAND_CALL( curandGenerateUniform( *( curandGenerator_t * )m_gen, m_randomBufGPU, size_t( m_localLogicSize ) ));
    launchSetInitialKernel( ( GPUQbit* )m_data, m_randomBufGPU, m_localLogicSize );
    m_needRecalcFitness = true;
}

//------------------------------------------------------------

void QGPUIndivid::evolve( const QBaseIndivid& bestInd )
{
    const QGPUIndivid& castedIndivid = ( const QGPUIndivid& )bestInd;
    computeThetaBuffer( castedIndivid );
    launchEvolveKernel( ( GPUQbit* )m_data, m_thetaBuf, m_localLogicSize );
}

//------------------------------------------------------------

bool QGPUIndivid::bcast( int root )
{
    if ( !QBaseIndivid::bcast( root ) )
        return false;

    int rank = 0;
    CHECK( MPI_Comm_rank( m_context.rowComm, &rank ) );

    if ( rank = root )
        SAFE_CALL( cudaMemcpy( m_cpuBuf, m_data, size_t( m_localLogicSize * sizeof( QBit ) ), cudaMemcpyDeviceToHost ) );

    CHECK( MPI_Bcast( m_cpuBuf, int( m_localLogicSize ), MPI_QBIT, root, m_context.rowComm ) );
    SAFE_CALL( cudaMemcpy( m_data, m_cpuBuf, size_t( m_localLogicSize * sizeof( QBit ) ), cudaMemcpyHostToDevice ) );
    return true;
}

//-----------------------------------------------------------

QBaseIndivid& QGPUIndivid::operator=( const QBaseIndivid& rInd )
{
    // CRAP
    const QGPUIndivid& castedIndivid = ( const QGPUIndivid& )rInd;
    if ( m_localLogicSize != castedIndivid.m_localLogicSize )
        throw std::string( "Trying to assing individs with different topologies. " ).append( __FUNCTION__ );

    switch ( rInd.getType() )
    {
        case INDIVID_TYPE_CPU:
        {
            SAFE_CALL( cudaMemcpy( m_data, castedIndivid.m_data, 
                size_t( m_localLogicSize * sizeof( m_data[0] ) ), cudaMemcpyHostToDevice ) );
            break;
        }

        case INDIVID_TYPE_GPU:
        {
            SAFE_CALL( cudaMemcpy( m_data, castedIndivid.m_data, 
                size_t( m_localLogicSize * sizeof( m_data[0] ) ), cudaMemcpyDeviceToDevice ) );
            break;
        }

        default:
        {
            throw std::string( "QCPUIndivid is trying to assign individ with unknown type. " ).append( __FUNCTION__ );  
            break;
        }
    }

    *m_observeState     = *castedIndivid.m_observeState;
    m_fitness           = castedIndivid.m_fitness;
    m_needRecalcFitness = castedIndivid.m_needRecalcFitness;

    return *this;
}

//-----------------------------------------------------------

void QGPUIndivid::runObserveKernel( bool* gpuData ) const
{
    SAFE_CURAND_CALL( curandGenerateUniform( *( curandGenerator_t * )m_gen, m_randomBufGPU, size_t( m_localLogicSize ) ));
    launchObserveKernel( ( GPUQbit* )m_data, m_randomBufGPU, gpuData, m_localLogicSize );
}

//-----------------------------------------------------------

void QGPUIndivid::computeThetaBuffer( const QGPUIndivid& bestInd ) const
{
    for ( long long i = 0; i < m_localLogicSize; ++i )
        m_isBetterBuf[i] = m_fitness > bestInd.m_fitness;

    SAFE_CALL( cudaMemcpy( m_isBetterGPUBuf, m_isBetterBuf, size_t( m_localLogicSize ), cudaMemcpyHostToDevice ) );

    launchThetaBufferKernel( ( GPUQbit* )m_data, ( GPUQbit* )bestInd.m_data, 
        m_observeState->gpuBuffer(), bestInd.m_observeState->gpuBuffer(), m_gpuThetaFiled, 
        m_isBetterGPUBuf, m_thetaBuf, m_localLogicSize );
}

//------------------------------------------------------------

int QGPUIndivid::selectDevice()
{
   int deviceCount = 0;
   int suitableDevice = -1;

   SAFE_CALL( cudaGetDeviceCount( &deviceCount ) );
   if ( deviceCount <= 0 )
       throw std::string( "Not found CUDA devices. " ).append( __FUNCTION__ );

   for ( int device = 0; device < deviceCount; ++device )
   {
       cudaDeviceProp devProp;
       SAFE_CALL( cudaGetDeviceProperties ( &devProp, device ) );

       if( devProp.major >= 2 )
          suitableDevice = device ;
   }

   if ( suitableDevice <= -1 )
       throw std::string( "Not found suitable CUDA device. " ).append( __FUNCTION__ );

   return suitableDevice;
}

//------------------------------------------------------------
}

#endif
