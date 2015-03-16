#if !defined( QINDIVID_GPU_H ) && defined( GPU )
#define QINDIVID_GPU_H

#include "qindivid_base.h"

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

class QGPUIndivid: public QBaseIndivid
{
public:
    QGPUIndivid( long long size, MPI_Comm generalComm, MPI_Comm rowComm, int coords[2] );
    ~QGPUIndivid() OVERRIDE;

        // QBaseIndivid

    EIndividType getType() const OVERRIDE { return INDIVID_TYPE_GPU; };

    bool resize( long long newSize ) OVERRIDE;
    void setInitial() OVERRIDE;
    
    void evolve( const QBaseIndivid& bestInd ) OVERRIDE;

    bool bcast( int root ) OVERRIDE;

    QBaseIndivid& operator=( const QBaseIndivid& rInd ) OVERRIDE;

public:
    void runObserveKernel( bool* gpuData ) const;
    
private:
    int selectDevice();
    void computeThetaBuffer( const QGPUIndivid& bestInd ) const;

private:
    static int deviceId;

    QBit* m_cpuBuf;

    static BASETYPE* m_gpuThetaFiled;
    BASETYPE* m_thetaBuf;
    bool* m_isBetterGPUBuf;
    bool* m_isBetterBuf;

    void* m_gen;
    BASETYPE* m_randomBufGPU;
};

//------------------------------------------------------------
}

//------------------------------------------------------------
#endif
