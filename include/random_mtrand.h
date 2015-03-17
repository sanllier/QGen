#if defined( RANDOM_MTRAND_H ) && defined( MTRAND )
#define RANDOM_MTRAND_H

#include "random.h"
#include "mtrand.h"

//------------------------------------------------------------

class RandomMTRand: public IRandom
{
public:
    RandomMTRand();
    ~RandomMTRand();

    void setSeed( unsigned seed ) OVERRIDE
    {
        m_randomizer.seed(seed);
    }

    BASETYPE next() OVERRIDE
    {
        return BASETYPE( m_randomizer() );
    }

private:
    static MTRand_closed m_randomizer;
};

//------------------------------------------------------------
#endif
