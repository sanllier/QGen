#ifndef RANDOMIZER_H
#define RANDOMIZER_H

#include "defs.h"
#include "random.h"

//------------------------------------------------------------

class Randomizer
{
    Randomizer( IRandom* random, unsigned seed ) 
    {
        m_rand = random;
        m_rand->setSeed( seed );
    }

    ~Randomizer() 
    {
        delete m_rand;
    }

public:
    static Randomizer* Create( IRandom* random, unsigned seed )
    {
        if ( m_self != 0 )
            delete m_self;

        m_self = new Randomizer( random, seed );
        return m_self;
    }

    static void Dispose()
    {
        delete m_self;
        m_self = 0;
    }

    static IRandom* Get()
    {
        return m_self->m_rand;
    }

private:
    IRandom* m_rand;
    static Randomizer* m_self;
};

//------------------------------------------------------------

#endif
