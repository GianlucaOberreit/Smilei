#ifndef PARTCOMPTIME2D2ORDER_H
#define PARTCOMPTIME2D2ORDER_H

#include "PartCompTime.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PartCompTime2D2Order
//! Evaluation du temps de calcul des macro-particules en mode AM Order 2
//  --------------------------------------------------------------------------------------------------------------------
class PartCompTime2D2Order final : public PartCompTime
{
public:
    PartCompTime2D2Order();
    ~PartCompTime2D2Order() override final {};
    
    // -------------------------------------------------------------------------
    //! Evaluate the time (simple precision) to compute all particles
    //! in the current patch with vectorized operators
    //! @param count the numer of particles per cell
    //! @aram vecto_time time in vector mode
    //! @aram scalar_time time in scalar mode
    // -------------------------------------------------------------------------
    virtual void operator() (   const std::vector<int> &count,
                        float &vecto_time,
                        float &scalar_time ) override final;
    
    inline float __attribute__((always_inline)) getParticleComputationTimeVecto( const float log_particle_number );
    inline float __attribute__((always_inline)) getParticleComputationTimeScalar( const float log_particle_number );
    
private:

};//END class PartCompTime2D2Order

#endif
