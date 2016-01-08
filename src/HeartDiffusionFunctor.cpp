#include "HeartDiffusionFunctor.hpp"
#include <boost/shared_ptr.hpp>

using namespace LifeV;
using namespace ripple;
// ===================================================
// Set Methods
// ===================================================

HeartDiffusionFunctor::HeartDiffusionFunctor(  ):
    M_dataFile(*(boost::shared_ptr<GetPot>( new GetPot() ) ) ),
    M_stimulusMode(""),
    M_restPotential(-84.)
{}

HeartDiffusionFunctor::HeartDiffusionFunctor ( GetPot& dataFile ):
    M_dataFile(dataFile),
    M_stimulusMode(dataFile("stim/mode", "S1S2")),
    M_restPotential(dataFile ("electric/physics/u0", -84.))
{
    if(M_stimulusMode == "S1S2"){

      Int S1_count(dataFile("S1_count", 1));
      Real S1_start(dataFile("S1_start", 100));
      Real S1_duration(dataFile("S1_duration", 5));
      Real S1_interval(dataFile("S1_interval", 400));
      Real S1_distance(dataFile("S1_distance", 0.5));
      Real S1_current(dataFile("S1_current", 20.0));
      Real S2_start(dataFile("S2_start", 100));
      Real S2_duration(dataFile("S2_duration", 5));
      Real S2_distance(dataFile("S2_distance", 0.5));
      Real S2_current(dataFile("S2_current", 25.0));

      axisElecPtrType S1(new axisElecType());
      axisElecPtrType S2(new axisElecType());
/*
      S1->addAxis(0, Real(0.0), S1_distance);
      S2->addAxis(1, Real(0.0), S2_distance);
      for(Int i = 0; i < S1_count; i++){
        for(Real t = 0; t < S1_duration; t += 1.0 ){
          S1->setCurrent(
              static_cast<unsigned int>( t + S1_start + i * S1_interval),
              S1_current
          );
        }
      }
      for(Real t = 0; t < S2_duration; t += 1.0 ){
        S2->setCurrent(
            static_cast<unsigned int>( t + S2_start ),
            S2_current
        );
      }
      this->vecPtrStimElectrode.push_back(static_cast<elecPtrType>(S1));
      this->vecPtrStimElectrode.push_back(static_cast<elecPtrType>(S2));
*/
    } // S1S2 setting
}

Real
HeartDiffusionFunctor::setStimulus ( const Real& t, const Real& x, const Real& y, const Real& z, const ID&   id) const
{
    Real ret = 0.;
    /*
    for( std::vector<elecPtrType>::const_iterator it = this->vecPtrStimElectrode.begin();
        it != vecPtrStimElectrode.end(); it++)
    {
      //ret += (*it)->getCurrent(static_cast<unsigned int>(t), x, y, z);
    }
    */
    return ret;
}

Real
HeartDiffusionFunctor::setInitialScalar (
    const Real& /*t*/,
    const Real& /*x*/,
    const Real& /*y*/,
    const Real& /*z*/,
    const ID&   /*id*/ )
{
    return M_restPotential;
}



Real HeartDiffusionFunctor::setZeroScalar (
    const Real& /*t*/,
    const Real& /*x*/,
    const Real& /*y*/,
    const Real& /*z*/,
    const ID&   /*id*/ )
{
    return 0.;
}

// ===================================================
// Get Methods
// ===================================================

HeartDiffusionFunctor::funcType
HeartDiffusionFunctor::stimulus()
{
    funcType f;
    f =  boost::bind (&HeartDiffusionFunctor::setStimulus, this, _1, _2, _3, _4, _5);
    return f;
}

HeartDiffusionFunctor::funcType
HeartDiffusionFunctor::initialScalar()
{
    funcType f;
    f =  boost::bind (&HeartDiffusionFunctor::setInitialScalar, this, _1, _2, _3, _4, _5);
     return f;
}

HeartDiffusionFunctor::funcType
HeartDiffusionFunctor::zeroScalar()
{
    funcType f;
    f =  boost::bind (&HeartDiffusionFunctor::setZeroScalar, this, _1, _2, _3, _4, _5);
    return f;
}
