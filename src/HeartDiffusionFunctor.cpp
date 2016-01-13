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
    M_stimCurrent(20.0),
    M_restPotential(-84.)
{}

HeartDiffusionFunctor::HeartDiffusionFunctor ( GetPot& dataFile ):
    M_dataFile(dataFile),
    M_stimulusMode(dataFile("stim/mode", "S1S2")),
    M_stimCurrent(dataFile("stim/array/current", 20.0)),
    M_restPotential(dataFile ("electric/physics/u0", -84.))
{
    // Array electordes
    Int elec_numX(dataFile("stim/array/numX", 3));
    Int elec_numY(dataFile("stim/array/numY", 3));
    Real elec_off(dataFile("stim/array/offset", 0.1));
    Real elec_int(dataFile("stim/array/interval", 0.1));
    Real elec_size(dataFile("stim/array/size", 0.05));
    
    for(int i=0; i<elec_numY;i++){
      for(int j=0; j<elec_numX; j++){
        axisElecPtrType elec(new axisElecType());
        elec->addAxis(0, Real(elec_off + j * elec_int), elec_size);
        elec->addAxis(1, Real(elec_off + i * elec_int), elec_size);
        this->vecPtrStimElectrode.push_back(static_cast<elecPtrType>(elec));
      }
    }
    
    if(M_stimulusMode == "S1S2"){

      Int S1_count(dataFile("stim/S1S2/S1_count", 1));
      Real S1_start(dataFile("stim/S1S2/S1_start", 100));
      Real S1_duration(dataFile("stim/S1S2/S1_duration", 5));
      Real S1_interval(dataFile("stim/S1S2/S1_interval", 400));
      Real S1_distance(dataFile("stim/S1S2/S1_distance", 0.5));
      Real S1_current(dataFile("Sstim/S1S2/S1_current", 20.0));
      Real S2_start(dataFile("stim/S1S2/S2_start", 100));
      Real S2_duration(dataFile("stim/S1S2/S2_duration", 5));
      Real S2_distance(dataFile("stim/S1S2/S2_distance", 0.5));
      Real S2_current(dataFile("stim/S1S2/S2_current", 25.0));

      axisElecPtrType S1(new axisElecType());
      axisElecPtrType S2(new axisElecType());

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
    } // S1S2 setting

    
}

void 
HeartDiffusionFunctor::stimulate(unsigned int time, unsigned int idElectrode)
{
  //Real current(20.0);
  if( idElectrode < vecPtrStimElectrode.size()){
    this->vecPtrStimElectrode[idElectrode]->setCurrent(
        time, 
        //current
        M_stimCurrent
    );
  }
}

Real
HeartDiffusionFunctor::setStimulus ( const Real& t, const Real& x, const Real& y, const Real& z, const ID&   id) const
{
    Real ret = 0.;
    for( std::vector<elecPtrType>::const_iterator it = this->vecPtrStimElectrode.begin();
        it != vecPtrStimElectrode.end(); it++)
    {
      ret += (*it)->getCurrent(static_cast<unsigned int>(t), x, y, z);
    }
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
