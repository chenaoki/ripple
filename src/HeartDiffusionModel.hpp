#ifndef _HEARTMODELPROPERTY_H_
#define  _HEARTMODELPROPERTY_H_
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/fem/TimeData.hpp>

#include "HeartDiffusionFunctor.hpp"

namespace LifeV{

class HeartDiffusionModel : 
  public MeshData, 
  public TimeData
{

  private:
    
    bool         M_hasFibers;
   
    Real        M_diffusivity;
    Real        M_longitudinalConductivity;
    Real        M_transversalConductivity;
    Real        M_membraneCapacitance;
    Real        M_volumeSurfaceRatio;
    Real        M_conductivityRatio;

    std::string M_fibersDirectory;
    std::string M_fibersFile;
    std::string M_uOrder;

public:
    HeartDiffusionModel();
    HeartDiffusionModel( const GetPot& dataFile );
    HeartDiffusionModel ( const HeartDiffusionModel& model );
    HeartDiffusionModel& operator= ( const HeartDiffusionModel& model );
    virtual    ~HeartDiffusionModel() {};

public:
    void                      setup ( const GetPot& dataFile );

public:
    bool                      hasFibers()                           const{return M_hasFibers;};
    const Real&         volumeSurfaceRatio()           const{return M_volumeSurfaceRatio;};
    const Real&         longitudinalConductivity()      const{return M_longitudinalConductivity;};
    const Real&         transversalConductivity()     const{return M_transversalConductivity;};
    const Real&         membraneCapacitance()      const{return M_membraneCapacitance;};
    const Real&          diffusivity()                            const{return M_diffusivity;};
    const Real&          conductivityRatio()               const{return M_conductivityRatio;};
    std::string             fibersFile()                             const{return M_fibersFile;};
    std::string             fibersDirectory()                             const{return M_fibersDirectory;};
    std::string             uOrder()                                const{return M_uOrder;};
};

} // namespace LifeV
#endif
