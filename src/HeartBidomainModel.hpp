#ifndef _DATABIDOMAIN_H_
#define _DATABIDOMAIN_H_
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/fem/TimeData.hpp>

#include "HeartDiffusionFunctor.hpp"

namespace LifeV{

class HeartBidomainModel:
    public MeshData,
    public TimeData
{

private:

  bool        M_fibersFormat;
  bool        M_hasFibers;


  Int         M_heartDiffusionFactor;

  Real        M_externalDiffusivity;
  Real        M_internalDiffusivity;
  Real        M_longitudinalExternalConductivity;
  Real        M_longitudinalInternalConductivity;
  Real        M_membraneCapacitance;
  Real        M_transversalExternalConductivity;
  Real        M_transversalInternalConductivity;
  Real        M_volumeSurfaceRatio;

  std::string M_fibersDirectory;
  std::string M_fibersFile;
  std::string M_uOrder;
  
public:

  HeartBidomainModel();
  HeartBidomainModel ( const GetPot& dataFile );
  HeartBidomainModel ( const HeartBidomainModel& dataBidomain );
  virtual ~HeartBidomainModel() {};
  HeartBidomainModel& operator= ( const HeartBidomainModel& dataBidomain );

public:
    
    void showMe ( std::ostream& output = std::cout );
    void setup ( const GetPot& dataFile );
    
    //---------------
    // Getters
    //---------------
    
    //! FE space order
    std::string             uOrder()                         const
    {
        return M_uOrder;
    }

    //! Chi
    const Real&             volumeSurfaceRatio()             const
    {
        return M_volumeSurfaceRatio;
    }
    //! fiber File
    std::string             fibersFile()                     const
    {
        return M_fibersFile;
    }

    const Int&              heartDiffusionFactor()           const
    {
        return M_heartDiffusionFactor;
    }

    const bool&             hasFibers()                      const
    {
        return M_hasFibers;
    }

    //! format vct
    const bool&             fibersFormat()                   const
    {
        return M_fibersFormat;
    }

    //! sigma_l
    const Real&             longitudinalInternalConductivity() const
    {
        return M_longitudinalInternalConductivity;
    }
    const Real&             longitudinalExternalConductivity() const
    {
        return M_longitudinalExternalConductivity;
    }

    //! sigma_t
    const Real&             transversalInternalConductivity() const
    {
        return M_transversalInternalConductivity;
    }
    const Real&             transversalExternalConductivity() const
    {
        return M_transversalExternalConductivity;
    }

    //! Cm
    const Real&             membraneCapacitance()             const
    {
        return M_membraneCapacitance;
    }
    //! D
    const Real&             internalDiffusivity()             const
    {
        return M_internalDiffusivity;
    }
    //! Post_dir
    const Real&             externalDiffusivity()             const
    {
        return M_externalDiffusivity;
    }

};

}
#endif
