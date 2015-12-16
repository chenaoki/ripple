//@HEADER
/*
*******************************************************************************
Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University
This file is part of LifeV.
LifeV is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
LifeV is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License
along with LifeV.  If not, see <http://www.gnu.org/licenses/>.
*******************************************************************************
*/
//@HEADER
/*!
@file
@brief File containing a class for handling Monodomain data with GetPot
@date 11−2007
@author Lucia Mirabella <lucia.mirabella@gmail.com>, Mauro Perego <perego.mauro@gmail.com>
@contributor Simone Rossi <simone.rossi@epfl.ch>, Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
@mantainer Simone Rossi <simone.rossi@epfl.ch>
*/
#ifndef _HEARTMODELPROPERTY_H_
#define  _HEARTMODELPROPERTY_H_
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/fem/TimeData.hpp>

#include "HeartDiffusionFunctor.hpp"
namespace LifeV{
/*!
\class DataMonodomain
Base class which holds usual data for the Monodomain model solvers
*/
class HeartDiffusionModel : public MeshData, public TimeData
{
public:
    enum EnumIonConductivity{
        RM = 1,
        LR
    };

private:
    // Electrical properties
    Real        M_diffusivity;
    Real        M_longitudinalConductivity;
    Real        M_transversalConductivity;
    Real        M_membraneCapacitance;
    Real        M_volumeSurfaceRatio;
    Real        M_conductivityRatio;

    // Finite element properties
    std::string M_uOrder;

    // Fiber properties
    bool         M_hasFibers;
    std::string M_fibersDirectory;
    std::string M_fibersFile;

public:
    HeartDiffusionModel();
    HeartDiffusionModel( const GetPot& dataFile );
    HeartDiffusionModel ( const HeartDiffusionModel& model );
    HeartDiffusionModel& operator= ( const HeartDiffusionModel& model );
    virtual    ~HeartDiffusionModel() {}

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
