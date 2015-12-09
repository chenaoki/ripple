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
  @brief Cardiac Electrophysiology Test
  @author Lucia Mirabella <lucia.mirabella@mail.polimi.it> and Mauro Perego <mauro.perego@polimi.it>
  @date 11-2007
  @contributors Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
  @last update 11-2010
 */

#ifndef __HEART_H
#define __HEART_H

#define MONODOMAIN

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>

#ifdef MONODOMAIN
#include <lifev/heart/solver/HeartMonodomainSolver.hpp>
#else
#include <lifev/heart/solver/HeartBidomainSolver.hpp>
#endif
#include <lifev/heart/solver/HeartIonicSolver.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <stdexcept>

namespace LifeV
{
/*!
  \class Heart

  3D Action potential propagation class


 */

class HeartException : public std::domain_error {
public:
    HeartException(const std::string& cause) :
    std::domain_error("Heart : " + cause){}
};

class Heart
{
public:

    typedef RegionMesh<LinearTetra>                                mesh_Type;
    typedef HeartIonicSolver< mesh_Type >                       ionicSolver_Type;
    typedef FESpace< mesh_Type, MapEpetra >               feSpace_Type;
#ifdef MONODOMAIN
    typedef HeartMonodomainData                                      elecData_Type;
    typedef HeartMonodomainSolver< mesh_Type >       elecSolver_Type;
#elif BIDOMAIN
    typedef HeartBidomainData                                             elecData_Type;
    typedef HeartBidomainSolver< mesh_Type >              elecSolver_Type;
#endif
    typedef elecSolver_Type::vector_Type                          vector_Type;
    typedef elecSolver_Type::matrix_Type                          matrix_Type;

    // pointer types
    typedef boost::shared_ptr< elecSolver_Type >          elecSolverPtr_Type;
    typedef boost::shared_ptr< vector_Type >                  vectorPtr_Type;
    typedef boost::shared_ptr< matrix_Type >                  matrixPtr_Type;
    typedef boost::shared_ptr< feSpace_Type >              feSpacePtr_Type;
    typedef boost::shared_ptr< ionicSolver_Type >         ionicSolverPtr_Type;

private:

    /*Object members*/
    elecData_Type M_data;
    HeartIonicData M_dataIonic;
    MeshData M_meshData;
    BCFunctionBase M_uZero;
    BCHandler M_bcH;

    /*Pointer members*/
    boost::shared_ptr< HeartFunctors > M_heartFctPtr;
    boost::shared_ptr< Exporter<mesh_Type > > M_exporterPtr;
    elecSolverPtr_Type M_electricSolverPtr;
    ionicSolverPtr_Type M_ionicModelPtr;
    feSpacePtr_Type M_uFESpacePtr;
#if BIDOMAIN
    feSpacePtr_Type M_FESpacePtr;
#endif

    /* current buffers */
    vectorPtr_Type M_Uptr;
    vectorPtr_Type M_Ueptr;
    vectorPtr_Type M_Fptr;

public:
    Heart ( GetPot& dataFile );
    virtual ~Heart();
    void step();
    void computeRhs ( vector_Type& rhs );

private:
    static Real zero_scalar_function ( const Real& /* t */,
                     const Real& /* x */,
                     const Real& /* y */,
                     const Real& /* z */,
                     const ID& /* i */ ){ return 0.; };


};

}
#endif /* __HEART_H */
