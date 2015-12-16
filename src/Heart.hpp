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

#ifndef __HEART_H
#define __HEART_H

//#define BIDOMAIN

#include <stdexcept>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>

#include "HeartDiffusionSolver.hpp"
#include "HeartIonicSolver.hpp"

#ifdef BIDOMAIN
#include <lifev/heart/solver/HeartBidomainSolver.hpp>
#endif

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
    typedef HeartDiffusionModel                                          elecModel_Type;

    // typedef LifeV::LinearSolver AlgoSolver;
    typedef LifeV::SolverAztecOO AlgoSolver;
    typedef FESpace< mesh_Type, MapEpetra >               feSpace_Type;
    typedef HeartIonicSolver< mesh_Type, AlgoSolver>    ionicSolver_Type;
    typedef HeartDiffusionSolver< elecModel_Type, mesh_Type, AlgoSolver>  elecSolver_Type;

    //typedef HeartMonodomainSolver< RegionMesh<LinearTetra> >::vector_Type   vector_Type;
    //typedef HeartMonodomainSolver<RegionMesh<LinearTetra> >::matrix_Type        matrix_Type;
    typedef elecSolver_Type::vector_Type                         vector_Type;
    typedef elecSolver_Type::matrix_Type                         matrix_Type;

    // pointer types

#ifdef BIDOMAIN
    typedef HeartBidomainData                                            elecData_Type;
    typedef HeartBidomainSolver< mesh_Type >                elecSolver_Type;
#endif

private:
    /*Object members*/
    elecModel_Type     M_elecModel;
    HeartIonicData       M_ionData;
    MeshData               M_meshData;
    BCFunctionBase    M_uZero;
    BCHandler              M_bcH;

    /* Shared pointer managed members*/
    boost::shared_ptr< Epetra_Comm >                 M_commPtr;
    boost::shared_ptr< HeartDiffusionFunctor >     M_elecFctrPtr;
    boost::shared_ptr< elecSolver_Type >             M_elecSolverPtr;
    boost::shared_ptr< ionicSolver_Type >            M_ionicSolverPtr;
    boost::shared_ptr< feSpace_Type >                 M_uFESpacePtr;
    boost::shared_ptr< vector_Type >                 M_Uptr;
    boost::shared_ptr< vector_Type >                 M_Fptr;
    boost::shared_ptr< Exporter<mesh_Type > >  M_exporterPtr;

#ifdef BIDOMAIN
    boost::shared_ptr< feSpace_Type>                  M_FESpacePtr;
    boost::shared_ptr< vector_Type >                    M_Ueptr;
#endif


public:
    Heart ();
    Heart ( GetPot& dataFile );
    virtual ~Heart();

public:
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
