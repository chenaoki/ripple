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
#include "HeartBidomainSolver.hpp"
#include "HeartBidomainModel.hpp"
#endif

namespace LifeV
{

class HeartException : public std::domain_error {
public:
    HeartException(const std::string& cause) :
    std::domain_error("Heart : " + cause){}
};

class Heart
{
public:

    typedef RegionMesh<LinearTetra>                                        mesh_Type;
    // typedef LifeV::LinearSolver                                            AlgoSolver;
    typedef LifeV::SolverAztecOO                                           AlgoSolver;
    typedef FESpace< mesh_Type, MapEpetra >                                feSpace_Type;
    typedef HeartIonicSolver< mesh_Type, AlgoSolver>                       ionicSolver_Type;
#ifdef BIDOMAIN
    typedef HeartBidomainModel                                             elecModel_Type;
    typedef HeartBidomainSolver< mesh_Type, AlgoSolver >                   elecSolver_Type;
#else
    typedef HeartDiffusionModel                                            elecModel_Type;
    typedef HeartDiffusionSolver< elecModel_Type, mesh_Type, AlgoSolver>   elecSolver_Type;
#endif
    typedef elecSolver_Type::vector_Type                                   vector_Type;
    typedef elecSolver_Type::matrix_Type                                   matrix_Type;

    // pointer types

private:
    /*Object members*/
    elecModel_Type       M_elecModel;
    HeartIonicData       M_ionData;
    MeshData             M_meshData;
    BCFunctionBase       M_uZero;
    BCHandler            M_bcH;

    /* Shared pointer managed members*/
    boost::shared_ptr< Epetra_Comm >                 M_commPtr;
    boost::shared_ptr< HeartDiffusionFunctor >       M_elecFctrPtr;
    boost::shared_ptr< elecSolver_Type >             M_elecSolverPtr;
    boost::shared_ptr< ionicSolver_Type >            M_ionicSolverPtr;
    boost::shared_ptr< feSpace_Type >                M_uFESpacePtr;

public:
    boost::shared_ptr< vector_Type >                 M_Uptr;
    boost::shared_ptr< vector_Type >                 M_Fptr;
    boost::shared_ptr< Exporter<mesh_Type > >  M_exporterPtr;

#ifdef BIDOMAIN
    boost::shared_ptr< feSpace_Type>                  M_FESpacePtr;
    boost::shared_ptr< vector_Type >                  M_Ueptr;
#endif


public:
    Heart ();
    Heart ( GetPot& dataFile );
    virtual ~Heart();

public:
    void step();
    void stimulate(unsigned int idElectrode);
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
