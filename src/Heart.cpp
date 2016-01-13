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

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <opencv/cv.hpp>
#include <opencv/highgui.h>

#include "Heart.hpp"

using namespace LifeV;

//! Identifiers for heart boundaries
const Int EPICARDIUM    = 40;
const Int ENDOCARDIUM   = 60;
const Int TRUNC_SEC     = 50;


// ===================================================
//! Constructors
// ===================================================

Heart::Heart ( GetPot& dataFile ):
M_uZero ( Heart::zero_scalar_function ),
M_bcH(),
M_meshData(dataFile, "electric/space_discretization"),
M_elecModel(dataFile),
M_ionData(dataFile),
M_commPtr(new Epetra_MpiComm ( MPI_COMM_WORLD )),
M_elecFctrPtr ( new HeartDiffusionFunctor ( dataFile )),
M_elecSolverPtr(NULL),
M_ionicSolverPtr(NULL),
M_uFESpacePtr ( NULL ),
M_Uptr(NULL),
M_Fptr(NULL)
#ifdef BIDOMAIN
, M_FESpacePtr ( NULL ),
M_Ueptr(NULL),
#endif
{

    // check params
   std::string uOrder =  dataFile ( "electric/space_discretization/u_order", "P1");
    if ( uOrder.compare ("P1") != 0 ){
        throw HeartException("Heart exception");
    }
    bool verbose = (M_commPtr->MyPID() == 0);

    //! Boundary conditions handler and function
    M_bcH.addBC ( "Endo",      ENDOCARDIUM,    Natural,    Full,   M_uZero,  1 );
    M_bcH.addBC ( "Epi",       EPICARDIUM,     Natural,    Full,   M_uZero,  1 );
    M_bcH.addBC ( "Trunc",     TRUNC_SEC,      Natural,    Full,   M_uZero,  1 );

    //--------------------
    // Object constructions
    //--------------------

    if (!M_commPtr->MyPID() ){
        std::cout << "My PID = " << M_commPtr->MyPID() << std::endl;
    }

    //! Local Mesh
    boost::shared_ptr<mesh_Type> localMeshPtr;
    boost::shared_ptr<mesh_Type > fullMeshPtr ( new mesh_Type ( M_commPtr ) );
    readMesh (*fullMeshPtr, M_meshData);
    MeshPartitioner< mesh_Type >   meshPart (fullMeshPtr, M_commPtr);
    localMeshPtr = meshPart.meshPartition();

     // FE space
    M_uFESpacePtr.reset(new feSpace_Type (
        localMeshPtr,
        feTetraP1,
        quadRuleTetra15pt,
        quadRuleTria3pt,
        1,
        M_commPtr)
    );

     // Electric Solver
    M_elecSolverPtr.reset( new HeartDiffusionSolver< HeartDiffusionModel, RegionMesh<LinearTetra > >(
        M_elecModel,
        *M_uFESpacePtr,
        M_bcH,
        M_commPtr)
    );

    switch ( dataFile ("electric/physics/ion_solver", 1)  ){
    case 1:
        M_ionicSolverPtr.reset (new RogersMcCulloch< mesh_Type > (
            M_ionData,
            *localMeshPtr,
            *M_uFESpacePtr,
            *M_commPtr) );
        break;
    case 2:
        M_ionicSolverPtr.reset (new LuoRudy< mesh_Type > (
            M_ionData,
            *localMeshPtr,
            *M_uFESpacePtr,
            *M_commPtr) );
        break;
    case 3:
        M_ionicSolverPtr.reset (new MitchellSchaeffer< mesh_Type > (
            M_ionData,
            *localMeshPtr,
            *M_uFESpacePtr,
            *M_commPtr) );
        break;
    default:
        throw HeartException("electric/physics/ion_model value not supported");
        break;
    }

    // Current state buffers construction
    M_Uptr.reset( new vector_Type (M_elecSolverPtr->solutionTransmembranePotential(), Repeated ));
#ifdef BIDOMAIN
    M_Ueptr.reset( new vector_Type (M_elecSolverPtr->solutionExtraPotential(), Repeated ));
#endif
    M_Fptr.reset( new vector_Type (M_elecSolverPtr->fiberVector(), Repeated ));

    /*-----------------------*/
    //  Initialization
    /*-----------------------*/

    // Elec solver
    M_elecSolverPtr->setup ( dataFile );
    M_elecSolverPtr->initialize ( M_elecFctrPtr->initialScalar() );
#ifdef BIDOMAIN
    M_elecSolverPtr->initialize ( M_elecFctrPtr->initialScalar(), M_elecFctrPtr->zeroScalar() );
#endif
    M_elecSolverPtr->buildSystem( );
    M_ionicSolverPtr->initialize( );

    //! Building time-independent part of the system
    MPI_Barrier (MPI_COMM_WORLD);
    M_elecModel.setTime(0);
    M_elecSolverPtr->resetPreconditioner();

    //! Exporter
    if ( dataFile ( "exporter/type", "ensight").compare ("hdf5") == 0) {
        M_exporterPtr.reset ( new ExporterHDF5<mesh_Type > ( dataFile, "heart" ) );
        M_exporterPtr->setPostDir ( "./" ); // This is a test to see if M_post_dir is working
        M_exporterPtr->setMeshProcId ( localMeshPtr, M_commPtr->MyPID() );
    }
    M_exporterPtr->addVariable ( ExporterData<mesh_Type>::ScalarField,  "potential", M_uFESpacePtr,
                            M_Uptr, UInt (0) );
#ifdef BIDOMAIN
    M_exporterPtr->addVariable ( ExporterData<mesh_Type>::ScalarField,  "potential_e", M_FESpacePtr,
                            M_Ueptr, UInt (0) );
#endif
    if (M_elecModel.hasFibers() )
        M_exporterPtr->addVariable ( ExporterData<mesh_Type>::VectorField,
                                "fibers",
                                M_uFESpacePtr,
                                M_Fptr,
                                UInt (0) );
    M_exporterPtr->postProcess ( 0 );

    MPI_Barrier (MPI_COMM_WORLD);

}

Heart::~Heart(void){
    //cvDestroyWindow("Heart");
}

// ===================================================
//! Methods
// ===================================================

void 
Heart::stimulate(unsigned int idElectrode)
{
  this->M_elecFctrPtr->stimulate(
    static_cast<unsigned int>(M_elecModel.time())+1,
    idElectrode
  );
}

void
Heart::step()
{
    static const Real dt     = M_elecModel.timeStep();
    static const Real tFinal = M_elecModel.endTime();
    static Real time(0.0);
    static Int iter = 0;

    vector_Type rhs ( M_elecSolverPtr->getMap() );
    Real normu;
    Real meanu;
    Real minu;

    iter++;
    time += dt;
    if( iter == 1){
        std::cout << " step , time(s), norm u, mean u, max u,  profile(s) " << std::endl;
    }
    if(time > tFinal + dt / 2. ){
        throw HeartException("finish");
    }

    MPI_Barrier (MPI_COMM_WORLD);

    // solve
    M_elecModel.setTime (time);
    M_ionicSolverPtr->solveIonicModel ( M_elecSolverPtr->solutionTransmembranePotential(), M_elecModel.timeStep() );
    rhs *= 0;
    computeRhs ( rhs );
    M_elecSolverPtr->updatePDESystem ( rhs );
    M_elecSolverPtr->PDEiterate ( M_bcH );

    // export
    *M_Uptr = M_elecSolverPtr->solutionTransmembranePotential();
#ifdef BIDOMAIN
    *M_Ueptr = M_elecSolverPtr->solutionExtraPotential();
#endif
    M_exporterPtr->postProcess ( time );

    // analyze
    normu = M_elecSolverPtr->solutionTransmembranePotential().norm2();
    M_elecSolverPtr->solutionTransmembranePotential().epetraVector().MeanValue (&meanu);
    M_elecSolverPtr->solutionTransmembranePotential().epetraVector().MaxValue (&minu);

    MPI_Barrier (MPI_COMM_WORLD);

    std::cout << iter << ",";
    std::cout << M_elecModel.time() << ",";
    std::cout << normu << ",";
    std::cout << meanu << ",";
    std::cout << minu << ",";
    std::cout << std::endl;

}


void
Heart::computeRhs ( vector_Type& rhs)
{
    bool verbose = (M_commPtr->MyPID() == 0);
    if (verbose)
    {
        std::cout << "  f-  Computing Rhs ...        " << "\n" << std::flush;
    }
    LifeChrono chrono;
    chrono.start();

    //! u, w with repeated map
    vector_Type uVecRep (M_elecSolverPtr->solutionTransmembranePotential(), Repeated);
    M_ionicSolverPtr->updateRepeated();
    VectorElemental elvec_Iapp ( M_elecSolverPtr->potentialFESpace().fe().nbFEDof(), 2 ),
                    elvec_u ( M_elecSolverPtr->potentialFESpace().fe().nbFEDof(), 1 ),
                    elvec_Iion ( M_elecSolverPtr->potentialFESpace().fe().nbFEDof(), 1 );

    for (UInt iVol = 0; iVol < M_elecSolverPtr->potentialFESpace().mesh()->numVolumes(); ++iVol)
    {
        M_elecSolverPtr->potentialFESpace().fe().updateJacQuadPt ( M_elecSolverPtr->potentialFESpace().mesh()->volumeList ( iVol ) );
        elvec_Iapp.zero();
        elvec_u.zero();
        elvec_Iion.zero();

        UInt eleIDu = M_elecSolverPtr->potentialFESpace().fe().currentLocalId();
        UInt nbNode = ( UInt ) M_elecSolverPtr->potentialFESpace().fe().nbFEDof();

        //! Filling local elvec_u with potential values in the nodes
        for ( UInt iNode = 0 ; iNode < nbNode ; iNode++ )
        {
            Int  ig = M_elecSolverPtr->potentialFESpace().dof().localToGlobalMap ( eleIDu, iNode );
            elvec_u.vec() [ iNode ] = uVecRep[ig];
        }

        M_ionicSolverPtr->updateElementSolution (eleIDu);
        M_ionicSolverPtr->computeIonicCurrent (M_elecModel.membraneCapacitance(), elvec_Iion, elvec_u, M_elecSolverPtr->potentialFESpace() );

        //! Computing the current source of the righthand side, repeated
        AssemblyElemental::source (M_elecFctrPtr->stimulus(),
                                   elvec_Iapp,
                                   M_elecSolverPtr->potentialFESpace().fe(),
                                   M_elecModel.time(),
                                   0);
        AssemblyElemental::source (M_elecFctrPtr->stimulus(),
                                   elvec_Iapp,
                                   M_elecSolverPtr->potentialFESpace().fe(),
                                   M_elecModel.time(),
                                   1);

        //! Assembling the righthand side
        for ( UInt i = 0 ; i < M_elecSolverPtr->potentialFESpace().fe().nbFEDof() ; i++ )
        {
            Int  ig = M_elecSolverPtr->potentialFESpace().dof().localToGlobalMap ( eleIDu, i );
            rhs.sumIntoGlobalValues (ig, (M_elecModel.conductivityRatio() * elvec_Iapp.vec() [i] +
                                          elvec_Iapp.vec() [i + nbNode]) /
                                     (1 + M_elecModel.conductivityRatio() ) + M_elecModel.volumeSurfaceRatio() * elvec_Iion.vec() [i] );
        }
    }
    rhs.globalAssemble();
    Real coeff = M_elecModel.volumeSurfaceRatio() * M_elecModel.membraneCapacitance() / M_elecModel.timeStep();
    vector_Type tmpvec (M_elecSolverPtr->solutionTransmembranePotential() );
    tmpvec *= coeff;
    rhs += M_elecSolverPtr->massMatrix() * tmpvec;
    MPI_Barrier (MPI_COMM_WORLD);
    chrono.stop();
    if (verbose)
    {
        std::cout << "done in " << chrono.diff() << " s." << std::endl;
    }
}

#ifdef BIDOMAIN
void
Heart::computeRhs ( vector_Type& rhs )
{
    bool verbose = (M_commPtr->MyPID() == 0);
    if (verbose)
    {
        std::cout << "  f-  Computing Rhs ...        " << "\n" << std::flush;
    }
    LifeChrono chrono;
    chrono.start();

    //! u, w with repeated map
    vector_Type uVecRep (M_elecSolverPtr->solutionTransmembranePotential(), Repeated);
    M_ionicSolverPtr->updateRepeated();

    VectorElemental elvec_Iapp ( M_elecSolverPtr->potentialFESpace().fe().nbFEDof(), 2 ),
                    elvec_u ( M_elecSolverPtr->potentialFESpace().fe().nbFEDof(), 1 ),
                    elvec_Iion ( M_elecSolverPtr->potentialFESpace().fe().nbFEDof(), 1 );
    for (UInt iVol = 0; iVol < M_elecSolverPtr->potentialFESpace().mesh()->numVolumes(); ++iVol)
    {
        M_elecSolverPtr->potentialFESpace().fe().updateJacQuadPt ( M_elecSolverPtr->potentialFESpace().mesh()->volumeList ( iVol ) );
        elvec_u.zero();
        elvec_Iion.zero();
        elvec_Iapp.zero();

        UInt eleIDu = M_elecSolverPtr->potentialFESpace().fe().currentLocalId();
        UInt nbNode = ( UInt ) M_elecSolverPtr->potentialFESpace().fe().nbFEDof();
        for ( UInt iNode = 0 ; iNode < nbNode ; iNode++ )
        {
            Int ig = M_elecSolverPtr->potentialFESpace().dof().localToGlobalMap ( eleIDu, iNode );
            elvec_u.vec() [ iNode ] = uVecRep[ig];
        }

        UInt eleID = M_elecSolverPtr->potentialFESpace().fe().currentLocalId();
        M_ionicSolverPtr->updateElementSolution (eleID);
        M_ionicSolverPtr->computeIonicCurrent (M_elecModel.membraneCapacitance(), elvec_Iion, elvec_u, M_elecSolverPtr->potentialFESpace() );

        //! Computing Iapp
        AssemblyElemental::source (M_elecFctrPtr->stimulus(),
                                   elvec_Iapp,
                                   M_elecSolverPtr->potentialFESpace().fe(),
                                   M_elecModel.time(), 0);
        AssemblyElemental::source (M_elecFctrPtr->stimulus(),
                                   elvec_Iapp,
                                   M_elecSolverPtr->potentialFESpace().fe(),
                                   M_elecModel.time(),
                                   1);
        UInt totalUDof  = M_elecSolverPtr->potentialFESpace().map().map (Unique)->NumGlobalElements();

        for ( UInt iNode = 0 ; iNode < nbNode ; iNode++ )
        {
            Int ig = M_elecSolverPtr->potentialFESpace().dof().localToGlobalMap ( eleIDu, iNode );
            rhs.sumIntoGlobalValues (ig, elvec_Iapp.vec() [iNode] +
                                     M_elecModel.volumeSurfaceRatio() * elvec_Iion.vec() [iNode] );
            rhs.sumIntoGlobalValues (ig + totalUDof,
                                     -elvec_Iapp.vec() [iNode + nbNode] -
                                     M_elecModel.volumeSurfaceRatio() * elvec_Iion.vec() [iNode] );
        }
    }
    rhs.globalAssemble();

    rhs += M_elecSolverPtr->matrMass() * M_elecModel.volumeSurfaceRatio() *
           M_elecModel.membraneCapacitance() * M_elecSolverPtr->BDFIntraExtraPotential().time_der (M_elecModel.timeStep() );

    MPI_Barrier (MPI_COMM_WORLD);

    chrono.stop();
    if (verbose)
    {
        std::cout << "done in " << chrono.diff() << " s." << std::endl;
    }
}
#endif
