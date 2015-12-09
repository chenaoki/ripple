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
  @contributors Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>, Simone Rossi <simone.rossi@epfl.ch>
  @last update 11-2010
 */


#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <opencv/cv.hpp>
#include <opencv/highgui.h>

#include "heart.hpp"

using namespace LifeV;

//! Identifiers for heart boundaries
const Int EPICARDIUM    = 40;
const Int ENDOCARDIUM   = 60;
const Int TRUNC_SEC     = 50;


// ===================================================
//! Constructors
// ===================================================

Heart::Heart ( GetPot& dataFile ):
M_heartFctPtr ( NULL ),
M_ionicModelPtr(NULL),
M_uFESpacePtr ( NULL ),
#ifdef BIDOMAIN
M_FESpacePtr ( NULL ),
#endif
M_Uptr(NULL),
M_Ueptr(NULL),
M_Fptr(NULL),
M_uZero ( Heart::zero_scalar_function ),
M_bcH(),
M_meshData(),
M_data(),
M_dataIonic(dataFile)
{

    /*-- Initialization --*/
    LifeChrono chronoinitialsettings;
    chronoinitialsettings.start();

    M_heartFctPtr.reset (new HeartFunctors ( dataFile) );
    //! Pointer to access functors
    M_heartFctPtr->M_comm.reset (new Epetra_MpiComm ( MPI_COMM_WORLD ) );
    if (!M_heartFctPtr->M_comm->MyPID() )
    {
        std::cout << "My PID = " << M_heartFctPtr->M_comm->MyPID() << std::endl;
    }

    // check params
   std::string uOrder =  M_heartFctPtr->M_dataFile ( "electric/space_discretization/u_order", "P1");
    if ( uOrder.compare ("P1") != 0 ){
        throw HeartException("Heart exception");
    }
    bool verbose = (M_heartFctPtr->M_comm->MyPID() == 0);


    //--------------------
    // Constructions
    //--------------------
    //! Boundary conditions handler and function
    M_bcH.addBC ( "Endo",      ENDOCARDIUM,    Natural,    Full,   M_uZero,  1 );
    M_bcH.addBC ( "Epi",       EPICARDIUM,     Natural,    Full,   M_uZero,  1 );
    M_bcH.addBC ( "Trunc",     TRUNC_SEC,      Natural,    Full,   M_uZero,  1 );

    // Setup electric model data
    M_data.setup(dataFile, M_heartFctPtr);

    //! Local Mesh construction
    boost::shared_ptr<mesh_Type> localMeshPtr;
    M_meshData.setup (dataFile, "electric/space_discretization");
    boost::shared_ptr<mesh_Type > fullMeshPtr ( new mesh_Type ( M_heartFctPtr->M_comm ) );
    readMesh (*fullMeshPtr, M_meshData);
    MeshPartitioner< mesh_Type >   meshPart (fullMeshPtr, M_heartFctPtr->M_comm);
    localMeshPtr = meshPart.meshPartition();

     // Finie element space construction
    M_uFESpacePtr.reset(new feSpace_Type ( localMeshPtr, feTetraP1, quadRuleTetra15pt, quadRuleTria3pt,1,M_heartFctPtr->M_comm));
#if BIDOMAIN
    M_FESpacePtr.reset(new feSpace_Type (localMeshPtr, feTetraP1,quadRuleTetra15pt,quadRuleTria3pt,2,M_heartFctPtr->M_comm));
#endif

     // Solver construction
#ifdef MONODOMAIN
    M_electricSolverPtr.reset( new elecSolver_Type(M_data, *M_uFESpacePtr, M_bcH, M_heartFctPtr->M_comm));
#elif BIDOMAIN
    M_electricSolverPtr.reset(new elecSolver_Type(M_data, *M_FESpacePtr, *M_uFESpacePtr, M_bcH, M_heartFctPtr->M_comm));
#endif

    switch ( dataFile ("electric/physics/ion_model", 1)  ){
    case 1:
        M_ionicModelPtr.reset (new RogersMcCulloch< mesh_Type > (M_dataIonic,
                                                            *localMeshPtr,
                                                            *M_uFESpacePtr,
                                                            *M_heartFctPtr->M_comm) );
        break;
    case 2:
        M_ionicModelPtr.reset (new LuoRudy< mesh_Type > (M_dataIonic,
                                                    *localMeshPtr,
                                                    *M_uFESpacePtr,
                                                    *M_heartFctPtr->M_comm) );
        break;
    case 3:
        M_ionicModelPtr.reset (new MitchellSchaeffer< mesh_Type > (M_dataIonic,
                                                              *localMeshPtr,
                                                              *M_uFESpacePtr,
                                                              *M_heartFctPtr->M_comm) );
        break;
    default:
        throw HeartException("electric/physics/ion_model value not supported");
        break;
    }

    // Current state buffers construction
    M_Uptr.reset( new vector_Type (M_electricSolverPtr->solutionTransmembranePotential(), Repeated ));
#ifdef BIDOMAIN
    M_Ueptr.reset( new vector_Type (M_electricSolverPtr->solutionExtraPotential(), Repeated ));
#endif
    M_Fptr.reset( new vector_Type (M_electricSolverPtr->fiberVector(), Repeated ));

    /*-----------------------*/
    //  Setup
    /*-----------------------*/

    M_electricSolverPtr->setup ( M_heartFctPtr->M_dataFile );
#ifdef MONODOMAIN
    M_electricSolverPtr->initialize ( M_heartFctPtr->initialScalar() );
#elif BIDOMAIN
    M_electricSolverPtr->initialize ( M_heartFctPtr->initialScalar(), M_heartFctPtr->zeroScalar() );
#endif
    M_ionicModelPtr->initialize( );

    //! Building time-independent part of the system
    M_electricSolverPtr->buildSystem( );
    MPI_Barrier (MPI_COMM_WORLD);
    M_data.setTime(0);
    M_electricSolverPtr->resetPreconditioner();

    //! Setting generic Exporter postprocessing
    if (M_heartFctPtr->M_dataFile ( "exporter/type", "ensight").compare ("hdf5") == 0) {
        M_exporterPtr.reset ( new ExporterHDF5<mesh_Type > ( M_heartFctPtr->M_dataFile, "heart" ) );
        M_exporterPtr->setPostDir ( "./" ); // This is a test to see if M_post_dir is working
        M_exporterPtr->setMeshProcId ( localMeshPtr, M_heartFctPtr->M_comm->MyPID() );
    }
    M_exporterPtr->addVariable ( ExporterData<mesh_Type>::ScalarField,  "potential", M_uFESpacePtr,
                            M_Uptr, UInt (0) );
#ifdef BIDOMAIN
    M_exporterPtr->addVariable ( ExporterData<mesh_Type>::ScalarField,  "potential_e", M_FESpacePtr,
                            M_Ueptr, UInt (0) );
#endif

    if (M_data.hasFibers() )
        M_exporterPtr->addVariable ( ExporterData<mesh_Type>::VectorField,
                                "fibers",
                                M_uFESpacePtr,
                                M_Fptr,
                                UInt (0) );
    M_exporterPtr->postProcess ( 0 );

    MPI_Barrier (MPI_COMM_WORLD);
    chronoinitialsettings.stop();

    //cv::namedWindow("Heart", cv::WINDOW_AUTOSIZE);

}

Heart::~Heart(void){
    //cvDestroyWindow("Heart");
}

// ===================================================
//! Methods
// ===================================================

void
Heart::step()
{

    /*
    cv::Mat img = cv::imread("lena.png", CV_LOAD_IMAGE_COLOR);
    cv::imshow("Heart", img);
    cv::waitKey(3);
    */

    static const Real dt     = M_data.timeStep();
    static const Real tFinal = M_data.endTime();
    static Real time(0.0);
    static Int iter = 0;

    LifeChrono chrono;
    vector_Type rhs ( M_electricSolverPtr->getMap() );
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

    chrono.start();
    MPI_Barrier (MPI_COMM_WORLD);

    // solve
    M_data.setTime (time);
    M_ionicModelPtr->solveIonicModel ( M_electricSolverPtr->solutionTransmembranePotential(), M_data.timeStep() );
    rhs *= 0;
    computeRhs ( rhs );
    M_electricSolverPtr->updatePDESystem ( rhs );
    M_electricSolverPtr->PDEiterate ( M_bcH );

    // export
    *M_Uptr = M_electricSolverPtr->solutionTransmembranePotential();
#ifdef BIDOMAIN
    *M_Ueptr = M_electricSolverPtr->solutionExtraPotential();
#endif
    M_exporterPtr->postProcess ( time );

    // analyze
    normu = M_electricSolverPtr->solutionTransmembranePotential().norm2();
    M_electricSolverPtr->solutionTransmembranePotential().epetraVector().MeanValue (&meanu);
    M_electricSolverPtr->solutionTransmembranePotential().epetraVector().MaxValue (&minu);

    MPI_Barrier (MPI_COMM_WORLD);
    chrono.stop();

    std::cout << iter << ",";
    std::cout << M_data.time() << ",";
    std::cout << normu << ",";
    std::cout << meanu << ",";
    std::cout << minu << ",";
    std::cout << chrono.diff();
    std::cout << std::endl;

}


#ifdef MONODOMAIN
void
Heart::computeRhs ( vector_Type& rhs)
{
    bool verbose = (M_heartFctPtr->M_comm->MyPID() == 0);
    if (verbose)
    {
        std::cout << "  f-  Computing Rhs ...        " << "\n" << std::flush;
    }
    LifeChrono chrono;
    chrono.start();

    //! u, w with repeated map
    vector_Type uVecRep (M_electricSolverPtr->solutionTransmembranePotential(), Repeated);
    M_ionicModelPtr->updateRepeated();
    VectorElemental elvec_Iapp ( M_electricSolverPtr->potentialFESpace().fe().nbFEDof(), 2 ),
                    elvec_u ( M_electricSolverPtr->potentialFESpace().fe().nbFEDof(), 1 ),
                    elvec_Iion ( M_electricSolverPtr->potentialFESpace().fe().nbFEDof(), 1 );

    for (UInt iVol = 0; iVol < M_electricSolverPtr->potentialFESpace().mesh()->numVolumes(); ++iVol)
    {
        M_electricSolverPtr->potentialFESpace().fe().updateJacQuadPt ( M_electricSolverPtr->potentialFESpace().mesh()->volumeList ( iVol ) );
        elvec_Iapp.zero();
        elvec_u.zero();
        elvec_Iion.zero();

        UInt eleIDu = M_electricSolverPtr->potentialFESpace().fe().currentLocalId();
        UInt nbNode = ( UInt ) M_electricSolverPtr->potentialFESpace().fe().nbFEDof();

        //! Filling local elvec_u with potential values in the nodes
        for ( UInt iNode = 0 ; iNode < nbNode ; iNode++ )
        {
            Int  ig = M_electricSolverPtr->potentialFESpace().dof().localToGlobalMap ( eleIDu, iNode );
            elvec_u.vec() [ iNode ] = uVecRep[ig];
        }

        M_ionicModelPtr->updateElementSolution (eleIDu);
        M_ionicModelPtr->computeIonicCurrent (M_data.membraneCapacitance(), elvec_Iion, elvec_u, M_electricSolverPtr->potentialFESpace() );

        //! Computing the current source of the righthand side, repeated
        AssemblyElemental::source (M_heartFctPtr->stimulus(),
                                   elvec_Iapp,
                                   M_electricSolverPtr->potentialFESpace().fe(),
                                   M_data.time(),
                                   0);
        AssemblyElemental::source (M_heartFctPtr->stimulus(),
                                   elvec_Iapp,
                                   M_electricSolverPtr->potentialFESpace().fe(),
                                   M_data.time(),
                                   1);

        //! Assembling the righthand side
        for ( UInt i = 0 ; i < M_electricSolverPtr->potentialFESpace().fe().nbFEDof() ; i++ )
        {
            Int  ig = M_electricSolverPtr->potentialFESpace().dof().localToGlobalMap ( eleIDu, i );
            rhs.sumIntoGlobalValues (ig, (M_data.conductivityRatio() * elvec_Iapp.vec() [i] +
                                          elvec_Iapp.vec() [i + nbNode]) /
                                     (1 + M_data.conductivityRatio() ) + M_data.volumeSurfaceRatio() * elvec_Iion.vec() [i] );
        }
    }
    rhs.globalAssemble();
    Real coeff = M_data.volumeSurfaceRatio() * M_data.membraneCapacitance() / M_data.timeStep();
    vector_Type tmpvec (M_electricSolverPtr->solutionTransmembranePotential() );
    tmpvec *= coeff;
    rhs += M_electricSolverPtr->massMatrix() * tmpvec;
    MPI_Barrier (MPI_COMM_WORLD);
    chrono.stop();
    if (verbose)
    {
        std::cout << "done in " << chrono.diff() << " s." << std::endl;
    }
}
#else
void
Heart::computeRhs ( vector_Type& rhs )
{
    bool verbose = (M_heartFctPtr->M_comm->MyPID() == 0);
    if (verbose)
    {
        std::cout << "  f-  Computing Rhs ...        " << "\n" << std::flush;
    }
    LifeChrono chrono;
    chrono.start();

    //! u, w with repeated map
    vector_Type uVecRep (M_electricSolverPtr->solutionTransmembranePotential(), Repeated);
    M_ionicModelPtr->updateRepeated();

    VectorElemental elvec_Iapp ( M_electricSolverPtr->potentialFESpace().fe().nbFEDof(), 2 ),
                    elvec_u ( M_electricSolverPtr->potentialFESpace().fe().nbFEDof(), 1 ),
                    elvec_Iion ( M_electricSolverPtr->potentialFESpace().fe().nbFEDof(), 1 );
    for (UInt iVol = 0; iVol < M_electricSolverPtr->potentialFESpace().mesh()->numVolumes(); ++iVol)
    {
        M_electricSolverPtr->potentialFESpace().fe().updateJacQuadPt ( M_electricSolverPtr->potentialFESpace().mesh()->volumeList ( iVol ) );
        elvec_u.zero();
        elvec_Iion.zero();
        elvec_Iapp.zero();

        UInt eleIDu = M_electricSolverPtr->potentialFESpace().fe().currentLocalId();
        UInt nbNode = ( UInt ) M_electricSolverPtr->potentialFESpace().fe().nbFEDof();
        for ( UInt iNode = 0 ; iNode < nbNode ; iNode++ )
        {
            Int ig = M_electricSolverPtr->potentialFESpace().dof().localToGlobalMap ( eleIDu, iNode );
            elvec_u.vec() [ iNode ] = uVecRep[ig];
        }

        UInt eleID = M_electricSolverPtr->potentialFESpace().fe().currentLocalId();
        M_ionicModelPtr->updateElementSolution (eleID);
        M_ionicModelPtr->computeIonicCurrent (M_data.membraneCapacitance(), elvec_Iion, elvec_u, M_electricSolverPtr->potentialFESpace() );

        //! Computing Iapp
        AssemblyElemental::source (M_heartFctPtr->stimulus(),
                                   elvec_Iapp,
                                   M_electricSolverPtr->potentialFESpace().fe(),
                                   M_data.time(), 0);
        AssemblyElemental::source (M_heartFctPtr->stimulus(),
                                   elvec_Iapp,
                                   M_electricSolverPtr->potentialFESpace().fe(),
                                   M_data.time(),
                                   1);
        UInt totalUDof  = M_electricSolverPtr->potentialFESpace().map().map (Unique)->NumGlobalElements();

        for ( UInt iNode = 0 ; iNode < nbNode ; iNode++ )
        {
            Int ig = M_electricSolverPtr->potentialFESpace().dof().localToGlobalMap ( eleIDu, iNode );
            rhs.sumIntoGlobalValues (ig, elvec_Iapp.vec() [iNode] +
                                     M_data.volumeSurfaceRatio() * elvec_Iion.vec() [iNode] );
            rhs.sumIntoGlobalValues (ig + totalUDof,
                                     -elvec_Iapp.vec() [iNode + nbNode] -
                                     M_data.volumeSurfaceRatio() * elvec_Iion.vec() [iNode] );
        }
    }
    rhs.globalAssemble();

    rhs += M_electricSolverPtr->matrMass() * M_data.volumeSurfaceRatio() *
           M_data.membraneCapacitance() * M_electricSolverPtr->BDFIntraExtraPotential().time_der (M_data.timeStep() );

    MPI_Barrier (MPI_COMM_WORLD);

    chrono.stop();
    if (verbose)
    {
        std::cout << "done in " << chrono.diff() << " s." << std::endl;
    }
}
#endif
