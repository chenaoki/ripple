#include "HeartDiffusionSolver.hpp"

using namespace std;

//
// IMPLEMENTATION
//

template<typename Model, typename Mesh>
LifeV::HeartDiffusionSolver<Model, Mesh>
::HeartDiffusionSolver (
    const Model&                                         model,
    LifeV::FESpace<Mesh, LifeV::MapEpetra>&                uFESpace,
    LifeV::BCHandler&                                                 BCh_u,
    boost::shared_ptr<Epetra_Comm>&       comm
) :
    M_model                   ( model ),
    M_uFESpace               ( uFESpace ),
    M_comm                   ( comm ),
    M_me                     ( M_comm->MyPID() ),
    M_BChandlerElectric      ( &BCh_u ),
    M_setBC                  ( true ),
    M_localMap               ( M_uFESpace.map() ),
    M_localMapVector         (M_localMap + M_localMap + M_localMap),
    M_massMatrix             ( ),
    M_stiffnessMatrix        ( ),
    M_matrNoBC               ( ),
    M_rhsNoBC                ( M_localMap ),
    M_solutionTransmembranePotential      ( M_localMap ),
    M_fiberVector                  ( M_localMapVector, Repeated ),
    M_residual               ( M_localMap ),
    M_linearSolver           ( ),
    M_preconditioner         ( ),
    M_updated                ( false ),
    M_reusePreconditioner    ( true ),
    M_resetPreconditioner    ( true ),
    M_maxIteration           ( -1 ),
    M_recomputeMatrix        ( false ),
    M_stiffnessElementaryMatrix ( M_uFESpace.fe().nbFEDof(), 1, 1 ),
    M_massElementaryMatrix   ( M_uFESpace.fe().nbFEDof(), 1, 1 )
{
    if (M_model.hasFibers() )
    {
        std::ifstream fibers (M_model.fibersFile().c_str() );
        std::cout << "fiber_file: " <<  M_model.fibersFile().c_str() << std::endl;
        UInt NumGlobalElements = M_localMapVector.map (Repeated)->NumGlobalElements();
        std::vector<Real> fiber_global_vector (NumGlobalElements);
        for ( UInt i = 0; i < NumGlobalElements; ++i)
        {
            fibers >> fiber_global_vector[i];
        }
        UInt NumMyElements = M_localMapVector.map (Repeated)->NumMyElements();
        for (UInt j = 0; j < NumMyElements; ++j)
        {
            UInt ig = M_localMapVector.map (Repeated)->MyGlobalElements() [j];
            M_fiberVector[ig] = fiber_global_vector[ig];
        }
        std::cout << std::endl;
        fiber_global_vector.clear();
    }
}

template<typename Model, typename Mesh>
void LifeV::HeartDiffusionSolver<Model, Mesh>
::setup ( const GetPot& dataFile )
{
    M_diagonalize = dataFile ( "electric/space_discretization/diagonalize",  1. );
    M_reusePreconditioner   = dataFile ( "electric/prec/reuse", true);
    M_linearSolver.setCommunicator (M_comm);
    M_linearSolver.setDataFromGetPot ( dataFile, "electric/solver" );
    M_maxIteration = dataFile ( "electric/solver/max_iter", -1);
    std::string precType = dataFile ( "electric/prec/prectype", "Ifpack");
    M_preconditioner.reset ( PRECFactory::instance().createObject ( precType ) );
    ASSERT (M_preconditioner.get() != 0, "monodomainSolver : Preconditioner not set");
    M_preconditioner->setDataFromGetPot ( dataFile, "electric/prec" );
}

template<typename Model, typename Mesh>
void LifeV::HeartDiffusionSolver<Model, Mesh>
::buildSystem()
{
    M_massMatrix.reset  ( new matrix_Type (M_localMap) );
    M_stiffnessMatrix.reset ( new matrix_Type (M_localMap) );
    M_comm->Barrier();

    for ( UInt iVol = 0; iVol < M_uFESpace.mesh()->numVolumes(); iVol++ ){

        M_stiffnessElementaryMatrix.zero();
        M_massElementaryMatrix.zero();
        M_uFESpace.fe().updateFirstDeriv ( M_uFESpace.mesh()->volumeList ( iVol ) );
        if (M_model.hasFibers() ){
            stiff ( M_model.longitudinalConductivity(),
                    M_model.transversalConductivity(),
                    M_fiberVector,
                    M_stiffnessElementaryMatrix,
                    M_uFESpace.fe(),
                    M_uFESpace.dof(),
                    0,
                    0);
        }else{
            AssemblyElemental::stiff (
                M_model.diffusivity(),
                M_stiffnessElementaryMatrix,
                M_uFESpace.fe(),
                0, 0 );
        }

        AssemblyElemental::mass (
            1.,
            M_massElementaryMatrix,
            M_uFESpace.fe(),
            0, 0 );

        assembleMatrix (
            *M_stiffnessMatrix,
             M_stiffnessElementaryMatrix,
             M_uFESpace.fe(),
             M_uFESpace.fe(),
             M_uFESpace.dof(),
             M_uFESpace.dof(),
             0, 0,
             0, 0);

        assembleMatrix (
            *M_massMatrix,
             M_massElementaryMatrix,
             M_uFESpace.fe(),
             M_uFESpace.fe(),
             M_uFESpace.dof(),
             M_uFESpace.dof(),
             0, 0,
             0, 0);
    }

    massCoefficient = M_model.volumeSurfaceRatio() * M_model.membraneCapacitance() / M_model.timeStep();
    M_comm->Barrier();

    M_stiffnessMatrix->globalAssemble();
    M_massMatrix->globalAssemble();
    M_comm->Barrier();

    M_matrNoBC.reset (new matrix_Type (M_localMap, M_stiffnessMatrix->meanNumEntries() ) );
    *M_matrNoBC += *M_stiffnessMatrix;
    *M_matrNoBC += *M_massMatrix * massCoefficient;
    M_matrNoBC->globalAssemble();

    return;
}

template<typename Model, typename Mesh>
void LifeV::HeartDiffusionSolver<Model, Mesh>
::initialize ( const HeartDiffusionFunctor::funcType& u0 )
{
    vector_Type u (M_uFESpace.map() );
    M_uFESpace.interpolate (u0, u, 0.);
    this->initialize (u);
}

template<typename Model, typename Mesh>
void LifeV::HeartDiffusionSolver<Model, Mesh>
::initialize ( const vector_Type& u0 )
{
    M_solutionTransmembranePotential = u0;
}

template<typename Model, typename Mesh>
void LifeV::HeartDiffusionSolver<Model, Mesh>
::PDEiterate ( bcHandlerRaw_Type& bch )
{
    matrixPtr_Type matrFull ( new matrix_Type (*M_matrNoBC) );
    vector_Type    rhsFull = M_rhsNoBC;

    // boundary conditions update
    M_comm->Barrier();
    applyBoundaryConditions ( *matrFull, rhsFull, bch);
    M_comm->Barrier();
    solveSystem ( matrFull, rhsFull );
    return;
}

template<typename Model, typename Mesh>
void LifeV::HeartDiffusionSolver<Model, Mesh>
::updatePDESystem (
    Real alpha,
    vector_Type& sourceVec )
{
    M_rhsNoBC = sourceVec;
    M_rhsNoBC.globalAssemble();

    M_updated = false;
    if (M_recomputeMatrix) buildSystem();

    M_matrNoBC.reset (new matrix_Type (M_localMap, M_stiffnessMatrix->meanNumEntries() ) );
    *M_matrNoBC += *M_stiffnessMatrix;
    *M_matrNoBC += *M_massMatrix * alpha;

    M_updated = true;
    M_matrNoBC->globalAssemble();

    return;
}

template<typename Model, typename Mesh>
void LifeV::HeartDiffusionSolver<Model, Mesh>
::updatePDESystem (vector_Type& sourceVec )
{
    M_rhsNoBC = sourceVec;
    M_rhsNoBC.globalAssemble();
    return;
}

//*******************************************************************************//

template<typename Model, typename Mesh>
void LifeV::HeartDiffusionSolver<Model, Mesh>
::applyBoundaryConditions (
    matrix_Type&        matrix,
    vector_Type&        rhs,
    bcHandlerRaw_Type&  BCh )
{
    // BC manage for the PDE
    if ( !BCh.bcUpdateDone() )
    {
        BCh.bcUpdate ( *M_uFESpace.mesh(), M_uFESpace.feBd(), M_uFESpace.dof() );
    }
    vector_Type rhsFull (M_rhsNoBC, Repeated, Zero);
    bcManage ( matrix, rhs, *M_uFESpace.mesh(), M_uFESpace.dof(),
               BCh, M_uFESpace.feBd(), 1., M_model.time() );
    rhs = rhsFull;
    if ( BCh.hasOnlyEssential() && M_diagonalize )
    {
        matrix.diagonalize ( 1 * dim_u(),
                             M_diagonalize,
                             rhs,
                             0.);
    }
    return;
}

template<typename Model, typename Mesh>
void LifeV::HeartDiffusionSolver<Model, Mesh>
::solveSystem ( matrixPtr_Type  matrFull, vector_Type&    rhsFull )
{
    LifeChrono chrono;

    M_linearSolver.setMatrix (*matrFull);

    if ( !M_reusePreconditioner || M_resetPreconditioner )
    {
        M_preconditioner->buildPreconditioner (matrFull);
        Real condest = M_preconditioner->condest();
        M_linearSolver.setPreconditioner (M_preconditioner);
        M_resetPreconditioner = false;
    }

    Int numIter = M_linearSolver.solve (M_solutionTransmembranePotential, rhsFull);
    if (numIter > M_maxIteration) M_resetPreconditioner = true;
    M_comm->Barrier();
    return;
}
