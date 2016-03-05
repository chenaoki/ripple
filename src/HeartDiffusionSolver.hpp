#ifndef _MONODOMAINSOLVER_H_
#define _MONODOMAINSOLVER_H_

#include <boost/shared_ptr.hpp>
#include <lifev/core/array/MatrixElemental.hpp>
#include <lifev/core/array/VectorElemental.hpp>
#include <lifev/core/fem/AssemblyElemental.hpp>
#include <lifev/core/fem/Assembly.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/core/algorithm/SolverAztecOO.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/fem/SobolevNorms.hpp>
#include <lifev/core/fem/GeometricMap.hpp>
#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/TimeAdvanceBDF.hpp>

#include "HeartDiffusionModel.hpp"
#include "HeartDiffusionFunctor.hpp"
#include "HeartStiffnessFibers.hpp"

namespace LifeV
{

//! monodomainSolver - Class featuring the usual solver for monodomain equations
template <
    typename Model,
    typename Mesh,
    typename AlgoSolver = LifeV::SolverAztecOO>
class HeartDiffusionSolver
{

public:
    typedef BCHandler                                                        bcHandlerRaw_Type;
    typedef boost::shared_ptr<bcHandlerRaw_Type>       bcHandler_Type;
    typedef typename AlgoSolver::matrix_type                 matrix_Type;
    typedef typename boost::shared_ptr<matrix_Type>                     matrixPtr_Type;
    typedef typename AlgoSolver::vector_type                 vector_Type;
    typedef typename AlgoSolver::prec_raw_type           preconditionerRaw_Type;
    typedef typename AlgoSolver::prec_type                    preconditioner_Type;

protected:
    const Model&                                                  M_model;
    const boost::shared_ptr<Epetra_Comm>        M_comm;
    FESpace<Mesh, MapEpetra>&                        M_uFESpace;
    Int                                                                     M_me;
    BCHandler*                                                     M_BChandlerElectric;
    bool                                                                   M_setBC;
    MapEpetra                                                       M_localMap;
    MapEpetra                                                       M_localMapVector;
    matrixPtr_Type                                                M_massMatrix;
    matrixPtr_Type                                                M_stiffnessMatrix;
    matrixPtr_Type                                                M_matrNoBC;
    vector_Type                                                     M_rhsNoBC;
    vector_Type                                                     M_solutionTransmembranePotential;
    vector_Type                                                     M_fiberVector;
    vector_Type                                                     M_residual;
    AlgoSolver                                                       M_linearSolver;
    preconditioner_Type                                        M_preconditioner;
    Real                                                                   M_diagonalize;
    bool                                                                   M_updated;
    bool                                                                   M_reusePreconditioner;
    bool                                                                   M_resetPreconditioner;
    Int                                                                       M_maxIteration;
    bool                                                                   M_recomputeMatrix;
    MatrixElemental                                              M_stiffnessElementaryMatrix;
    MatrixElemental                                              M_massElementaryMatrix;
    Real                                                                   massCoefficient;

public:
    HeartDiffusionSolver (
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
    };
    virtual ~HeartDiffusionSolver() {};

public:

    // Getters
    const vector_Type& solutionTransmembranePotential() const{ return M_solutionTransmembranePotential; };
    const vector_Type& fiberVector() const{ return M_fiberVector; };
    const vector_Type& residual() const{ return M_residual; };
    FESpace<Mesh, MapEpetra>& potentialFESpace(){ return M_uFESpace; };
    Epetra_Map const& getRepeatedMapEpetra() const{ return *M_localMap.map (Repeated); };
    Epetra_Map const& getRepeatedMapEpetraVec() const{ return *M_localMapVector.map (Repeated); };
    MapEpetra const& getMap() const{ return M_localMap; };
    matrix_Type& massMatrix(){ return *M_massMatrix; };
    UInt dim_u() const{ return M_uFESpace.dim(); };

    // Setters
    void setBC ( BCHandler& BCh_u )
    {
        M_BChandlerElectric = &BCh_u;
        M_setBC = true;
        return;
    };
    void postProcessing (bool _writeMesh = false);
    void resetPreconditioner(){ M_resetPreconditioner = true; };
    void recomputeMatrix (bool const recomp){ M_recomputeMatrix = recomp; };

public:

    void setup ( const GetPot& dataFile )
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
    };

    void buildSystem()
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
    };

    void initialize ( const HeartDiffusionFunctor::funcType& u0 )
    {
        vector_Type u (M_uFESpace.map() );
        M_uFESpace.interpolate (u0, u, 0.);
        this->initialize (u);
    };

    void initialize ( const vector_Type& u0 )
    {
        M_solutionTransmembranePotential = u0;
    };

    void updatePDESystem (Real alpha, vector_Type&  sourceVec)
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
    };

    void updatePDESystem ( vector_Type& sourceVec )
    {
        M_rhsNoBC = sourceVec;
        M_rhsNoBC.globalAssemble();
        return;
    };

    void PDEiterate ( bcHandlerRaw_Type& bch )
    {
        matrixPtr_Type matrFull ( new matrix_Type (*M_matrNoBC) );
        vector_Type    rhsFull = M_rhsNoBC;

        // boundary conditions update
        M_comm->Barrier();
        applyBoundaryConditions ( *matrFull, rhsFull, bch);
        M_comm->Barrier();
        solveSystem ( matrFull, rhsFull );
        return;
    };

private:

    void solveSystem (  matrixPtr_Type matrFull, vector_Type&   rhsFull )
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
    };

    void applyBoundaryConditions (
        matrix_Type& matrix,
        vector_Type& rhs,
        bcHandlerRaw_Type& BCh)
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
    };

}; // class MonodomainSolver


} // namespace LifeV
#endif //_MONODOMAINSOLVER_H_
