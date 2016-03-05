#ifndef _BIDOMAINSOLVER_H_
#define _BIDOMAINSOLVER_H_

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

#include "HeartBidomainModel.hpp"
#include "HeartDiffusionFunctor.hpp"
#include "HeartStiffnessFibers.hpp"

namespace LifeV
{

const UInt nbComp = 2;

template < typename Mesh,
         typename SolverType = LifeV::SolverAztecOO >
class HeartBidomainSolver
{

public:

    typedef HeartBidomainModel data_type;
    typedef Real ( *Function ) ( const Real&, const Real&, const Real&,
                                 const Real&, const ID& );
    typedef boost::function < Real ( const Real&, const Real&, const Real&,
                                     const Real&, const ID& ) > source_Type;
    typedef Mesh mesh_Type;
    typedef BCHandler                             bchandlerRaw_Type;
    typedef typename SolverType::matrix_type      matrix_Type;
    typedef boost::shared_ptr<matrix_Type>        matrixPtr_Type;
    typedef typename SolverType::vector_type      vector_Type;
    typedef typename SolverType::prec_raw_type    precRaw_Type;
    typedef typename SolverType::prec_type        prec_Type;

public:

    HeartBidomainSolver ( const data_type&          dataType,
                          FESpace<Mesh, MapEpetra>& pFESpace,
                          FESpace<Mesh, MapEpetra>& uFESpace,
                          BCHandler&                bcHandler,
                          boost::shared_ptr<Epetra_Comm>& comm );

    virtual ~HeartBidomainSolver();

public:

    //! Updates sources, bc treatment and solve the bidomain system
    virtual void PDEiterate ( bchandlerRaw_Type& bch );

    //! Sets up the system solver
    virtual void setup        ( const GetPot& dataFile );

    //! Builds time independent parts of PDE system
    virtual void buildSystem();

    //! Updates time dependent parts of PDE system
    virtual void updatePDESystem (Real alpha,
                                  vector_Type& sourceVec);

    //! Updates time dependent parts of PDE system
    virtual void updatePDESystem ( vector_Type& sourceVec );

    //! Initialize
    void initialize ( const source_Type&, const source_Type&  );
    void initialize ( const Function&, const Function&  );
    void initialize ( const vector_Type&, const vector_Type& );

    //! Returns the local solution vector
    const vector_Type& solutionTransmembranePotential() const
    {
        return M_solutionTransmembranePotential;
    }
    const vector_Type& solutionExtraPotential() const
    {
        return M_solutionExtraPotential;
    }
    const vector_Type& solutionIntraExtraPotential() const
    {
        return M_solutionIntraExtraPotential;
    }
    const vector_Type& solutionTransmembranePotentialExtrapolated() const
    {
        return M_solutionTransmembranePotentialExtrapolated;
    }
    const vector_Type& fiberVector() const
    {
        return M_fiberVector;
    }

    //! Returns the local residual vector
    const vector_Type& residual() const
    {
        return M_residual;
    }

    //! Returns u FE space
    FESpace<Mesh, MapEpetra>& potentialFESpace()
    {
        return M_uFESpace;
    }

    //! Sets Bidomain BCs
    void setBC (BCHandler& BCh_u)
    {
        M_BCh_electric = &BCh_u;
        M_setBC = true;
    }

    void resetPreconditioner()
    {
        M_resetPreconditioner = true;
    }

    //! Return maps getMap
    Epetra_Map const& getRepeatedMapEpetra() const
    {
        return *M_localMap.map (Repeated);
    }

    Epetra_Map const& getRepeatedMapEpetraVec() const
    {
        return *M_localMapVec.map (Repeated);
    }

    MapEpetra const& getMap() const
    {
        return M_localMap;
    }

    void recomputeMatrix (bool const recomp)
    {
        M_recomputeMatrix = recomp;
    }

    matrix_Type& matrMass()
    {
        return *M_matrMass;
    }

    TimeAdvanceBDF<vector_Type>& BDFIntraExtraPotential()
    {
        return M_BDFIntraExtraPotential;
    }

    //@}

protected:

    //! Solves PDE system
    void solveSystem            (  matrixPtr_Type matrFull,
                                   vector_Type&   rhsFull );

    //! Apply BC
    void applyBoundaryConditions (  matrix_Type&        matrix,
                                    vector_Type&        rhs,
                                    bchandlerRaw_Type& BCh);

    //! compute mean of vector x
    Real computeMean ( vector_Type& x );

    //! removes a scalar from each entry of vector x
    void removeValue ( vector_Type& x, Real& value );

    //! Data
    const data_type&               M_data;

    //! u FE space
    FESpace<Mesh, MapEpetra>&      M_pFESpace;
    FESpace<Mesh, MapEpetra>&      M_uFESpace;

    //! MPI communicator
    const boost::shared_ptr<Epetra_Comm> M_comm;
    Int                            M_me;

    //! Bidomain BC
    BCHandler*                     M_BCh_electric;
    bool                           M_setBC;

    //! Map
    MapEpetra                      M_localMap;
    MapEpetra                      M_localMap_u;
    MapEpetra                      M_localMapVec;

    //! mass matrix
    matrixPtr_Type                 M_matrMass;

    //! Stiff matrix: D*stiff
    matrixPtr_Type                 M_matrStiff;

    matrixPtr_Type                 M_matrNoBC;

    //! Right hand side for the PDE
    vector_Type                    M_rhsNoBC;

    //! Global solution _u
    vector_Type                    M_solutionIntraExtraPotential;

    vector_Type                    M_solutionTransmembranePotential;

    vector_Type                    M_solutionExtraPotential;

    vector_Type                    M_solutionIntraExtraPotentialExtrapolated;

    vector_Type                    M_solutionTransmembranePotentialExtrapolated;

    //! Local fibers vector
    vector_Type                    M_fiberVector;

    //! residual
    vector_Type                    M_residual;

    //! Solver
    SolverType                     M_linearSolver;

    prec_Type                      M_prec;

    Real                           M_diagonalize;

    //! Boolean that indicates if output is sent to cout
    bool                           M_verbose;

    //! Boolean that indicates if the matrix is updated for the current iteration
    bool                           M_updated;

    //! Boolean that indicates if the precond has to be recomputed
    bool                           M_reusePrec;
    bool                           M_resetPreconditioner;

    //! Integer storing the max number of solver iteration with prec recomputing
    Int                            M_maxIterSolver;

    //! Boolean that indicates if the matrix has to be recomputed
    bool                           M_recomputeMatrix;

    TimeAdvanceBDF<vector_Type>            M_BDFIntraExtraPotential;
private:

    //! Elementary matrices
    MatrixElemental                        M_elmatStiff;
    MatrixElemental                        M_elmatMass;
    Real                           massCoeff;
    UInt potentialFESpaceDimension() const
    {
        return M_uFESpace.dim();
    }
}; // class BidomainSolver



//
// IMPLEMENTATION
//

//! Constructors
template<typename Mesh, typename SolverType>
HeartBidomainSolver<Mesh, SolverType>::
HeartBidomainSolver ( const data_type&          dataType,
                      FESpace<Mesh, MapEpetra>& pFESpace,
                      FESpace<Mesh, MapEpetra>& uFESpace,
                      BCHandler&                BCh_u,
                      boost::shared_ptr<Epetra_Comm>&  comm ) :
    M_data                   ( dataType ),
    M_pFESpace               ( pFESpace ),
    M_uFESpace               ( uFESpace ),
    M_comm                   ( comm ),
    M_me                     ( M_comm->MyPID() ),
    M_BCh_electric           ( &BCh_u ),
    M_setBC                  ( true ),
    M_localMap               ( pFESpace.map() ),
    M_localMap_u             ( uFESpace.map() ),
    M_localMapVec            (M_localMap_u + M_localMap_u + M_localMap_u),
    M_matrMass               ( ),
    M_matrStiff              ( ),
    M_matrNoBC               ( ),
    M_rhsNoBC                ( M_localMap ),
    M_solutionIntraExtraPotential               ( M_localMap ),
    M_solutionTransmembranePotential            ( M_localMap_u ),
    M_solutionExtraPotential                    ( M_localMap_u ),
    M_solutionIntraExtraPotentialExtrapolated   ( M_localMap ),
    M_solutionTransmembranePotentialExtrapolated ( M_localMap_u ),
    M_fiberVector           ( M_localMapVec, Repeated ),
    M_residual               ( M_localMap ),
    M_linearSolver           ( ),
    M_prec                   ( ),
    M_verbose                ( M_me == 0),
    M_updated                ( false ),
    M_reusePrec              ( true ),
    M_resetPreconditioner              ( true ),
    M_maxIterSolver          ( -1 ),
    M_recomputeMatrix        ( false ),
    M_BDFIntraExtraPotential ( ),
    M_elmatStiff             ( M_pFESpace.fe().nbFEDof(), 2, 2 ),
    M_elmatMass              ( M_pFESpace.fe().nbFEDof(), 2, 2 )
{
    if (M_data.hasFibers() )
    {
        std::stringstream MyPID;
        std::ifstream fibers (M_data.fibersFile().c_str() );
        char line[1024];
        UInt numPoints = M_localMap_u.map (Unique)->NumGlobalElements();
        UInt numGlobalElements = M_localMapVec.map (Unique)->NumGlobalElements();
        std::vector<Real> fiberGlobalVector (numGlobalElements);
        if (M_data.fibersFormat() )
        {
            fibers.getline (line, 1024);
            for ( UInt i = 0; i < numPoints; ++i)
                for (UInt k = 0; k < 3; ++k)
                {
                    fibers >> fiberGlobalVector[i + k * numPoints];
                }
        }
        else
        {
            for (UInt i = 0 ; i < numGlobalElements; ++i)
            {
                fibers >> fiberGlobalVector[i];
            }
        }
        UInt numMyElements = M_localMapVec.map (Repeated)->NumMyElements();
        for (UInt j = 0; j < numMyElements; ++j)
        {
            UInt ig = M_localMapVec.map (Repeated)->MyGlobalElements() [j];
            M_fiberVector[ig] = fiberGlobalVector[ig];
        }
        std::cout << std::endl;
        fiberGlobalVector.clear();
    }
};

template<typename Mesh, typename SolverType>
HeartBidomainSolver<Mesh, SolverType>::
~HeartBidomainSolver()
{
}

template<typename Mesh, typename SolverType>
void HeartBidomainSolver<Mesh, SolverType>::setup ( const GetPot& dataFile )
{
    M_diagonalize = dataFile ( "electric/space_discretization/diagonalize",  1. );
    M_reusePrec   = dataFile ( "electric/prec/reuse", true);
    M_linearSolver.setCommunicator (M_comm);
    M_linearSolver.setDataFromGetPot ( dataFile, "electric/solver" );
    M_maxIterSolver = dataFile ( "electric/solver/max_iter", -1);
    std::string precType = dataFile ( "electric/prec/prectype", "Ifpack");
    M_prec.reset ( PRECFactory::instance().createObject ( precType ) );
    ASSERT (M_prec.get() != 0, "bidomainSolver : Preconditioner not set");
    M_prec->setDataFromGetPot ( dataFile, "electric/prec" );
    M_BDFIntraExtraPotential.setup(1,1);
}

template<typename Mesh, typename SolverType>
void HeartBidomainSolver<Mesh, SolverType>::buildSystem()
{
    M_matrMass.reset ( new matrix_Type (M_localMap) );
    M_matrStiff.reset ( new matrix_Type (M_localMap) );

    if (M_verbose)
    {
        std::cout << "  f-  Computing constant matrices ...        ";
    }

    LifeChrono chrono;
    LifeChrono chronoDer;
    LifeChrono chronoStiff;
    LifeChrono chronoMass;
    LifeChrono chronoStiffAssemble;
    LifeChrono chronoMassAssemble;
    LifeChrono chronoZero;

    M_comm->Barrier();
    chrono.start();

    //! Elementary computation and matrix assembling
    //! Loop on elements

    for ( UInt iVol = 0; iVol < M_pFESpace.mesh()->numVolumes(); iVol++ )
    {
        chronoZero.start();
        M_elmatStiff.zero();
        M_elmatMass.zero();
        chronoZero.stop();

        chronoStiff.start();
        switch (M_data.heartDiffusionFactor() )
        {
            case 0:
                chronoDer.start();
                M_pFESpace.fe().updateFirstDeriv ( M_pFESpace.mesh()->volumeList ( iVol ) );
                chronoDer.stop();
                if (M_data.hasFibers() )
                {
                    stiff ( M_data.longitudinalInternalConductivity(),
                            M_data.transversalInternalConductivity(),
                            M_fiberVector,
                            M_elmatStiff,
                            M_pFESpace.fe(),
                            M_pFESpace.dof(), 0, 0);
                    stiff ( M_data.longitudinalExternalConductivity(),
                            M_data.transversalExternalConductivity(),
                            M_fiberVector,
                            M_elmatStiff,
                            M_pFESpace.fe(),
                            M_pFESpace.dof(), 1, 1);
                }
                else
                {
                    AssemblyElemental::stiff (
                        M_data.internalDiffusivity(), 
                        M_elmatStiff, 
                        M_pFESpace.fe(), 0, 0 );
                    AssemblyElemental::stiff ( 
                        M_data.externalDiffusivity(), 
                        M_elmatStiff,  
                        M_pFESpace.fe(), 1, 1 );
                }
                break;
        }
        chronoStiff.stop();

        chronoMass.start();
        AssemblyElemental::mass ( 1., M_elmatMass, M_pFESpace.fe(), 0, 0 );
        AssemblyElemental::mass ( -1., M_elmatMass, M_pFESpace.fe(), 0, 1 );
        AssemblyElemental::mass ( -1., M_elmatMass, M_pFESpace.fe(), 1, 0 );
        AssemblyElemental::mass ( 1., M_elmatMass, M_pFESpace.fe(), 1, 1 );
        chronoMass.stop();

        chronoStiffAssemble.start();
        for ( UInt iComp = 0; iComp < nbComp; iComp++ )
        {
            assembleMatrix ( *M_matrStiff,
                             M_elmatStiff,
                             M_pFESpace.fe(),
                             M_pFESpace.fe(),
                             M_pFESpace.dof(),
                             M_pFESpace.dof(),
                             iComp, iComp,
                             iComp * potentialFESpaceDimension(), iComp * potentialFESpaceDimension() );
        }
        chronoStiffAssemble.stop();

        chronoMassAssemble.start();
        for ( UInt iComp = 0; iComp < nbComp; iComp++ )
        {
            for ( UInt jComp = 0; jComp < nbComp; jComp++ )
            {
                assembleMatrix ( *M_matrMass,
                                 M_elmatMass,
                                 M_pFESpace.fe(),
                                 M_pFESpace.fe(),
                                 M_pFESpace.dof(),
                                 M_pFESpace.dof(),
                                 iComp, jComp,
                                 iComp * potentialFESpaceDimension(), jComp * potentialFESpaceDimension() );
            }
        }
        chronoMassAssemble.stop();
    }

    massCoeff = M_data.volumeSurfaceRatio() * M_data.membraneCapacitance() *
                M_BDFIntraExtraPotential.coefficientExtrapolationFirstDerivative(0) / M_data.timeStep();
     

    M_comm->Barrier();
    chrono.stop();

    if (M_verbose)
    {
        std::cout << "done in " << chrono.diff() << " s.\n" << std::flush;
    }
    if (M_verbose)
    {
        std::cout << "  f-  Finalizing the matrices     ...        " << std::flush;
    }

    chrono.start();
    M_matrStiff->globalAssemble();
    M_matrMass->globalAssemble();
    M_comm->Barrier();

    M_matrNoBC.reset (new matrix_Type (M_localMap, M_matrStiff->meanNumEntries() ) );

    //! Computing 1.0/dt * M + K
    *M_matrNoBC += *M_matrStiff;
    *M_matrNoBC += *M_matrMass * massCoeff;
    M_matrNoBC->globalAssemble();
    chrono.stop();
    if (M_verbose)
    {
        std::cout << "done in " << chrono.diff() << " s." << std::endl;
    }

}

template<typename Mesh, typename SolverType>
void HeartBidomainSolver<Mesh, SolverType>::
initialize ( const source_Type& ui0, const source_Type& ue0 )
{

    vector_Type intracellularPotential (M_uFESpace.map() );
    vector_Type extracellularPotential (M_uFESpace.map() );

    M_uFESpace.interpolate (ui0, intracellularPotential, 0.);
    M_uFESpace.interpolate (ue0, extracellularPotential, 0.);
    initialize (intracellularPotential, extracellularPotential);
}

template<typename Mesh, typename SolverType>
void HeartBidomainSolver<Mesh, SolverType>::
initialize ( const Function& ui0, const Function& ue0 )
{

    vector_Type intracellularPotential (M_uFESpace.map() );
    vector_Type extracellularPotential (M_uFESpace.map() );
    M_uFESpace.interpolate (ui0, intracellularPotential, 0.);
    M_uFESpace.interpolate (ue0, extracellularPotential, 0.);
    initialize (intracellularPotential, extracellularPotential);
}

template<typename Mesh, typename SolverType>
void HeartBidomainSolver<Mesh, SolverType>::
initialize ( const vector_Type& ui0, const vector_Type& ue0 )
{
    M_solutionIntraExtraPotential = ui0;
    M_solutionIntraExtraPotential.add (ue0, M_uFESpace.dof().numTotalDof() );

    for ( Int i = 0 ; i < M_solutionTransmembranePotential.epetraVector().MyLength() ; i++ )
    {
        Int ig = M_solutionTransmembranePotential.blockMap().MyGlobalElements() [i];
        M_solutionTransmembranePotential[ig] = M_solutionIntraExtraPotential[ig] -
                                               M_solutionIntraExtraPotential[ig + potentialFESpaceDimension()];
        M_solutionExtraPotential[ig] = M_solutionIntraExtraPotential[ig +
                                                                     potentialFESpaceDimension()];
    }
    M_solutionIntraExtraPotentialExtrapolated = M_solutionIntraExtraPotential;
    M_solutionTransmembranePotentialExtrapolated = M_solutionTransmembranePotential;
    M_BDFIntraExtraPotential.setInitialCondition (M_solutionIntraExtraPotential);
    //M_BDFIntraExtraPotential.showMe();
}

template<typename Mesh, typename SolverType>
void HeartBidomainSolver<Mesh, SolverType>::
updatePDESystem (Real alpha, vector_Type& sourceVec)
{

    LifeChrono chrono;
    if (M_verbose)
        std::cout << "  f-  Updating mass term and right hand side... "
                  << std::flush;

    chrono.start();

    M_rhsNoBC = sourceVec;
    M_rhsNoBC.globalAssemble();

    chrono.stop();

    if (M_verbose)
    {
        std::cout << "done in " << chrono.diff() << " s.\n"  << std::flush;
    }

    M_updated = false;

    if (M_recomputeMatrix)
    {
        buildSystem();
    }

    if (M_verbose)
        std::cout << "  f-  Copying the matrices ...                 "
                  << std::flush;

    chrono.start();

    M_matrNoBC.reset (new matrix_Type (M_localMap, M_matrStiff->meanNumEntries() ) );

    *M_matrNoBC += *M_matrStiff;

    *M_matrNoBC += *M_matrMass * alpha;

    chrono.stop();
    if (M_verbose) std::cout << "done in " << chrono.diff() << " s.\n"
                                 << std::flush;

    M_updated = true;
    M_matrNoBC->globalAssemble();
}

template<typename Mesh, typename SolverType>
void HeartBidomainSolver<Mesh, SolverType>::
updatePDESystem (vector_Type& sourceVec )
{

    LifeChrono chrono;

    if (M_verbose)
        std::cout << "  f-  Updating right hand side... "
                  << std::flush;

    chrono.start();
    M_rhsNoBC = sourceVec;
    M_rhsNoBC.globalAssemble();
    chrono.stop();

    if (M_verbose)
    {
        std::cout << "done in " << chrono.diff() << " s.\n"  << std::flush;
    }
}

template<typename Mesh, typename SolverType>
void HeartBidomainSolver<Mesh, SolverType>::PDEiterate ( bchandlerRaw_Type& bch )
{

    LifeChrono chrono;
    chrono.start();

    matrixPtr_Type matrFull ( new matrix_Type (*M_matrNoBC) );
    vector_Type    rhsFull = M_rhsNoBC;

    chrono.stop();

    if (M_verbose)
        std::cout << "done in " << chrono.diff() << " s.\n"
                  << std::flush;

    // boundary conditions update
    M_comm->Barrier();
    if (M_verbose) std::cout << "  f-  Applying boundary conditions ...         "
                                 << std::flush;

    chrono.start();
    applyBoundaryConditions ( *matrFull, rhsFull, bch);

    chrono.stop();

    M_comm->Barrier();

    if (M_verbose)
    {
        std::cout << "done in " << chrono.diff() << " s.\n" << std::flush;
    }

    //! Solving the system
    solveSystem ( matrFull, rhsFull );

    //    M_residual  = M_rhsNoBC;
    //    M_residual -= *M_matrNoBC*M_solutionIntraExtraPotential;

} // iterate()


template<typename Mesh, typename SolverType>
void HeartBidomainSolver<Mesh, SolverType>::solveSystem ( matrixPtr_Type  matrFull,
                                                          vector_Type&    rhsFull )
{
    LifeChrono chrono;

    if (M_verbose)
    {
        std::cout << "  f-  Setting up the solver ...                ";
    }

    chrono.start();
    M_linearSolver.setMatrix (*matrFull);
    chrono.stop();

    if (M_verbose)
        std::cout << "done in " << chrono.diff() << " s.\n"
                  << std::flush;

    if ( !M_reusePrec || M_resetPreconditioner )
    {
        chrono.start();

        if (M_verbose)
        {
            std::cout << "  f-  Computing the precond ...                ";
        }
        M_prec->buildPreconditioner (matrFull);

        Real condest = M_prec->condest();

        M_linearSolver.setPreconditioner (M_prec);

        chrono.stop();
        if (M_verbose)
        {
            std::cout << "done in " << chrono.diff() << " s.\n";
            std::cout << "         Estimated condition number = " << condest << "\n" <<  std::flush;
        }

        M_resetPreconditioner = false;
    }
    else
    {
        if (M_verbose)
        {
            std::cout << "  f-  Reusing  precond ...                \n" <<  std::flush;
        }
    }

    chrono.start();

    if (M_verbose)
    {
        std::cout << "  f-  Solving system ...                                ";
    }

    Int numIter = M_linearSolver.solve (M_solutionIntraExtraPotential, rhsFull);

    chrono.stop();

    if (M_verbose)
    {
        std::cout << "\ndone in " << chrono.diff()
                  << " s. ( " << numIter << "  iterations. ) \n"
                  << std::flush;
    }

    if (numIter > M_maxIterSolver)
    {
        M_resetPreconditioner = true;
    }

    for ( Int i = 0 ; i < M_solutionTransmembranePotential.epetraVector().MyLength() ; i++ )
    {
        UInt ig = M_solutionTransmembranePotential.blockMap().MyGlobalElements() [i];
        M_solutionExtraPotential[ig] = M_solutionIntraExtraPotential[ig + potentialFESpaceDimension()];
    }
    Real meanExtraPotential = computeMean (M_solutionExtraPotential);
    removeValue (M_solutionIntraExtraPotential, meanExtraPotential);

    for ( Int i = 0 ; i < M_solutionTransmembranePotential.epetraVector().MyLength() ; i++ )
    {
        UInt ig = M_solutionTransmembranePotential.blockMap().MyGlobalElements() [i];
        M_solutionTransmembranePotential[ig] = M_solutionIntraExtraPotential[ig] - M_solutionIntraExtraPotential[ig + potentialFESpaceDimension()];
        M_solutionExtraPotential[ig] = M_solutionIntraExtraPotential[ig + potentialFESpaceDimension()];
    }

    M_BDFIntraExtraPotential.shiftRight (M_solutionIntraExtraPotential);
    M_BDFIntraExtraPotential.extrapolation(M_solutionIntraExtraPotentialExtrapolated);

    for ( Int i = 0 ; i < M_solutionTransmembranePotential.epetraVector().MyLength() ; i++ )
    {
        UInt ig = M_solutionTransmembranePotential.blockMap().MyGlobalElements() [i];
        M_solutionTransmembranePotentialExtrapolated[ig] = M_solutionIntraExtraPotentialExtrapolated[ig] - M_solutionIntraExtraPotentialExtrapolated[ig + potentialFESpaceDimension()];
    }
}

template<typename Mesh, typename SolverType>
void HeartBidomainSolver<Mesh, SolverType>::applyBoundaryConditions (matrix_Type& matrix,
                                                                     vector_Type& rhs,
                                                                     bchandlerRaw_Type& BCh )
{
    // BC manage for the PDE
    if ( !BCh.bcUpdateDone() )
    {
        BCh.bcUpdate ( *M_pFESpace.mesh(), M_pFESpace.feBd(), M_pFESpace.dof() );
    }

    vector_Type rhsFull (M_rhsNoBC, Repeated, Zero);

    //    rhsFull.Import(M_rhsNoBC, Zero); // ignoring non-local entries, Otherwise they are summed up lately

    bcManage ( matrix, rhs, *M_pFESpace.mesh(), M_pFESpace.dof(), BCh, M_pFESpace.feBd(), 1.,
               M_data.time() );
    rhs = rhsFull;

    if ( BCh.hasOnlyEssential() && M_diagonalize )
    {
        matrix.diagonalize ( 1 * potentialFESpaceDimension(),
                             M_diagonalize,
                             rhs,
                             0.);
    }

} // applyBoundaryCondition

template<typename Mesh, typename SolverType>
Real HeartBidomainSolver<Mesh, SolverType>::computeMean ( vector_Type& x )
{
    Real meanExtraPotential (0.);
    x.epetraVector().MeanValue (&meanExtraPotential);
    return meanExtraPotential;
} // computeMean()

template<typename Mesh, typename SolverType>
void HeartBidomainSolver<Mesh, SolverType>::removeValue ( vector_Type& x, Real& value )
{
    for ( Int i = 0 ; i < x.epetraVector().MyLength() ; i++ )
    {
        Int ig = x.blockMap().MyGlobalElements() [i];
        x[ig] -= value;
    }

} // removeMean()

} // namespace LifeV


#endif //_BIDOMAINSOLVER_H_
