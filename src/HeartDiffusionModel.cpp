#include "HeartDiffusionModel.hpp"
#include "HeartIonicSolver.hpp"

namespace LifeV
{

HeartDiffusionModel::HeartDiffusionModel(const GetPot& dataFile ) :
    MeshData ( dataFile, "electric/space_discretization" ),
    TimeData ( dataFile, "electric/time_discretization" )
{
    setup( dataFile);
}

HeartDiffusionModel::HeartDiffusionModel() :
    MeshData                        ( ),
    TimeData                        ( ),
    M_hasFibers                     ( ),
    M_diffusivity                   ( ),
    M_longitudinalConductivity      ( ),
    M_membraneCapacitance           ( ),
    M_transversalConductivity       ( ),
    M_volumeSurfaceRatio            ( ),
    M_conductivityRatio             ( ),
    M_fibersDirectory               ( ),
    M_fibersFile                    ( ),
    M_uOrder                        ( )
{
}

HeartDiffusionModel::HeartDiffusionModel ( const HeartDiffusionModel& model ) :
    MeshData                        ( model ),
    TimeData                        ( model ),
    M_hasFibers                     ( model.M_hasFibers ),
    M_diffusivity                   ( model.M_diffusivity ),
    M_longitudinalConductivity      ( model.M_longitudinalConductivity ),
    M_membraneCapacitance           ( model.M_membraneCapacitance ),
    M_transversalConductivity       ( model.M_transversalConductivity ),
    M_volumeSurfaceRatio            ( model.M_volumeSurfaceRatio ),
    M_conductivityRatio             ( model.M_conductivityRatio ),
    M_fibersDirectory               ( model.M_fibersDirectory ),
    M_fibersFile                    ( model.M_fibersFile ),
    M_uOrder                        ( model.M_uOrder )
{
}


// ===================================================
// Methods
// ===================================================
HeartDiffusionModel&
HeartDiffusionModel::operator= ( const HeartDiffusionModel& model )
{
    if ( this != &model )
    {
        M_hasFibers                                    = model.M_hasFibers;
        M_diffusivity                                     = model.M_diffusivity;
        M_longitudinalConductivity          = model.M_longitudinalConductivity;
        M_membraneCapacitance           = model.M_membraneCapacitance;
        M_transversalConductivity           = model.M_transversalConductivity;
        M_volumeSurfaceRatio                = model.M_volumeSurfaceRatio;
        M_conductivityRatio                      = model.M_conductivityRatio;
        M_fibersDirectory                           = model.M_fibersDirectory;
        M_fibersFile                                     = model.M_fibersFile;
        M_uOrder                                         = model.M_uOrder;
    }
    return *this;
}

void
HeartDiffusionModel::setup (  const GetPot& dataFile )
{

    dynamic_cast< MeshData* >(this) -> setup( dataFile, "electric/space_discretization" );
    dynamic_cast< TimeData* >(this) -> setup( dataFile, "electric/time_discretization" );

    switch(dataFile ("electric/physics/ion_model", 1))
    {
    case EnumHeartIonicSolverType::RM:
        M_diffusivity              = dataFile ("electric/physics/D" , 0.0156);  // 0.0156   [1/Ohm/cm]    L^2/T*D,  L=0.099 cm, T=0.63 ms D=1,  //RogersMcCulloch1994
        M_longitudinalConductivity = dataFile ("electric/physics/sigmal", 0.0328);   // 0.0328   [1/Ohm/cm]   sigmal_LR * D_RM/D_LR
        M_transversalConductivity  = dataFile ("electric/physics/sigmat", 0.00699);   // 0.00699  [1/Ohm/cm]   sigmat_LR * D_RM/D_LR
        break;
    case EnumHeartIonicSolverType::LR:
    case EnumHeartIonicSolverType::MS:
    case EnumHeartIonicSolverType::OR:
        M_diffusivity = dataFile ("electric/physics/D" , 5.7e-4); // 5.7e-4 [1/Ohm/cm]              sigmal/3 + sigmat*2/3
        M_longitudinalConductivity = dataFile ("electric/physics/sigmal", 1.2e-3); // 1.2e-3  [1/Ohm/cm]   sigmal_i*sigmal_e/(sigmal_i+sigmal_e)    ColliPavarinoTaccardi2005
        M_transversalConductivity  = dataFile ("electric/physics/sigmat", 2.56e-4); // 2.56e-4 [1/Ohm/cm]   sigmat_i*sigmat_e/(sigmat_i+sigmat_e)    ColliPavarinoTaccardi2005
        break;
    default:
        std::cout << "error : invalid ion model type" << std::endl;
        return;
        break;
    }
    M_volumeSurfaceRatio           = dataFile ("electric/physics/Chi", 1e3);  // [1e-3 1/cm]    ColliPavarinoTaccardi2005
    M_membraneCapacitance      = dataFile ("electric/physics/Cm", 1e-3);  // [1e-3 mF/cm2]  ColliPavarinoTaccardi2005
    M_conductivityRatio               =  dataFile ("electric/physics/lambda", 0.66667); // 0.66667 [adim]       sigmal_e/sigmal_i

    M_uOrder                               = dataFile ( "electric/space_discretization/u_order", "P1");

    M_hasFibers                          = dataFile ( "electric/space_discretization/has_fibers", 0);
    if ( M_hasFibers ){
        std::string fibersDirectory = dataFile ( "electric/space_discretization/fibers_dir", this->meshDir().c_str() );
        std::string fibersFile = this -> meshFile();
        fibersFile.replace ( fibersFile.find (".mesh"), 5, "fibers" );
        M_fibersFile = fibersDirectory + dataFile ( "electric/space_discretization/fibers_file", fibersFile.c_str() );
        std::cout << "Fibers File: " << M_fibersFile << std::endl;
    }else{
        M_fibersFile = "";
        std::cout << "Fibers not included!" << std::endl;
    }
}

} // namespace LifeV
