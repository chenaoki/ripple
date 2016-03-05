#include "HeartBidomainModel.hpp"
#include "HeartIonicSolver.hpp"

namespace LifeV
{

HeartBidomainModel::HeartBidomainModel (const GetPot& dataFile ) :
    MeshData ( dataFile, "electric/space_discretization" ),
    TimeData ( dataFile, "electric/time_discretization" )
{
    setup ( dataFile);
}

HeartBidomainModel::HeartBidomainModel() :
    MeshData                            ( ),
    TimeData                            ( ),
    M_fibersFormat                      ( ),
    M_hasFibers                         ( ),
    M_heartDiffusionFactor              ( ),
    M_externalDiffusivity               ( ),
    M_internalDiffusivity               ( ),
    M_longitudinalExternalConductivity  ( ),
    M_longitudinalInternalConductivity  ( ),
    M_membraneCapacitance               ( ),
    M_transversalExternalConductivity   ( ),
    M_transversalInternalConductivity   ( ),
    M_volumeSurfaceRatio                ( ),
    M_fibersDirectory                   ( ),
    M_fibersFile                        ( ),
    M_uOrder                            ( )
{
}

HeartBidomainModel::HeartBidomainModel ( const HeartBidomainModel& dataBidomain ) :
    MeshData                            ( dataBidomain ),
    TimeData                            ( dataBidomain ),
    M_fibersFormat                      ( dataBidomain.M_fibersFormat ),
    M_hasFibers                         ( dataBidomain.M_hasFibers ),
    M_heartDiffusionFactor              ( dataBidomain.M_heartDiffusionFactor ),
    M_externalDiffusivity               ( dataBidomain.M_externalDiffusivity ),
    M_internalDiffusivity               ( dataBidomain.M_internalDiffusivity ),
    M_longitudinalExternalConductivity  ( dataBidomain.M_longitudinalExternalConductivity ),
    M_longitudinalInternalConductivity  ( dataBidomain.M_longitudinalInternalConductivity ),
    M_membraneCapacitance               ( dataBidomain.M_membraneCapacitance ),
    M_transversalExternalConductivity   ( dataBidomain.M_transversalExternalConductivity ),
    M_transversalInternalConductivity   ( dataBidomain.M_transversalInternalConductivity ),
    M_volumeSurfaceRatio                ( dataBidomain.M_volumeSurfaceRatio ),
    M_fibersDirectory                   ( dataBidomain.M_fibersDirectory ),
    M_fibersFile                        ( dataBidomain.M_fibersFile ),
    M_uOrder                            ( dataBidomain.M_uOrder )
{
}


// ===================================================
// Methods
// ===================================================
HeartBidomainModel&
HeartBidomainModel::operator= ( const HeartBidomainModel& dataBidomain )
{
    if ( this != &dataBidomain )
    {

        M_fibersFormat                      = dataBidomain.M_fibersFormat;
        M_hasFibers                         = dataBidomain.M_hasFibers;
        M_heartDiffusionFactor              = dataBidomain.M_heartDiffusionFactor;
        M_externalDiffusivity               = dataBidomain.M_externalDiffusivity;
        M_internalDiffusivity               = dataBidomain.M_internalDiffusivity;
        M_longitudinalExternalConductivity  = dataBidomain.M_longitudinalExternalConductivity;
        M_longitudinalInternalConductivity  = dataBidomain.M_longitudinalInternalConductivity;
        M_membraneCapacitance               = dataBidomain.M_membraneCapacitance;
        M_transversalExternalConductivity   = dataBidomain.M_transversalExternalConductivity;
        M_transversalInternalConductivity   = dataBidomain.M_transversalInternalConductivity;
        M_volumeSurfaceRatio                = dataBidomain.M_volumeSurfaceRatio;
        M_fibersDirectory                   = dataBidomain.M_fibersDirectory;
        M_fibersFile                        = dataBidomain.M_fibersFile;
        M_uOrder                            = dataBidomain.M_uOrder;
    }

    return *this;
}


void
HeartBidomainModel::setup ( const GetPot& dataFile )
{
    M_volumeSurfaceRatio                   = dataFile ("electric/physics/Chi", 1e3); // [1e-3 1/cm]   ColliPavarinoTaccardi2005
    M_membraneCapacitance                  = dataFile ("electric/physics/Cm", 1e-3);    // [1e-3 mF/cm2]   ColliPavarinoTaccardi2005
    if ( dataFile ("electric/physics/ion_model", 1) == 1)
    {
        M_internalDiffusivity              = dataFile ("electric/physics/D_i" , 3.3e-2);     // 3.3e-2  [1/Ohm/cm]   D_i_LR * D_RM/D_LR  see dataMonodomain
        M_externalDiffusivity              = dataFile ("electric/physics/D_e" , 4.29e-2);    // 4.29e-2 [1/Ohm/cm]  D_e_LR * D_RM/D_LR  see dataMonodomain
        M_longitudinalInternalConductivity = dataFile ("electric/physics/sigmal_i", 8.19e-2); // 8.19e-2 [1/Ohm/cm]   sigmal_i_LR * D_RM/D_LR
        M_transversalInternalConductivity  = dataFile ("electric/physics/sigmat_i", 8.6e-3);    // 8.6e-3  [1/Ohm/cm]   sigmat_i_LR * D_RM/D_LR
        M_longitudinalExternalConductivity = dataFile ("electric/physics/sigmal_e", 5.46e-2); // 5.46e-2 [1/Ohm/cm]  sigmal_e_LR * D_RM/D_LR
        M_transversalExternalConductivity  = dataFile ("electric/physics/sigmat_e", 3.69e-2); // 3.69e-2 [1/Ohm/cm]   sigmat_e_LR * D_RM/D_LR
    }
    else if ( dataFile ("electric/physics/ion_model", 1) == 2)
    {
        M_internalDiffusivity              = dataFile ("electric/physics/D_i" , 1.21e-3);       // sigmal_i/3 + sigmat_i*2/3
        M_externalDiffusivity              = dataFile ("electric/physics/D_e" , 1.57e-3);       // sigmal_e/3 + sigmat_e*2/3
        M_longitudinalInternalConductivity = dataFile ("electric/physics/sigmal_i", 3e-3);      // 3e-3      [1/Ohm/cm]   ColliPavarinoTaccardi2005
        M_transversalInternalConductivity  = dataFile ("electric/physics/sigmat_i", 3.1525e-4); // 3.1525e-4 [1/Ohm/cm]   ColliPavarinoTaccardi2005
        M_longitudinalExternalConductivity = dataFile ("electric/physics/sigmal_e", 2e-3);      // 2e-3      [1/Ohm/cm]   ColliPavarinoTaccardi2005
        M_transversalExternalConductivity  = dataFile ("electric/physics/sigmat_e", 1.3514e-3); // 1.3514e-3 [1/Ohm/cm]   ColliPavarinoTaccardi2005
    }
    else if ( dataFile ("electric/physics/ion_model", 1) == 3)
    {
        M_internalDiffusivity              = dataFile ("electric/physics/D_i" , 3.3e-2);     // 3.3e-2  [1/Ohm/cm]   D_i_LR * D_RM/D_LR  see dataMonodomain
        M_externalDiffusivity              = dataFile ("electric/physics/D_e" , 4.29e-2);    // 4.29e-2 [1/Ohm/cm]      D_e_LR * D_RM/D_LR  see dataMonodomain
        M_longitudinalInternalConductivity = dataFile ("electric/physics/sigmal_i", 8.19e-2); // 8.19e-2 [1/Ohm/cm]   sigmal_i_LR * D_RM/D_LR
        M_transversalInternalConductivity  = dataFile ("electric/physics/sigmat_i", 8.6e-3);     // 8.6e-3  [1/Ohm/cm]   sigmat_i_LR * D_RM/D_LR
        M_longitudinalExternalConductivity = dataFile ("electric/physics/sigmal_e", 5.46e-2); // 5.46e-2 [1/Ohm/cm]        sigmal_e_LR * D_RM/D_LR
        M_transversalExternalConductivity  = dataFile ("electric/physics/sigmat_e", 3.69e-2); // 3.69e-2 [1/Ohm/cm]        sigmat_e_LR * D_RM/D_LR
    }
    M_heartDiffusionFactor                 = dataFile ("electric/physics/heart_diff_fct", 0);
    M_uOrder                               = dataFile ( "electric/space_discretization/u_order", "P1");
    M_fibersFormat                         = dataFile ("electric/space_discretization/fibers_format", 0);
    M_hasFibers                            = dataFile ( "electric/space_discretization/has_fibers", 0);
    if ( M_hasFibers )
    {
        std::string fibersDirectory = dataFile ( "electric/space_discretization/fibers_dir", this->meshDir().c_str() );
        std::string fibersFile = this -> meshFile();
        fibersFile.replace ( fibersFile.find (".mesh"), 5, "fibers");
        M_fibersFile = fibersDirectory + dataFile ( "electric/space_discretization/fibers_file", fibersFile.c_str() );
        std::cout << "Fibers File: " << M_fibersFile << std::endl;
    }
    else
    {
        M_fibersFile = "";
        std::cout << "Fibers not included!" << std::endl;
    }
}


void
HeartBidomainModel::showMe ( std::ostream& output )
{
    output << "\n*** Values for data [fluid/physics]\n\n";
    output << "endtime   = " << endTime() << std::endl;
    output << "\n*** Values for data [fluid/miscellaneous]\n\n";
}

}
