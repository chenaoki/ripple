/*--
 OHaraRudy model
--*/

#ifndef _OHARA_RUDY_MODEL_
#define _OHARA_RUDY_MODEL_

#include "HeartIonicSolver.hpp"

namespace LifeV{

    template < typename Mesh,
    typename SolverType = LifeV::SolverAztecOO >
    class OHaraRudy : public virtual HeartIonicSolver<Mesh, SolverType>
    {
    public:
        typedef typename HeartIonicSolver<Mesh, SolverType>::data_Type data_Type;
        typedef typename HeartIonicSolver<Mesh, SolverType>::vector_Type vector_Type;
        typedef typename HeartIonicSolver<Mesh, SolverType>::function_Type function_Type;

    protected:

        vector_Type                    M_ionicCurrent;
        vector_Type                    M_ionicCurrentRepeated;
        vector_Type                    M_APD_flag;
        vector_Type                    M_vo;
        vector_Type                    M_vdot;

        VectorElemental              M_elemVecIonicCurrent;

        /*--------------------------------------------------*/
        /* ORd Human Model Basic Parameters     */
        /*--------------------------------------------------*/
        // Reversal Potentials
        vector_Type                    M_ENa;
        vector_Type                    M_EK;
        vector_Type                    M_EKs;
        // Currents
        vector_Type                    M_INa;
        vector_Type                    M_INaL;
        vector_Type                    M_Ito;
        vector_Type                    M_ICaL;
        vector_Type                    M_ICaNa;
        vector_Type                    M_ICaK;
        vector_Type                    M_IKr;
        vector_Type                    M_IKs;
        vector_Type                    M_IK1;
        vector_Type                    M_INaCa_i;
        vector_Type                    M_INaCa_ss;
        vector_Type                    M_INaCa;
        vector_Type                    M_INaK;
        vector_Type                    M_IKb;
        vector_Type                    M_INab;
        vector_Type                    M_ICab;
        vector_Type                    M_IpCa;
        // Fluxes
        vector_Type                    M_Jdiff; // diffusion current
        vector_Type                    M_JdiffNa;
        vector_Type                    M_JdiffK;
        vector_Type                    M_Jleak;
        vector_Type                    M_Jrel; //total Ca2+ release, via ryanodine receptors, from jsr to myoplasm
        vector_Type                    M_Jrelnp;
        vector_Type                    M_Jrelp;
        vector_Type                    M_Jup; //total Ca2+ uptake, via SERCA pump, from myoplasm to nsr
        vector_Type                    M_Jtr; //Ca2+ translocation from nsr to jsr
        // concentration of ions
        vector_Type                    M_nai;             // [Na+] _i
        vector_Type                    M_nass;          // [Na+]_ss
        vector_Type                    M_ki;               // [K+]_i
        vector_Type                    M_kss;            // [K+]_ss
        vector_Type                    M_cai;              // [Ca2+]_i
        vector_Type                    M_cass;           // [Ca2+]_ss
        vector_Type                    M_cansr;          // [Ca2+]_network SR
        vector_Type                    M_cajsr;           // [Ca2+]_junctional SR
        vector_Type                    M_CaMKa; // fraction of active CAMK binding sites
        vector_Type                    M_CaMKb; // fraction of CaMK bindings sites bound to Ca2+ / Calmodulin
        vector_Type                    M_CaMKt;  // fraction of autonomous CAMK binding sites with trapped calmodulin
        //------------------------------//
        // Time Dependent Gate
        //------------------------------//
        // Calmodulin
        vector_Type                    M_hf;// fast development of inactivation for CaMK phosphorylated fast INa
        vector_Type                    M_hs;// slow development of inactivation for CaMK phosphorylated fast INa
        vector_Type                    M_j;//recovery from inactivation for CaMK phosphorylated fast INa
        vector_Type                    M_hsp;
        vector_Type                    M_jp;
        // INa
        vector_Type                    M_m; // activation of  fast Na
        vector_Type                    M_mL; // activation for late INa
        vector_Type                    M_hL;  // inactivation for late INa
        vector_Type                    M_hLp; //
        // Ito
        vector_Type                    M_a; // activation for Ito
        vector_Type                    M_iF; // fast inactivation for Ito
        vector_Type                    M_iS; // slow inactivation for Ito
        vector_Type                    M_ap;
        vector_Type                    M_iFp;
        vector_Type                    M_iSp;
        // ICaL
        vector_Type                    M_d;
        vector_Type                    M_ff;
        vector_Type                    M_fs;
        vector_Type                    M_fcaf;
        vector_Type                    M_fcas;
        vector_Type                    M_jca;
        vector_Type                    M_nca;
        vector_Type                    M_ffp;
        vector_Type                    M_fcafp;
        // IKr
        vector_Type                    M_xrf;
        vector_Type                    M_xrs;
        vector_Type                    M_xs1;
        vector_Type                    M_xs2;
        vector_Type                    M_xk1;

        // Temporary variables
        Real _dt, _ENa,_EK,_EKs,_INa,_INaL,_Ito,_ICaL,_ICaNa,_ICaK,_IKr,_IKs,_IK1,_INaCa_i,_INaCa_ss,_INaCa,_INaK,_IKb,_INab,_ICab,_IpCa,_Jdiff,_JdiffNa,_JdiffK,_Jleak,_Jrel,_Jrelnp,_Jrelp,_Jup,_Jtr,_nai,_nass,_ki,_kss,_cai,_cass,_cansr,_cajsr,_CaMKa,_CaMKb,_CaMKt,_hf,_hs,_j,_hsp,_jp,_m,_mL,_hL,_hLp,_a,_iF,_iS,_ap,_iFp,_iSp,_d,_ff,_fs,_fcaf,_fcas,_jca,_nca,_ffp,_fcafp,_xrf,_xrs,_xs1,_xs2,_xk1, _ionicCurrent, _APD_flag, _vo, _vdot;

        // const parameters
        const int k_celltype;
        const Real k_nao;
        const Real k_cao;
        const Real k_ko;
        //buffer paramaters
        const Real k_BSRmax;
        const Real k_KmBSR;
        const Real k_BSLmax;
        const Real k_KmBSL;
        const Real k_cmdnmax;
        const Real k_kmcmdn;
        const Real k_trpnmax;
        const Real k_kmtrpn;
        const Real k_csqnmax;
        const Real k_kmcsqn;
        //CaMK paramaters
        const Real k_aCaMK;
        const Real k_bCaMK;
        const Real k_CaMKo;
        const Real k_KmCaM;
        const Real k_KmCaMK;
        //physical constants;
        const Real k_R;
        const Real k_T;
        const Real k_F;
        //cell geometry
        const Real k_L;
        const Real k_rad;
        const Real k_vcell;
        const Real k_Ageo;
        const Real k_Acap;
        const Real k_vmyo;
        const Real k_vmito;
        const Real k_vsr;
        const Real k_vnsr;
        const Real k_vjsr;
        const Real k_vss;
        //
        const Real k_vrest;

    private:
        void getElementValue(Int ig);
        void setElementValue(Int ig);
        void globalAssembleAll();
        void calculateIonicCurrent(Real v);
        void calculateBuffers(Real v);

    public:
        OHaraRudy(
                  const data_Type&          dataType,
                  const Mesh&               mesh,
                  FESpace<Mesh, MapEpetra>& uFEspace,
                  Epetra_Comm&              comm );
        virtual ~OHaraRudy(){};
        void updateRepeated();
        void updateElementSolution ( UInt eleID);
        void solveIonicModel ( const vector_Type& u, const Real timeStep );
        void computeIonicCurrent ( Real Capacitance,
                                  VectorElemental& elvec,
                                  VectorElemental& elvec_u,
                                  FESpace<Mesh, MapEpetra>& uFESpace );
        void initialize();

    };

        //! Constructor
    template<typename Mesh, typename SolverType>
    OHaraRudy<Mesh, SolverType>::
    OHaraRudy ( const data_Type& dataType,
               const Mesh& mesh,
               FESpace<Mesh, MapEpetra>& uFEspace,
               Epetra_Comm& comm ):
    HeartIonicSolver<Mesh, SolverType> ( dataType, mesh, uFEspace, comm),
    M_APD_flag(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_vo(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_vdot(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_ionicCurrent(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_ionicCurrentRepeated(M_ionicCurrent, Repeated),
    M_elemVecIonicCurrent ( HeartIonicSolver<Mesh, SolverType>::M_uFESpace.fe().nbFEDof(), 1 ),
    M_ENa(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_EK(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_EKs(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_INa(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_INaL(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_Ito(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_ICaL(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_ICaNa(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_ICaK(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_IKr(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_IKs(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_IK1(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_INaCa_i(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_INaCa_ss(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_INaCa(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_INaK(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_IKb(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_INab(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_ICab(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_IpCa(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_Jdiff(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_JdiffNa(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_JdiffK(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_Jleak(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_Jrel(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_Jrelnp(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_Jrelp(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_Jup(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_Jtr(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_nai(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_nass(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_ki(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_kss(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_cai(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_cass(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_cansr(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_cajsr(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_CaMKa(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_CaMKb(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_CaMKt(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_hf(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_hs(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_j(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_hsp(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_jp(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_m(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_mL(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_hL(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_hLp(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_a(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_iF(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_iS(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_ap(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_iFp(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_iSp(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_d(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_ff(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_fs(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_fcaf(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_fcas(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_jca(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_nca(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_ffp(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_fcafp(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_xrf(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_xrs(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_xs1(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_xs2(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    M_xk1(HeartIonicSolver<Mesh, SolverType>::M_localMap),
    // timeStep
    _dt(0.005),
    //initial values for state variables, there are 41 of them
    _ENa(0),
    _EK(0),
    _EKs(0),
    _INa(0),
    _INaL(0),
    _Ito(0),
    _ICaL(0),
    _ICaNa(0),
    _ICaK(0),
    _IKr(0),
    _IKs(0),
    _IK1(0),
    _INaCa_i(0),
    _INaCa_ss(0),
    _INaCa(0),
    _INaK(0),
    _IKb(0),
    _INab(0),
    _ICab(0),
    _IpCa(0),
    _Jdiff(0),
    _JdiffNa(0),
    _JdiffK(0),
    _Jleak(0),
    _Jrel(0),
    _Jrelnp(0),
    _Jrelp(0),
    _Jup(0),
    _Jtr(0),
    _nai(7),
    _nass(_nai),
    _ki(145),
    _kss(_ki),
    _cai(1.0e-4),
    _cass(_cai),
    _cansr(1.2),
    _cajsr(_cansr),
    _CaMKa(0),
    _CaMKb(0),
    _CaMKt(0),
    _hf(1),
    _hs(1),
    _j(1),
    _hsp(1),
    _jp(1),
    _m(0),
    _mL(0),
    _hL(1),
    _hLp(1),
    _a(0),
    _iF(1),
    _iS(1),
    _ap(0),
    _iFp(1),
    _iSp(1),
    _d(0),
    _ff(1),
    _fs(1),
    _fcaf(1),
    _fcas(1),
    _jca(1),
    _nca(0),
    _ffp(1),
    _fcafp(1),
    _xrf(0),
    _xrs(0),
    _xs1(0),
    _xs2(0),
    _xk1(1),
    _ionicCurrent(0.0),
    _vo(-83.9),
    _vdot(0.0),
    _APD_flag(0.0),
    // const values
    k_celltype(0),
    k_nao(140.0),
    k_cao(1.8),
    k_ko(5.4),
    k_BSRmax(0.047),
    k_KmBSR(0.00087),
    k_BSLmax(1.124),
    k_KmBSL(0.0087),
    k_cmdnmax(0.05),
    k_kmcmdn(0.00238),
    k_trpnmax(0.07),
    k_kmtrpn(0.0005),
    k_csqnmax(10.0),
    k_kmcsqn(0.8),
    k_aCaMK(0.05),
    k_bCaMK(0.00068),
    k_CaMKo(0.05),
    k_KmCaM(0.0015),
    k_KmCaMK(0.15),
    k_R(8314.0),
    k_T(310.0),
    k_F(96485.0),
    k_L(0.01),
    k_rad(0.0011),
    k_vcell(1000 * 3.14 * k_rad * k_rad * k_L),
    k_Ageo(2 * 3.14 * k_rad * k_rad + 2 * 3.14 * k_rad * k_L),
    k_Acap(2 * k_Ageo),
    k_vmyo(0.68 * k_vcell),
    k_vmito(0.26 * k_vcell),
    k_vsr(0.06 * k_vcell),
    k_vnsr(0.0552 * k_vcell),
    k_vjsr(0.0048 * k_vcell),
    k_vss(0.02 * k_vcell),
    k_vrest(-87.5)
    {}

    template<typename Mesh, typename SolverType>
    void OHaraRudy<Mesh, SolverType>::
    getElementValue(Int ig){
        _nai = M_nai[ig];
        _nass = M_nass[ig];
        _ki = M_ki[ig];
        _kss = M_kss[ig];
        _cai = M_cai[ig];
        _cass = M_cass[ig];
        _cansr = M_cansr[ig];
        _cajsr = M_cajsr[ig];
        _m = M_m[ig];
        _hf = M_hf[ig];
        _hs = M_hs[ig];
        _j = M_j[ig];
        _hsp = M_hsp[ig];
        _jp = M_jp[ig];
        _mL = M_mL[ig];
        _hL = M_hL[ig];
        _hLp = M_hLp[ig];
        _a = M_a[ig];
        _iF = M_iF[ig];
        _iS = M_iS[ig];
        _ap = M_ap[ig];
        _iFp = M_iFp[ig];
        _iSp = M_iSp[ig];
        _d = M_d[ig];
        _ff = M_ff[ig];
        _fs = M_fs[ig];
        _fcaf = M_fcaf[ig];
        _fcas = M_fcas[ig];
        _jca = M_jca[ig];
        _nca = M_nca[ig];
        _ffp = M_ffp[ig];
        _fcafp = M_fcafp[ig];
        _xrf = M_xrf[ig];
        _xrs = M_xrs[ig];
        _xs1 = M_xs1[ig];
        _xs2 = M_xs2[ig];
        _xk1 = M_xk1[ig];
        _Jrelnp = M_Jrelnp[ig];
        _Jrelp = M_Jrelp[ig];
        _CaMKt = M_CaMKt[ig];
        _ionicCurrent = M_ionicCurrent[ig];
         _APD_flag = M_APD_flag[ig];
         _vo = M_vo[ig];
         _vdot = M_vdot[ig];

        return;
    }

    template<typename Mesh, typename SolverType>
    void OHaraRudy<Mesh, SolverType>::
    setElementValue(Int ig){
        M_ENa.epetraVector().ReplaceGlobalValue( ig, 0, _ENa);
        M_EK.epetraVector().ReplaceGlobalValue( ig, 0, _EK);
        M_EKs.epetraVector().ReplaceGlobalValue( ig, 0, _EKs);
        M_INa.epetraVector().ReplaceGlobalValue( ig, 0, _INa);
        M_INaL.epetraVector().ReplaceGlobalValue( ig, 0, _INaL);
        M_Ito.epetraVector().ReplaceGlobalValue( ig, 0, _Ito);
        M_ICaL.epetraVector().ReplaceGlobalValue( ig, 0, _ICaL);
        M_ICaNa.epetraVector().ReplaceGlobalValue( ig, 0, _ICaNa);
        M_ICaK.epetraVector().ReplaceGlobalValue( ig, 0, _ICaK);
        M_IKr.epetraVector().ReplaceGlobalValue( ig, 0, _IKr);
        M_IKs.epetraVector().ReplaceGlobalValue( ig, 0, _IKs);
        M_IK1.epetraVector().ReplaceGlobalValue( ig, 0, _IK1);
        M_INaCa_i.epetraVector().ReplaceGlobalValue( ig, 0, _INaCa_i);
        M_INaCa_ss.epetraVector().ReplaceGlobalValue( ig, 0, _INaCa_ss);
        M_INaCa.epetraVector().ReplaceGlobalValue( ig, 0, _INaCa);
        M_INaK.epetraVector().ReplaceGlobalValue( ig, 0, _INaK);
        M_IKb.epetraVector().ReplaceGlobalValue( ig, 0, _IKb);
        M_INab.epetraVector().ReplaceGlobalValue( ig, 0, _INab);
        M_ICab.epetraVector().ReplaceGlobalValue( ig, 0, _ICab);
        M_IpCa.epetraVector().ReplaceGlobalValue( ig, 0, _IpCa);
        M_Jdiff.epetraVector().ReplaceGlobalValue( ig, 0, _Jdiff);
        M_JdiffNa.epetraVector().ReplaceGlobalValue( ig, 0, _JdiffNa);
        M_JdiffK.epetraVector().ReplaceGlobalValue( ig, 0, _JdiffK);
        M_Jleak.epetraVector().ReplaceGlobalValue( ig, 0, _Jleak);
        M_Jrel.epetraVector().ReplaceGlobalValue( ig, 0, _Jrel);
        M_Jrelnp.epetraVector().ReplaceGlobalValue( ig, 0, _Jrelnp);
        M_Jrelp.epetraVector().ReplaceGlobalValue( ig, 0, _Jrelp);
        M_Jup.epetraVector().ReplaceGlobalValue( ig, 0, _Jup);
        M_Jtr.epetraVector().ReplaceGlobalValue( ig, 0, _Jtr);
        M_nai.epetraVector().ReplaceGlobalValue( ig, 0, _nai);
        M_nass.epetraVector().ReplaceGlobalValue( ig, 0, _nass);
        M_ki.epetraVector().ReplaceGlobalValue( ig, 0, _ki);
        M_kss.epetraVector().ReplaceGlobalValue( ig, 0, _kss);
        M_cai.epetraVector().ReplaceGlobalValue( ig, 0, _cai);
        M_cass.epetraVector().ReplaceGlobalValue( ig, 0, _cass);
        M_cansr.epetraVector().ReplaceGlobalValue( ig, 0, _cansr);
        M_cajsr.epetraVector().ReplaceGlobalValue( ig, 0, _cajsr);
        M_CaMKa.epetraVector().ReplaceGlobalValue( ig, 0, _CaMKa);
        M_CaMKb.epetraVector().ReplaceGlobalValue( ig, 0, _CaMKb);
        M_CaMKt.epetraVector().ReplaceGlobalValue( ig, 0, _CaMKt);
        M_hf.epetraVector().ReplaceGlobalValue( ig, 0, _hf);
        M_hs.epetraVector().ReplaceGlobalValue( ig, 0, _hs);
        M_j.epetraVector().ReplaceGlobalValue( ig, 0, _j);
        M_hsp.epetraVector().ReplaceGlobalValue( ig, 0, _hsp);
        M_jp.epetraVector().ReplaceGlobalValue( ig, 0, _jp);
        M_m.epetraVector().ReplaceGlobalValue( ig, 0, _m);
        M_mL.epetraVector().ReplaceGlobalValue( ig, 0, _mL);
        M_hL.epetraVector().ReplaceGlobalValue( ig, 0, _hL);
        M_hLp.epetraVector().ReplaceGlobalValue( ig, 0, _hLp);
        M_a.epetraVector().ReplaceGlobalValue( ig, 0, _a);
        M_iF.epetraVector().ReplaceGlobalValue( ig, 0, _iF);
        M_iS.epetraVector().ReplaceGlobalValue( ig, 0, _iS);
        M_ap.epetraVector().ReplaceGlobalValue( ig, 0, _ap);
        M_iFp.epetraVector().ReplaceGlobalValue( ig, 0, _iFp);
        M_iSp.epetraVector().ReplaceGlobalValue( ig, 0, _iSp);
        M_d.epetraVector().ReplaceGlobalValue( ig, 0, _d);
        M_ff.epetraVector().ReplaceGlobalValue( ig, 0, _ff);
        M_fs.epetraVector().ReplaceGlobalValue( ig, 0, _fs);
        M_fcaf.epetraVector().ReplaceGlobalValue( ig, 0, _fcaf);
        M_fcas.epetraVector().ReplaceGlobalValue( ig, 0, _fcas);
        M_jca.epetraVector().ReplaceGlobalValue( ig, 0, _jca);
        M_nca.epetraVector().ReplaceGlobalValue( ig, 0, _nca);
        M_ffp.epetraVector().ReplaceGlobalValue( ig, 0, _ffp);
        M_fcafp.epetraVector().ReplaceGlobalValue( ig, 0, _fcafp);
        M_xrf.epetraVector().ReplaceGlobalValue( ig, 0, _xrf);
        M_xrs.epetraVector().ReplaceGlobalValue( ig, 0, _xrs);
        M_xs1.epetraVector().ReplaceGlobalValue( ig, 0, _xs1);
        M_xs2.epetraVector().ReplaceGlobalValue( ig, 0, _xs2);
        M_xk1.epetraVector().ReplaceGlobalValue( ig, 0, _xk1);
        M_ionicCurrent.epetraVector().ReplaceGlobalValue( ig, 0, _ionicCurrent );
        M_APD_flag.epetraVector().ReplaceGlobalValue( ig, 0, _APD_flag );
        M_vo.epetraVector().ReplaceGlobalValue( ig, 0, _vo );
        M_vdot.epetraVector().ReplaceGlobalValue( ig, 0, _vdot );
        return;
    }

    template<typename Mesh, typename SolverType>
    void OHaraRudy<Mesh, SolverType>::
    initialize( )
    {
        M_EK.epetraVector().PutScalar( _EK );
        M_EKs.epetraVector().PutScalar( _EKs );
        M_INa.epetraVector().PutScalar( _INa );
        M_INaL.epetraVector().PutScalar( _INaL );
        M_Ito.epetraVector().PutScalar( _Ito );
        M_ICaL.epetraVector().PutScalar( _ICaL );
        M_ICaNa.epetraVector().PutScalar( _ICaNa );
        M_ICaK.epetraVector().PutScalar( _ICaK );
        M_IKr.epetraVector().PutScalar( _IKr );
        M_IKs.epetraVector().PutScalar( _IKs );
        M_IK1.epetraVector().PutScalar( _IK1 );
        M_INaCa_i.epetraVector().PutScalar( _INaCa_i );
        M_INaCa_ss.epetraVector().PutScalar( _INaCa_ss );
        M_INaCa.epetraVector().PutScalar( _INaCa );
        M_INaK.epetraVector().PutScalar( _INaK );
        M_IKb.epetraVector().PutScalar( _IKb );
        M_INab.epetraVector().PutScalar( _INab );
        M_ICab.epetraVector().PutScalar( _ICab );
        M_IpCa.epetraVector().PutScalar( _IpCa );
        M_Jdiff.epetraVector().PutScalar( _Jdiff );
        M_JdiffNa.epetraVector().PutScalar( _JdiffNa );
        M_JdiffK.epetraVector().PutScalar( _JdiffK );
        M_Jleak.epetraVector().PutScalar( _Jleak );
        M_Jrel.epetraVector().PutScalar( _Jrel );
        M_Jrelnp.epetraVector().PutScalar( _Jrelnp );
        M_Jrelp.epetraVector().PutScalar( _Jrelp );
        M_Jup.epetraVector().PutScalar( _Jup );
        M_Jtr.epetraVector().PutScalar( _Jtr );
        M_nai.epetraVector().PutScalar( _nai );
        M_nass.epetraVector().PutScalar( _nass );
        M_ki.epetraVector().PutScalar( _ki );
        M_kss.epetraVector().PutScalar( _kss );
        M_cai.epetraVector().PutScalar( _cai );
        M_cass.epetraVector().PutScalar( _cass );
        M_cansr.epetraVector().PutScalar( _cansr );
        M_cajsr.epetraVector().PutScalar( _cajsr );
        M_CaMKa.epetraVector().PutScalar( _CaMKa );
        M_CaMKb.epetraVector().PutScalar( _CaMKb );
        M_CaMKt.epetraVector().PutScalar( _CaMKt );
        M_hf.epetraVector().PutScalar( _hf );
        M_hs.epetraVector().PutScalar( _hs );
        M_j.epetraVector().PutScalar( _j );
        M_hsp.epetraVector().PutScalar( _hsp );
        M_jp.epetraVector().PutScalar( _jp );
        M_m.epetraVector().PutScalar( _m );
        M_mL.epetraVector().PutScalar( _mL );
        M_hL.epetraVector().PutScalar( _hL );
        M_hLp.epetraVector().PutScalar( _hLp );
        M_a.epetraVector().PutScalar( _a );
        M_iF.epetraVector().PutScalar( _iF );
        M_iS.epetraVector().PutScalar( _iS );
        M_ap.epetraVector().PutScalar( _ap );
        M_iFp.epetraVector().PutScalar( _iFp );
        M_iSp.epetraVector().PutScalar( _iSp );
        M_d.epetraVector().PutScalar( _d );
        M_ff.epetraVector().PutScalar( _ff );
        M_fs.epetraVector().PutScalar( _fs );
        M_fcaf.epetraVector().PutScalar( _fcaf );
        M_fcas.epetraVector().PutScalar( _fcas );
        M_jca.epetraVector().PutScalar( _jca );
        M_nca.epetraVector().PutScalar( _nca );
        M_ffp.epetraVector().PutScalar( _ffp );
        M_fcafp.epetraVector().PutScalar( _fcafp );
        M_xrf.epetraVector().PutScalar( _xrf );
        M_xrs.epetraVector().PutScalar( _xrs );
        M_xs1.epetraVector().PutScalar( _xs1 );
        M_xs2.epetraVector().PutScalar( _xs2 );
        M_xk1.epetraVector().PutScalar( _xk1 );
        M_ionicCurrent.epetraVector().PutScalar( _ionicCurrent );
        M_APD_flag.epetraVector().PutScalar( _APD_flag );
        M_vo.epetraVector().PutScalar( _vo );
        M_vdot.epetraVector().PutScalar( _vdot );
        globalAssembleAll();

    }

    template<typename Mesh, typename SolverType>
    void OHaraRudy<Mesh, SolverType>::
    globalAssembleAll(){
        M_ENa.globalAssemble();
        M_EK.globalAssemble();
        M_EKs.globalAssemble();
        M_INa.globalAssemble();
        M_INaL.globalAssemble();
        M_Ito.globalAssemble();
        M_ICaL.globalAssemble();
        M_ICaNa.globalAssemble();
        M_ICaK.globalAssemble();
        M_IKr.globalAssemble();
        M_IKs.globalAssemble();
        M_IK1.globalAssemble();
        M_INaCa_i.globalAssemble();
        M_INaCa_ss.globalAssemble();
        M_INaCa.globalAssemble();
        M_INaK.globalAssemble();
        M_IKb.globalAssemble();
        M_INab.globalAssemble();
        M_ICab.globalAssemble();
        M_IpCa.globalAssemble();
        M_Jdiff.globalAssemble();
        M_JdiffNa.globalAssemble();
        M_JdiffK.globalAssemble();
        M_Jleak.globalAssemble();
        M_Jrel.globalAssemble();
        M_Jrelnp.globalAssemble();
        M_Jrelp.globalAssemble();
        M_Jup.globalAssemble();
        M_Jtr.globalAssemble();
        M_nai.globalAssemble();
        M_nass.globalAssemble();
        M_ki.globalAssemble();
        M_kss.globalAssemble();
        M_cai.globalAssemble();
        M_cass.globalAssemble();
        M_cansr.globalAssemble();
        M_cajsr.globalAssemble();
        M_CaMKa.globalAssemble();
        M_CaMKb.globalAssemble();
        M_CaMKt.globalAssemble();
        M_hf.globalAssemble();
        M_hs.globalAssemble();
        M_j.globalAssemble();
        M_hsp.globalAssemble();
        M_jp.globalAssemble();
        M_m.globalAssemble();
        M_mL.globalAssemble();
        M_hL.globalAssemble();
        M_hLp.globalAssemble();
        M_a.globalAssemble();
        M_iF.globalAssemble();
        M_iS.globalAssemble();
        M_ap.globalAssemble();
        M_iFp.globalAssemble();
        M_iSp.globalAssemble();
        M_d.globalAssemble();
        M_ff.globalAssemble();
        M_fs.globalAssemble();
        M_fcaf.globalAssemble();
        M_fcas.globalAssemble();
        M_jca.globalAssemble();
        M_nca.globalAssemble();
        M_ffp.globalAssemble();
        M_fcafp.globalAssemble();
        M_xrf.globalAssemble();
        M_xrs.globalAssemble();
        M_xs1.globalAssemble();
        M_xs2.globalAssemble();
        M_xk1.globalAssemble();
        M_ionicCurrent.globalAssemble();
        M_APD_flag.globalAssemble();
        M_vo.globalAssemble();
        M_vdot.globalAssemble();
        return;
    }

    template<typename Mesh, typename SolverType>
    void OHaraRudy<Mesh, SolverType>::
    calculateIonicCurrent(Real v)
    {

        /*----------------------------------*/
        /*         Reverse Potential       */
        /*----------------------------------*/

        _ENa=(k_R * k_T/k_F) * log(k_nao/_nai);
        _EK=(k_R * k_T/k_F) * log(k_ko/_ki);
        _EKs=(k_R * k_T/k_F) * log((k_ko + 0.01833 * k_nao)/(_ki + 0.01833 * _nai));

        _CaMKb=k_CaMKo * (1.0 - _CaMKt)/(1.0 + k_KmCaM/_cass);
        _CaMKa=_CaMKb + _CaMKt;

        /*----------------------------------*/
        /*      Ionic Gate & Current      */
        /*----------------------------------*/

        Real vffrt=v * k_F * k_F/(k_R * k_T);
        Real vfrt=v * k_F/(k_R * k_T);
        Real mss=1.0/(1.0 + exp(( - (v + 39.57))/9.871));
        Real tm=1.0/(6.765 * exp((v + 11.64)/34.77) + 8.552 * exp( - (v + 77.42)/5.955));
        _m=mss - (mss - _m) * exp( - _dt/tm);

        Real hss=1.0/(1 + exp((v + 82.90)/6.086));
        Real thf=1.0/(1.432e-5 * exp( - (v + 1.196)/6.285) + 6.149 * exp((v + 0.5096)/20.27));
        Real ths=1.0/(0.009794 * exp( - (v + 17.95)/28.05) + 0.3343 * exp((v + 5.730)/56.66));
        Real Ahf=0.99;
        Real Ahs=1.0 - Ahf;
        _hf=hss - (hss - _hf) * exp( - _dt/thf);
        _hs=hss - (hss - _hs) * exp( - _dt/ths);

        Real h=Ahf * _hf + Ahs * _hs;
        Real jss=hss;
        Real tj=2.038 + 1.0/(0.02136 * exp( - (v + 100.6)/8.281) + 0.3052 * exp((v + 0.9941)/38.45));
        _j=jss - (jss - _j) * exp( - _dt/tj);

        Real hssp=1.0/(1 + exp((v + 89.1)/6.086));
        Real thsp=3.0 * ths;
        _hsp=hssp - (hssp - _hsp) * exp( - _dt/thsp);

        Real hp=Ahf * _hf + Ahs * _hsp;
        Real tjp=1.46 * tj;
        _jp=jss - (jss - _jp) * exp( - _dt/tjp);

        Real GNa=75;
        Real fINap=(1.0/(1.0 + k_KmCaMK/_CaMKa));
        _INa=GNa * (v - _ENa) * _m * _m * _m * ((1.0 - fINap) * h*_j + fINap * hp * _jp);

        Real mLss=1.0/(1.0 + exp(( - (v + 42.85))/5.264));
        Real tmL=tm;
        _mL=mLss - (mLss - _mL) * exp( - _dt/tmL);

        Real hLss=1.0/(1.0 + exp((v + 87.61)/7.488));
        Real thL=200.0;
        _hL=hLss - (hLss - _hL) * exp( - _dt/thL);

        Real hLssp=1.0/(1.0 + exp((v + 93.81)/7.488));
        Real thLp=3.0 * thL;
        _hLp=hLssp - (hLssp - _hLp) * exp( - _dt/thLp);

        Real GNaL=0.0075;
        if (k_celltype==1)
            GNaL *= 0.6;
        Real fINaLp=(1.0/(1.0 + k_KmCaMK/_CaMKa));
        _INaL=GNaL * (v - _ENa) * _mL * ((1.0 - fINaLp) * _hL + fINaLp * _hLp);

        Real ass=1.0/(1.0 + exp(( - (v - 14.34))/14.82));
        Real ta=1.0515/(1.0/(1.2089 * (1.0 + exp( - (v - 18.4099)/29.3814))) + 3.5/(1.0 + exp((v + 100.0)/29.3814)));
        _a=ass - (ass - _a) * exp( - _dt/ta);

        Real iss=1.0/(1.0 + exp((v + 43.94)/5.711));
        Real delta_epi;
        if (k_celltype==1){
            delta_epi=1.0 - (0.95/(1.0 + exp((v + 70.0)/5.0)));
        }else{
            delta_epi=1.0;
        }
        Real tiF=4.562 + 1/(0.3933 * exp(( - (v + 100.0))/100.0) + 0.08004 * exp((v + 50.0)/16.59));
        Real tiS=23.62 + 1/(0.001416 * exp(( - (v + 96.52))/59.05) + 1.780e-8 * exp((v + 114.1)/8.079));
        tiF *= delta_epi;
        tiS *= delta_epi;
        Real AiF=1.0/(1.0 + exp((v - 213.6)/151.2));
        Real AiS=1.0 - AiF;
        _iF=iss - (iss - _iF) * exp( - _dt/tiF);
        _iS=iss - (iss - _iS) * exp( - _dt/tiS);

        Real i=AiF * _iF + AiS * _iS;
        Real assp=1.0/(1.0 + exp(( - (v - 24.34))/14.82));
        _ap=assp - (assp - _ap) * exp( - _dt/ta);

        Real dti_develop=1.354 + 1.0e-4/(exp((v - 167.4)/15.89) + exp( - (v - 12.23)/0.2154));
        Real dti_recover=1.0 - 0.5/(1.0 + exp((v + 70.0)/20.0));
        Real tiFp=dti_develop * dti_recover * tiF;
        Real tiSp=dti_develop * dti_recover * tiS;
        _iFp=iss - (iss - _iFp) * exp( - _dt/tiFp);
        _iSp=iss - (iss - _iSp) * exp( - _dt/tiSp);

        Real ip=AiF * _iFp + AiS * _iSp;
        Real Gto=0.02;
        if (k_celltype==1)
            Gto *= 4.0;
        if (k_celltype==2)
            Gto *= 4.0;
        Real fItop=(1.0/(1.0 + k_KmCaMK/_CaMKa));
        _Ito=Gto * (v - _EK) * ((1.0 - fItop) * _a * i + fItop * _ap * ip);

        Real dss=1.0/(1.0 + exp(( - (v + 3.940))/4.230));
        Real td=0.6 + 1.0/(exp( - 0.05 * (v + 6.0)) + exp(0.09 * (v + 14.0)));
        _d=dss - (dss - _d) * exp( - _dt/td);

        Real fss=1.0/(1.0 + exp((v + 19.58)/3.696));
        Real tff=7.0 + 1.0/(0.0045 * exp( - (v + 20.0)/10.0) + 0.0045 * exp((v + 20.0)/10.0));
        Real tfs=1000.0 + 1.0/(0.000035 * exp( - (v + 5.0)/4.0) + 0.000035 * exp((v + 5.0)/6.0));
        Real Aff=0.6;
        Real Afs=1.0 - Aff;
        _ff=fss - (fss - _ff) * exp( - _dt/tff);
        _fs=fss - (fss - _fs) * exp( - _dt/tfs);

        Real f=Aff * _ff + Afs * _fs;
        Real fcass=fss;
        Real tfcaf=7.0 + 1.0/(0.04 * exp( - (v - 4.0)/7.0) + 0.04 * exp((v - 4.0)/7.0));
        Real tfcas=100.0 + 1.0/(0.00012 * exp( - v/3.0) + 0.00012 * exp(v/7.0));
        Real Afcaf=0.3 + 0.6/(1.0 + exp((v - 10.0)/10.0));
        Real Afcas=1.0 - Afcaf;
        _fcaf=fcass - (fcass - _fcaf) * exp( - _dt/tfcaf);
        _fcas=fcass - (fcass - _fcas) * exp( - _dt/tfcas);

        Real fca=Afcaf * _fcaf + Afcas * _fcas;
        Real tjca=75.0;
        _jca=fcass - (fcass - _jca) * exp( - _dt/tjca);

        Real tffp=2.5 * tff;
        _ffp=fss - (fss - _ffp) * exp( - _dt/tffp);

        Real fp=Aff * _ffp + Afs * _fs;
        Real tfcafp=2.5 * tfcaf;
        _fcafp=fcass - (fcass - _fcafp) * exp( - _dt/tfcafp);

        Real fcap=Afcaf * _fcafp + Afcas * _fcas;
        Real Kmn=0.002;
        Real k2n=1000.0;
        Real km2n=_jca * 1.0;
        Real anca=1.0/(k2n/km2n + pow(1.0 + Kmn/_cass,4.0));
        _nca=anca * k2n/km2n - (anca * k2n/km2n - _nca) * exp( - km2n * _dt);

        Real PhiCaL=4.0 * vffrt * (_cass * exp(2.0 * vfrt) - 0.341 * k_cao)/(exp(2.0 * vfrt) - 1.0);
        Real PhiCaNa=1.0 * vffrt * (0.75 * _nass * exp(1.0 * vfrt) - 0.75 * k_nao)/(exp(1.0 * vfrt) - 1.0);
        Real PhiCaK=1.0 * vffrt * (0.75 * _kss * exp(1.0 * vfrt) - 0.75 * k_ko)/(exp(1.0 * vfrt) - 1.0);
        Real zca=2.0;
        Real PCa=0.0001;
        if (k_celltype==1)
            PCa *= 1.2;
        if (k_celltype==2)
            PCa *= 2.5;
        Real PCap=1.1 * PCa;
        Real PCaNa=0.00125 * PCa;
        Real PCaK=3.574e-4 * PCa;
        Real PCaNap=0.00125 * PCap;
        Real PCaKp=3.574e-4 * PCap;
        Real fICaLp=(1.0/(1.0 + k_KmCaMK/_CaMKa));
        _ICaL=(1.0 - fICaLp) * PCa * PhiCaL * _d * (f * (1.0 - _nca) + _jca * fca * _nca) + fICaLp * PCap * PhiCaL * _d * (fp * (1.0 - _nca) + _jca * fcap * _nca);
        _ICaNa=(1.0 - fICaLp) * PCaNa * PhiCaNa * _d * (f * (1.0 - _nca) + _jca * fca * _nca) + fICaLp * PCaNap * PhiCaNa * _d * (fp * (1.0 - _nca) + _jca * fcap * _nca);
        _ICaK=(1.0 - fICaLp) * PCaK * PhiCaK * _d * (f * (1.0 - _nca) + _jca * fca * _nca) + fICaLp * PCaKp * PhiCaK * _d * (fp * (1.0 - _nca) + _jca * fcap * _nca);

        Real xrss=1.0/(1.0 + exp(( - (v + 8.337))/6.789));
        Real txrf=12.98 + 1.0/(0.3652 * exp((v - 31.66)/3.869) + 4.123e-5 * exp(( - (v - 47.78))/20.38));
        Real txrs=1.865 + 1.0/(0.06629 * exp((v - 34.70)/7.355) + 1.128e-5 * exp(( - (v - 29.74))/25.94));
        Real Axrf=1.0/(1.0 + exp((v + 54.81)/38.21));
        Real Axrs=1.0 - Axrf;
        _xrf=xrss - (xrss - _xrf) * exp( - _dt/txrf);
        _xrs=xrss - (xrss - _xrs) * exp( - _dt/txrs);

        Real xr=Axrf * _xrf + Axrs * _xrs;
        Real rkr=1.0/(1.0 + exp((v + 55.0)/75.0)) * 1.0/(1.0 + exp((v - 10.0)/30.0));
        Real GKr=0.046;
        if (k_celltype==1)
            GKr *= 1.3;
        if (k_celltype==2)
            GKr *= 0.8;
        _IKr=GKr * sqrt(k_ko/5.4) * xr * rkr * (v - _EK);

        Real xs1ss=1.0/(1.0 + exp(( - (v + 11.60))/8.932));
        Real txs1=817.3 + 1.0/(2.326e-4 * exp((v + 48.28)/17.80) + 0.001292 * exp(( - (v + 210.0))/230.0));
        _xs1=xs1ss - (xs1ss - _xs1) * exp( - _dt/txs1);

        Real xs2ss=xs1ss;
        Real txs2=1.0/(0.01 * exp((v - 50.0)/20.0) + 0.0193 * exp(( - (v + 66.54))/31.0));
        _xs2=xs2ss - (xs2ss - _xs2) * exp( - _dt/txs2);

        Real KsCa=1.0 + 0.6/(1.0 + pow(3.8e-5/_cai,1.4));
        Real GKs=0.0034;
        if (k_celltype==1)
            GKs *= 1.4;
        _IKs=GKs * KsCa * _xs1 * _xs2 * (v - _EKs);

        Real xk1ss=1.0/(1.0 + exp( - (v + 2.5538 * k_ko + 144.59)/(1.5692 * k_ko + 3.8115)));
        Real txk1=122.2/(exp(( - (v + 127.2))/20.36) + exp((v + 236.8)/69.33));
        _xk1=xk1ss - (xk1ss - _xk1) * exp( - _dt/txk1);

        Real rk1=1.0/(1.0 + exp((v + 105.8 - 2.6 * k_ko)/9.493));
        Real GK1=0.1908;
        if (k_celltype==1)
            GK1 *= 1.2;
        if (k_celltype==2)
            GK1 *= 1.3;
        _IK1=GK1 * sqrt(k_ko) * rk1 * _xk1 * (v - _EK);

        Real kna1=15.0;
        Real kna2=5.0;
        Real kna3=88.12;
        Real kasymm=12.5;
        Real wna=6.0e4;
        Real wca=6.0e4;
        Real wnaca=5.0e3;
        Real kcaon=1.5e6;
        Real kcaoff=5.0e3;
        Real qna=0.5224;
        Real qca=0.1670;
        Real hca=exp((qca * v * k_F)/(k_R * k_T));
        Real hna=exp((qna * v * k_F)/(k_R * k_T));
        Real h1=1 + _nai/kna3 * (1 + hna);
        Real h2=(_nai * hna)/(kna3 * h1);
        Real h3=1.0/h1;
        Real h4=1.0 + _nai/kna1 * (1 + _nai/kna2);
        Real h5=_nai * _nai/(h4 * kna1 * kna2);
        Real h6=1.0/h4;
        Real h7=1.0 + k_nao/kna3 * (1.0 + 1.0/hna);
        Real h8=k_nao/(kna3 * hna * h7);
        Real h9=1.0/h7;
        Real h10=kasymm + 1.0 + k_nao/kna1 * (1.0 + k_nao/kna2);
        Real h11=k_nao * k_nao/(h10 * kna1 * kna2);
        Real h12=1.0/h10;
        Real k1=h12 * k_cao * kcaon;
        Real k2=kcaoff;
        Real k3p=h9 * wca;
        Real k3pp=h8 * wnaca;
        Real k3=k3p + k3pp;
        Real k4p=h3 * wca/hca;
        Real k4pp=h2 * wnaca;
        Real k4=k4p + k4pp;
        Real k5=kcaoff;
        Real k6=h6 * _cai * kcaon;
        Real k7=h5 * h2 * wna;
        Real k8=h8 * h11 * wna;
        Real x1=k2 * k4 * (k7 + k6) + k5 * k7 * (k2 + k3);
        Real x2=k1 * k7 * (k4 + k5) + k4 * k6 * (k1 + k8);
        Real x3=k1 * k3 * (k7 + k6) + k8 * k6 * (k2 + k3);
        Real x4=k2 * k8 * (k4 + k5) + k3 * k5 * (k1 + k8);
        Real E1=x1/(x1 + x2 + x3 + x4);
        Real E2=x2/(x1 + x2 + x3 + x4);
        Real E3=x3/(x1 + x2 + x3 + x4);
        Real E4=x4/(x1 + x2 + x3 + x4);
        Real KmCaAct=150.0e-6;
        Real allo=1.0/(1.0 + pow(KmCaAct/_cai,2.0));
        Real zna=1.0;
        Real JncxNa=3.0 * (E4 * k7 - E1 * k8) + E3 * k4pp - E2 * k3pp;
        Real JncxCa=E2 * k2 - E1 * k1;
        Real Gncx=0.0008;
        if (k_celltype==1)
            Gncx *= 1.1;
        if (k_celltype==2)
            Gncx *= 1.4;
        _INaCa_i=0.8 * Gncx * allo * (zna * JncxNa + zca * JncxCa);

        h1=1 + _nass/kna3 * (1 + hna);
        h2=(_nass * hna)/(kna3 * h1);
        h3=1.0/h1;
        h4=1.0 + _nass/kna1 * (1 + _nass/kna2);
        h5=_nass * _nass/(h4 * kna1 * kna2);
        h6=1.0/h4;
        h7=1.0 + k_nao/kna3 * (1.0 + 1.0/hna);
        h8=k_nao/(kna3 * hna * h7);
        h9=1.0/h7;
        h10=kasymm + 1.0 + k_nao/kna1 * (1 + k_nao/kna2);
        h11=k_nao * k_nao/(h10 * kna1 * kna2);
        h12=1.0/h10;
        k1=h12 * k_cao * kcaon;
        k2=kcaoff;
        k3p=h9 * wca;
        k3pp=h8 * wnaca;
        k3=k3p + k3pp;
        k4p=h3 * wca/hca;
        k4pp=h2 * wnaca;
        k4=k4p + k4pp;
        k5=kcaoff;
        k6=h6 * _cass * kcaon;
        k7=h5 * h2 * wna;
        k8=h8 * h11 * wna;
        x1=k2 * k4 * (k7 + k6) + k5 * k7 * (k2 + k3);
        x2=k1 * k7 * (k4 + k5) + k4 * k6 * (k1 + k8);
        x3=k1 * k3 * (k7 + k6) + k8 * k6 * (k2 + k3);
        x4=k2 * k8 * (k4 + k5) + k3 * k5 * (k1 + k8);
        E1=x1/(x1 + x2 + x3 + x4);
        E2=x2/(x1 + x2 + x3 + x4);
        E3=x3/(x1 + x2 + x3 + x4);
        E4=x4/(x1 + x2 + x3 + x4);
        KmCaAct=150.0e-6;
        allo=1.0/(1.0 + pow(KmCaAct/_cass,2.0));
        JncxNa=3.0 * (E4 * k7 - E1 * k8) + E3 * k4pp - E2 * k3pp;
        JncxCa=E2 * k2 - E1 * k1;
        _INaCa_ss=0.2 * Gncx * allo * (zna * JncxNa + zca * JncxCa);
        _INaCa=_INaCa_i + _INaCa_ss;

        Real k1p=949.5;
        Real k1m=182.4;
        Real k2p=687.2;
        Real k2m=39.4;
        k3p=1899.0;
        Real k3m=79300.0;
        k4p=639.0;
        Real k4m=40.0;
        Real Knai0=9.073;
        Real Knao0=27.78;
        Real delta= - 0.1550;
        Real Knai=Knai0 * exp((delta * v * k_F)/(3.0 * k_R * k_T));
        Real Knao=Knao0 * exp(((1.0 - delta) * v * k_F)/(3.0 * k_R * k_T));
        Real Kki=0.5;
        Real Kko=0.3582;
        Real MgADP=0.05;
        Real MgATP=9.8;
        Real Kmgatp=1.698e-7;
        Real H=1.0e-7;
        Real eP=4.2;
        Real Khp=1.698e-7;
        Real Knap=224.0;
        Real Kxkur=292.0;
        Real P=eP/(1.0 + H/Khp + _nai/Knap + _ki/Kxkur);
        Real a1=(k1p * pow(_nai/Knai,3.0))/(pow(1.0 + _nai/Knai,3.0) + pow(1.0 + _ki/Kki,2.0) - 1.0);
        Real b1=k1m * MgADP;
        Real a2=k2p;
        Real b2=(k2m * pow(k_nao/Knao,3.0))/(pow(1.0 + k_nao/Knao,3.0) + pow(1.0 + k_ko/Kko,2.0) - 1.0);
        Real a3=(k3p * pow(k_ko/Kko,2.0))/(pow(1.0 + k_nao/Knao,3.0) + pow(1.0 + k_ko/Kko,2.0) - 1.0);
        Real b3=(k3m * P*H)/(1.0 + MgATP/Kmgatp);
        Real a4=(k4p * MgATP/Kmgatp)/(1.0 + MgATP/Kmgatp);
        Real b4=(k4m * pow(_ki/Kki,2.0))/(pow(1.0 + _nai/Knai,3.0) + pow(1.0 + _ki/Kki,2.0) - 1.0);
        x1=a4 * a1 * a2 + b2 * b4 * b3 + a2 * b4 * b3 + b3 * a1 * a2;
        x2=b2 * b1 * b4 + a1 * a2 * a3 + a3 * b1 * b4 + a2 * a3 * b4;
        x3=a2 * a3 * a4 + b3 * b2 * b1 + b2 * b1 * a4 + a3 * a4 * b1;
        x4=b4 * b3 * b2 + a3 * a4 * a1 + b2 * a4 * a1 + b3 * b2 * a1;
        E1=x1/(x1 + x2 + x3 + x4);
        E2=x2/(x1 + x2 + x3 + x4);
        E3=x3/(x1 + x2 + x3 + x4);
        E4=x4/(x1 + x2 + x3 + x4);
        Real zk=1.0;
        Real JnakNa=3.0 * (E1 * a3 - E2 * b3);
        Real JnakK=2.0 * (E4 * b1 - E3 * a1);
        Real Pnak=30;
        if (k_celltype==1)
            Pnak *= 0.9;
        if (k_celltype==2)
            Pnak *= 0.7;
        _INaK=Pnak * (zna * JnakNa + zk * JnakK);

        Real xkb=1.0/(1.0 + exp( - (v - 14.48)/18.34));
        Real GKb=0.003;
        _IKb=GKb * xkb * (v - _EK);

        Real PNab=3.75e-10;
        _INab=PNab * vffrt * (_nai * exp(vfrt) - k_nao)/(exp(vfrt) - 1.0);

        Real PCab=2.5e-8;
        _ICab=PCab * 4.0 * vffrt * (_cai * exp(2.0 * vfrt) - 0.341 * k_cao)/(exp(2.0 * vfrt) - 1.0);

        Real GpCa=0.0005;
        _IpCa=GpCa * _cai/(0.0005 + _cai);

        // Summation of all ionic currents
        _ionicCurrent = _INa+_INaL+_Ito+_ICaL+_ICaNa+_ICaK+_IKr+_IKs+_IK1+_INaCa+_INaK+_INab+_IKb+_IpCa+_ICab;
    }

    template<typename Mesh, typename SolverType>
    void OHaraRudy<Mesh, SolverType>::
    calculateBuffers(Real v)
    {

        /*----------------------------------*/
        /*                 Fluxes                 */
        /*----------------------------------*/

        Real _CaMKb=k_CaMKo * (1.0 - _CaMKt)/(1.0 + k_KmCaM/_cass);
        _CaMKa=_CaMKb + _CaMKt;
        _CaMKt +=_dt * (k_aCaMK * _CaMKb * (_CaMKb + _CaMKt) - k_bCaMK * _CaMKt);

        // Fluxes
        _JdiffNa=(_nass - _nai)/2.0;
        _JdiffK=(_kss - _ki)/2.0;
        _Jdiff=(_cass - _cai)/0.2;

        Real bt=4.75;
        Real a_rel=0.5 * bt;
        Real Jrel_inf=a_rel * ( - _ICaL)/(1.0 + pow(1.5/_cajsr,8.0));
        if (k_celltype==2)
            Jrel_inf *= 1.7;
        Real tau_rel=bt/(1.0 + 0.0123/_cajsr);
        if (tau_rel<0.005)
            tau_rel=0.005;
        _Jrelnp=Jrel_inf - (Jrel_inf - _Jrelnp) * exp( - _dt/tau_rel);

        Real btp=1.25 * bt;
        Real a_relp=0.5 * btp;
        Real Jrel_infp=a_relp * ( - _ICaL)/(1.0 + pow(1.5/_cajsr,8.0));
        if (k_celltype==2)
            Jrel_infp *= 1.7;
        Real tau_relp=btp/(1.0 + 0.0123/_cajsr);
        if (tau_relp<0.005)
            tau_relp=0.005;
        _Jrelp=Jrel_infp - (Jrel_infp - _Jrelp) * exp( - _dt/tau_relp);

        Real fJrelp=(1.0/(1.0 + k_KmCaMK/_CaMKa));
        _Jrel=(1.0 - fJrelp) * _Jrelnp + fJrelp * _Jrelp;

        Real Jupnp=0.004375 * _cai/(_cai + 0.00092);
        Real Jupp=2.75 * 0.004375 * _cai/(_cai + 0.00092 - 0.00017);
        if (k_celltype==1){
            Jupnp *= 1.3;
            Jupp *= 1.3;
        }
        Real fJupp=(1.0/(1.0 + k_KmCaMK/_CaMKa));
        _Jleak=0.0039375 * _cansr/15.0;
        _Jup=(1.0 - fJupp) * Jupnp + fJupp * Jupp - _Jleak;

        _Jtr=(_cansr - _cajsr)/100.0;

        /*----------------------------------*/
        /*           Concentrations         */
        /*----------------------------------*/

        _nai +=_dt * ( - (_INa + _INaL + 3.0 * _INaCa_i + 3.0 * _INaK + _INab) * k_Acap/(k_F * k_vmyo) + _JdiffNa * k_vss/k_vmyo);
        _nass +=_dt * ( - (_ICaNa + 3.0 * _INaCa_ss) * k_Acap/(k_F * k_vss) - _JdiffNa);
        _ki +=_dt * ( - (_Ito + _IKr + _IKs + _IK1 + _IKb - 2.0 * _INaK) * k_Acap/(k_F * k_vmyo) + _JdiffK * k_vss/k_vmyo);
        _kss +=_dt * ( - (_ICaK) * k_Acap/(k_F * k_vss) - _JdiffK);
        Real Bcai;
        if (k_celltype==1){
            Bcai=1.0/(1.0 + 1.3 * k_cmdnmax * k_kmcmdn/pow(k_kmcmdn + _cai,2.0) + k_trpnmax * k_kmtrpn/pow(k_kmtrpn + _cai,2.0));
        }else{
            Bcai=1.0/(1.0 + k_cmdnmax * k_kmcmdn/pow(k_kmcmdn + _cai,2.0) + k_trpnmax * k_kmtrpn/pow(k_kmtrpn + _cai,2.0));
        }
        _cai += _dt * (Bcai * ( - (_IpCa + _ICab - 2.0 * _INaCa_i) * k_Acap/(2.0 * k_F * k_vmyo) - _Jup * k_vnsr/k_vmyo + _Jdiff * k_vss/k_vmyo));

        Real Bcass=1.0/(1.0 + k_BSRmax * k_KmBSR/pow(k_KmBSR + _cass,2.0) + k_BSLmax * k_KmBSL/pow(k_KmBSL + _cass,2.0));
        _cass +=_dt * (Bcass * ( - (_ICaL - 2.0 * _INaCa_ss) * k_Acap/(2.0 * k_F * k_vss) + _Jrel * k_vjsr/k_vss - _Jdiff));

        _cansr +=_dt * (_Jup - _Jtr * k_vjsr/k_vnsr);
        Real Bcajsr=1.0/(1.0 + k_csqnmax * k_kmcsqn/pow(k_kmcsqn + _cajsr,2.0));
        _cajsr +=_dt * (Bcajsr * (_Jtr - _Jrel));

    }

    template<typename Mesh, typename SolverType>
    void OHaraRudy<Mesh, SolverType>::
    updateRepeated( )
    {
        M_ionicCurrentRepeated = M_ionicCurrent;
    }

    template<typename Mesh, typename SolverType>
    void OHaraRudy<Mesh, SolverType>::
    updateElementSolution ( UInt eleID)
    {
        M_elemVecIonicCurrent.zero();
        UInt ig;
        //! Filling local elvec with recovery variable values in the nodes
        for ( UInt iNode = 0 ; iNode < HeartIonicSolver<Mesh, SolverType>::M_uFESpace.fe().nbFEDof() ; iNode++ )
        {
            ig = HeartIonicSolver<Mesh, SolverType>::M_uFESpace.dof().localToGlobalMap ( eleID, iNode );
            M_elemVecIonicCurrent.vec() [ iNode ] = M_ionicCurrentRepeated[ig];
        }

    }

    template<typename Mesh, typename SolverType>
    void OHaraRudy<Mesh, SolverType>::
    solveIonicModel ( const vector_Type& u, const Real timeStep )
    {

        //! Solving dw/dt=eta2 (u/vp -  eta3 w)
        LifeChrono chronoionmodelsolve;
        chronoionmodelsolve.start();
        HeartIonicSolver<Mesh, SolverType>::M_comm->Barrier();

        _dt = timeStep;
        Real t;

        for ( Int i = 0 ; i < u.epetraVector().MyLength() ; i++ ){
            Int ig = u.blockMap().MyGlobalElements()[i];
            Real v = u[ig];
            getElementValue( ig );

            //std::cerr <<  "index " << ig ;

            /*----------------------------------*/
            /*               APD flags              */
            /*----------------------------------*/
            Real vdot_old=_vdot;
            _vdot = ( v - _vo ) / _dt;
            if( _APD_flag < 0. && v>-40 && _vdot < vdot_old ){
                _APD_flag=1.0;
            }
            if( _APD_flag > 0. && v <  0.9 * k_vrest ){
                _APD_flag=-1.0;
            }
            //std::cerr << "v: " << v << ", vo:" << _vo << ",";

            t = 0.0; // [sec]
            while( t < timeStep / 1000 ){

                // define dt
                if ( _APD_flag > 0.0 && v < 0.7* k_vrest ){
                    _dt = 0.005; // [sec]
                }else if (fabs(v-_vo)<0.2){
                    _dt=fabs(0.8/_vdot);
                    if (_dt>1.0) _dt=1.0;
                }else if (fabs(v-_vo)>0.8){
                    _dt = fabs(0.2/_vdot); // [msec]
                }
                if( _dt < timeStep / 100000 ) _dt = timeStep / 100000;

                t += _dt;
                //std::cerr <<  t << " < " << timeStep / 1000 << " @ " << _dt << std::endl;
                getElementValue( ig );
                calculateIonicCurrent( v );
                calculateBuffers( v );
                _vo = v;
                setElementValue( ig );
            }
            //std::cerr <<  "done " <<std::endl;

        }
        globalAssembleAll();
        return;
    }

    template<typename Mesh, typename SolverType>
    void OHaraRudy<Mesh, SolverType>::
    computeIonicCurrent (
                         Real Capacitance,
                         VectorElemental& elvec,
                         VectorElemental& /*elvec_u*/,
                         FESpace<Mesh, MapEpetra>& uFESpace )
    {
        Real Iion_ig;
        for ( UInt ig = 0; ig < uFESpace.fe().nbQuadPt(); ig++ )
        {
            Iion_ig = 0.;
            for ( UInt i = 0; i < uFESpace.fe().nbFEDof(); i++ )
            {
                Iion_ig += M_elemVecIonicCurrent ( i ) * uFESpace.fe().phi ( i, ig ); //
            }
            for ( UInt i = 0; i < uFESpace.fe().nbFEDof(); i++ )
            {
                // divide by 1000 to convert microA in mA
                elvec ( i ) -= Iion_ig * Capacitance * uFESpace.fe().phi ( i, ig ) * uFESpace.fe().weightDet ( ig );
            }
        }
    }

} // namespace LifeV

#endif //_OHARA_RUDY_MODEL_
