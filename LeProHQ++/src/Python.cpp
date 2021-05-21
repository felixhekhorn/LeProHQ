#include <boost/python.hpp>

#include "InclusiveLeptoProduction.h"
#include "FullyDiffLeptoProduction.h"
#include "Projection.hpp"
#include "Flags.hpp"
#include "FullyDiff/histT.hpp"

using namespace boost::python;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(FullyDiffLeptoProduction_activateHistogram_overloads, activateHistogram, 3, 5)

BOOST_PYTHON_MODULE(LeProHQpp)
{
    enum_<Projection>("Projection", "projection type")
        .value("F2", F2)
        .value("FL", FL)
        .value("x2g1", x2g1)
        .value("xF3", xF3)
        .value("g4", g4)
        .value("gL", gL)
    ;
    def("isParityConservingProj", &isParityConservingProj, "does projection conserve parity?");
    def("strProj", &strProj, "cast projection to string");
    
    class_<InclusiveLeptoProduction>("InclusiveLeptoProduction", "application class for full inclusive lepto production", init<uint,double,double>())
        // global getter and setter
        .def("setProjection", &InclusiveLeptoProduction::setProjection, "sets projection")
        .def("setNumberOfLightFlavours", &InclusiveLeptoProduction::setNumberOfLightFlavours, "sets number of light flavours")
        .def("setQ2", &InclusiveLeptoProduction::setQ2, "sets transferred momentum Q2 > 0")
        .def("setM2", &InclusiveLeptoProduction::setM2, "sets mass of heavy quark; unsets partonic center of mass energy")
        .def("setDelta", &InclusiveLeptoProduction::setDelta, "sets phase space slice")
        .def("getIntegrationConfig", &InclusiveLeptoProduction::getIntegrationConfig, "get matching IntegrationConfig", return_value_policy<reference_existing_object>())
        .def("getIntegrationOutput", &InclusiveLeptoProduction::getIntegrationOutput, "get current/latest IntegrationOutput")
        // partonic setter
        .def("setPartonicEta", &InclusiveLeptoProduction::setPartonicEta, "sets partonic eta")
        .def("setPartonicS", &InclusiveLeptoProduction::setPartonicS, "sets partonic s = (k1+q)²")
        // partonic coefficient functions
        .def("cg0_VV", &InclusiveLeptoProduction::cg0_VV)
        .def("cg0_VA", &InclusiveLeptoProduction::cg0_VA)
        .def("cg0_AA", &InclusiveLeptoProduction::cg0_AA)
        .def("dq1_VV", &InclusiveLeptoProduction::dq1_VV)
        .def("dq1_VA", &InclusiveLeptoProduction::dq1_VA)
        .def("dq1_AA", &InclusiveLeptoProduction::dq1_AA)
        .def("cq1_VV", &InclusiveLeptoProduction::cq1_VV)
        .def("cq1_VA", &InclusiveLeptoProduction::cq1_VA)
        .def("cq1_AA", &InclusiveLeptoProduction::cq1_AA)
        .def("cqBarF1_VV", &InclusiveLeptoProduction::cqBarF1_VV)
        .def("cqBarF1_VA", &InclusiveLeptoProduction::cqBarF1_VA)
        .def("cqBarF1_AA", &InclusiveLeptoProduction::cqBarF1_AA)
        .def("cg1_VV", &InclusiveLeptoProduction::cg1_VV)
        .def("cg1_VA", &InclusiveLeptoProduction::cg1_VA)
        .def("cg1_AA", &InclusiveLeptoProduction::cg1_AA)
        .def("cgBarF1_VV", &InclusiveLeptoProduction::cgBarF1_VV)
        .def("cgBarF1_VA", &InclusiveLeptoProduction::cgBarF1_VA)
        .def("cgBarF1_AA", &InclusiveLeptoProduction::cgBarF1_AA)
        .def("cgBarR1_VV", &InclusiveLeptoProduction::cgBarR1_VV)
        .def("cgBarR1_VA", &InclusiveLeptoProduction::cgBarR1_VA)
        .def("cgBarR1_AA", &InclusiveLeptoProduction::cgBarR1_AA)
        .def("cgBar1_VV", &InclusiveLeptoProduction::cgBar1_VV)
        .def("cgBar1_VA", &InclusiveLeptoProduction::cgBar1_VA)
        .def("cgBar1_AA", &InclusiveLeptoProduction::cgBar1_AA)
        // hadronic setter
        .def("setPdf", &InclusiveLeptoProduction::setPdf, "sets pdf")
        .def("setMuR2", &InclusiveLeptoProduction::setMuR2, "sets factors for renormalisation scale")
        .def("setMuF2", &InclusiveLeptoProduction::setMuF2, "sets factors for factorisation scale")
        .def("setMu2", &InclusiveLeptoProduction::setMu2, "sets factors for common scale")
        .def("setLambdaQCD", &InclusiveLeptoProduction::setLambdaQCD, "sets AlphaS_Analytic and its lambda_{QCD,f}")
        .def("setAlphaSByLHAPDF", &InclusiveLeptoProduction::setAlphaSByLHAPDF, "sets Alpha_S by one LHAPDF")
        .def("setXBjorken", &InclusiveLeptoProduction::setXBjorken, "sets Bjorken x")
        .def("setHadronicS", &InclusiveLeptoProduction::setHadronicS, "sets hadronic Sh = (p+q)²")
        .def("flags", &InclusiveLeptoProduction::flags, "manipulate controlling flags", return_value_policy<reference_existing_object>())
        // hadronic structure functions
        .def("F",   &InclusiveLeptoProduction::F)
        .def("dF_dHAQTransverseMomentum",     &InclusiveLeptoProduction::dF_dHAQTransverseMomentum)
        .def("dF_dHAQTransverseMomentumScaling",     &InclusiveLeptoProduction::dF_dHAQTransverseMomentumScaling)
        .def("dF_dHAQRapidity",     &InclusiveLeptoProduction::dF_dHAQRapidity)
        .def("dF_dHAQFeynmanX",     &InclusiveLeptoProduction::dF_dHAQFeynmanX)
        // leptonic setter
        .def("setAlphaEM", &InclusiveLeptoProduction::setAlphaEM, "sets electro-magnetic coupling constant")
        .def("setPolarizeBeams", &InclusiveLeptoProduction::setPolarizeBeams, "use polarized beams?")
        .def("setLeptonicS", &InclusiveLeptoProduction::setLeptonicS, "sets leptonic Sl = (p+l1)²")
        .def("setQ2min", &InclusiveLeptoProduction::setQ2min, "sets constant cut on lower Q2")
        .def("setQ2minByHVQDIS", &InclusiveLeptoProduction::setQ2minByHVQDIS, "cut on lower Q2 as done by HVQDIS")
        // leptonic cross sections
        .def("sigma",   &InclusiveLeptoProduction::sigma, "computes the leptonic cross section")
    ;
    
    class_<FullyDiffLeptoProduction>("FullyDiffLeptoProduction", "application class for fully differential lepto production", init<uint,double,double,double,double,double>())
        // global getter and setter
        .def("setProjection", &FullyDiffLeptoProduction::setProjection, "sets projection")
        .def("setNumberOfLightFlavours", &FullyDiffLeptoProduction::setNumberOfLightFlavours, "sets number of light flavours")
        .def("setQ2", &FullyDiffLeptoProduction::setQ2, "sets transferred momentum Q2 > 0")
        .def("setM2", &FullyDiffLeptoProduction::setM2, "sets mass of heavy quark; unsets partonic center of mass energy")
        .def("setXTilde", &FullyDiffLeptoProduction::setXTilde, "sets xTilde")
        .def("setOmega", &FullyDiffLeptoProduction::setOmega, "sets omega")
        .def("setDeltax", &FullyDiffLeptoProduction::setDeltax, "sets deltax")
        .def("setDeltay", &FullyDiffLeptoProduction::setDeltay, "sets deltay")
        .def("getIntegrationConfig", &FullyDiffLeptoProduction::getIntegrationConfig, "get matching IntegrationConfig", return_value_policy<reference_existing_object>())
        .def("getIntegrationOutput", &FullyDiffLeptoProduction::getIntegrationOutput, "get current/latest IntegrationOutput")
        // partonic setter
        .def("setPartonicEta", &FullyDiffLeptoProduction::setPartonicEta, "sets partonic eta")
        .def("setPartonicS", &FullyDiffLeptoProduction::setPartonicS, "sets partonic s = (k1+q)²")
        // partonic coefficient functions
        .def("cg0_VV", &FullyDiffLeptoProduction::cg0_VV)
        .def("cg0_VA", &FullyDiffLeptoProduction::cg0_VA)
        .def("cg0_AA", &FullyDiffLeptoProduction::cg0_AA)
        .def("dq1_VV", &FullyDiffLeptoProduction::dq1_VV)
        .def("dq1_VA", &FullyDiffLeptoProduction::dq1_VA)
        .def("dq1_AA", &FullyDiffLeptoProduction::dq1_AA)
        .def("cq1_VV", &FullyDiffLeptoProduction::cq1_VV)
        .def("cq1_VA", &FullyDiffLeptoProduction::cq1_VA)
        .def("cq1_AA", &FullyDiffLeptoProduction::cq1_AA)
        .def("cqBarF1_VV", &FullyDiffLeptoProduction::cqBarF1_VV)
        .def("cqBarF1_VA", &FullyDiffLeptoProduction::cqBarF1_VA)
        .def("cqBarF1_AA", &FullyDiffLeptoProduction::cqBarF1_AA)
        .def("cg1_VV", &FullyDiffLeptoProduction::cg1_VV)
        .def("cg1_VA", &FullyDiffLeptoProduction::cg1_VA)
        .def("cg1_AA", &FullyDiffLeptoProduction::cg1_AA)
        .def("cgBarF1_VV", &FullyDiffLeptoProduction::cgBarF1_VV)
        .def("cgBarF1_VA", &FullyDiffLeptoProduction::cgBarF1_VA)
        .def("cgBarF1_AA", &FullyDiffLeptoProduction::cgBarF1_AA)
        .def("cgBarR1_VV", &FullyDiffLeptoProduction::cgBarR1_VV)
        .def("cgBarR1_VA", &FullyDiffLeptoProduction::cgBarR1_VA)
        .def("cgBarR1_AA", &FullyDiffLeptoProduction::cgBarR1_AA)
        .def("cgBar1_VV", &FullyDiffLeptoProduction::cgBar1_VV)
        .def("cgBar1_VA", &FullyDiffLeptoProduction::cgBar1_VA)
        .def("cgBar1_AA", &FullyDiffLeptoProduction::cgBar1_AA)
        // hadronic setter
        .def("setPdf", &FullyDiffLeptoProduction::setPdf, "sets pdf")
        .def("setMuR2", &FullyDiffLeptoProduction::setMuR2, "sets factors for renormalisation scale")
        .def("setMuF2", &FullyDiffLeptoProduction::setMuF2, "sets factors for factorisation scale")
        .def("setMu2", &FullyDiffLeptoProduction::setMu2, "sets factors for common scale")
        .def("setLambdaQCD", &FullyDiffLeptoProduction::setLambdaQCD, "sets AlphaS_Analytic and its lambda_{QCD,f}")
        .def("setAlphaSByLHAPDF", &FullyDiffLeptoProduction::setAlphaSByLHAPDF, "sets Alpha_S by one LHAPDF")
        .def("setXBjorken", &FullyDiffLeptoProduction::setXBjorken, "sets Bjorken x")
        .def("setHadronicS", &FullyDiffLeptoProduction::setHadronicS, "sets hadronic Sh = (p+q)²")
        .def("flags", &FullyDiffLeptoProduction::flags, "manipulate controlling flags", return_value_policy<reference_existing_object>())
        // hadronic structure functions
        .def("F",   &FullyDiffLeptoProduction::F)
        .def("activateHistogram", &FullyDiffLeptoProduction::activateHistogram, FullyDiffLeptoProduction_activateHistogram_overloads())
        // leptonic setter
        .def("setAlphaEM", &FullyDiffLeptoProduction::setAlphaEM, "sets electro-magnetic coupling constant")
        .def("setPolarizeBeams", &FullyDiffLeptoProduction::setPolarizeBeams, "use polarized beams?")
        .def("setLeptonicS", &FullyDiffLeptoProduction::setLeptonicS, "sets leptonic Sl = (p+l1)²")
        .def("setQ2min", &FullyDiffLeptoProduction::setQ2min, "sets constant cut on lower Q2")
        .def("setQ2minByHVQDIS", &FullyDiffLeptoProduction::setQ2minByHVQDIS, "cut on lower Q2 as done by HVQDIS")
        // leptonic cross sections
        .def("sigma",   &FullyDiffLeptoProduction::sigma, "computes the leptonic cross section")
    ;
    
    struct Flags_pickle_suite : boost::python::pickle_suite {
        static boost::python::tuple getinitargs(Flags const& f) {
            return boost::python::make_tuple(f.useGluonicChannel,f.useQuarkChannel,
                f.useLeadingOrder, f.useNextToLeadingOrder,
                f.usePhoton, f.usePhotonZ, f.useZ);
        }
    };
    class_<Flags>("Flags", "controls active channels, orders and bosons", init<bool,bool,bool,bool,bool,bool,bool>())
        .def_readwrite("useGluonicChannel", &Flags::useGluonicChannel, "use gluonic initial state?")
        .def_readwrite("useQuarkChannel", &Flags::useQuarkChannel, "use light quark initial state?")
        .def_readwrite("useLeadingOrder", &Flags::useLeadingOrder, "use leading order calculations?")
        .def_readwrite("useNextToLeadingOrder", &Flags::useNextToLeadingOrder, "use pure(!) next-to-leading order corrections?")
        .def_readwrite("usePhoton", &Flags::usePhoton, "use photon exchange?")
        .def_readwrite("usePhotonZ", &Flags::usePhotonZ, "use interference between photon and Z exchange?")
        .def_readwrite("useZ", &Flags::useZ, "use Z exchange?")
        .def("toYAML", &Common::IntegrationOutput::toYAML, "return YAML representation as string")
        .def_pickle(Flags_pickle_suite())
    ;
    
    struct DynamicScaleFactors_pickle_suite : boost::python::pickle_suite {
        static boost::python::tuple getinitargs(DynamicScaleFactors const& i) {
            return boost::python::make_tuple(i.cM2, i.cQ2, i.cHAQTransverseMomentum, i.cHQPairTransverseMomentum);
        }
    };
    class_<DynamicScaleFactors>("DynamicScaleFactors", "computes dynamic scales", init<double,double,double,double>())
        .def_readwrite("cM2", &DynamicScaleFactors::cM2, "factor to m2")
        .def_readwrite("cQ2", &DynamicScaleFactors::cQ2, "factor to Q2")
        .def_readwrite("cHAQTransverseMomentum", &DynamicScaleFactors::cHAQTransverseMomentum, "factor to p_{Qbar,T)^2")
        .def_readwrite("cHQPairTransverseMomentum", &DynamicScaleFactors::cHQPairTransverseMomentum, "factor to (p_{Q,T}+p_{Qbar,T})^2")
        .def(self * cdbl())
        .def(cdbl() * self)
        .def("toYAML", &DynamicScaleFactors::toYAML, "return YAML representation as string")
        .def_pickle(DynamicScaleFactors_pickle_suite())
    ;
    
    struct IntegrationOutput_pickle_suite : boost::python::pickle_suite {
        static boost::python::tuple getinitargs(Common::IntegrationOutput const& i) {
            return boost::python::make_tuple(i.result, i.error, i.MC_chi2, i.MC_chi2inter);
        }
    };
    class_<Common::IntegrationOutput>("IntegrationOutput", "stores meta data for a single integration", init<double,double,double,uint>())
        .def_readwrite("result", &Common::IntegrationOutput::result, "result")
        .def_readwrite("error", &Common::IntegrationOutput::error, "error")
        .def_readwrite("MC_chi2", &Common::IntegrationOutput::MC_chi2, "chi^2")
        .def_readwrite("MC_chi2inter", &Common::IntegrationOutput::MC_chi2inter, "number of iteration to converge chi^2")
        .def("toYAML", &Common::IntegrationOutput::toYAML, "return YAML representation as string")
        .def_pickle(IntegrationOutput_pickle_suite())
    ;
    
    enum_<FullyDiff::histT>("FullyDiffHistT", "historgram type")
        .value("log10z", FullyDiff::histT::log10z)
        .value("log10xi", FullyDiff::histT::log10xi)
        .value("Theta1", FullyDiff::histT::Theta1)
        
        .value("HQPairInvMass", FullyDiff::histT::HQPairInvMass)
        .value("HQPairDeltaPhi", FullyDiff::histT::HQPairDeltaPhi)
        .value("HQPairTransverseMomentum", FullyDiff::histT::HQPairTransverseMomentum)
        .value("HQPairConeSizeVariable", FullyDiff::histT::HQPairConeSizeVariable)
        
        .value("HAQRapidity", FullyDiff::histT::HAQRapidity)
        .value("HAQTransverseMomentum", FullyDiff::histT::HAQTransverseMomentum)
        .value("HAQTransverseMomentumScaling", FullyDiff::histT::HAQTransverseMomentumScaling)
        .value("HAQFeynmanX", FullyDiff::histT::HAQFeynmanX)
        
        .value("x", FullyDiff::histT::x)
        .value("y", FullyDiff::histT::y)
        .value("Theta2", FullyDiff::histT::Theta2)
    ;
    
    scope intCfg = class_<Common::IntegrationConfig,Common::IntegrationConfig*>("IntegrationConfig", "configurates a single integration")
        .def_readwrite("method", &Common::IntegrationConfig::method, "used integration method")
        .def_readwrite("verbosity", &Common::IntegrationConfig::verbosity, "level of output")
        .def_readwrite("calls", &Common::IntegrationConfig::calls, "calls")
        .def_readwrite("MC_warmupCalls", &Common::IntegrationConfig::MC_warmupCalls, "calls for warmup")
        .def_readwrite("MC_iterations", &Common::IntegrationConfig::MC_iterations, "iterations")
        .def_readwrite("MC_warmupIterations", &Common::IntegrationConfig::MC_warmupIterations, "iterations for warmup")
        .def_readwrite("MC_adaptChi2", &Common::IntegrationConfig::MC_adaptChi2, "iterate until |chi2-1| < 0.5?")
        .def_readwrite("Dvegas_bins", &Common::IntegrationConfig::Dvegas_bins, "number of bins")
        .def_readwrite("GslQag_epsabs", &Common::IntegrationConfig::GslQag_epsabs, "absolute error in gsl_integrate_qag")
        .def_readwrite("GslQag_epsrel", &Common::IntegrationConfig::GslQag_epsrel, "relative error in gsl_integrate_qag")
        .def_readwrite("GslQag_key", &Common::IntegrationConfig::GslQag_key, "number of nodes in gsl_integrate_qag")
        .def("toYAML", &Common::IntegrationConfig::toYAML, "return YAML representation as string")
    ;
    intCfg.attr("GSL_INTEG_GAUSS15") = (int)GSL_INTEG_GAUSS15;
    intCfg.attr("GSL_INTEG_GAUSS21") = (int)GSL_INTEG_GAUSS21;
    intCfg.attr("GSL_INTEG_GAUSS31") = (int)GSL_INTEG_GAUSS31;
    intCfg.attr("GSL_INTEG_GAUSS41") = (int)GSL_INTEG_GAUSS41;
    intCfg.attr("GSL_INTEG_GAUSS51") = (int)GSL_INTEG_GAUSS51;
    intCfg.attr("GSL_INTEG_GAUSS61") = (int)GSL_INTEG_GAUSS61;
}
