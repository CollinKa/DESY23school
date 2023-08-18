// -*- C++ -*-


//just comment everything unless needed
// the cut doesn't match to what you have you should change it!
//this is a template
//identify what is in the paper(histogram,edit the plot file)
//reconstruct the jet of final state, consider only jet larger than specific Pt
//jet cluster momentum togeter. add the four P together that form a jet

//I need to use FastJets class. It is need to PP -collision.

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class CMS_2021_I1876550 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2021_I1876550);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      //const UnstableParticles ufs(Cuts::abseta < 4.9);//should I contiue to use this one? check paper!
      //const ChargedFinalState cfs;
      //declare(cfs,"FS"); 
      //declare(ufs,"UFS");
      declare(UnstableParticles(),"UFS");

      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      //FastJets jetfs(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      //declare(jetfs, "jets");

      // FinalState of direct photons and bare muons and electrons in the event
      //DirectFinalState photons(Cuts::abspid == PID::PHOTON);
      //DirectFinalState bare_leps(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);

      // Dress the bare direct leptons with direct photons within dR < 0.1,
      // and apply some fiducial cuts on the dressed leptons
      //Cut lepton_cuts = Cuts::abseta < 2.5 && Cuts::pT > 20*GeV;
      //DressedLeptons dressed_leps(photons, bare_leps, 0.1, lepton_cuts);
      //declare(dressed_leps, "leptons");

      // Missing momentum
      //declare(MissingMomentum(fs), "MET");

      // Book histograms
/*
      // specify custom binning
      book(_h["XXXX"], "myh1", 20, 0.0, 100.0);
      book(_h["YYYY"], "myh2", logspace(20, 1e-2, 1e3));
      book(_h["ZZZZ"], "myh3", {0.0, 1.0, 2.0, 4.0, 8.0, 16.0});
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      book(_h["AAAA"], 1, 1, 1);
      book(_p["BBBB"], 2, 1, 1);
      book(_c["CCCC"], 3, 1, 1);
*/
      book(_h["D*P"], 1, 1, 1);
      book(_h["D0P"], 1, 1, 2);
      book(_h["DP"], 1, 1, 3);

      book(_h["D*eta"], 2, 1, 1);
      book(_h["D0eta"], 2, 1, 2);
      book(_h["Deta"], 2, 1, 3);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {


      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      //particle code: 411 = D+  421 = "D0" 413 = "D*+"
      //we only interest about the "+" D meson
      for(const Particle& p : ufs.particles()) {
        if (p.abspid() == 411 ||  p.abspid() == 413 || p.abspid() == 421) {
          ConstGenVertexPtr gv = p.genParticle()->production_vertex();
          bool nonPromt = false;
          
          //if have vertex? bug
          /*
          if(gv) {
            for (ConstGenVertexPtr pi: HepMCUtils::particles(gv,Relatives::ANCESTORS)){ //check the particle history
              const PdgId pid2 = pi->pdg_id();
              if (PID::isHadron(pid2) && PID::hasBottom(pid2)) {  //hasBottom:bottom hardon or bottom quark?
                nonPromt = true;
                break;
              } 
            }
          }
          */
          
        
        double abseta = p.abseta();
        double pT = p.perp();
        double GeV = p.E();
        double eta = p.eta();
        if(abseta > 2.1 && !nonPromt) break;
        if(pT < 4. || pT > 100.) break;
        //D+
        if (p.abspid()==411){
          _h["DP"]->fill(pT/GeV);
          _h["Deta"]->fill(eta/GeV);
        }
          
        //D0
        if (p.abspid()==421){
          _h["D0P"]->fill(pT/GeV);
          _h["D0eta"]->fill(eta/GeV);
        }
          
        //D*
        if (p.abspid()==413) {
          _h["D*P"]->fill(pT/GeV);
          _h["D*eta"]->fill(eta/GeV);
        }
          
        }

      }

      
      // Retrieve dressed leptons, sorted by pT
      //Particles leptons = apply<FinalState>(event, "leptons").particles();
      
      //Should I apply a leptop cut? like this?
      //Particles quark = apply<FinalState>(event, "quarks").particles();

      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      //Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 30*GeV);

      //Q: fast Jets?
      //Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 4*GeV);
      

      // Remove all jets within dR < 0.2 of a dressed lepton
      //idiscardIfAnyDeltaRLess(jets, leptons, 0.2);

      // Select jets ghost-associated to B-hadrons with a certain fiducial selection
      //Jets bjets = filter_select(jets, hasBTag(Cuts::pT > 5*GeV && Cuts::abseta < 2.5));
      //Note: I am interest about the D_meson! with 4<pT<100GeV & |eta| < 2.1
      //Jets bjets = filter_select(jets, hasDTag(Cuts::pT < 100*GeV && Cuts::pT > 4*GeV && Cuts::abseta < 2.1));


      // Veto event if there are no b-jets
      //if (bjets.empty()) vetoEvent;
      //Note: check D jets
      //if (Djets.empty()) vetoEvent;
      

      // Apply a missing-momentum cut
      //if (apply<MissingMomentum>(event, "MET").missingPt() < 30*GeV)  vetoEvent;

      // Fill histogram with leading b-jet pT
      //_h["XXXX"]->fill(bjets[0].pT()/GeV);

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h["XXXX"]); // normalize to unity
      normalize(_h["YYYY"], crossSection()/picobarn); // normalize to generated cross-section in pb (no cuts)
      scale(_h["ZZZZ"], crossSection()/picobarn/sumW()); // norm to generated cross-section in pb (after cuts)

    }

    /// @}


    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMS_2021_I1876550);

}
