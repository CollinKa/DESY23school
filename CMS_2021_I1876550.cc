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


      //FastJets jets(veto, FastJets::ANTIKT, 0.6); //I probably need to change the number 0.6  see https://arxiv.org/pdf/1112.4432.pdf
      //declare(jets, "jets");
      declare(UnstableParticles(),"UFS");


      // Book histograms

      book(_h["D*P"], 1, 1, 1);
      book(_h["D0P"], 1, 1, 2);
      book(_h["DP"], 1, 1, 3);

      book(_h["D*eta"], 2, 1, 1);
      book(_h["D0eta"], 2, 1, 2);
      book(_h["Deta"], 2, 1, 3);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // get the jets
      /*
      Jets jets;
      for (const Jet& jet : apply<FastJets>(event, "jets").jetsByPt(4.0*GeV)) {
        if ( jet.abseta() < 2.1 ) jets.push_back(jet);
      }
      */

      //get the excited prompt D meson
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      //particle code: 411 = D+  421 = "D0" 413 = "D*+"
      //we only interest about the "+" D meson
      for(const Particle& p : ufs.particles()) {
        if (p.abspid() == 411 ||  p.abspid() == 413 || p.abspid() == 421) {
          ConstGenVertexPtr gv = p.genParticle()->production_vertex(); //gv represent the where the particle being created
          
          if (gv) {  // Check if gv is not a null pointer
            double x = gv->position().x();
            double y = gv->position().y();
            double z = gv->position().z();
            std::cout << x << " " << y << " " << z << std::endl;
          }
          
          //bool nonPromt = false;
          bool nonPromt = true;
          
          if(gv) {
            //for (ConstGenVertexPtr pi: HepMCUtils::particles(gv,Relatives::ANCESTORS)){ //check the particle history  issue:GenVertex has no member name pdg_id
            //for (const GenParticle* pi: HepMCUtils::particles(gv,Relatives::ANCESTORS)){   //‘GenParticle’ does not name a type; did you mean ‘Particle’?
            //for (const GenParticle& pi : HepMCUtils::particles(gv,Relatives::ANCESTORS)){
              //const PdgId pid2 = pi.pdg_id();   //I really don't know hot to use it
            for(const Particle& p : HepMCUtils::particles(gv,Relatives::ANCESTORS)){    //get the ancestor of particle at specific vertex
              //const PdgId pid2 = p.pdgId();
              const PdgId pid2 = p.pid();
              if (PID::isHadron(pid2) && PID::hasBottom(pid2)) {  //B Hadron -> non prompt
                nonPromt = true;
                break; // change another particle
              }
              
              //adding new cut
              if (PID::isHadron(pid2) && PID::hasCharm(pid2)) {  //B Hadron -> non prompt
                nonPromt = false;

        
                double abseta = p.abseta();
                double pT = p.perp();
                //double GeV = p.E();
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
          }         
        }

      }

      


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
