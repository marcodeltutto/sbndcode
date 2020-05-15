////////////////////////////////////////////////////////////////////////
// Class:       AnalyseMichels
// Plugin Type: analyzer (art v3_05_00)
// File:        AnalyseMichels_module.cc
//
// Generated at Fri May  1 11:14:14 2020 by Georgia Chisnall using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include <TTree.h>
#include <TH1D.h>
#include <vector>
#include <string>


namespace sbnd {
  class AnalyseMichels;
}


class sbnd::AnalyseMichels : public art::EDAnalyzer {
public:
  explicit AnalyseMichels(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  AnalyseMichels(AnalyseMichels const&) = delete;
  AnalyseMichels(AnalyseMichels&&) = delete;
  AnalyseMichels& operator=(AnalyseMichels const&) = delete;
  AnalyseMichels& operator=(AnalyseMichels&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  TTree *fTree;
  TH1D *fMichelEnergyHist;

  unsigned int fEventID;

  //  std::string fMCParticleLabel;
  unsigned int fNMCParticles;
  std::vector<int> *fMCParticlePDG;
  std::vector<float> *fEnergy;
  std::vector<int> *fMother;
};


sbnd::AnalyseMichels::AnalyseMichels(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}, fMCParticlePDG(nullptr), fEnergy(nullptr), fMother(nullptr)  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void sbnd::AnalyseMichels::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  fEventID = e.id().event();
  
  fMCParticlePDG->clear();  
  fEnergy->clear();
  fMother->clear();
  fNMCParticles = 0;
  
  art::Handle< std::vector<simb::MCParticle> > mcparticleHandle;
  std::vector< art::Ptr<simb::MCParticle> > mcparticleVect;       // Create vector of pointers to MCParticle objects

  if(e.getByLabel("largeant", mcparticleHandle))     // Make sure artHandle is from module largeant
    art::fill_ptr_vector(mcparticleVect, mcparticleHandle);  // Fill the vector with art::Ptr MCParticle objects
  
  for (auto const& mcp : mcparticleVect){     // Not sure why we have to use a reference to a pointer
    fMCParticlePDG->push_back(mcp->PdgCode());
    fEnergy->push_back(mcp->E());
    fMother->push_back(mcp->Mother());
    fNMCParticles++;

    if(mcp->PdgCode() == 11){      // && std::abs(mcp->Mother()) == 13){    // if particle is michel electron, add energy to hist
      float michelenergy = mcp->E();
      fMichelEnergyHist->Fill(michelenergy);
    }
  }

  fTree->Fill();
}

void sbnd::AnalyseMichels::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("ptree", "Analyser Output Tree");
  fTree->Branch("eventID", &fEventID, "eventID/i");
  fTree->Branch("nMCParticles", &fNMCParticles, "nMCParticles/i");
  fTree->Branch("PDG", &fMCParticlePDG);
  fTree->Branch("Energy", &fEnergy);
  fTree->Branch("Mother", &fMother);

  fMichelEnergyHist = tfs->make<TH1D>("michelEnergyHist", "Energy of michels;E_{michels} (GeV);Events", 600, 0, 0.06);
}

void sbnd::AnalyseMichels::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(sbnd::AnalyseMichels)
