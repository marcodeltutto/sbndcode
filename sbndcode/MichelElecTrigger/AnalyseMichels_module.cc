////////////////////////////////////////////////////////////////////////
// Class:       AnalyseMichels
// Plugin Type: analyzer (art v3_04_00)
// File:        AnalyseMichels_module.cc
//
// Generated at Wed Feb 26 09:53:33 2020 by Georgia Chisnall using cetskelgen
// from cetlib version v3_09_00.
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
#include "lardataobj/RecoBase/PFParticle.h"

#include <TTree.h>
#include <vector>
#include <string>


namespace sbnd {
  class AnalyseMichels;
}


// The class AnalyseMichels is a derived class that inherits publicly from the base class art::EDAnalyzer
// * The AnalyseMichels class is a subclass of art::EDAnalyzer. Includes required functions we must write and define variables that will persist across events

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
  // Output tree declaration
  TTree *fTree;

  // Variables to fill output tree
  unsigned int fEventID;

  // Variables to access and analyse
  unsigned int fNPFParticles;
  unsigned int fNPrimaries;
  int fNPrimaryDaughters;
  std::string fPFParticleLabel;

};


sbnd::AnalyseMichels::AnalyseMichels(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  // Get the produced name for the PFParticles from the configuration fhicl file analysisConfig.fcl
  fPFParticleLabel = p.get<std::string>("PFParticleLabel");

}

void sbnd::AnalyseMichels::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  // Define event ID variable
  fEventID = e.id().event();

  // Initialise counters for this event
  fNPFParticles = 0;
  fNPrimaries   = 0;
  
  // Access the PFParticles from Pandora
  art::Handle< std::vector<recob::PFParticle> > pfparticleHandle;
  std::vector< art::Ptr<recob::PFParticle>  > pfparticleVect;
  if(e.getByLabel(fPFParticleLabel, pfparticleHandle)) // Make sure the handle is valid
    art::fill_ptr_vector(pfparticleVect, pfparticleHandle);

  if(!pfparticleVect.size()) return;    // If there are no reconstructed particles, skip the event
  fNPFParticles = pfparticleVect.size();

  // Initiate muon ID to be non-physical so we can check we've found it later
  size_t muonID = 99999;

  // if we aren't looking at a primary muon particle, move on to next particle in list
  for(const art::Ptr<recob::PFParticle> &pfp : pfparticleVect){
    if(!(pfp->IsPrimary() && std::abs()pfp->PdgCode() == 13)) continue;
    
    muonID = pfp->Self();
    fNPrimaryDaughters = pfp->NumDaughters();
    fNPrimaries++;
  }

  if(muonID == 99999) return; //If we haven't found a muon, skip event

  // Fill the output TTree with all the relevant variables
  fTree->Fill();

}

void sbnd::AnalyseMichels::beginJob()
{
  // Implementation of optional member function here.
  // TFile service defines the TTree and writes it to the output file
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("tree", "Analyser Output Tree");

  // Add branches to TTree
  fTree->Branch("eventID", &fEventID, "eventID/i");
  fTree->Branch("nPFParticles", &fNPFParticles, "nPFParticles/i");
  fTree->Branch("nPrimaries", &fNPrimaries, "nPrimaries/i");
  fTree->Branch("nPrimaryDaughters", &fNPrimaryDaughters, "nPrimaryDaughters/I");
}

void sbnd::AnalyseMichels::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(sbnd::AnalyseMichels)
