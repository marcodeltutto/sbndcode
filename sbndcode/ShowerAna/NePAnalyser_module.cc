////////////////////////////////////////////////////////////////////////
// Class:       NePAnalyser
// Plugin Type: analyzer (art v3_04_00)
// File:        NePAnalyser_module.cc
//
// Generated at Tue Mar 31 10:43:53 2020 by Ala Zglam using cetskelgen
// from cetlib version v3_09_00.
////////////////////////////////////////////////////////////////////////

//Larsoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

// Additional framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Persistency/Common/FindManyP.h"

// Root includes
#include <TTree.h>
#include <TLorentzVector.h>
#include "TProfile.h"

// C++ includes
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <fstream>

namespace test {
  class NePAnalyser;
}


class test::NePAnalyser : public art::EDAnalyzer {
public:
  explicit NePAnalyser(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NePAnalyser(NePAnalyser const&) = delete;
  NePAnalyser(NePAnalyser&&) = delete;
  NePAnalyser& operator=(NePAnalyser const&) = delete;
  NePAnalyser& operator=(NePAnalyser&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  // Output tree declaration
  TTree *fTree;
  // Variables to fill the output tree with
  unsigned int fEventID;
  //unsigned int fNPFParticles;
  //unsigned int fNPrimaries;
  //int fNPrimaryDaughters;
  std::string fPFParticleLabel;
  std::string fTrackLabel;

};


test::NePAnalyser::NePAnalyser(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  // Get the producer name for the PFParticles from the configuration fhicle file
  //  NePanalysisConfig.fcl
  //
  fPFParticleLabel  = p.get<std::string>("PFParticleLabel");
  fTrackLabel       = p.get<std::string>("TrackLabel");
  //fTrackLengthHist  = tfs->make<TH1D>("trackLengthHist", "Reconstructed track lengths",70,0,340);
  //fCalorimetryLabel =p.get<std::string>("CalorimetryLabel");
}

void test::NePAnalyser::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  // Declaration and Initiate variables
  int NoProton =0;
  int NoDaughter =0;
  // Define our event ID variable
  fEventID = e.id().event();
  //fDaughterTrackLengths->clear();
  //fDaughterTrackdEdx->clear();
  //fDaughterTrackResidualRange->clear();
  //******************************************
  // Initialise our counters for this event
  //fNPFParticles = 0;
  //fNPrimaries   = 0;
  // Access the PFParticles from Pandora
  art::Handle< std::vector<recob::PFParticle> > pfparticleHandle;
  std::vector< art::Ptr<recob::PFParticle> > pfparticleVect;
  // Make sure the handle is valid
  if(e.getByLabel(fPFParticleLabel, pfparticleHandle))
    //Fill the vector with art::Ptr PFParticles
    art::fill_ptr_vector(pfparticleVect, pfparticleHandle);
  if(!pfparticleVect.size()) return; // if there are no
  //fNPFParticles = pfparticleVect.size();

  art::Handle< std::vector<recob::Track> > trackHandle;
  std::vector<art::Ptr<recob::Track> > trackVect;  
  if(e.getByLabel(fTrackLabel, trackHandle)) // Make sure the handle is valid
      art::fill_ptr_vector(trackVect, trackHandle); //fill the vector with art::Ptr Tracks

  //Create a map of trackIDs to MCParticles
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;
  std::map<int, const simb::MCParticle*> particleMap;
  const sim::ParticleList& particles = particleInventory->ParticleList();
  for(auto const& particleIter : particles){
    int trackID = particleIter.first;
    std::cout << "Track ID = " << trackID << std::endl;
    if(particles.IsPrimary(trackID)) continue;
    //    for(sim::ParticleList::const_iterator particleIter = particles.begin(); particleIter != particles.end(); ++particleIter){
      const simb::MCParticle *particle = particleIter.second;
      auto ProcessPart = particle->Process();
      if(ProcessPart =="primary") continue;
      //std::cout<<" Process = "<< ProcessPart << std::endl;
      particleMap[particle->TrackId()] = particle;
    }
  ////////////////////////////////////////////////////////////////////////////////////////////////
  //                          Looping inside the map of MCParticles                            //
  //////////////////////////////////////////////////////////////////////////////////////////////
  for(auto const& mapiter : particleMap){
    //Declaration and Initiate variables
    const simb::MCParticle *particle = mapiter.second;
    int DaughterPdgcode = particle->PdgCode();
    NoDaughter = particle->NumberDaughters();
    //int particleId = particle->TrackId();
    if(DaughterPdgcode == 2212){
      ++NoProton;
    }
  }
    //
    std::cout << " " <<"Event Id"<< " | "<< "No Proton"<< " | "<<"NoDaughter" << std::endl;
    std::cout << "----------------------------------------------------------------"<< std::endl;
    std::cout << " " << fEventID << " | " << NoProton << " | " << NoDaughter << std::endl;
    std::cout << "----------------------------------------------------------------"<< std::endl;
  //----------------------------------------------------------------------------------------------------
  // Fill the output TTree with all relevant variables
  fTree->Fill();
}

void test::NePAnalyser::beginJob()
{
  // Implementation of optional member function here.
  // The TFileService is used to define the TTree and writing it to the output file
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("tree", "Analyser Output Tree");
  // Add branches to our TTree
  fTree->Branch("eventID", &fEventID, "eventID/i");
  //fTree->Branch("nPFParticles",&fNPFParticles,"nPFParticles/i");
  //fTree->Branch("nPrimaries",&fNPrimaries,"nPrimaries/i");
  //fTree->Branch("nPrimaryDaughters",&fNPrimaryDaughters,"nPrimaryDaughters/I");
  //fTree->Branch("daughterTrackLengths",&fDaughterTrackLengths);
  //fTree->Branch("daughterTrackdEdx",&fDaughterTrackdEdx);
  //fTree->Branch("daughterTrackResidualRange", &fDaughterTrackResidualRange);

}

void test::NePAnalyser::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(test::NePAnalyser)
