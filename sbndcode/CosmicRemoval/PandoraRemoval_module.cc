////////////////////////////////////////////////////////////////////////
// Class:       PandoraRemoval
// Module Type: analyzer
// File:        PandoraRemoval_module.cc
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbndcode/CosmicRemoval/CosmicRemovalUtils/CosmicRemovalUtils.h"

// LArSoft includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larcorealg/CoreUtils/NumericUtils.h" // util::absDiff()
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Utilities/Exception.h"
#include "larsim/MCCheater/BackTrackerService.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

#include "Pandora/PdgTable.h"

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <https://root.cern.ch/doc/master/annotated.html>
#include "TVector3.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TFile.h"

// C++ includes
#include <map>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

namespace {
  // Local namespace for local functions
  // Declare here, define later

}

namespace sbnd {

  class PandoraRemoval : public art::EDAnalyzer {
  public:

    // Describes configuration parameters of the module
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
 
      // One Atom for each parameter
      fhicl::Atom<art::InputTag> SimModuleLabel {
        Name("SimModuleLabel"),
        Comment("tag of detector simulation data product")
      };

      fhicl::Atom<art::InputTag> TpcTrackModuleLabel {
        Name("TpcTrackModuleLabel"),
        Comment("tag of TPC track producer data product")
      };

      fhicl::Atom<art::InputTag> ShowerModuleLabel {
        Name("ShowerModuleLabel"),
        Comment("tag of shower producer data product")
      };
      
      fhicl::Atom<art::InputTag> PandoraLabel {
        Name("PandoraLabel"),
        Comment("tag of pandora data product")
      };

      fhicl::Atom<bool> Verbose {
        Name("Verbose"),
        Comment("Print information about what's going on")
      };

    }; // Config
 
    using Parameters = art::EDAnalyzer::Table<Config>;
 
    // Constructor: configures module
    explicit PandoraRemoval(Parameters const& config);
 
    // Called once, at start of the job
    virtual void beginJob() override;
 
    // Called once per event
    virtual void analyze(const art::Event& event) override;

    // Called once, at end of the job
    virtual void endJob() override;

    typedef art::Handle< std::vector<recob::PFParticle> > PFParticleHandle;
    typedef std::map< size_t, art::Ptr<recob::PFParticle> > PFParticleIdMap;
    typedef std::vector< art::Ptr<recob::PFParticle> > PFParticleVector;
    typedef std::vector< art::Ptr<recob::Track> > TrackVector;
    typedef std::vector< art::Ptr<recob::Shower> > ShowerVector;

  private:

    // fcl file parameters
    art::InputTag fSimModuleLabel;      ///< name of detsim producer
    art::InputTag fTpcTrackModuleLabel; ///< name of TPC track producer
    art::InputTag fShowerModuleLabel;
    art::InputTag fPandoraLabel;
    bool          fVerbose;             ///< print information about what's going on

    // histograms

    int nNu = 0;
    int nNuIsNu = 0;
    int nNuIsCr = 0;
    int nMu = 0;
    int nNuIsMu = 0;
    int nCrIsMu = 0;
    int nCr = 0;
    int nCrIsNu = 0;
    int nCrIsCr = 0;

    void GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap);

    void PrintOutScores(const art::Event &evt, const PFParticleHandle &pfParticleHandle) const;

    void GetFinalStatePFParticleVectors(const PFParticleIdMap &pfParticleMap, PFParticleVector &crParticles, PFParticleVector &nuParticles);

    void CollectTracksAndShowers(const PFParticleVector &particles, const PFParticleHandle &pfParticleHandle, const art::Event &evt, TrackVector &tracks, ShowerVector &showers);

  }; // class PandoraRemoval


  // Constructor
  PandoraRemoval::PandoraRemoval(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel       (config().SimModuleLabel())
    , fTpcTrackModuleLabel  (config().TpcTrackModuleLabel())
    , fShowerModuleLabel    (config().ShowerModuleLabel())
    , fPandoraLabel         (config().PandoraLabel())
    , fVerbose              (config().Verbose())
  {

  } // PandoraRemoval()


  void PandoraRemoval::beginJob()
  {
    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;

    // Initial output
    if(fVerbose) std::cout<<"----------------- Pandora Removal Ana Module -------------------"<<std::endl;

  }// PandoraRemoval::beginJob()


  void PandoraRemoval::analyze(const art::Event& event)
  {

    // Fetch basic event info
    if(fVerbose){
      std::cout<<"============================================"<<std::endl
               <<"Run = "<<event.run()<<", SubRun = "<<event.subRun()<<", Event = "<<event.id().event()<<std::endl
               <<"============================================"<<std::endl;
    }

    PFParticleHandle pfParticleHandle;
    event.getByLabel(fPandoraLabel, pfParticleHandle);

    if( !pfParticleHandle.isValid() ){
      std::cout<<"Failed to find the PFParticles."<<std::endl;
      return;
    }

    PFParticleIdMap pfParticleMap;
    this->GetPFParticleIdMap(pfParticleHandle, pfParticleMap);

    if(fVerbose) this->PrintOutScores(event, pfParticleHandle);

    std::vector<art::Ptr<recob::PFParticle>> crParticles;
    std::vector<art::Ptr<recob::PFParticle>> nuParticles;

    this->GetFinalStatePFParticleVectors(pfParticleMap, crParticles, nuParticles);

    std::vector<art::Ptr<recob::Track>> tracks;
    std::vector<art::Ptr<recob::Shower>> showers;
    this->CollectTracksAndShowers(nuParticles, pfParticleHandle, event, tracks, showers);

    std::vector<art::Ptr<recob::Track>> crTracks;
    std::vector<art::Ptr<recob::Shower>> crShowers;
    this->CollectTracksAndShowers(crParticles, pfParticleHandle, event, crTracks, crShowers);

    std::cout << "Consolidated event summary:" << std::endl;
    std::cout << "  - Number of primary cosmic-ray PFParticles   : " << crParticles.size() << std::endl;
    std::cout << "    ... of which are track-like   : " << crTracks.size() << std::endl;
    std::cout << "    ... of which are showers-like : " << crShowers.size() << std::endl;
    std::cout << "  - Number of neutrino final-state PFParticles : " << nuParticles.size() << std::endl;
    std::cout << "    ... of which are track-like   : " << tracks.size() << std::endl;
    std::cout << "    ... of which are showers-like : " << showers.size() << std::endl;


    //----------------------------------------------------------------------------------------------------------
    //                                          TRUTH MATCHING
    //----------------------------------------------------------------------------------------------------------
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    // Retrieve all the truth info in the events
    auto particleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fSimModuleLabel);

    // Loop over true particles in readout window
    std::map<int, simb::MCParticle> particles;
    //std::vector<simb::MCParticle> parts;
    std::vector<int> nuParticleIds;
    std::vector<int> lepParticleIds;

    for (auto const& particle: (*particleHandle)){

      int partId = particle.TrackId();
      particles[partId] = particle;
      //parts.push_back(particle);
      art::Ptr<simb::MCTruth> truth = pi_serv->TrackIdToMCTruth_P(partId);
      if(truth->Origin() == simb::kBeamNeutrino){
        geo::Point_t vtx;
        vtx.SetX(truth->GetNeutrino().Nu().Vx()); vtx.SetY(truth->GetNeutrino().Nu().Vy()); vtx.SetZ(truth->GetNeutrino().Nu().Vz());
        double time = particle.T() * 1e-3; // [us]
        if(fVerbose && particle.Mother()==0) 
          std::cout<<"Nu VTX = "<<vtx<<" ID = "<<partId<<" pdg = "<<particle.PdgCode()<<" time = "<<time<<" length = "<<particle.Trajectory().TotalLength()
                   <<" start = ("<<particle.Vx()<<", "<<particle.Vy()<<", "<<particle.Vz()<<") end = ("<<particle.EndX()<<", "<<particle.EndY()<<", "<<particle.EndZ()<<")\n";
        if(!CosmicRemovalUtils::InFiducial(vtx, 0, 0)) continue;
        if(std::abs(particle.PdgCode())==13 && particle.Mother()==0){ 
          lepParticleIds.push_back(partId);
          nMu++;
        }
        nuParticleIds.push_back(partId);
      }
    }

    if(fVerbose) std::cout<<"Number of true particles = "<<particles.size()<<"\n";

    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTpcTrackModuleLabel);
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTpcTrackModuleLabel);

    for(auto const& track : tracks){
      nNu++;
      std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(track->ID());
      int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
      if (particles.find(trueId) == particles.end()){ 
        if (fVerbose) std::cout<<"No valid true track!\n\n"; 
        continue; 
      }
      if(std::find(lepParticleIds.begin(), lepParticleIds.end(), trueId) != lepParticleIds.end()){
        nNuIsMu++;
      }
      if(std::find(nuParticleIds.begin(), nuParticleIds.end(), trueId) != nuParticleIds.end()){
        if(fVerbose) std::cout<<"From neutrino!\n";
        nNuIsNu++;
      }
      else{
        nNuIsCr++;
      }
    }

    for(auto const& crTrack : crTracks){
      nCr++;
      std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(crTrack->ID());
      int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
      if (particles.find(trueId) == particles.end()){ 
        if (fVerbose) std::cout<<"No valid true track!\n\n"; 
        continue; 
      }
      if(std::find(lepParticleIds.begin(), lepParticleIds.end(), trueId) != lepParticleIds.end()){
        nCrIsMu++;
      }
      if(std::find(nuParticleIds.begin(), nuParticleIds.end(), trueId) != nuParticleIds.end()){
        if(fVerbose) std::cout<<"From neutrino!\n";
        nCrIsNu++;
      }
      else{
        nCrIsCr++;
      }
    }
      

    //----------------------------------------------------------------------------------------------------------
    //                                           CRT COSMIC ID
    //----------------------------------------------------------------------------------------------------------

    //----------------------------------------------------------------------------------------------------------
    //                                          NUMU SELECTION
    //----------------------------------------------------------------------------------------------------------

  } // PandoraRemoval::analyze()


  void PandoraRemoval::endJob(){

    std::cout<<"Total nu tracks             ="<<nNu<<"\n"
             <<"Total mu tracks             ="<<nMu<<"\n"
             <<"Nu is nu                    ="<<nNuIsNu<<"\n"
             <<"Nu is mu                    ="<<nNuIsMu<<"\n"
             <<"Nu is cr                    ="<<nNuIsCr<<"\n"
             <<"-----------------------------------------------------------\n"
             <<"Total cos tracks            ="<<nCr<<"\n"
             <<"Cr is cr                    ="<<nCrIsCr<<"\n"
             <<"Cr is nu                    ="<<nCrIsNu<<"\n"
             <<"Cr is mu                    ="<<nCrIsMu<<"\n";

  } // PandoraRemoval::endJob()

  void PandoraRemoval::GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap){
      for (unsigned int i = 0; i < pfParticleHandle->size(); ++i){
          const art::Ptr<recob::PFParticle> pParticle(pfParticleHandle, i);
          if (!pfParticleMap.insert(PFParticleIdMap::value_type(pParticle->Self(), pParticle)).second){
              std::cout << "  Unable to get PFParticle ID map, the input PFParticle collection has repeat IDs!" <<"\n";
          }
      }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void PandoraRemoval::PrintOutScores(const art::Event &evt, const PFParticleHandle &pfParticleHandle) const{

      // Get the associations between PFParticles and larpandoraobj::PFParticleMetadata
      art::FindManyP< larpandoraobj::PFParticleMetadata > pfPartToMetadataAssoc(pfParticleHandle, evt, fPandoraLabel);

      for (unsigned int i = 0; i < pfParticleHandle->size(); ++i){
          const std::vector< art::Ptr<larpandoraobj::PFParticleMetadata> > &pfParticleMetadataList(pfPartToMetadataAssoc.at(i));
          if (!pfParticleMetadataList.empty()){
              const art::Ptr<recob::PFParticle> pParticle(pfParticleHandle, i);
              for (unsigned int j=0; j<pfParticleMetadataList.size(); ++j){
                  const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(pfParticleMetadataList.at(j));
                  const larpandoraobj::PFParticleMetadata::PropertiesMap &pfParticlePropertiesMap(pfParticleMetadata->GetPropertiesMap());
                  if (!pfParticlePropertiesMap.empty())
                      std::cout << " Found PFParticle " << pParticle->Self() << " with: " << std::endl;
                  for (larpandoraobj::PFParticleMetadata::PropertiesMap::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it)
                      std::cout << "  - " << it->first << " = " << it->second << std::endl;
              }
          }
      }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------
      
  void PandoraRemoval::GetFinalStatePFParticleVectors(const PFParticleIdMap &pfParticleMap, PFParticleVector &crParticles, PFParticleVector &nuParticles){

      for (PFParticleIdMap::const_iterator it = pfParticleMap.begin(); it != pfParticleMap.end(); ++it){
          const art::Ptr<recob::PFParticle> pParticle(it->second);

          // Only look for primary particles
          if (!pParticle->IsPrimary()) continue;

          // Check if this particle is identified as the neutrino
          const int pdg(pParticle->PdgCode());
          const bool isNeutrino(std::abs(pdg) == pandora::NU_E || std::abs(pdg) == pandora::NU_MU || std::abs(pdg) == pandora::NU_TAU);

          // All non-neutrino primary particles are reconstructed under the cosmic hypothesis
          if (!isNeutrino){
              crParticles.push_back(pParticle);
              continue;
          }

          // ATTN. We are filling nuParticles under the assumption that there is only one reconstructed neutrino identified per event.
          //       If this is not the case please handle accordingly
          if (!nuParticles.empty()){
              std::cout << "  This event contains multiple reconstructed neutrinos!" << "\n";
          }

          // Add the daughters of the neutrino PFParticle to the nuPFParticles vector
          for (const size_t daughterId : pParticle->Daughters()){
              if (pfParticleMap.find(daughterId) == pfParticleMap.end())
                  std::cout << "  Invalid PFParticle collection!" <<"\n";

              nuParticles.push_back(pfParticleMap.at(daughterId));
          }
      }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------
      
  void PandoraRemoval::CollectTracksAndShowers(const PFParticleVector &particles, const PFParticleHandle &pfParticleHandle, const art::Event &evt, TrackVector &tracks, ShowerVector &showers)
  {
      // Get the associations between PFParticles and tracks/showers from the event
      art::FindManyP< recob::Track > pfPartToTrackAssoc(pfParticleHandle, evt, fTpcTrackModuleLabel);
      art::FindManyP< recob::Shower > pfPartToShowerAssoc(pfParticleHandle, evt, fShowerModuleLabel);
     
      for (const art::Ptr<recob::PFParticle> &pParticle : particles){
          const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pParticle.key()));
          const std::vector< art::Ptr<recob::Shower> > associatedShowers(pfPartToShowerAssoc.at(pParticle.key()));
          const unsigned int nTracks(associatedTracks.size());
          const unsigned int nShowers(associatedShowers.size());

          // Check if the PFParticle has no associated tracks or showers
          if (nTracks == 0 && nShowers == 0){
              std::cout << "  No tracks or showers were associated to PFParticle " << pParticle->Self() << std::endl;
              continue;
          }

          // Check if there is an associated track
          if (nTracks == 1 && nShowers == 0){
              tracks.push_back(associatedTracks.front());
              continue;
          }

          // Check if there is an associated shower
          if (nTracks == 0 && nShowers == 1){
              showers.push_back(associatedShowers.front());
              continue;
          }

          std::cout << "  There were " << nTracks << " tracks and " << nShowers << " showers associated with PFParticle " << pParticle->Self() << "\n";
      }
  }


  DEFINE_ART_MODULE(PandoraRemoval)
} // namespace sbnd

