////////////////////////////////////////////////////////////////////////
// Class:       PandoraSelection
// Module Type: analyzer
// File:        PandoraSelection_module.cc
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbndcode/CosmicRemoval/CosmicRemovalUtils/CosmicRemovalUtils.h"
#include "sbndcode/CRT/CRTUtils/CRTAnaUtils.h"
#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/CRTTruthRecoAlg.h"

// LArSoft includes
#include "lardataobj/AnalysisBase/T0.h"
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
#include "TCanvas.h"
#include "TPolyLine3D.h"

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

  struct NuTrackCosId{
    recob::Track track;
    int trueId;
    bool inFiducial;
    bool startStops;
    bool endStops;
    bool startInFidStop;
    bool endInFidStop;
    bool diffTpc;
    bool matchesCrtTrack;
    bool isStitched;
    crt::CRTHit closestCrtHit;
    double closestCrtHitDistance;
    bool crossesApa;
  };

  class PandoraSelection : public art::EDAnalyzer {
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

      fhicl::Atom<art::InputTag> CrtHitModuleLabel {
        Name("CrtHitModuleLabel"),
        Comment("tag of CRT hit producer data product")
      };

      fhicl::Atom<art::InputTag> TpcTrackModuleLabel {
        Name("TpcTrackModuleLabel"),
        Comment("tag of TPC track producer data product")
      };

      fhicl::Atom<art::InputTag> CaloModuleLabel {
        Name("CaloModuleLabel"),
        Comment("tag of calorimetry producer data product")
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

      fhicl::Atom<double> Fiducial {
        Name("Fiducial"),
        Comment("Fiducial volume cut")
      };

      fhicl::Atom<double> FiducialTop {
        Name("FiducialTop"),
        Comment("Fiducial volume cut for top of TPC")
      };

      fhicl::Atom<double> FiducialStop {
        Name("FiducialStop"),
        Comment("Fiducial volume cut for stopping tracks")
      };

      fhicl::Atom<double> CpaStitchDistance {
        Name("CpaStitchDistance"),
        Comment("")
      };
      
      fhicl::Atom<double> CpaStitchAngle {
        Name("CpaStitchAngle"),
        Comment("")
      };
      
      fhicl::Atom<double> CpaXDifference {
        Name("CpaXDifference"),
        Comment("")
      };
      
      fhicl::Atom<double> ApaDistance {
        Name("ApaDistance"),
        Comment("")
      };
      
      fhicl::Atom<double> ResRgMin {
        Name("ResRgMin"),
        Comment("")
      };
      
      fhicl::Atom<double> ResRgMax {
        Name("ResRgMax"),
        Comment("")
      };
      
      fhicl::Atom<double> DEdxMax {
        Name("DEdxMax"),
        Comment("")
      };
      
      fhicl::Atom<double> StoppingChi2Limit {
        Name("StoppingChi2Limit"),
        Comment("")
      };
      
      fhicl::Atom<double> MinTrackLength {
        Name("MinTrackLength"),
        Comment("")
      };
      
      fhicl::Atom<double> TrackDirectionFrac {
        Name("TrackDirectionFrac"),
        Comment("")
      };
      
      fhicl::Atom<double> DistanceLimit {
        Name("DistanceLimit"),
        Comment("")
      };
      
      fhicl::Atom<double> MaxAngleDiff {
        Name("MaxAngleDiff"),
        Comment("")
      };
      
      fhicl::Atom<double> MaxDistance {
        Name("MaxDistance"),
        Comment("")
      };

      fhicl::Atom<double> BeamTimeLimit {
        Name("BeamTimeLimit"),
        Comment("")
      };


    }; // Config
 
    using Parameters = art::EDAnalyzer::Table<Config>;
 
    // Constructor: configures module
    explicit PandoraSelection(Parameters const& config);
 
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
    art::InputTag fCrtHitModuleLabel; ///< name of CRT track producer
    art::InputTag fTpcTrackModuleLabel; ///< name of TPC track producer
    art::InputTag fCaloModuleLabel; ///< name of calorimetry producer
    art::InputTag fShowerModuleLabel;
    art::InputTag fPandoraLabel;
    bool          fVerbose;             ///< print information about what's going on
    double        fFiducial;
    double        fFiducialTop;
    double        fFiducialStop;
    double        fCpaStitchDistance;
    double        fCpaStitchAngle;
    double        fCpaXDifference;
    double        fApaDistance;
    double        fResRgMin;
    double        fResRgMax;
    double        fDEdxMax;
    double        fStoppingChi2Limit;
    double        fMinTrackLength;
    double        fTrackDirectionFrac;
    double        fDistanceLimit;
    double        fMaxAngleDiff;
    double        fMaxDistance;
    double        fBeamTimeLimit;

    geo::GeometryCore const* fGeometryService;                 ///< pointer to Geometry provider
    CRTTruthRecoAlg truthAlg;
    // histograms
    int nPfParts = 0;
    int nPfPartsRemain = 0;
    int nNuPfp = 0;
    int nNuMuPfp = 0;
    int nCrPfp = 0;
    int nNuPfpFid = 0;
    int nNuMuPfpFid = 0;
    int nCrPfpFid = 0;
    int nNuPfpRemain = 0;
    int nNuMuPfpRemain = 0;
    int nCrPfpRemain = 0;

    void GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap);

    void PrintOutScores(const art::Event &evt, const PFParticleHandle &pfParticleHandle) const;

    void GetFinalStatePFParticleVectors(const PFParticleIdMap &pfParticleMap, PFParticleVector &crParticles, PFParticleVector &nuParticles);

    void CollectTracksAndShowers(const PFParticleVector &particles, const PFParticleHandle &pfParticleHandle, const art::Event &evt, TrackVector &tracks, ShowerVector &showers);

    void DrawTrueTracks(std::vector<std::vector<NuTrackCosId>> pfps, std::map<int, simb::MCParticle> particles, bool truth, bool tpcTracks);

  }; // class PandoraSelection


  // Constructor
  PandoraSelection::PandoraSelection(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel       (config().SimModuleLabel())
    , fCrtHitModuleLabel    (config().CrtHitModuleLabel())
    , fTpcTrackModuleLabel  (config().TpcTrackModuleLabel())
    , fCaloModuleLabel      (config().CaloModuleLabel())
    , fShowerModuleLabel    (config().ShowerModuleLabel())
    , fPandoraLabel         (config().PandoraLabel())
    , fVerbose              (config().Verbose())
    , fFiducial             (config().Fiducial())
    , fFiducialTop          (config().FiducialTop())
    , fFiducialStop         (config().FiducialStop())
    , fCpaStitchDistance    (config().CpaStitchDistance())
    , fCpaStitchAngle       (config().CpaStitchAngle())
    , fCpaXDifference       (config().CpaXDifference())
    , fApaDistance          (config().ApaDistance())
    , fResRgMin             (config().ResRgMin())
    , fResRgMax             (config().ResRgMax())
    , fDEdxMax              (config().DEdxMax())
    , fStoppingChi2Limit    (config().StoppingChi2Limit())
    , fMinTrackLength       (config().MinTrackLength())
    , fTrackDirectionFrac   (config().TrackDirectionFrac())
    , fDistanceLimit        (config().DistanceLimit())
    , fMaxAngleDiff         (config().MaxAngleDiff())
    , fMaxDistance          (config().MaxDistance())
    , fBeamTimeLimit        (config().BeamTimeLimit())
  {
    fGeometryService    = lar::providerFrom<geo::Geometry>();

  } // PandoraSelection()


  void PandoraSelection::beginJob()
  {
    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;

    // Initial output
    if(fVerbose) std::cout<<"----------------- Pandora Selection Ana Module -------------------"<<std::endl;

  }// PandoraSelection::beginJob()


  void PandoraSelection::analyze(const art::Event& event)
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

    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTpcTrackModuleLabel);
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTpcTrackModuleLabel);
    art::FindManyP<anab::Calorimetry> findManyCalo(tpcTrackHandle, event, fCaloModuleLabel);
    art::FindManyP< recob::Track > pfPartToTrackAssoc(pfParticleHandle, event, fTpcTrackModuleLabel);

    //----------------------------------------------------------------------------------------------------------
    //                                          TRUTH MATCHING
    //----------------------------------------------------------------------------------------------------------
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    // Retrieve all the truth info in the events
    auto particleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fSimModuleLabel);

    // Loop over true particles in readout window
    std::map<int, simb::MCParticle> particles;
    std::vector<int> nuParticleIds;
    std::vector<int> lepParticleIds;

    for (auto const& particle: (*particleHandle)){

      int partId = particle.TrackId();
      particles[partId] = particle;
      art::Ptr<simb::MCTruth> truth = pi_serv->TrackIdToMCTruth_P(partId);
      if(truth->Origin() == simb::kBeamNeutrino){
        geo::Point_t vtx;
        vtx.SetX(truth->GetNeutrino().Nu().Vx()); vtx.SetY(truth->GetNeutrino().Nu().Vy()); vtx.SetZ(truth->GetNeutrino().Nu().Vz());
        if(!CosmicRemovalUtils::InFiducial(vtx, 0, 0)) continue;
        if(std::abs(particle.PdgCode())==13 && particle.Mother()==0){ 
          lepParticleIds.push_back(partId);
        }
        nuParticleIds.push_back(partId);
      }
    }

    if(fVerbose) std::cout<<"Number of true particles = "<<particles.size()<<"\n";

    int nPfPart = -1;
    //Loop over the pfparticle map
    for (PFParticleIdMap::const_iterator it = pfParticleMap.begin(); it != pfParticleMap.end(); ++it){
      nPfPart++;

      const art::Ptr<recob::PFParticle> pParticle(it->second);
      // Only look for primary particles
      if (!pParticle->IsPrimary()) continue;

      // Check if this particle is identified as the neutrino
      const int pdg(pParticle->PdgCode());
      const bool isNeutrino(std::abs(pdg) == pandora::NU_E || std::abs(pdg) == pandora::NU_MU || std::abs(pdg) == pandora::NU_TAU);
      //Find neutrino pfparticle
      if(!isNeutrino) continue;
      if(fVerbose) std::cout<<"\nNeutrino PfPart number = "<<nPfPart<<"\n";

      //Get the tracks associated with thedaughters
      int nDaught = 0;
      int nDaughtTracks = 0;
      bool isTrueNu = false;
      bool isTrueNuMu = false;
      for (const size_t daughterId : pParticle->Daughters()){
        if(fVerbose) std::cout<<"Daughter "<<nDaught<<":\n";
        nDaught++;
        art::Ptr<recob::PFParticle> pParticle = pfParticleMap.at(daughterId);

        const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pParticle.key()));
        if(fVerbose) std::cout<<associatedTracks.size()<<" tracks associated with it\n";
        if(associatedTracks.size() != 1) continue;

        nDaughtTracks++;

        //If cosmic ID then  veto pfparticle
        recob::Track tpcTrack = *associatedTracks.front();
        std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
        int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
        if (particles.find(trueId) == particles.end()){ 
          if(fVerbose) std::cout<<"No valid true track!\n\n"; 
          continue;
        }
        if(std::find(nuParticleIds.begin(), nuParticleIds.end(), trueId) != nuParticleIds.end()) {if(fVerbose) std::cout<<"From neutrino!\n"; isTrueNu = true;}
        else{ if(fVerbose) std::cout<<"From cosmic\n";}
        if(std::find(lepParticleIds.begin(), lepParticleIds.end(), trueId) != lepParticleIds.end()){ if(fVerbose) std::cout<<"Primary muon\n"; isTrueNuMu = true;}
      }

      if(isTrueNu) nNuPfp++;
      else nCrPfp++;
      if(isTrueNuMu) nNuMuPfp++;
    }

    //----------------------------------------------------------------------------------------------------------
    //                                          CRT RECONSTRUCTION
    //----------------------------------------------------------------------------------------------------------
    // Retrieve list of CRT tracks
    art::Handle< std::vector<crt::CRTHit> > crtHitHandle;
    std::vector<art::Ptr<crt::CRTHit> > crtHitsPtr;
    if (event.getByLabel(fCrtHitModuleLabel, crtHitHandle))
      art::fill_ptr_vector(crtHitsPtr, crtHitHandle);

    // Do track reconstruction from the hits
    std::vector<crt::CRTHit> crtHits;
    std::vector<crt::CRTHit> crtHitsNoTop;
    for(auto const& crtHit : (*crtHitHandle)){
      crtHits.push_back(crtHit);
      if(crtHit.tagger != "volTaggerTopHigh_0"){
        crtHitsNoTop.push_back(crtHit);
      }
    }
    std::vector<crt::CRTTrack> crtTracks = CRTAnaUtils::CreateCRTTracks(crtHitsPtr, 0.2, 30., true, 25.);

    std::vector<recob::Track> tpcTracksTPC0;
    std::vector<recob::Track> tpcTracksTPC1;
    std::vector<double> trackT0s;
    // Loop over the tpc tracks
    for(auto const& track : (*tpcTrackHandle)){
      std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(track.ID());
      int tpc = CosmicRemovalUtils::DetectedInTPC(hits);
      double startX = track.Start().X();
      double endX = track.End().X();
      if(tpc == 0 && !(startX>0 || endX>0)) tpcTracksTPC0.push_back(track);
      else if(tpc == 1 && !(startX<0 || endX<0)) tpcTracksTPC1.push_back(track);

      int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
      if (particles.find(trueId) == particles.end()) continue;
      if (track.Length() > 10.) trackT0s.push_back(particles[trueId].T() * 1e-3);
    }

    // Get t0s returned from pandora
    std::map<int, std::vector<double>> t0Map;
    art::FindManyP<anab::T0> findManyT0(pfParticleHandle, event, fPandoraLabel);
    for (unsigned int i = 0; i < pfParticleHandle->size(); ++i){
      const art::Ptr<recob::PFParticle> pParticle(pfParticleHandle, i);
      int pfpKey = pParticle->Self();
      const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pfpKey));
      const std::vector< art::Ptr<anab::T0> > associatedT0s(findManyT0.at(pfpKey));
      for(size_t j = 0; j < associatedTracks.size(); j++){
        int trackId = associatedTracks[j]->ID();
        for(size_t k = 0; k < associatedT0s.size(); k++){
          t0Map[trackId].push_back(associatedT0s[k]->Time()*1e-3);
        }
      }
    }


    if(fVerbose) std::cout<<"\nCosmic removal:\n";
    std::vector<std::vector<NuTrackCosId>> pfpNuTracks;
    nPfPart = -1;
    //Loop over the pfparticle map
    for (PFParticleIdMap::const_iterator it = pfParticleMap.begin(); it != pfParticleMap.end(); ++it){
      nPfPart++;

      const art::Ptr<recob::PFParticle> pParticle(it->second);
      // Only look for primary particles
      if (!pParticle->IsPrimary()) continue;

      // Check if this particle is identified as the neutrino
      const int pdg(pParticle->PdgCode());
      const bool isNeutrino(std::abs(pdg) == pandora::NU_E || std::abs(pdg) == pandora::NU_MU || std::abs(pdg) == pandora::NU_TAU);
      //Find neutrino pfparticle
      if(!isNeutrino) continue;
      if(fVerbose) std::cout<<"\nNeutrino PfPart number = "<<nPfPart<<"\n";

      nPfParts++;

      //Get the tracks associated with thedaughters
      int nDaught = 0;
      bool isTrueNu = false;
      bool isTrueNuMu = false;
      std::vector<NuTrackCosId> nuTracks;
      for (const size_t daughterId : pParticle->Daughters()){

        art::Ptr<recob::PFParticle> pParticle = pfParticleMap.at(daughterId);

        const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pParticle.key()));
        if(associatedTracks.size() != 1) continue;

        //If cosmic ID then  veto pfparticle
        recob::Track tpcTrack = *associatedTracks.front();
        std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
        std::vector<art::Ptr<anab::Calorimetry>> calos = findManyCalo.at(tpcTrack.ID());

        int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
        if(std::find(nuParticleIds.begin(), nuParticleIds.end(), trueId) != nuParticleIds.end()) isTrueNu = true;
        if(std::find(lepParticleIds.begin(), lepParticleIds.end(), trueId) != lepParticleIds.end()) isTrueNuMu = true;

        NuTrackCosId nuTrack;
        nuTrack.track = tpcTrack;
        nuTrack.trueId = trueId;

        if(fVerbose) std::cout<<"Daughter "<<nDaught<<":\n";
        nDaught++;

        //---------------------------------  FIDUCIAL VOLUME CUT --------------------------------------
        // Remove any tracks that enter and exit the fiducial volume
        bool startInFiducial = CosmicRemovalUtils::InFiducial(tpcTrack.Vertex(), fFiducial, fFiducialTop);
        bool endInFiducial = CosmicRemovalUtils::InFiducial(tpcTrack.End(), fFiducial, fFiducialTop);
        if(!startInFiducial && !endInFiducial)  nuTrack.inFiducial = false;
        else nuTrack.inFiducial = true;

        //---------------------------------  STOPPING PARTICLE CUT --------------------------------------
        nuTrack.startInFidStop =  CosmicRemovalUtils::InFiducial(tpcTrack.Vertex(), fFiducialStop, fFiducialStop);
        nuTrack.endInFidStop = CosmicRemovalUtils::InFiducial(tpcTrack.End(), fFiducialStop, fFiducialStop);
        nuTrack.startStops = CosmicRemovalUtils::StoppingEnd(calos, tpcTrack.Vertex(), fResRgMin, fResRgMax, 
                                                          fDEdxMax, fStoppingChi2Limit);
        nuTrack.endStops = CosmicRemovalUtils::StoppingEnd(calos, tpcTrack.End(), fResRgMin, fResRgMax,
                                                        fDEdxMax, fStoppingChi2Limit);

         //---------------------------------  DIFFERENT TPC CUT --------------------------------------
        // Remove any tracks that are detected in one TPC and reconstructed in another
        int tpc = CosmicRemovalUtils::DetectedInTPC(hits);
        double startX = tpcTrack.Start().X();
        double endX = tpcTrack.End().X();
        // Check if track is stitched
        // If it is check the start/end points are in same TPC
        if(tpc == 0 && (startX>0 || endX>0)) nuTrack.diffTpc = true;
        else if(tpc == 1 && (startX<0 || endX<0)) nuTrack.diffTpc = true;
        else nuTrack.diffTpc = false;

        //--------------------------------- CPA STITCHING CUT --------------------------------------
        double stitchTime = -99999;
        bool stitchExit = false;
        // Try to match tracks from CPA crossers
        if(tpc == 0){
          std::pair<double, bool> stitchResults = CosmicRemovalUtils::T0FromCpaStitching(tpcTrack, 
                                                                                         tpcTracksTPC1, 
                                                                                         fCpaStitchDistance, 
                                                                                         fCpaStitchAngle, 
                                                                                         fCpaXDifference, 
                                                                                         fFiducial, 
                                                                                         fFiducialTop);
          stitchTime = stitchResults.first;
          stitchExit = stitchResults.second;
        }
        else if(tpc == 1){
          std::pair<double, bool> stitchResults = CosmicRemovalUtils::T0FromCpaStitching(tpcTrack, 
                                                                                         tpcTracksTPC0, 
                                                                                         fCpaStitchDistance, 
                                                                                         fCpaStitchAngle, 
                                                                                         fCpaXDifference, 
                                                                                         fFiducial, 
                                                                                         fFiducialTop);
          stitchTime = stitchResults.first;
          stitchExit = stitchResults.second;
        }
        if(stitchTime == -99999 && t0Map[tpcTrack.ID()].size()>0) stitchTime = t0Map[tpcTrack.ID()][0];
        // If tracks are stitched, get time and remove any outside of beam window
        if(stitchTime != -99999 && (std::abs(stitchTime) > fBeamTimeLimit || stitchExit)) nuTrack.isStitched = true;
        else nuTrack.isStitched = false;
       
        //--------------------------------- CRTTRACK MATCHING --------------------------------------
        // Try to get T0 from CRTTracks, if there's a match then remove the track
        double trackTime = CRTAnaUtils::T0FromCRTTracks(tpcTrack, crtTracks, tpc, fMaxAngleDiff, fMaxDistance);
        if(trackTime != -99999) nuTrack.matchesCrtTrack = true;
        else nuTrack.matchesCrtTrack = false;
       
        //---------------------------------- CRTHIT MATCHING ---------------------------------------
        // Try to get T0 from CRTHits, remove any matched outside the beam
        std::pair<crt::CRTHit, double> closestHit = CRTAnaUtils::ClosestCRTHit(tpcTrack, crtHitsNoTop, tpc, fTrackDirectionFrac);
        nuTrack.closestCrtHit = closestHit.first;
        nuTrack.closestCrtHitDistance = closestHit.second;

        //----------------- APA CROSS MATCHING FOR THROUGH GOING PARTICLES -------------------------
        // Match APA crossers with times from CRT tracks that cross the APA
        double crossTimeThrough = CosmicRemovalUtils::T0FromApaCross(tpcTrack, trackT0s, tpc, fFiducial, 2.);
        nuTrack.crossesApa = false;
        if(crossTimeThrough != -99999 && std::abs(crossTimeThrough) > fBeamTimeLimit){ 
          // Check that the end that doesn't cross the APA exits the TPC
          if(tpc == 0){
            if(tpcTrack.Vertex().X() < tpcTrack.End().X() && (!endInFiducial || nuTrack.endStops)) nuTrack.crossesApa = true;
            else if(!startInFiducial || nuTrack.startStops) nuTrack.crossesApa = true;
          }
          else if(tpc == 1){
            if(tpcTrack.Vertex().X() > tpcTrack.End().X() && (!endInFiducial || nuTrack.endStops)) nuTrack.crossesApa = true;
            else if(!startInFiducial || nuTrack.startStops) nuTrack.crossesApa = true;
          }
        }

        nuTracks.push_back(nuTrack);
       
      }

      if(nuTracks.size() == 0){
        if(fVerbose) std::cout<<"No tracks associated with pfparticle!\n";
        continue;
      }

      std::sort(nuTracks.begin(), nuTracks.end(), [](auto& left, auto& right){
              return left.track.Length() > right.track.Length();});
      NuTrackCosId nuTrack = nuTracks[0];
      if(!nuTrack.inFiducial){
        if(fVerbose) std::cout<<"Track not in fiducial volume\n";
        if(isTrueNu) nNuPfpFid++;
        else nCrPfpFid++;
        if(isTrueNuMu) nNuMuPfpFid++;
        continue;
      }
      if(nuTrack.diffTpc){
        if(fVerbose) std::cout<<"Different TPC\n";
        continue;
      }
      if(nuTrack.matchesCrtTrack){
        if(fVerbose) std::cout<<"Matches CRT track\n";
        continue;
      }
      if(nuTrack.isStitched){
        if(fVerbose) std::cout<<"Track stitched across CPA\n";
        continue;
      }

      bool singleTrack = false;
      if(nuTracks.size() > 1){
        // More than one track associated
        if(fVerbose) std::cout<<"More than one track associated with pfparticle!\n";
        // Look at the two longest tracks
        bool cutVeto = false;
        recob::Track longestTrack = nuTrack.track;
        size_t secTrack = 99999;
        TVector3 start = longestTrack.Vertex<TVector3>();
        TVector3 end = longestTrack.End<TVector3>();
        TVector3 diff1;
        TVector3 diff2;
        bool t1StartEnd = 0;
        bool t2StartEnd = 0;
        for(size_t i = 1; i < nuTracks.size(); i++){
          recob::Track secondTrack = nuTracks[i].track;
          TVector3 start2 = secondTrack.Vertex<TVector3>();
          TVector3 end2 = secondTrack.End<TVector3>();
          double minDist = 5;
          if((start-start2).Mag() < minDist){
            minDist = (start-start2).Mag();
            secTrack = i;
            diff1 = end - start;
            diff2 = end2 - start2;
            t1StartEnd = 1;
            t2StartEnd = 1;
          }
          if((start-end2).Mag() < minDist){
            minDist = (start-end2).Mag();
            secTrack = i;
            diff1 = end - start;
            diff2 = start2 - end2;
            t1StartEnd = 1;
          }
          if((end-start2).Mag() < minDist){
            minDist = (end-start2).Mag();
            secTrack = i;
            diff1 = start - end;
            diff2 = end2 - start2;
            t2StartEnd = 1;
          }
          if((end-end2).Mag() < minDist){
            secTrack = i;
            diff1 = start - end;
            diff2 = start2 - end2;
          }
          if(secTrack != 99999) continue;
        }

        if(secTrack != 99999){
          double angle = diff1.Angle(diff2);
          // Are the angles compatible with them being split tracks, near 180
          if(angle < 2.5) cutVeto = true;
        }

        if(!cutVeto && secTrack != 99999){
          NuTrackCosId nuTrack2 = nuTracks[secTrack];
          //Do the fiducial volume cut assuming tracks are the same
          TVector3 mergeStart = nuTrack.track.Vertex<TVector3>();
          //bool mergeStartStops = nuTrack.startStops;
          bool mergeStartInFidStop = nuTrack.startInFidStop;
          if(t1StartEnd){ 
            mergeStart = nuTrack.track.End<TVector3>();
            //mergeStartStops = nuTrack.endStops;
            mergeStartInFidStop = nuTrack.endInFidStop;
          }
          TVector3 mergeEnd = nuTrack2.track.Vertex<TVector3>();
          //bool mergeEndStops = nuTrack2.startStops;
          bool mergeEndInFidStop = nuTrack2.startInFidStop;
          if(t2StartEnd){ 
            mergeEnd = nuTrack2.track.End<TVector3>();
            //mergeEndStops = nuTrack2.endStops;
            mergeEndInFidStop = nuTrack2.endInFidStop;
          }
          if(!mergeStartInFidStop && !mergeEndInFidStop){
            if(fVerbose) std::cout<<"Track not in fiducial volume\n";
            if(isTrueNu) nNuPfpFid++;
            else nCrPfpFid++;
            if(isTrueNuMu) nNuMuPfpFid++;
            continue;
          }
          //Do the Diff tpc, crt track and cpa stitch cut for the second track
          //Find the vertex
          //Do the stopping cut for the longest track and assuming tracks are split
          if((nuTrack.startStops && !nuTrack.endInFidStop && nuTrack.track.End().Y() > 0) 
              || (nuTrack.endStops && !nuTrack.startInFidStop && nuTrack.track.Vertex().Y() > 0)){
              //|| (mergeStartStops && !mergeEndInFidStop && mergeEnd.Y() > 0)
              //|| (mergeEndStops && !mergeStartInFidStop && mergeStart.Y() > 0)){
            if(fVerbose) std::cout<<"Stopping particle\n";
            continue;
          }
          //Do crt hit cut for both tracks
          double crtHitTime = std::abs(((double)(int)nuTrack.closestCrtHit.ts1_ns)*1e-4);
          double crtHitTime2 = std::abs(((double)(int)nuTrack2.closestCrtHit.ts1_ns)*1e-4);
          if((nuTrack.closestCrtHitDistance < fDistanceLimit && crtHitTime > fBeamTimeLimit)
              || (nuTrack2.track.Length()>20 && nuTrack2.closestCrtHitDistance < fDistanceLimit && crtHitTime2 > fBeamTimeLimit)){
            if(fVerbose) std::cout<<"Matches CRT hit outside of beam\n";
            continue;
          }
          //Do apa cross cut for both tracks
          if(nuTrack.crossesApa || (nuTrack2.track.Length()>10 && nuTrack2.crossesApa)){
            if(fVerbose) std::cout<<"Crosses APA\n";
            continue;
          }
        }

        //Behaviour if second track is not valid same as if only one track
        if(secTrack == 99999){
          singleTrack = true;
        }
      }

      if(nuTracks.size() == 1 || singleTrack){
        // Only one track associated
        if(fVerbose) std::cout<<"One track associated with pfparticle!\n";
        if((nuTrack.startStops && !nuTrack.endInFidStop && nuTrack.track.End().Y() > 0) 
            || (nuTrack.endStops && !nuTrack.startInFidStop && nuTrack.track.Vertex().Y() > 0)){
          if(fVerbose) std::cout<<"Stopping particle\n";
          continue;
        }
        double crtHitTime = ((double)(int)nuTrack.closestCrtHit.ts1_ns)*1e-4;
        if(nuTrack.closestCrtHitDistance < fDistanceLimit && std::abs(crtHitTime) > fBeamTimeLimit){
          if(fVerbose) std::cout<<"Matches CRT hit outside of beam\n";
          continue;
        }
        if(nuTrack.crossesApa){
          if(fVerbose) std::cout<<"Crosses APA\n";
          continue;
        }
      }
      
      pfpNuTracks.push_back(nuTracks);

      //If particle passes then select
      nPfPartsRemain++;

      if(isTrueNu){ std::cout<<"NU EVENT REMAINING!\n"; nNuPfpRemain++; }
      else{ std::cout<<"COSMIC REMAINING!\n"; nCrPfpRemain++; }
      if(isTrueNuMu){ std::cout<<"NU MU EVENT REMAINING!\n"; nNuMuPfpRemain++; }
    }

    //DrawTrueTracks(pfpNuTracks, particles, true, true);

    
  } // PandoraSelection::analyze()


  void PandoraSelection::endJob(){

    std::cout<<"Original number of slices        = "<<nPfParts<<"\n"
             <<"From neutrinos                   = "<<nNuPfp<<"\n"
             <<"From cosmic rays                 = "<<nCrPfp<<"\n"
             <<"From nu mu                       = "<<nNuMuPfp<<"\n"
             <<"---------------------------------------------\n"
             <<"Neutrinos removed by fid vol     = "<<nNuPfpFid<<"\n"
             <<"Cosmic rays removed by fid vol   = "<<nCrPfpFid<<"\n"
             <<"Nu mu removed by fid vol         = "<<nNuMuPfpFid<<"\n"
             <<"---------------------------------------------\n"
             <<"Number remaining after cosmic ID = "<<nPfPartsRemain<<"\n"
             <<"From neutrinos                   = "<<nNuPfpRemain<<"\n"
             <<"From cosmic rays                 = "<<nCrPfpRemain<<"\n"
             <<"From nu mu                       = "<<nNuMuPfpRemain<<"\n";

  } // PandoraSelection::endJob()

  void PandoraSelection::GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap){
      for (unsigned int i = 0; i < pfParticleHandle->size(); ++i){
          const art::Ptr<recob::PFParticle> pParticle(pfParticleHandle, i);
          if (!pfParticleMap.insert(PFParticleIdMap::value_type(pParticle->Self(), pParticle)).second){
              std::cout << "  Unable to get PFParticle ID map, the input PFParticle collection has repeat IDs!" <<"\n";
          }
      }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void PandoraSelection::PrintOutScores(const art::Event &evt, const PFParticleHandle &pfParticleHandle) const{

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
      
  void PandoraSelection::GetFinalStatePFParticleVectors(const PFParticleIdMap &pfParticleMap, PFParticleVector &crParticles, PFParticleVector &nuParticles){

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
      
  void PandoraSelection::CollectTracksAndShowers(const PFParticleVector &particles, const PFParticleHandle &pfParticleHandle, const art::Event &evt, TrackVector &tracks, ShowerVector &showers)
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


  // Function to draw true and reco tracks
  void PandoraSelection::DrawTrueTracks(std::vector<std::vector<NuTrackCosId>> pfps, std::map<int, simb::MCParticle> particles, bool truth, bool tpcTracks){

    // Create a canvas 
    TCanvas *c1 = new TCanvas("c2","",700,700);
    
    double xmin = -200; //-2.0 * fGeometryService->DetHalfWidth();
    double xmax = 200; //2.0 * fGeometryService->DetHalfWidth();
    double ymin = -fGeometryService->DetHalfHeight();
    double ymax = fGeometryService->DetHalfHeight();
    double zmin = 0.;
    double zmax = fGeometryService->DetLength();
    double rmin[3] = {xmin, ymin, zmin};
    double rmax[3] = {0, ymax, zmax};
    truthAlg.DrawCube(c1, rmin, rmax, 1);
    double rmin1[3] = {0, ymin, zmin};
    double rmax1[3] = {xmax, ymax, zmax};
    truthAlg.DrawCube(c1, rmin1, rmax1, 1);

    // Draw the true particles
    TPolyLine3D *trajectories[100];
    TPolyLine3D *tpctrack[100];

    // Plot the tracks
    int ntrk = 0;
    for(size_t i = 0; i < pfps.size(); i++){
      std::vector<NuTrackCosId> pfp = pfps[i];
      for(size_t k = 0; k < pfp.size(); k++){
        // Get the start and end points
        recob::Track tr = pfp[k].track;
        size_t npts = tr.NumberTrajectoryPoints();
        tpctrack[ntrk] = new TPolyLine3D(npts);
        for(size_t j = 0; j < npts; j++){
          auto& pos = tr.LocationAtPoint(j);
          if(pos.X()==-999 || (pos.X()==0&&pos.Y()==0)) continue;
          tpctrack[ntrk]->SetPoint(j, pos.X(), pos.Y(), pos.Z());
        }
        // Draw a line between them
        tpctrack[ntrk]->SetLineColor(3);
        tpctrack[ntrk]->SetLineWidth(2);
        if(tpcTracks) tpctrack[ntrk]->Draw();
        ntrk++;
        int trueId = pfp[k].trueId;
        if (particles.find(trueId) == particles.end()) continue;
        int nTraj = particles[trueId].NumberTrajectoryPoints();
        trajectories[ntrk] = new TPolyLine3D(nTraj);
        int ipt = 0;
        for(int j = 0; j < nTraj; j++){
          double px = particles[trueId].Vx(j);
          double py = particles[trueId].Vy(j);
          double pz = particles[trueId].Vz(j);
          if(abs(px) < 250 && py < 250 && py > -250 && pz < 550 && pz > -50){
            trajectories[ntrk]->SetPoint(ipt, px, py, pz);
            ipt++;
          }
        }
        trajectories[ntrk]->SetLineColor(4);
        trajectories[ntrk]->SetLineWidth(2);
        if(truth) trajectories[ntrk]->Draw();
      }
    }
    c1->SaveAs("aPandoraPlot.root");

  } // TPCCosmicRemoval::DrawTrueTracks()


  DEFINE_ART_MODULE(PandoraSelection)
} // namespace sbnd

