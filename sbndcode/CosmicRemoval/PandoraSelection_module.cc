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
#include "sbndcode/CRT/CRTUtils/CRTT0MatchAlg.h"
#include "sbndcode/CosmicRemoval/CosmicRemovalAlgs/FiducialVolumeCosmicTagAlg.h"
#include "sbndcode/CosmicRemoval/CosmicRemovalAlgs/StoppingParticleCosmicTagAlg.h"
#include "sbndcode/CosmicRemoval/CosmicRemovalAlgs/GeometryCosmicTagAlg.h"
#include "sbndcode/CosmicRemoval/CosmicRemovalAlgs/CpaCrossCosmicTagAlg.h"
#include "sbndcode/CosmicRemoval/CosmicRemovalAlgs/ApaCrossCosmicTagAlg.h"
#include "sbndcode/CosmicRemoval/CosmicRemovalAlgs/CrtTrackCosmicTagAlg.h"

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
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larcorealg/CoreUtils/NumericUtils.h" // util::absDiff()
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

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
    //crt::CRTHit closestCrtHit;
    //double closestCrtHitDistance;
    double crtHitTime;
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

      fhicl::Atom<art::InputTag> CrtTrackModuleLabel {
        Name("CrtTrackModuleLabel"),
        Comment("tag of CRT track producer data product")
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
      
      fhicl::Atom<double> DistanceLimit {
        Name("DistanceLimit"),
        Comment("")
      };

      fhicl::Atom<double> BeamTimeLimit {
        Name("BeamTimeLimit"),
        Comment("")
      };

      fhicl::Table<FiducialVolumeCosmicTagAlg::Config> FVTagAlg {
        Name("FVTagAlg"),
      };

      fhicl::Table<StoppingParticleCosmicTagAlg::Config> SPTagAlg {
        Name("SPTagAlg"),
      };

      fhicl::Table<CpaCrossCosmicTagAlg::Config> CCTagAlg {
        Name("CCTagAlg"),
      };

      fhicl::Table<ApaCrossCosmicTagAlg::Config> ACTagAlg {
        Name("ACTagAlg"),
      };

      fhicl::Table<CRTT0MatchAlg::Config> CRTT0Alg {
        Name("CRTT0Alg"),
      };

      fhicl::Table<CrtTrackCosmicTagAlg::Config> CTTagAlg {
        Name("CTTagAlg"),
      };

      fhicl::Table<trkf::TrajectoryMCSFitter::Config> fitter {
        Name("fitter"),
      };

    }; // Inputs

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
    art::InputTag fCrtTrackModuleLabel; ///< name of CRT track producer
    art::InputTag fTpcTrackModuleLabel; ///< name of TPC track producer
    art::InputTag fCaloModuleLabel; ///< name of calorimetry producer
    art::InputTag fShowerModuleLabel;
    art::InputTag fPandoraLabel;
    bool          fVerbose;             ///< print information about what's going on
    double        fFiducial;
    double        fFiducialTop;
    double        fFiducialStop;
    double        fDistanceLimit;
    double        fBeamTimeLimit;

    geo::GeometryCore const* fGeometryService;                 ///< pointer to Geometry provider
    detinfo::DetectorProperties const* fDetectorProperties;    ///< pointer to detector properties provider
    detinfo::DetectorClocks const* fDetectorClocks;            ///< pointer to detector clocks provider
    CRTTruthRecoAlg truthAlg;
    FiducialVolumeCosmicTagAlg fvTag;
    StoppingParticleCosmicTagAlg spTag;
    GeometryCosmicTagAlg geoTag;
    CpaCrossCosmicTagAlg ccTag;
    ApaCrossCosmicTagAlg acTag;
    CRTT0MatchAlg crtT0Alg;
    CrtTrackCosmicTagAlg ctTag;
    // Momentum fitters
    trkf::TrajectoryMCSFitter     fMcsFitter; 
    trkf::TrackMomentumCalculator fRangeFitter;
    
    // histograms
    TH1D *hNuMuMomTruth;

    TH1D *hNuMuMom;
    TH1D *hDirtMuMom;
    TH1D *hCrMuMom;
    TH1D *hCrMuMomHigh;
    TH1D *hOtherMuMom;

    TH1D *hNuMuMomRemain;
    TH1D *hOtherMuMomRemain;
    TH1D *hDirtMuMomRemain;
    TH1D *hCrMuMomRemain;
    TH1D *hCrMuMomRemainHigh;

    TH1D *hNuMuMomSelect;
    TH1D *hDirtMuMomSelect;
    TH1D *hOtherMuMomSelect;
    TH1D *hCrMuMomSelect;
    TH1D *hCrMuMomSelectHigh;

    TH1D *hNuMuRecoMom;
    TH1D *hNuRecoMom;
    TH1D *hCrRecoMom;
    TH1D *hDirtRecoMom;
    TH1D *hOtherRecoMom;

    TH1D *hNuMuRecoMomRemain;
    TH1D *hNuRecoMomRemain;
    TH1D *hCrRecoMomRemain;
    TH1D *hDirtRecoMomRemain;
    TH1D *hOtherRecoMomRemain;

    TH1D *hNuMuRecoMomSelect;
    TH1D *hNuRecoMomSelect;
    TH1D *hCrRecoMomSelect;
    TH1D *hDirtRecoMomSelect;
    TH1D *hOtherRecoMomSelect;

    int nPfParts = 0;
    int nPfPartsRemain = 0;
    int nPfPartsSelect = 0;
    int nNuPfp = 0;
    int nNuMuPfp = 0;
    int nCrPfp = 0;
    int nDirtPfp = 0;
    int nOtherPfp = 0;
    int nShowerPfp = 0;

    int nNuPfpFid = 0;
    int nNuMuPfpFid = 0;
    int nCrPfpFid = 0;
    int nDirtPfpFid = 0;
    int nOtherPfpFid = 0;

    int nNuPfpRemain = 0;
    int nNuMuPfpRemain = 0;
    int nCrPfpRemain = 0;
    int nDirtPfpRemain = 0;
    int nOtherPfpRemain = 0;

    int nNuPfpSelect = 0;
    int nNuMuPfpSelect = 0;
    int nCrPfpSelect = 0;
    int nDirtPfpSelect = 0;
    int nOtherPfpSelect = 0;

    void GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap);

    void DrawTrueTracks(std::vector<std::vector<NuTrackCosId>> pfps, std::map<int, simb::MCParticle> particles, bool truth, bool tpcTracks);

  }; // class PandoraSelection


  // Constructor
  PandoraSelection::PandoraSelection(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel       (config().SimModuleLabel())
    , fCrtHitModuleLabel    (config().CrtHitModuleLabel())
    , fCrtTrackModuleLabel  (config().CrtTrackModuleLabel())
    , fTpcTrackModuleLabel  (config().TpcTrackModuleLabel())
    , fCaloModuleLabel      (config().CaloModuleLabel())
    , fShowerModuleLabel    (config().ShowerModuleLabel())
    , fPandoraLabel         (config().PandoraLabel())
    , fVerbose              (config().Verbose())
    , fFiducial             (config().Fiducial())
    , fFiducialTop          (config().FiducialTop())
    , fFiducialStop         (config().FiducialStop())
    , fDistanceLimit        (config().DistanceLimit())
    , fBeamTimeLimit        (config().BeamTimeLimit())
    , fvTag                 (config().FVTagAlg())
    , spTag                 (config().SPTagAlg())
    , ccTag                 (config().CCTagAlg())
    , acTag                 (config().ACTagAlg())
    , crtT0Alg              (config().CRTT0Alg())
    , ctTag                 (config().CTTagAlg())
    , fMcsFitter            (config().fitter)
  {
    fGeometryService    = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 
    fDetectorClocks = lar::providerFrom<detinfo::DetectorClocksService>(); 

  } // PandoraSelection()


  void PandoraSelection::beginJob()
  {
    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;

    hNuMuMomTruth      = tfs->make<TH1D>("hNuMuMomTruth",      "", 20, 0, 2);

    hNuMuMom           = tfs->make<TH1D>("hNuMuMom",           "", 20, 0, 2);
    hDirtMuMom         = tfs->make<TH1D>("hDirtMuMom",         "", 20, 0, 2);
    hOtherMuMom        = tfs->make<TH1D>("hOtherMuMom",        "", 20, 0, 2);
    hCrMuMom           = tfs->make<TH1D>("hCrMuMom",           "", 20, 0, 2);
    hCrMuMomHigh       = tfs->make<TH1D>("hCrMuMomHigh",       "", 20, 0, 10);

    hNuMuMomRemain     = tfs->make<TH1D>("hNuMuMomRemain",     "", 20, 0, 2);
    hOtherMuMomRemain  = tfs->make<TH1D>("hOtherMuMomRemain",  "", 20, 0, 2);
    hDirtMuMomRemain   = tfs->make<TH1D>("hDirtMuMomRemain",   "", 20, 0, 2);
    hCrMuMomRemain     = tfs->make<TH1D>("hCrMuMomRemain",     "", 20, 0, 2);
    hCrMuMomRemainHigh = tfs->make<TH1D>("hCrMuMomRemainHigh", "", 20, 0, 10);

    hNuMuMomSelect     = tfs->make<TH1D>("hNuMuMomSelect",     "", 20, 0, 2);
    hDirtMuMomSelect   = tfs->make<TH1D>("hDirtMuMomSelect",   "", 20, 0, 2);
    hOtherMuMomSelect  = tfs->make<TH1D>("hOtherMuMomSelect",  "", 20, 0, 2);
    hCrMuMomSelect     = tfs->make<TH1D>("hCrMuMomSelect",     "", 20, 0, 2);
    hCrMuMomSelectHigh = tfs->make<TH1D>("hCrMuMomSelectHigh", "", 20, 0, 10);

    hNuMuRecoMom       = tfs->make<TH1D>("hNuMuRecoMom",       "", 20, 0, 2);
    hNuRecoMom         = tfs->make<TH1D>("hNuRecoMom",         "", 20, 0, 2);
    hOtherRecoMom      = tfs->make<TH1D>("hOtherRecoMom",      "", 20, 0, 2);
    hDirtRecoMom       = tfs->make<TH1D>("hDirtRecoMom",       "", 20, 0, 2);
    hCrRecoMom         = tfs->make<TH1D>("hCrRecoMom",         "", 20, 0, 2);

    hNuMuRecoMomRemain  = tfs->make<TH1D>("hNuMuRecoMomRemain",  "", 20, 0, 2);
    hNuRecoMomRemain    = tfs->make<TH1D>("hNuRecoMomRemain",    "", 20, 0, 2);
    hOtherRecoMomRemain = tfs->make<TH1D>("hOtherRecoMomRemain", "", 20, 0, 2);
    hDirtRecoMomRemain  = tfs->make<TH1D>("hDirtRecoMomRemain",  "", 20, 0, 2);
    hCrRecoMomRemain    = tfs->make<TH1D>("hCrRecoMomRemain",    "", 20, 0, 2);

    hNuMuRecoMomSelect  = tfs->make<TH1D>("hNuMuRecoMomSelect",  "", 20, 0, 2);
    hNuRecoMomSelect    = tfs->make<TH1D>("hNuRecoMomSelect",    "", 20, 0, 2);
    hOtherRecoMomSelect = tfs->make<TH1D>("hOtherRecoMomSelect", "", 20, 0, 2);
    hDirtRecoMomSelect  = tfs->make<TH1D>("hDirtRecoMomSelect",  "", 20, 0, 2);
    hCrRecoMomSelect    = tfs->make<TH1D>("hCrRecoMomSelect",    "", 20, 0, 2);

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

    //----------------------------------------------------------------------------------------------------------
    //                                          GETTING PRODUCTS
    //----------------------------------------------------------------------------------------------------------

    // Get truth info and matching
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    // Retrieve all the truth info in the events
    auto particleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fSimModuleLabel);

    // Get PFParticles from pandora
    PFParticleHandle pfParticleHandle;
    event.getByLabel(fPandoraLabel, pfParticleHandle);
    if( !pfParticleHandle.isValid() ){
      if(fVerbose) std::cout<<"Failed to find the PFParticles."<<std::endl;
      return;
    }
    PFParticleIdMap pfParticleMap;
    this->GetPFParticleIdMap(pfParticleHandle, pfParticleMap);
    // Get PFParticle to track associations
    art::FindManyP< recob::Track > pfPartToTrackAssoc(pfParticleHandle, event, fTpcTrackModuleLabel);
    // Get t0s returned from pandora
    art::FindManyP<anab::T0> findManyT0(pfParticleHandle, event, fPandoraLabel);
    
    // Get track to hit and colorimetry associations
    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTpcTrackModuleLabel);
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTpcTrackModuleLabel);
    art::FindManyP<anab::Calorimetry> findManyCalo(tpcTrackHandle, event, fCaloModuleLabel);

    // Retrieve list of CRT hits
    auto crtHitHandle = event.getValidHandle<std::vector<crt::CRTHit>>(fCrtHitModuleLabel);
    std::vector<crt::CRTHit> crtHitsNoTop;
    for(auto const& crtHit : (*crtHitHandle)){
      if(crtHit.tagger != "volTaggerTopHigh_0"){
        crtHitsNoTop.push_back(crtHit);
      }
    }

    // Retrieve list of CRT tracks
    auto crtTrackHandle = event.getValidHandle<std::vector<crt::CRTTrack>>(fCrtTrackModuleLabel);
    std::vector<crt::CRTTrack> crtTracks = (*crtTrackHandle);

    //----------------------------------------------------------------------------------------------------------
    //                                          TRUTH MATCHING
    //----------------------------------------------------------------------------------------------------------

    // Record all true particles and sort by type
    std::map<int, simb::MCParticle> particles;
    std::vector<simb::MCParticle> parts;
    std::vector<int> nuParticleIds;
    std::vector<int> lepParticleIds;
    std::vector<int> dirtParticleIds;
    std::vector<int> crParticleIds;

    // Loop over all true particles
    for (auto const& particle: (*particleHandle)){
      // Store particle
      int partId = particle.TrackId();
      particles[partId] = particle;
      parts.push_back(particle);
      // Get MCTruth
      art::Ptr<simb::MCTruth> truth = pi_serv->TrackIdToMCTruth_P(partId);
      int pdg = std::abs(particle.PdgCode());

      // If origin is a neutrino
      if(truth->Origin() == simb::kBeamNeutrino){
        geo::Point_t vtx;
        vtx.SetX(truth->GetNeutrino().Nu().Vx()); vtx.SetY(truth->GetNeutrino().Nu().Vy()); vtx.SetZ(truth->GetNeutrino().Nu().Vz());
        // If neutrino vertex is not inside the TPC then call it a dirt particle
        if(!CosmicRemovalUtils::InFiducial(vtx, 0, 0)){ 
          dirtParticleIds.push_back(partId);
        }
        // If it's a primary muon
        else if(pdg==13 && particle.Mother()==0){ 
          lepParticleIds.push_back(partId);
          hNuMuMomTruth->Fill(particle.P());
        }
        // Other nu particles
        else{
          nuParticleIds.push_back(partId);
        }
      }

      // If origin is a cosmic ray
      else if(truth->Origin() == simb::kCosmicRay){
        crParticleIds.push_back(partId);
      }
    }

    //----------------------------------------------------------------------------------------------------------
    //                                      FAKE PDS RECONSTRUCTION
    //----------------------------------------------------------------------------------------------------------

    // Create fake flashes in each tpc
    std::pair<std::vector<double>, std::vector<double>> fakeFlashes = CosmicRemovalUtils::FakeTpcFlashes(parts);
    std::vector<double> fakeTpc0Flashes = fakeFlashes.first;
    std::vector<double> fakeTpc1Flashes = fakeFlashes.second;
    bool tpc0BeamFlash = CosmicRemovalUtils::BeamFlash(fakeTpc0Flashes, 4.);
    bool tpc1BeamFlash = CosmicRemovalUtils::BeamFlash(fakeTpc1Flashes, 4.);

    // If there are no flashes in time with the beam then ignore the event
    if(!tpc0BeamFlash && !tpc1BeamFlash) return;
 
    //----------------------------------------------------------------------------------------------------------
    //                                          COSMIC ID - CALCULATING CUTS
    //----------------------------------------------------------------------------------------------------------

    if(fVerbose) std::cout<<"\nCosmic removal:\n";

    std::vector<std::vector<NuTrackCosId>> pfpNuTracks;
    //Loop over the pfparticle map
    for (PFParticleIdMap::const_iterator it = pfParticleMap.begin(); it != pfParticleMap.end(); ++it){

      const art::Ptr<recob::PFParticle> pParticle(it->second);

      // Only look for primary particles
      if (!pParticle->IsPrimary()) continue;
      // Check if this particle is identified as the neutrino
      const int pdg(pParticle->PdgCode());
      const bool isNeutrino(std::abs(pdg) == pandora::NU_E || std::abs(pdg) == pandora::NU_MU || std::abs(pdg) == pandora::NU_TAU);
      //Find neutrino pfparticle
      if(!isNeutrino) continue;

      if(fVerbose) std::cout<<"\nNeutrino PfPart number = "<<pParticle->Self()<<"\n";
      nPfParts++;

      // Truth matching for pfparticle
      bool isTrueNu = false;
      bool isTrueNuMu = false;
      bool isDirt = false;
      bool isCosmic = false;
      bool isOther = false;

      // Create vector of track cosmic ID flags for each pfparticle
      std::vector<NuTrackCosId> nuTracks;
      // Loop over daughters of pfparticle
      for (const size_t daughterId : pParticle->Daughters()){

        // Get tracks associated with daughter
        art::Ptr<recob::PFParticle> pDaughter = pfParticleMap.at(daughterId);
        const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pDaughter.key()));
        if(associatedTracks.size() != 1) continue;

        // Get the first associated track FIXME deal with multiple
        recob::Track tpcTrack = *associatedTracks.front();
        // Get associated hits and calorimetry
        std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
        std::vector<art::Ptr<anab::Calorimetry>> calos = findManyCalo.at(tpcTrack.ID());

        // Get associated t0s from pandora
        const std::vector< art::Ptr<anab::T0> > associatedT0s(findManyT0.at(pDaughter.key()));

        if(fVerbose) std::cout<<"Daughter track "<<daughterId<<":\n";

        // Truth matching
        int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
        // Is track primary numu muon?
        if(std::find(lepParticleIds.begin(), lepParticleIds.end(), trueId) != lepParticleIds.end()){ 
          if(fVerbose) std::cout<<"->NU MU!\n";
          isTrueNuMu = true;
          hNuMuMom->Fill(particles[trueId].P());
        }
        // Is track from another nu particle in TPC?
        else if(std::find(nuParticleIds.begin(), nuParticleIds.end(), trueId) != nuParticleIds.end()){ 
          if(fVerbose) std::cout<<"->NU!\n";
          isTrueNu = true;
        }
        // Is track from nu interaction outside TPC?
        else if(std::find(dirtParticleIds.begin(), dirtParticleIds.end(), trueId) != dirtParticleIds.end()){ 
          if(fVerbose) std::cout<<"->DIRT!\n";
          isDirt = true;
          if(std::abs(particles[trueId].PdgCode())==13){ 
            hDirtMuMom->Fill(particles[trueId].P());
          }
        }
        // Is track from cosmic ray interaction?
        else if(std::find(crParticleIds.begin(), crParticleIds.end(), trueId) != crParticleIds.end()){
          if(fVerbose) std::cout<<"->COSMIC!\n";
          isCosmic = true;
          if(std::abs(particles[trueId].PdgCode())==13){ 
            hCrMuMom->Fill(particles[trueId].P());
            hCrMuMomHigh->Fill(particles[trueId].P());
          }
        }
        // Did truth matching fail or track from mystery origin?
        else{
          if(fVerbose) std::cout<<"->OTHER!\n";
          isOther = true;
          if(std::abs(particles[trueId].PdgCode())==13){ 
            hOtherMuMom->Fill(particles[trueId].P());
          }
        }

        // Record track and true ID
        NuTrackCosId nuTrack;
        nuTrack.track = tpcTrack;
        nuTrack.trueId = trueId;

        //---------------------------------  FIDUCIAL VOLUME CUT --------------------------------------
        // Remove any tracks that enter and exit the fiducial volume
        nuTrack.inFiducial = !fvTag.FiducialVolumeCosmicTag(tpcTrack);

        //---------------------------------  STOPPING PARTICLE CUT --------------------------------------
        nuTrack.startInFidStop =  CosmicRemovalUtils::InFiducial(tpcTrack.Vertex(), fFiducialStop, fFiducialStop);
        nuTrack.endInFidStop = CosmicRemovalUtils::InFiducial(tpcTrack.End(), fFiducialStop, fFiducialStop);

        nuTrack.startStops = spTag.StoppingEnd(tpcTrack.Vertex(), calos);
        nuTrack.endStops = spTag.StoppingEnd(tpcTrack.End(), calos);

         //---------------------------------  DIFFERENT TPC CUT --------------------------------------
        // Remove any tracks that are detected in one TPC and reconstructed in another
        int tpc = CosmicRemovalUtils::DetectedInTPC(hits);

        nuTrack.diffTpc = geoTag.GeometryCosmicTag(tpcTrack, hits, tpc0BeamFlash, tpc1BeamFlash);

        //--------------------------------- CPA STITCHING CUT --------------------------------------
        nuTrack.isStitched = ccTag.CpaCrossCosmicTag(tpcTrack, *tpcTrackHandle, findManyHits);

        if(associatedT0s.size() > 0){
          double pandoraTime = associatedT0s[0]->Time()*1e-3;
          if(pandoraTime < 0 || pandoraTime > fBeamTimeLimit) nuTrack.isStitched = true;
        }
    
        //--------------------------------- CRTTRACK MATCHING --------------------------------------
        // Try to get T0 from CRTTracks, if there's a match then remove the track
        nuTrack.matchesCrtTrack = ctTag.CrtTrackCosmicTag(tpcTrack, crtTracks, tpc);
       
        //---------------------------------- CRTHIT MATCHING ---------------------------------------
        // Try to get T0 from CRTHits, remove any matched outside the beam
        /*std::pair<crt::CRTHit, double> closestHit = crtT0Alg.ClosestCRTHit(tpcTrack, crtHitsNoTop, tpc);
        nuTrack.closestCrtHit = closestHit.first;
        nuTrack.closestCrtHitDistance = closestHit.second;*/
        nuTrack.crtHitTime = crtT0Alg.T0FromCRTHits(tpcTrack, crtHitsNoTop, tpc);

        //----------------- APA CROSS MATCHING FOR THROUGH GOING PARTICLES -------------------------
        // Match APA crossers with times from CRT tracks that cross the APA
        nuTrack.crossesApa = acTag.ApaCrossCosmicTag(tpcTrack, hits, fakeTpc0Flashes, fakeTpc1Flashes);

        nuTracks.push_back(nuTrack);
       
      }
      // Do some accounting FIXME deal with mixed events
      if(isTrueNuMu) nNuMuPfp++;
      else if(isTrueNu) nNuPfp++;
      else if(isDirt) nDirtPfp++;
      else if(isCosmic) nCrPfp++;
      else if(isOther) nOtherPfp++;
      else nShowerPfp++;

      // Don't look at pfps with only showers associated
      if(nuTracks.size() == 0){
        if(fVerbose) std::cout<<"No tracks associated with pfparticle!\n";
        continue;
      }

      // Sort tracks by length (longest chosen as muon)
      std::sort(nuTracks.begin(), nuTracks.end(), [](auto& left, auto& right){
              return left.track.Length() > right.track.Length();});
      NuTrackCosId nuTrack = nuTracks[0];

      // Calculate the momentum of the longest (muon) track
      double recoMuMomentum = 0.;
      const recob::Track& trk = nuTrack.track;
      // Exits if end is within 5cm of TPC boundaries
      bool exits = nuTrack.endInFidStop;
      double length = trk.Length();
      //Momentum calculation
      if(!exits){
        recoMuMomentum = fRangeFitter.GetTrackMomentum(length, 13);
      }
      else{
        recob::MCSFitResult mcsResult = fMcsFitter.fitMcs(trk);
        recoMuMomentum = mcsResult.bestMomentum();
      }

      // Record reco momentum before cosmic ID
      int trueId = nuTrack.trueId;
      //Is it from numu muon
      if(std::find(lepParticleIds.begin(), lepParticleIds.end(), trueId) != lepParticleIds.end()) {
        hNuMuRecoMom->Fill(recoMuMomentum);
      }
      //Is is from other nu particle
      else if(std::find(nuParticleIds.begin(), nuParticleIds.end(), trueId) != nuParticleIds.end()) {
        hNuRecoMom->Fill(recoMuMomentum);
      }
      // From dirt event
      else if(std::find(dirtParticleIds.begin(), dirtParticleIds.end(), trueId) != dirtParticleIds.end()) {
        hDirtRecoMom->Fill(recoMuMomentum);
      }
      // From cosmic ray
      else if(std::find(crParticleIds.begin(), crParticleIds.end(), trueId) != crParticleIds.end()) {
        hCrRecoMom->Fill(recoMuMomentum);
      }
      // Truth matching failed
      else{
        hOtherRecoMom->Fill(recoMuMomentum);
      }

      //----------------------------------------------------------------------------------------------------------
      //                                          COSMIC ID - APPLYING CUTS
      //----------------------------------------------------------------------------------------------------------

      if(!nuTrack.inFiducial){
        if(fVerbose) std::cout<<"Track not in fiducial volume\n";
        // Accounting
        if(isTrueNuMu) nNuMuPfpFid++;
        else if(isTrueNu) nNuPfpFid++;
        else if(isDirt) nDirtPfpFid++;
        else if(isCosmic) nCrPfpFid++;
        else if(isOther) nOtherPfpFid++;
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

      // Apply cut veto if there is a secondary track with angle between tracks consistent with neutrino interaction
      bool singleTrack = false;
      if(nuTracks.size() > 1){
        // More than one track associated
        if(fVerbose) std::cout<<"More than one track associated with pfparticle!\n";

        // Look at the two longest tracks
        size_t secTrack = 99999;
        TVector3 start = trk.Vertex<TVector3>();
        TVector3 end = trk.End<TVector3>();
        double angle = 99999;

        // Loop over the secondary tracks
        for(size_t i = 1; i < nuTracks.size(); i++){
          recob::Track secondTrack = nuTracks[i].track;
          // Only consider secondary tracks longer than 5 cm (try to exclude michel electrons)
          if(secondTrack.Length() < 5) continue;
          // If any of the secondary tracks meet the requirements continue
          if(secTrack != 99999 && angle < 2.5) continue;
          TVector3 start2 = secondTrack.Vertex<TVector3>();
          TVector3 end2 = secondTrack.End<TVector3>();
          double minDist = 5;
          // Do they share the same vertex? (no delta rays)
          if((start-start2).Mag() < minDist){
            minDist = (start-start2).Mag();
            secTrack = i;
            angle = (end - start).Angle(end2 - start2);
          }
        }

        // If angles compatible with them being split tracks (near 180) apply other cuts
        if(secTrack != 99999 && angle > 2.5){
          NuTrackCosId nuTrack2 = nuTracks[secTrack];

          // Do the fiducial volume cut assuming tracks are just one split track
          bool mergeStartInFidStop = nuTrack.endInFidStop;
          bool mergeEndInFidStop = nuTrack2.endInFidStop;
          if(!mergeStartInFidStop && !mergeEndInFidStop){
            if(fVerbose) std::cout<<"Track not in fiducial volume\n";
            // Accounting
            if(isTrueNuMu) nNuMuPfpFid++;
            else if(isTrueNu) nNuPfpFid++;
            else if(isDirt) nDirtPfpFid++;
            else if(isCosmic) nCrPfpFid++;
            else if(isOther) nOtherPfpFid++;
            continue;
          }

          // Do the Diff tpc, crt track and cpa stitch cut for the second track
          // Do the stopping cut for the longest track and assuming tracks are split
          TVector3 mergeStart = nuTrack.track.End<TVector3>();
          TVector3 mergeEnd = nuTrack2.track.End<TVector3>();
          bool mergeStartStops = nuTrack.endStops;
          bool mergeEndStops = nuTrack2.endStops;
          if((nuTrack.startStops && !nuTrack.endInFidStop && nuTrack.track.End().Y() > 0) 
              || (nuTrack.endStops && !nuTrack.startInFidStop && nuTrack.track.Vertex().Y() > 0)
              || (mergeStartStops && !mergeEndInFidStop && mergeEnd.Y() > 0)
              || (mergeEndStops && !mergeStartInFidStop && mergeStart.Y() > 0)){
            if(fVerbose) std::cout<<"Stopping particle\n";
            continue;
          }

          // Do crt hit cut for both tracks
          double crtHitTime = nuTrack.crtHitTime;
          double crtHitTime2 = nuTrack2.crtHitTime;
          if((crtHitTime != -99999 && (crtHitTime < 0 || crtHitTime > fBeamTimeLimit))
             || (crtHitTime2 != -99999 && (crtHitTime2 < 0 || crtHitTime2 > fBeamTimeLimit))){
            if(fVerbose) std::cout<<"Matches CRT hit outside of beam\n";
            continue;
          }

          // Do apa cross cut for both tracks
          if(nuTrack.crossesApa || (/*nuTrack2.track.Length()>10 && */nuTrack2.crossesApa)){
            if(fVerbose) std::cout<<"Crosses APA\n";
            continue;
          }
        }

        // Behaviour if second track is not valid same as if only one track
        if(secTrack == 99999){
          singleTrack = true;
        }
      }

      //If there's only one track then apply secondary cuts
      if(nuTracks.size() == 1 || singleTrack){

        if(fVerbose) std::cout<<"One track associated with pfparticle!\n";

        // Apply stopping particle cut
        if((nuTrack.startStops && !nuTrack.endInFidStop && nuTrack.track.End().Y() > 0) 
            || (nuTrack.endStops && !nuTrack.startInFidStop && nuTrack.track.Vertex().Y() > 0)){
          if(fVerbose) std::cout<<"Stopping particle\n";
          continue;
        }

        // Apply CRT hit matching cut
        double crtHitTime = nuTrack.crtHitTime;
        if(crtHitTime != -99999 && (crtHitTime < 0 || crtHitTime > fBeamTimeLimit)){
          if(fVerbose) std::cout<<"Matches CRT hit outside of beam\n";
          continue;
        }

        // Apply APA crossing cut
        if(nuTrack.crossesApa){
          if(fVerbose) std::cout<<"Crosses APA\n";
          continue;
        }
      }
      
      // For plotting
      pfpNuTracks.push_back(nuTracks);

      // If particle passes then select
      nPfPartsRemain++;
      // Sort by truth matching
      if(isTrueNuMu) { if(fVerbose) std::cout<<"NU MU EVENT REMAINING!\n"; nNuMuPfpRemain++; }
      else if(isTrueNu) { if(fVerbose) std::cout<<"NU EVENT REMAINING!\n"; nNuPfpRemain++; }
      else if(isDirt) { if(fVerbose) std::cout<<"DIRT EVENT REMAINING!\n"; nDirtPfpRemain++; }
      else if(isCosmic) { if(fVerbose) std::cout<<"COSMIC REMAINING!\n"; nCrPfpRemain++; }
      else if(isOther){ if(fVerbose) std::cout<<"OTHER REMAINING!\n"; nOtherPfpRemain++; }

      // Record the true muon momenta of the remaining PFParticles
      for(size_t i = 0; i < nuTracks.size(); i++){
        //Get the true ID
        int trueId2 = nuTracks[i].trueId;
        // Check it is truth matched
        if(particles.find(trueId2) == particles.end()) continue;
        //Check it matches a true muon
        if(std::abs(particles[trueId2].PdgCode()) != 13) continue;

        // Numu muon
        if(std::find(lepParticleIds.begin(), lepParticleIds.end(), trueId2) != lepParticleIds.end()) {
          hNuMuMomRemain->Fill(particles[trueId2].P());
        }
        // Dirt muon
        else if(std::find(dirtParticleIds.begin(), dirtParticleIds.end(), trueId2) != dirtParticleIds.end()) {
          hDirtMuMomRemain->Fill(particles[trueId2].P());
        }
        // Cosmic muon
        else if(std::find(crParticleIds.begin(), crParticleIds.end(), trueId2) != crParticleIds.end()) {
          hCrMuMomRemain->Fill(particles[trueId2].P());
          hCrMuMomRemainHigh->Fill(particles[trueId2].P());
        }
        // Other muon
        else{
          hOtherMuMomRemain->Fill(particles[trueId2].P());
        }
      }

      // Record reco momenta of longest tracks after cosmic ID
      // Numu muon
      if(std::find(lepParticleIds.begin(), lepParticleIds.end(), trueId) != lepParticleIds.end()) {
        hNuMuRecoMomRemain->Fill(recoMuMomentum);
      }
      // Other nu
      else if(std::find(nuParticleIds.begin(), nuParticleIds.end(), trueId) != nuParticleIds.end()) {
        hNuRecoMomRemain->Fill(recoMuMomentum);
      }
      // Dirt
      else if(std::find(dirtParticleIds.begin(), dirtParticleIds.end(), trueId) != dirtParticleIds.end()) {
        hDirtRecoMomRemain->Fill(recoMuMomentum);
      }
      // Cosmic
      else if(std::find(crParticleIds.begin(), crParticleIds.end(), trueId) != crParticleIds.end()) {
        hCrRecoMomRemain->Fill(recoMuMomentum);
      }
      // Truth matching failed
      else{
        hOtherRecoMomRemain->Fill(recoMuMomentum);
      }

      //----------------------------------------------------------------------------------------------------------
      //                                          NUMU SELECTION
      //----------------------------------------------------------------------------------------------------------

      // Get primary vertex from pandora
      TVector3 nuVtx = trk.Vertex<TVector3>();
      // Fiducial volume cut of vertex (from Gray/proposal)
      if(nuVtx.X() < -184.5 || nuVtx.X() > 184.15 || nuVtx.Y() < -185 || nuVtx.Y() > 185 || nuVtx.Z() < 15. || nuVtx.Z() > 420) continue;
      // If longest track is contained, check pid chi2, select mu/pi, select length > 50 cm
      if(!exits && length < 50.) continue;
      // If longest track exits, select length > 100 cm
      if(exits && length < 100.) continue;

      // If it passes these cuts then it is selected as a numuCC interaction
      nPfPartsSelect++;

      // Accounting
      if(isTrueNuMu) { if(fVerbose) std::cout<<"NU MU EVENT SELECTED!\n"; nNuMuPfpSelect++; }
      else if(isTrueNu) { if(fVerbose) std::cout<<"NU EVENT SELECTED!\n"; nNuPfpSelect++; }
      else if(isDirt) { if(fVerbose) std::cout<<"DIRT EVENT SELECTED!\n"; nDirtPfpSelect++; }
      else if(isCosmic) { if(fVerbose) std::cout<<"COSMIC SELECTED!\n"; nCrPfpSelect++; }
      else if(isOther){ if(fVerbose) std::cout<<"OTHER SELECTED!\n"; nOtherPfpSelect++; }

      // Record the true muon momenta of the selected PFParticles
      for(size_t i = 0; i < nuTracks.size(); i++){
        // Get the true ID
        int trueId2 = nuTracks[i].trueId;
        if(particles.find(trueId2) == particles.end()) continue;
        // Check it matches a true muon
        if(std::abs(particles[trueId2].PdgCode()) != 13) continue;
        // Is it from numu muon
        if(std::find(lepParticleIds.begin(), lepParticleIds.end(), trueId2) != lepParticleIds.end()) {
          hNuMuMomSelect->Fill(particles[trueId2].P());
        }
        // Dirt muon
        else if(std::find(dirtParticleIds.begin(), dirtParticleIds.end(), trueId2) != dirtParticleIds.end()) {
          hDirtMuMomSelect->Fill(particles[trueId2].P());
        }
        // Cosmic muon
        else if(std::find(crParticleIds.begin(), crParticleIds.end(), trueId2) != crParticleIds.end()) {
          hCrMuMomSelect->Fill(particles[trueId2].P());
          hCrMuMomSelectHigh->Fill(particles[trueId2].P());
        }
        // Other muon
        else{
          hOtherMuMomSelect->Fill(particles[trueId2].P());
        }
      }

      // Record reco momenta of longest tracks after selection
      // Numu muon
      if(std::find(lepParticleIds.begin(), lepParticleIds.end(), trueId) != lepParticleIds.end()) {
        hNuMuRecoMomSelect->Fill(recoMuMomentum);
      }
      // Other nu particle
      else if(std::find(nuParticleIds.begin(), nuParticleIds.end(), trueId) != nuParticleIds.end()) {
        hNuRecoMomSelect->Fill(recoMuMomentum);
      }
      // Dirt particle
      else if(std::find(dirtParticleIds.begin(), dirtParticleIds.end(), trueId) != dirtParticleIds.end()) {
        hDirtRecoMomSelect->Fill(recoMuMomentum);
      }
      // Cosmic ray particle
      else if(std::find(crParticleIds.begin(), crParticleIds.end(), trueId) != crParticleIds.end()) {
        hCrRecoMomSelect->Fill(recoMuMomentum);
      }
      // Truth matching failed
      else{
        hOtherRecoMomSelect->Fill(recoMuMomentum);
      }

    }

    //DrawTrueTracks(pfpNuTracks, particles, true, true);

    
  } // PandoraSelection::analyze()


  void PandoraSelection::endJob(){

    std::cout<<"---------------------------------------------\n"
             <<"Original number of slices        = "<<nPfParts<<"\n"
             <<"From numu CC                     = "<<nNuMuPfp<<"\n"
             <<"From other nu                    = "<<nNuPfp<<"\n"
             <<"From cosmic rays                 = "<<nCrPfp<<"\n"
             <<"From dirt interactions           = "<<nDirtPfp<<"\n"
             <<"From other                       = "<<nOtherPfp<<"\n"
             <<"Shower only                      = "<<nShowerPfp<<"\n"
             <<"---------------------------------------------\n"
             <<"Numu CC removed by fid vol       = "<<nNuMuPfpFid<<"\n"
             <<"Other nu removed by fid vol      = "<<nNuPfpFid<<"\n"
             <<"Cosmics removed by fid vol       = "<<nCrPfpFid<<"\n"
             <<"Dirt removed by fid vol          = "<<nDirtPfpFid<<"\n"
             <<"Other removed by fid vol         = "<<nOtherPfpFid<<"\n"
             <<"---------------------------------------------\n"
             <<"Number remaining after cosmic ID = "<<nPfPartsRemain<<"\n"
             <<"From numu CC                     = "<<nNuMuPfpRemain<<"\n"
             <<"From other nu                    = "<<nNuPfpRemain<<"\n"
             <<"From cosmic rays                 = "<<nCrPfpRemain<<"\n"
             <<"From dirt                        = "<<nDirtPfpRemain<<"\n"
             <<"From other                       = "<<nOtherPfpRemain<<"\n"
             <<"---------------------------------------------\n"
             <<"Number remaining after selection = "<<nPfPartsSelect<<"\n"
             <<"From numu CC                     = "<<nNuMuPfpSelect<<"\n"
             <<"From neutrinos                   = "<<nNuPfpSelect<<"\n"
             <<"From cosmic rays                 = "<<nCrPfpSelect<<"\n"
             <<"From dirt                        = "<<nDirtPfpSelect<<"\n"
             <<"From other                       = "<<nOtherPfpSelect<<"\n";

    std::ofstream myfile;
    myfile.open("results.txt");
    myfile<<nPfParts<<","<<nNuMuPfp<<","<<nNuPfp<<","<<nCrPfp<<","<<nDirtPfp<<","<<nOtherPfp<<","<<nShowerPfp<<","
          <<nNuMuPfpFid<<","<<nNuPfpFid<<","<<nCrPfpFid<<","<<nDirtPfpFid<<","<<nOtherPfpFid<<","
          <<nPfPartsRemain<<","<<nNuMuPfpRemain<<","<<nNuPfpRemain<<","<<nCrPfpRemain<<","<<nDirtPfpRemain<<","<<nOtherPfpRemain<<","
          <<nPfPartsSelect<<","<<nNuMuPfpSelect<<","<<nNuPfpSelect<<","<<nCrPfpSelect<<","<<nDirtPfpSelect<<","<<nOtherPfpSelect;
    myfile.close();

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

