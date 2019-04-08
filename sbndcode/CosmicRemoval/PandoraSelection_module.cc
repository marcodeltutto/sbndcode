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
    crt::CRTHit closestCrtHit;
    double closestCrtHitDistance;
    bool crossesApa;
  };

  class PandoraSelection : public art::EDAnalyzer {
  public:

    // Describes configuration parameters of the module
    struct Inputs {
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

    }; // Inputs

    // Describes configuration parameters of the module
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
 
      // One Atom for each parameter
      fhicl::Table<PandoraSelection::Inputs> inputs {
        Name("inputs"),
      };
      
      fhicl::Table<trkf::TrajectoryMCSFitter::Config> fitter {
        Name("fitter"),
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
    detinfo::DetectorProperties const* fDetectorProperties;    ///< pointer to detector properties provider
    detinfo::DetectorClocks const* fDetectorClocks;            ///< pointer to detector clocks provider
    CRTTruthRecoAlg truthAlg;
    // Momentum fitters
    trkf::TrajectoryMCSFitter     fMcsFitter; 
    trkf::TrackMomentumCalculator fRangeFitter;
    
    // histograms
    TH1D *hNuMuMomTruth;
    TH1D *hNuMuMom;
    TH1D *hCrMuMom;
    TH1D *hCrMuMomHigh;
    TH1D *hNuMuMomRemain;
    TH1D *hCrMuMomRemain;
    TH1D *hCrMuMomRemainHigh;
    TH1D *hNuMuMomSelect;
    TH1D *hCrMuMomSelect;
    TH1D *hCrMuMomSelectHigh;
    TH1D *hNuMuRecoMom;
    TH1D *hNuRecoMom;
    TH1D *hCrRecoMom;

    int nPfParts = 0;
    int nPfPartsRemain = 0;
    int nPfPartsSelect = 0;
    int nNuPfp = 0;
    int nNuMuPfp = 0;
    int nCrPfp = 0;
    int nNuPfpFid = 0;
    int nNuMuPfpFid = 0;
    int nCrPfpFid = 0;
    int nNuPfpRemain = 0;
    int nNuMuPfpRemain = 0;
    int nCrPfpRemain = 0;
    int nNuPfpSelect = 0;
    int nNuMuPfpSelect = 0;
    int nCrPfpSelect = 0;

    void GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap);

    void PrintOutScores(const art::Event &evt, const PFParticleHandle &pfParticleHandle) const;

    void GetFinalStatePFParticleVectors(const PFParticleIdMap &pfParticleMap, PFParticleVector &crParticles, PFParticleVector &nuParticles);

    void CollectTracksAndShowers(const PFParticleVector &particles, const PFParticleHandle &pfParticleHandle, const art::Event &evt, TrackVector &tracks, ShowerVector &showers);

    void DrawTrueTracks(std::vector<std::vector<NuTrackCosId>> pfps, std::map<int, simb::MCParticle> particles, bool truth, bool tpcTracks);

  }; // class PandoraSelection


  // Constructor
  PandoraSelection::PandoraSelection(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel       (config().inputs().SimModuleLabel())
    , fCrtHitModuleLabel    (config().inputs().CrtHitModuleLabel())
    , fTpcTrackModuleLabel  (config().inputs().TpcTrackModuleLabel())
    , fCaloModuleLabel      (config().inputs().CaloModuleLabel())
    , fShowerModuleLabel    (config().inputs().ShowerModuleLabel())
    , fPandoraLabel         (config().inputs().PandoraLabel())
    , fVerbose              (config().inputs().Verbose())
    , fFiducial             (config().inputs().Fiducial())
    , fFiducialTop          (config().inputs().FiducialTop())
    , fFiducialStop         (config().inputs().FiducialStop())
    , fCpaStitchDistance    (config().inputs().CpaStitchDistance())
    , fCpaStitchAngle       (config().inputs().CpaStitchAngle())
    , fCpaXDifference       (config().inputs().CpaXDifference())
    , fApaDistance          (config().inputs().ApaDistance())
    , fResRgMin             (config().inputs().ResRgMin())
    , fResRgMax             (config().inputs().ResRgMax())
    , fDEdxMax              (config().inputs().DEdxMax())
    , fStoppingChi2Limit    (config().inputs().StoppingChi2Limit())
    , fMinTrackLength       (config().inputs().MinTrackLength())
    , fTrackDirectionFrac   (config().inputs().TrackDirectionFrac())
    , fDistanceLimit        (config().inputs().DistanceLimit())
    , fMaxAngleDiff         (config().inputs().MaxAngleDiff())
    , fMaxDistance          (config().inputs().MaxDistance())
    , fBeamTimeLimit        (config().inputs().BeamTimeLimit())
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
    hCrMuMom           = tfs->make<TH1D>("hCrMuMom",           "", 20, 0, 2);
    hCrMuMomHigh       = tfs->make<TH1D>("hCrMuMomHigh",       "", 20, 0, 10);
    hNuMuMomRemain     = tfs->make<TH1D>("hNuMuMomRemain",     "", 20, 0, 2);
    hCrMuMomRemain     = tfs->make<TH1D>("hCrMuMomRemain",     "", 20, 0, 2);
    hCrMuMomRemainHigh = tfs->make<TH1D>("hCrMuMomRemainHigh", "", 20, 0, 10);
    hNuMuMomSelect     = tfs->make<TH1D>("hNuMuMomSelect",     "", 20, 0, 2);
    hCrMuMomSelect     = tfs->make<TH1D>("hCrMuMomSelect",     "", 20, 0, 2);
    hCrMuMomSelectHigh = tfs->make<TH1D>("hCrMuMomSelectHigh", "", 20, 0, 10);
    hNuMuRecoMom       = tfs->make<TH1D>("hNuMuRecoMom",       "", 20, 0, 2);
    hNuRecoMom         = tfs->make<TH1D>("hNuRecoMom",         "", 20, 0, 2);
    hCrRecoMom         = tfs->make<TH1D>("hCrRecoMom",         "", 20, 0, 2);

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
      if(fVerbose) std::cout<<"Failed to find the PFParticles."<<std::endl;
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

    std::vector<double> fakeTpc0Flashes;
    std::vector<double> fakeTpc1Flashes;
    double readoutWindowMuS  = fDetectorClocks->TPCTick2Time((double)fDetectorProperties->ReadOutWindowSize()); // [us]
    double driftTimeMuS = (2.*fGeometryService->DetHalfWidth()+3.)/fDetectorProperties->DriftVelocity(); // [us]

    for (auto const& particle: (*particleHandle)){

      int partId = particle.TrackId();
      particles[partId] = particle;
      art::Ptr<simb::MCTruth> truth = pi_serv->TrackIdToMCTruth_P(partId);
      int pdg = std::abs(particle.PdgCode());
      double time = particle.T() * 1e-3;
      if(truth->Origin() == simb::kBeamNeutrino){
        geo::Point_t vtx;
        vtx.SetX(truth->GetNeutrino().Nu().Vx()); vtx.SetY(truth->GetNeutrino().Nu().Vy()); vtx.SetZ(truth->GetNeutrino().Nu().Vz());
        if(!CosmicRemovalUtils::InFiducial(vtx, 0, 0)) continue;
        if(pdg==13 && particle.Mother()==0){ 
          lepParticleIds.push_back(partId);
          hNuMuMomTruth->Fill(particle.P());
        }
        nuParticleIds.push_back(partId);
      }
      //Check if time is in reconstructible window
      if(time < -driftTimeMuS || time > readoutWindowMuS) continue; 
      //Check if particle is visible, electron, muon, proton, pion, kaon, photon
      if(!(pdg==13||pdg==11||pdg==22||pdg==2212||pdg==211||pdg==321||pdg==111)) continue;
      //Loop over the trajectory
      int npts = particle.NumberTrajectoryPoints();
      double TPC0Energy = 0;
      double TPC1Energy = 0;
      for(int i = 1; i < npts; i++){
        geo::Point_t pt;
        pt.SetX(particle.Vx(i)); pt.SetY(particle.Vy(i)); pt.SetZ(particle.Vz(i));
        if(!CosmicRemovalUtils::InFiducial(pt, 0, 0)) continue;
        if(pt.X() < 0) TPC0Energy += particle.E(i-1) - particle.E(i);
        else TPC1Energy += particle.E(i-1) - particle.E(i);
      }
      if(TPC0Energy > 0.05) fakeTpc0Flashes.push_back(time);
      else if(TPC1Energy > 0.05) fakeTpc1Flashes.push_back(time);
    }

    bool tpc0BeamFlash = false;
    std::sort(fakeTpc0Flashes.begin(), fakeTpc0Flashes.end());
    double previousTime = -99999;
    for(size_t i = 0; i < fakeTpc0Flashes.size(); i++){
      double time = fakeTpc0Flashes[i];
      if(time > 0 && time < fBeamTimeLimit) tpc0BeamFlash = true;
      if(std::abs(time-previousTime) < 0.1){
        fakeTpc0Flashes.erase(fakeTpc0Flashes.begin()+i);
      }
      else previousTime = time;
    }

    bool tpc1BeamFlash = false;
    std::sort(fakeTpc1Flashes.begin(), fakeTpc1Flashes.end());
    previousTime = -99999;
    for(size_t i = 0; i < fakeTpc1Flashes.size(); i++){
      double time = fakeTpc1Flashes[i];
      if(time > 0 && time < fBeamTimeLimit) tpc1BeamFlash = true;
      if(std::abs(time-previousTime) < 0.1){
        fakeTpc1Flashes.erase(fakeTpc1Flashes.begin()+i);
      }
      else previousTime = time;
    }

    if(!tpc0BeamFlash && !tpc1BeamFlash) return;

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
        else{ 
          if(fVerbose) std::cout<<"From cosmic\n";
          if(std::abs(particles[trueId].PdgCode())==13){ 
            hCrMuMom->Fill(particles[trueId].P());
            hCrMuMomHigh->Fill(particles[trueId].P());
          }
        }
        if(std::find(lepParticleIds.begin(), lepParticleIds.end(), trueId) != lepParticleIds.end()){ 
          if(fVerbose) std::cout<<"Primary muon\n"; 
          isTrueNuMu = true;
          hNuMuMom->Fill(particles[trueId].P());
        }
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
      //FIXME should be true length in TPC not reco length
      std::pair<TVector3, TVector3> se = truthAlg.TpcCrossPoints(particles[trueId]);
      if ((se.second-se.first).Mag() > 10.) trackT0s.push_back(particles[trueId].T() * 1e-3);
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

        if(tpc0BeamFlash && !tpc1BeamFlash && tpc == 1) nuTrack.diffTpc = true;
        if(tpc1BeamFlash && !tpc0BeamFlash && tpc == 0) nuTrack.diffTpc = true;

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
        if(stitchTime != -99999 && (stitchTime < 0 || stitchTime > fBeamTimeLimit || stitchExit)) nuTrack.isStitched = true;
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
        nuTrack.crossesApa = false;
        if(tpc == 0){
          double crossTimeThrough = CosmicRemovalUtils::T0FromApaCross(tpcTrack, fakeTpc0Flashes, tpc, fFiducial, 2.);
          if(crossTimeThrough != -99999 && (crossTimeThrough<0 || crossTimeThrough>fBeamTimeLimit)) nuTrack.crossesApa = true;
        }
        if(tpc == 1){
          double crossTimeThrough = CosmicRemovalUtils::T0FromApaCross(tpcTrack, fakeTpc1Flashes, tpc, fFiducial, 2.);
          if(crossTimeThrough != -99999 && (crossTimeThrough<0 || crossTimeThrough>fBeamTimeLimit)) nuTrack.crossesApa = true;
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
        double angle = 99999;
        //bool t1StartEnd = 0;
        //bool t2StartEnd = 0;
        for(size_t i = 1; i < nuTracks.size(); i++){
          recob::Track secondTrack = nuTracks[i].track;
          if(secondTrack.Length() < 5) continue;
          if(secTrack != 99999 && angle < 2.5) continue;
          TVector3 start2 = secondTrack.Vertex<TVector3>();
          TVector3 end2 = secondTrack.End<TVector3>();
          double minDist = 5;
          if((start-start2).Mag() < minDist){
            minDist = (start-start2).Mag();
            secTrack = i;
            angle = (end - start).Angle(end2 - start2);
            //t1StartEnd = 1;
            //t2StartEnd = 1;
          }
        }

        if(secTrack != 99999){
          // Are the angles compatible with them being split tracks, near 180
          if(angle < 2.5) cutVeto = true;
        }

        if(!cutVeto && secTrack != 99999){
          NuTrackCosId nuTrack2 = nuTracks[secTrack];
          //Do the fiducial volume cut assuming tracks are the same
          TVector3 mergeStart = nuTrack.track.End<TVector3>();
          bool mergeStartStops = nuTrack.endStops;
          bool mergeStartInFidStop = nuTrack.endInFidStop;
          TVector3 mergeEnd = nuTrack2.track.End<TVector3>();
          bool mergeEndStops = nuTrack2.endStops;
          bool mergeEndInFidStop = nuTrack2.endInFidStop;
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
              || (nuTrack.endStops && !nuTrack.startInFidStop && nuTrack.track.Vertex().Y() > 0)
              || (mergeStartStops && !mergeEndInFidStop && mergeEnd.Y() > 0)
              || (mergeEndStops && !mergeStartInFidStop && mergeStart.Y() > 0)){
            if(fVerbose) std::cout<<"Stopping particle\n";
            continue;
          }
          //Do crt hit cut for both tracks
          double crtHitTime = ((double)(int)nuTrack.closestCrtHit.ts1_ns)*1e-4;//FIXME
          double crtHitTime2 = ((double)(int)nuTrack2.closestCrtHit.ts1_ns)*1e-4;//FIXME
          if((nuTrack.closestCrtHitDistance < fDistanceLimit && (crtHitTime<0 || crtHitTime>fBeamTimeLimit))
              || (nuTrack2.track.Length()>20 && nuTrack2.closestCrtHitDistance < fDistanceLimit && (crtHitTime2<0 || crtHitTime2>fBeamTimeLimit))){
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
        double crtHitTime = ((double)(int)nuTrack.closestCrtHit.ts1_ns)*1e-4;//FIXME
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

      if(isTrueNu){ if(fVerbose) std::cout<<"NU EVENT REMAINING!\n"; nNuPfpRemain++; }
      else{ if(fVerbose) std::cout<<"COSMIC REMAINING!\n"; nCrPfpRemain++; }
      if(isTrueNuMu){ if(fVerbose) std::cout<<"NU MU EVENT REMAINING!\n"; nNuMuPfpRemain++; }

      //Loop over the nuTracks
      for(size_t i = 0; i < nuTracks.size(); i++){
        //Get the true ID
        std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(nuTracks[i].track.ID());
        int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
        if(particles.find(trueId) == particles.end()) continue;
        //Check it matches a true muon
        if(std::abs(particles[trueId].PdgCode()) != 13) continue;
        //Is it from neutrino
        if(std::find(lepParticleIds.begin(), lepParticleIds.end(), trueId) != lepParticleIds.end()) {
          hNuMuMomRemain->Fill(particles[trueId].P());
        }
        //Is is from cosmic
        if(std::find(nuParticleIds.begin(), nuParticleIds.end(), trueId) == nuParticleIds.end()) {
          hCrMuMomRemain->Fill(particles[trueId].P());
          hCrMuMomRemainHigh->Fill(particles[trueId].P());
        }
      }

      //----------------------------------------------------------------------------------------------------------
      //                                          NUMU SELECTION
      //----------------------------------------------------------------------------------------------------------

      //Get primary vertex from pandora
      TVector3 nuVtx = nuTracks[0].track.Vertex<TVector3>();
      //Fiducial volume cut of vertex (16.5 cm from anode, 30 cm from top/bottom (?), 95cm total from front+back
      if(nuVtx.X() < -184.5 || nuVtx.X() > 184.15 || nuVtx.Y() < -185 || nuVtx.Y() > 185 || nuVtx.Z() < 15. || nuVtx.Z() > 420) continue;
      //If longest track is contained, check pid chi2, select mu/pi, select length > 50 cm
      bool exits = nuTracks[0].endInFidStop;
      double length = nuTracks[0].track.Length();
      if(!exits && length < 50.) continue;
      //If longest track exits, select length > 100 cm
      if(exits && length < 100.) continue;

      nPfPartsSelect++;
      if(isTrueNu){ if(fVerbose) std::cout<<"NU EVENT SELECTED!\n"; nNuPfpSelect++; }
      else{ if(fVerbose) std::cout<<"COSMIC SELECTED!\n"; nCrPfpSelect++; }
      if(isTrueNuMu){ if(fVerbose) std::cout<<"NU MU EVENT SELECTED!\n"; nNuMuPfpSelect++; }

      //Loop over the nuTracks
      for(size_t i = 0; i < nuTracks.size(); i++){
        //Get the true ID
        std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(nuTracks[i].track.ID());
        int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
        if(particles.find(trueId) == particles.end()) continue;
        //Check it matches a true muon
        if(std::abs(particles[trueId].PdgCode()) != 13) continue;
        //Is it from neutrino
        if(std::find(lepParticleIds.begin(), lepParticleIds.end(), trueId) != lepParticleIds.end()) {
          hNuMuMomSelect->Fill(particles[trueId].P());
        }
        //Is is from cosmic
        if(std::find(nuParticleIds.begin(), nuParticleIds.end(), trueId) == nuParticleIds.end()) {
          hCrMuMomSelect->Fill(particles[trueId].P());
          hCrMuMomSelectHigh->Fill(particles[trueId].P());
        }
      }

      double recoMuMomentum = 0.;
      const recob::Track& trk = nuTracks[0].track;
      //Momentum calculation
      if(!exits){
        recoMuMomentum = fRangeFitter.GetTrackMomentum(length, 13);
      }
      else{
        recob::MCSFitResult mcsResult = fMcsFitter.fitMcs(trk);
        recoMuMomentum = mcsResult.bestMomentum();
      }
      if(recoMuMomentum <= 0) continue;

      std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(nuTracks[0].track.ID());
      int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
      if(particles.find(trueId) == particles.end()) continue;
      //Is it from neutrino
      if(std::find(lepParticleIds.begin(), lepParticleIds.end(), trueId) != lepParticleIds.end()) {
        hNuMuRecoMom->Fill(recoMuMomentum);
      }
      //Is is from cosmic
      else if(std::find(nuParticleIds.begin(), nuParticleIds.end(), trueId) != nuParticleIds.end()) {
        hNuRecoMom->Fill(recoMuMomentum);
      }
      else{
        hCrRecoMom->Fill(recoMuMomentum);
      }

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
             <<"From nu mu                       = "<<nNuMuPfpRemain<<"\n"
             <<"---------------------------------------------\n"
             <<"Number remaining after selection = "<<nPfPartsSelect<<"\n"
             <<"From neutrinos                   = "<<nNuPfpSelect<<"\n"
             <<"From cosmic rays                 = "<<nCrPfpSelect<<"\n"
             <<"From nu mu                       = "<<nNuMuPfpSelect<<"\n";

    std::ofstream myfile;
    myfile.open("results.txt");
    myfile<<nPfParts<<","<<nNuPfp<<","<<nCrPfp<<","<<nNuMuPfp<<","<<nNuPfpFid<<","<<nCrPfpFid<<","<<nNuMuPfpFid<<","<<nPfPartsRemain<<","<<nNuPfpRemain<<","<<nCrPfpRemain<<","<<nNuMuPfpRemain<<","<<nPfPartsSelect<<","<<nNuPfpSelect<<","<<nCrPfpSelect<<","<<nNuMuPfpSelect;
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

