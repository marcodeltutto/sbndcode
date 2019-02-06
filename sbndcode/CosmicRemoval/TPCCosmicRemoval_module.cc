////////////////////////////////////////////////////////////////////////
// Class:       TPCCosmicRemoval
// Module Type: analyzer
// File:        TPCCosmicRemovalAna_module.cc
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbndcode/CRT/CRTUtils/CRTAnaUtils.h"
#include "sbndcode/CRT/CRTUtils/CRTT0MatchAlg.h"
#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/CRTTruthRecoAlg.h"

// LArSoft includes
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
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

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <https://root.cern.ch/doc/master/annotated.html>
#include "TVector3.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TFile.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TPolyLine3D.h"
#include "TGraphAsymmErrors.h"

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

  struct TrackMatch{
    int trueID;
    double trueTime;
    bool isNu;
    bool isLep;
    bool isRemoved;
  };

  struct RecoTruth{
    std::vector<simb::MCParticle> particles;
    std::vector<recob::Track> tpcTracks;
    std::map<int, TrackMatch> matchingMap;
  };

  class TPCCosmicRemoval : public art::EDAnalyzer {
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

      fhicl::Atom<art::InputTag> PidModuleLabel {
        Name("PidModuleLabel"),
        Comment("tag of PID producer data product")
      };

      fhicl::Atom<bool> Verbose {
        Name("Verbose"),
        Comment("Print information about what's going on")
      };

      fhicl::Atom<int> TrackID {
        Name("TrackID"),
        Comment("Track ID")
      };
      
    }; // Config
 
    using Parameters = art::EDAnalyzer::Table<Config>;
 
    // Constructor: configures module
    explicit TPCCosmicRemoval(Parameters const& config);
 
    // Called once, at start of the job
    virtual void beginJob() override;
 
    // Called once per event
    virtual void analyze(const art::Event& event) override;

    // Called once, at end of the job
    virtual void endJob() override;

    bool EntersFiducial(recob::Track track, double fiducial, double fiducialTop);

    bool InFiducial(geo::Point_t point, double fiducial, double fiducialTop);

    int DetectedInTPC(std::vector<art::Ptr<recob::Hit>> hits);

    double T0FromCpaStitching(recob::Track track, std::vector<recob::Track> tracks, double stitchDist, double stitchAngle);

    double T0FromApaCross(recob::Track track, std::vector<double> t0s, int tpc);

    bool IsStoppingMuon(art::Ptr<anab::ParticleID> pid, art::Ptr<anab::Calorimetry> calo);

    double StoppingEnd(std::vector<art::Ptr<anab::Calorimetry>> calo, geo::Point_t end, int id);
    
    // Function to draw true and reco tracks
    void DrawTrueTracks(RecoTruth truthMatch, bool truth, bool tpctracks, int id);

  private:

    // fcl file parameters
    art::InputTag fSimModuleLabel;      ///< name of detsim producer
    art::InputTag fCrtHitModuleLabel; ///< name of CRT track producer
    art::InputTag fTpcTrackModuleLabel; ///< name of TPC track producer
    art::InputTag fCaloModuleLabel; ///< name of calorimetry producer
    art::InputTag fPidModuleLabel; ///< name of calorimetry producer
    bool          fVerbose;             ///< print information about what's going on
    int           fTrackID;

    // histograms
    // True time of neutrino particles to get beam window
    TH1D* hTrueNuMuTime;
    // Contained length of primary muon within fiducial volume
    TH1D* hTrueNuMuLength;

    TH1D* hStoppingChi2;
    TH1D* hOtherChi2;
    TGraph* gdEdxResRg;

    TH1D* hFidTotalMu;

    TH1D* hFidTotalCos;

    TGraphAsymmErrors* fFidExitMu;
    TH1D* hFidExitCutMu;

    TGraphAsymmErrors* fFidExitCos;
    TH1D* hFidExitCutCos;

    TGraphAsymmErrors* fStitchAngleMu;
    TH1D* hStitchAngleTotalMu;
    TH1D* hStitchAngleCutMu;

    TGraphAsymmErrors* fStitchAngleCos;
    TH1D* hStitchAngleTotalCos;
    TH1D* hStitchAngleCutCos;

    TGraphAsymmErrors* fStitchDistMu;
    TH1D* hStitchDistTotalMu;
    TH1D* hStitchDistCutMu;

    TGraphAsymmErrors* fStitchDistCos;
    TH1D* hStitchDistTotalCos;
    TH1D* hStitchDistCutCos;

    TGraphAsymmErrors* fTrackAngleMu;
    TH1D* hTrackAngleTotalMu;
    TH1D* hTrackAngleCutMu;

    TGraphAsymmErrors* fTrackAngleCos;
    TH1D* hTrackAngleTotalCos;
    TH1D* hTrackAngleCutCos;

    TGraphAsymmErrors* fTrackDistMu;
    TH1D* hTrackDistTotalMu;
    TH1D* hTrackDistCutMu;

    TGraphAsymmErrors* fTrackDistCos;
    TH1D* hTrackDistTotalCos;
    TH1D* hTrackDistCutCos;

    TGraphAsymmErrors* fHitDistMu;
    TH1D* hHitDistTotalMu;
    TH1D* hHitDistCutMu;

    TGraphAsymmErrors* fHitDistCos;
    TH1D* hHitDistTotalCos;
    TH1D* hHitDistCutCos;

    TGraphAsymmErrors* fHitLengthMu;
    TH1D* hHitLengthTotalMu;
    TH1D* hHitLengthCutMu;

    TGraphAsymmErrors* fHitLengthCos;
    TH1D* hHitLengthTotalCos;
    TH1D* hHitLengthCutCos;

    // Other variables shared between different methods.
    geo::GeometryCore const* fGeometryService;                 ///< pointer to Geometry provider
    detinfo::DetectorProperties const* fDetectorProperties;    ///< pointer to detector properties provider
    detinfo::DetectorClocks const* fDetectorClocks;            ///< pointer to detector clocks provider

    CRTTruthRecoAlg truthAlg;

    double timeLim = 20.;

    // Performance Counters
    // Tracks from primary cosmic particles
    int nCosmicTracks = 0;
    int nCosmicRemoved = 0;
    int nCosmicTracksFid = 0;
    int nCosmicRemovedExit = 0;
    int nCosmicRemovedStop = 0;
    int nCosmicRemovedDiffTPC = 0;
    int nCosmicRemovedStitch = 0;
    int nCosmicRemovedTrack = 0;
    int nCosmicRemovedHit = 0;
    // Tracks from neutrino interactions
    int nNuTracks = 0;
    int nNuRemoved = 0;
    int nNuTracksFid = 0;
    int nNuRemovedExit = 0;
    int nNuRemovedStop = 0;
    int nNuRemovedDiffTPC = 0;
    int nNuRemovedStitch = 0;
    int nNuRemovedTrack = 0;
    int nNuRemovedHit = 0;

    int nLepTracks = 0;
    int nLepRemoved = 0;
    int nLepTracksFid = 0;
    int nLepRemovedExit = 0;
    int nLepRemovedStop = 0;
    int nLepRemovedDiffTPC = 0;
    int nLepRemovedStitch = 0;
    int nLepRemovedTrack = 0;
    int nLepRemovedHit = 0;

    int nCuts = 10;
    double fidStart = 2.5;
    double fidStep = 5.;
    double saStart = 2.;
    double saStep = 4.;
    double sdStart = 5.;
    double sdStep = 10.;
    double taStart = 2.;
    double taStep = 4.;
    double tdStart = 8;
    double tdStep = 16;
    double hdStart = 5;
    double hdStep = 10;
    double hlStart = 5;
    double hlStep = 10;

  }; // class TPCCosmicRemoval


  // Constructor
  TPCCosmicRemoval::TPCCosmicRemoval(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel       (config().SimModuleLabel())
    , fCrtHitModuleLabel    (config().CrtHitModuleLabel())
    , fTpcTrackModuleLabel  (config().TpcTrackModuleLabel())
    , fCaloModuleLabel      (config().CaloModuleLabel())
    , fPidModuleLabel       (config().PidModuleLabel())
    , fVerbose              (config().Verbose())
    , fTrackID              (config().TrackID())
  {

    // Get a pointer to the geometry service provider
    fGeometryService    = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 
    fDetectorClocks     = lar::providerFrom<detinfo::DetectorClocksService>(); 

  } // TPCCosmicRemoval()


  void TPCCosmicRemoval::beginJob()
  {
    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;

    hTrueNuMuTime = tfs->make<TH1D>("hTrueNuMuTime", "", 100, -10, 10);
    hTrueNuMuLength = tfs->make<TH1D>("hTrueNuMuLength", "", 100, 0, 500);

    hStoppingChi2 = tfs->make<TH1D>("StoppingChi2", "", 100, -1.5, 1.5);
    hOtherChi2 = tfs->make<TH1D>("OtherChi2", "", 100, -1.5, 1.5);
    gdEdxResRg       = tfs->makeAndRegister<TGraph>("dEdxResRg",      ";Residual range (cm);dE/dx");

    fFidExitMu       = tfs->makeAndRegister<TGraphAsymmErrors>("FidExitMu",      ";Fiducial volume cut (cm);Percentage removed");
    fStitchAngleMu   = tfs->makeAndRegister<TGraphAsymmErrors>("StitchAngleMu",  ";Stitch angle (deg);Percentage removed"      );
    fStitchDistMu    = tfs->makeAndRegister<TGraphAsymmErrors>("StitchDistMu",   ";Stitch dist (cm);Percentage removed"        );
    fTrackAngleMu    = tfs->makeAndRegister<TGraphAsymmErrors>("TrackAngleMu",   ";Max angle (deg);Percentage removed"         );
    fTrackDistMu     = tfs->makeAndRegister<TGraphAsymmErrors>("TrackDistMu",    ";Max dist (cm);Percentage removed"           );
    fHitDistMu       = tfs->makeAndRegister<TGraphAsymmErrors>("HitDistMu",      ";Max dist (cm);Percentage removed"           );
    fHitLengthMu     = tfs->makeAndRegister<TGraphAsymmErrors>("HitLengthMu",    ";Min length (cm);Percentage removed"         );

    fFidExitCos       = tfs->makeAndRegister<TGraphAsymmErrors>("FidExitCos",      ";Fiducial volume cut (cm);Percentage removed");
    fStitchAngleCos   = tfs->makeAndRegister<TGraphAsymmErrors>("StitchAngleCos",  ";Stitch angle (deg);Percentage removed"      );
    fStitchDistCos    = tfs->makeAndRegister<TGraphAsymmErrors>("StitchDistCos",   ";Stitch dist (cm);Percentage removed"        );
    fTrackAngleCos    = tfs->makeAndRegister<TGraphAsymmErrors>("TrackAngleCos",   ";Max angle (deg);Percentage removed"         );
    fTrackDistCos     = tfs->makeAndRegister<TGraphAsymmErrors>("TrackDistCos",    ";Max dist (cm);Percentage removed"           );
    fHitDistCos       = tfs->makeAndRegister<TGraphAsymmErrors>("HitDistCos",      ";Max dist (cm);Percentage removed"           );
    fHitLengthCos     = tfs->makeAndRegister<TGraphAsymmErrors>("HitLengthCos",    ";Min length (cm);Percentage removed"         );

    hFidTotalMu         = tfs->make<TH1D>("FidTotalMu",         "", 10, 0, 50);
    hStitchAngleTotalMu = tfs->make<TH1D>("StitchAngleTotalMu", "", 10, 0, 40);
    hStitchDistTotalMu  = tfs->make<TH1D>("StitchDistTotalMu",  "", 10, 0, 100);
    hTrackAngleTotalMu  = tfs->make<TH1D>("TrackAngleTotalMu",  "", 10, 0, 40);
    hTrackDistTotalMu   = tfs->make<TH1D>("TrackDistTotalMu",   "", 10, 0, 160);
    hHitDistTotalMu     = tfs->make<TH1D>("HitDistTotalMu",     "", 10, 0, 100);
    hHitLengthTotalMu   = tfs->make<TH1D>("HitLengthTotalMu",   "", 10, 0, 100);

    hFidTotalCos         = tfs->make<TH1D>("FidTotalCos",         "", 10, 0, 50);
    hStitchAngleTotalCos = tfs->make<TH1D>("StitchAngleTotalCos", "", 10, 0, 40);
    hStitchDistTotalCos  = tfs->make<TH1D>("StitchDistTotalCos",  "", 10, 0, 100);
    hTrackAngleTotalCos  = tfs->make<TH1D>("TrackAngleTotalCos",  "", 10, 0, 40);
    hTrackDistTotalCos   = tfs->make<TH1D>("TrackDistTotalCos",   "", 10, 0, 160);
    hHitDistTotalCos     = tfs->make<TH1D>("HitDistTotalCos",     "", 10, 0, 100);
    hHitLengthTotalCos   = tfs->make<TH1D>("HitLengthTotalCos",   "", 10, 0, 100);

    hFidExitCutMu     = tfs->make<TH1D>("FidExitCutMu",     "", 10, 0, 50);
    hStitchAngleCutMu = tfs->make<TH1D>("StitchAngleCutMu", "", 10, 0, 40);
    hStitchDistCutMu  = tfs->make<TH1D>("StitchDistCutMu",  "", 10, 0, 100);
    hTrackAngleCutMu  = tfs->make<TH1D>("TrackAngleCutMu",  "", 10, 0, 40);
    hTrackDistCutMu   = tfs->make<TH1D>("TrackDistCutMu",   "", 10, 0, 160);
    hHitDistCutMu     = tfs->make<TH1D>("HitDistCutMu",     "", 10, 0, 100);
    hHitLengthCutMu   = tfs->make<TH1D>("HitLengthCutMu",   "", 10, 0, 100);

    hFidExitCutCos     = tfs->make<TH1D>("FidExitCutCos",     "", 10, 0, 50);
    hStitchAngleCutCos = tfs->make<TH1D>("StitchAngleCutCos", "", 10, 0, 40);
    hStitchDistCutCos  = tfs->make<TH1D>("StitchDistCutCos",  "", 10, 0, 100);
    hTrackAngleCutCos  = tfs->make<TH1D>("TrackAngleCutCos",  "", 10, 0, 40);
    hTrackDistCutCos   = tfs->make<TH1D>("TrackDistCutCos",   "", 10, 0, 160);
    hHitDistCutCos     = tfs->make<TH1D>("HitDistCutCos",     "", 10, 0, 100);
    hHitLengthCutCos   = tfs->make<TH1D>("HitLengthCutCos",   "", 10, 0, 100);

    // Initial output
    if(fVerbose) std::cout<<"----------------- TPC Cosmic Removal Ana Module -------------------"<<std::endl;

  }// TPCCosmicRemoval::beginJob()


  void TPCCosmicRemoval::analyze(const art::Event& event)
  {
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

    // Fetch basic event info
    if(fVerbose){
      std::cout<<"============================================"<<std::endl
               <<"Run = "<<event.run()<<", SubRun = "<<event.subRun()<<", Event = "<<event.id().event()<<std::endl
               <<"============================================"<<std::endl;
    }

    // Retrieve all the truth info in the events
    auto particleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fSimModuleLabel);

    // Retrieve list of CRT tracks
    art::Handle< std::vector<crt::CRTHit> > crtHitHandle;
    std::vector<art::Ptr<crt::CRTHit> > crtHitsPtr;
    if (event.getByLabel(fCrtHitModuleLabel, crtHitHandle))
      art::fill_ptr_vector(crtHitsPtr, crtHitHandle);

    // Retrieve the TPC tracks
    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTpcTrackModuleLabel);

    // Get track to hit associations
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTpcTrackModuleLabel);
    art::FindManyP<anab::Calorimetry> findManyCalo(tpcTrackHandle, event, fCaloModuleLabel);
    art::FindManyP<anab::ParticleID> findManyPid(tpcTrackHandle, event, fPidModuleLabel);

    if(fVerbose) std::cout<<"Number of CRT hits = "<<crtHitHandle->size()<<std::endl
                          <<"Number of TPC tracks = "<<tpcTrackHandle->size()<<std::endl;

    // Do track reconstruction from the hits
    std::vector<crt::CRTHit> crtHits;
    for(auto const& crtHit : (*crtHitHandle)){
      crtHits.push_back(crtHit);
    }
    std::vector<crt::CRTTrack> crtTracks = CRTAnaUtils::CreateCRTTracks(crtHitsPtr, 0.2, 30., true, 25.); 
    if(fVerbose) std::cout<<"Number of CRTTracks = "<<crtTracks.size()<<std::endl;

    std::vector<std::vector<art::Ptr<crt::CRTHit>>> crtT0Ptr = CRTAnaUtils::CreateCRTTzeros(crtHitsPtr, 2.);
    std::vector<double> crtT0;
    for(size_t i = 0; i < crtT0Ptr.size(); i++){
      double t0 = 0;
      double npts = 0;
      for(size_t j = 0; j < crtT0Ptr[i].size(); j++){
        if(crtT0Ptr[i][j]->tagger != "volTaggerTopHigh_0"){
          t0 += crtT0Ptr[i][0]->ts1_ns*1e-4;
          npts++;
        }
      }
      if(t0 != 0) crtT0.push_back(t0/npts);
    }

    // Loop over true particles in readout window
    std::map<int, simb::MCParticle> particles;
    std::vector<simb::MCParticle> parts;
    std::vector<int> nuParticleIds;
    std::vector<int> lepParticleIds;

    for (auto const& particle: (*particleHandle)){

      int partId = particle.TrackId();
      particles[partId] = particle;
      parts.push_back(particle);
      art::Ptr<simb::MCTruth> truth = pi_serv->TrackIdToMCTruth_P(partId);
      if(truth->Origin() == simb::kBeamNeutrino){
        geo::Point_t vtx;
        vtx.SetX(truth->GetNeutrino().Nu().Vx()); vtx.SetY(truth->GetNeutrino().Nu().Vy()); vtx.SetZ(truth->GetNeutrino().Nu().Vz());
        double time = particle.T() * 1e-3; // [us]
        if(fVerbose && particle.Mother()==0) 
          std::cout<<"Nu VTX = "<<vtx<<" ID = "<<partId<<" pdg = "<<particle.PdgCode()<<" time = "<<time<<" length = "<<particle.Trajectory().TotalLength()
                   <<" start = ("<<particle.Vx()<<", "<<particle.Vy()<<", "<<particle.Vz()<<") end = ("<<particle.EndX()<<", "<<particle.EndY()<<", "<<particle.EndZ()<<")\n";
        if(!InFiducial(vtx, 0, 0)) continue;
        if(std::abs(particle.PdgCode())==13 && particle.Mother()==0) lepParticleIds.push_back(partId);
        nuParticleIds.push_back(partId);

        if(std::abs(particle.PdgCode())==13){ 
          hTrueNuMuTime->Fill(time);
          hTrueNuMuLength->Fill(particle.Trajectory().TotalLength());
        }

      }

    }

    if(fVerbose) std::cout<<"Number of true particles = "<<particles.size()<<"\n\n";

    std::map<int, TrackMatch> matchingMap;
    std::vector<recob::Track> tpcTracks;
    std::vector<recob::Track> tpcTracksTPC0;
    std::vector<recob::Track> tpcTracksTPC1;

    // Loop over the tpc tracks
    for(auto const& tpcTrack : (*tpcTrackHandle)){
      tpcTracks.push_back(tpcTrack);

      std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
      int tpc = DetectedInTPC(hits);
      double startX = tpcTrack.Start().X();
      double endX = tpcTrack.End().X();

      if(tpc == 0 && !(startX>0 || endX>0)) tpcTracksTPC0.push_back(tpcTrack);
      else if(tpc == 1 && !(startX<0 || endX<0)) tpcTracksTPC1.push_back(tpcTrack);

    }

    for(auto const& tpcTrack : (*tpcTrackHandle)){

      if(fVerbose) std::cout<<"------>Track "<<tpcTrack.ID()<<":\n";

      TrackMatch trackMatch;
      trackMatch.trueTime = -99999;
      trackMatch.isNu = false;
      trackMatch.isLep = false;
      trackMatch.isRemoved = false;

      // Match to the true particle
      std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
      int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
      trackMatch.trueID = trueId;
      if (particles.find(trueId) == particles.end()){ 
        if (fVerbose) std::cout<<"No valid true track!\n\n"; 
        matchingMap[tpcTrack.ID()] = trackMatch;
        continue; 
      }
      // Get the true T0
      double trueTime = particles[trueId].T()*1e-3;
      trackMatch.trueTime = trueTime;

      std::vector<art::Ptr<anab::Calorimetry>> calos = findManyCalo.at(tpcTrack.ID());
      if(fVerbose) std::cout<<"Calo size = "<<(calos[0]->dEdx()).size()<<" "<<(calos[1]->dEdx()).size()<<" "<<(calos[2]->dEdx()).size()<<"\n";
      if(tpcTrack.ID()==fTrackID){
        size_t npts = (calos[2]->dEdx()).size();
        for(size_t i = 0; i < npts; i++){
          if(calos[2]->dEdx()[i] < 30){
            gdEdxResRg->SetPoint(i, calos[2]->ResidualRange()[i], calos[2]->dEdx()[i]);
          }
        }
        gdEdxResRg->Draw();
      }
      std::vector<art::Ptr<anab::ParticleID>> pids = findManyPid.at(tpcTrack.ID());

      if(fVerbose) std::cout<<"trueID = "<<trueId<<" pdg = "<<particles[trueId].PdgCode()<<" time = "<<trueTime<<" reco length = "<<tpcTrack.Length()
                            <<" reco start = "<<tpcTrack.Start()<<" end = "<<tpcTrack.End()<<"\n";

      // Is track from a neutrino interaction
      bool isNu = false;
      bool isLep = false;
      bool isCos = false;
      bool removed = false;
      if(std::find(lepParticleIds.begin(), lepParticleIds.end(), trueId) != lepParticleIds.end()){
        // Set neutrino flag
        isLep = true;
        nLepTracks++;
        if(fVerbose) std::cout<<"Primary lepton!\n";
        //Fill total histograms
        for(size_t i = 0; i < 10; i++){
          hFidTotalMu->Fill(fidStart + i*fidStep);
          hStitchAngleTotalMu->Fill(saStart + i*saStep);
          hStitchDistTotalMu->Fill(sdStart + i*sdStep);
          hTrackAngleTotalMu->Fill(taStart + i*taStep);
          hTrackDistTotalMu->Fill(tdStart + i*tdStep); 
          hHitDistTotalMu->Fill(hdStart + i*hdStep);    
          hHitLengthTotalMu->Fill(hlStart + i*hlStep);  
        }
      }
      if(std::find(nuParticleIds.begin(), nuParticleIds.end(), trueId) != nuParticleIds.end()){
        // Set neutrino flag
        isNu = true;
        nNuTracks++;
        if(fVerbose) std::cout<<"From neutrino!\n";
      }
      // Else is it from a primary cosmic (true length > 500)
      else if(particles[trueId].Trajectory().TotalLength() > 500.){
        // Set primary cosmic flag
        isCos = true;
        nCosmicTracks++;
        if(fVerbose) std::cout<<"Primary cosmic!\n";
        //Fill total histograms
        for(size_t i = 0; i < 10; i++){
          hFidTotalCos->Fill(fidStart + i*fidStep);
          hStitchAngleTotalCos->Fill(saStart + i*saStep);
          hStitchDistTotalCos->Fill(sdStart + i*sdStep);
          hTrackAngleTotalCos->Fill(taStart + i*taStep);
          hTrackDistTotalCos->Fill(tdStart + i*tdStep); 
          hHitDistTotalCos->Fill(hdStart + i*hdStep);    
          hHitLengthTotalCos->Fill(hlStart + i*hlStep);  
        }
      }
      //geo::Point_t start {particles[trueId].Vx(), particles[trueId].Vy(), particles[trueId].Vz()};
      //geo::Point_t end {particles[trueId].EndX(), particles[trueId].EndY(), particles[trueId].EndZ()};
      /*if(isCos){// && ((InFiducial(end, 0, 0) && !InFiducial(start, 0,0))||(!InFiducial(end, 0, 0) && InFiducial(start, 0,0)))){
        if(!InFiducial(tpcTrack.Start(), 5, 5) && InFiducial(tpcTrack.End(), 5, 5)){
          hStoppingChi2->Fill(StoppingEnd(calos[2], tpcTrack.End(), tpcTrack.ID()));
        }
        if(InFiducial(tpcTrack.Start(), 5, 5) && !InFiducial(tpcTrack.End(), 5, 5)){
          hStoppingChi2->Fill(StoppingEnd(calos[2], tpcTrack.Start(), tpcTrack.ID()));
        }
      }
      else if (isLep){
        if(!InFiducial(tpcTrack.Start(), 5, 5) && InFiducial(tpcTrack.End(), 5, 5)){
          hOtherChi2->Fill(StoppingEnd(calos[2], tpcTrack.End(), tpcTrack.ID()));
        }
        if(InFiducial(tpcTrack.Start(), 5, 5) && !InFiducial(tpcTrack.End(), 5, 5)){
          hOtherChi2->Fill(StoppingEnd(calos[2], tpcTrack.Start(), tpcTrack.ID()));
        }
      }*/


      // FIDUCIAL VOLUME CUT
      // Remove any tracks that enter and exit the fiducial volume
      if(!removed && !InFiducial(tpcTrack.Start(), 10, 10) && !InFiducial(tpcTrack.End(), 10, 10)){ 
        removed = true;
        if(fVerbose) std::cout<<"REMOVED! - exits fiducial\n";
        if(isNu) nNuRemovedExit++;
        if(isLep) nLepRemovedExit++;
        // If primary cosmic flag
        if(isCos) nCosmicRemovedExit++;
      }

      if(isNu && !removed) nNuTracksFid++;
      if(isLep && !removed) nLepTracksFid++;
      if(isCos && !removed) nCosmicTracksFid++;


      /*std::cout<<"pids size = "<<pids.size()<<" calos size = "<<calos.size()<<"\n";
      if(!InFiducial(tpcTrack.Start(), 5, 5)){
        std::cout<<"Start not in fiducial\n";
      }
      if(!InFiducial(tpcTrack.End(), 5, 5)){
        std::cout<<"End not in fiducial\n";
      }*/
      /*if((!removed && IsStoppingMuon(pids[2], calos[2])) && (!InFiducial(tpcTrack.Start(), 5, 5) || !InFiducial(tpcTrack.End(), 5, 5))){
        removed = true;
        if(fVerbose) std::cout<<"REMOVED! - stopping muon\n";
        if(isNu) nNuRemoved++;
        if(isLep) nLepRemoved++;
        if(isCos) nCosmicRemoved++;
      }
      if(IsStoppingMuon(pids[2], calos[2]) && (!InFiducial(tpcTrack.Start(), 5, 5) || !InFiducial(tpcTrack.End(), 5, 5))){
        if(isNu) nNuRemovedStop++;
        if(isLep) nLepRemovedStop++;
        if(isCos) nCosmicRemovedStop++;
      }*/

      if(!removed && !InFiducial(tpcTrack.Start(), 5, 5) && InFiducial(tpcTrack.End(), 5, 5)){
        if(StoppingEnd(calos, tpcTrack.End(), tpcTrack.ID())){
          removed = true;
          if(fVerbose) std::cout<<"REMOVED! - stopping muon\n";
          if(isNu) nNuRemoved++;
          if(isLep) nLepRemoved++;
          if(isCos) nCosmicRemoved++;
        }
      }
      if(!InFiducial(tpcTrack.Start(), 5, 5) && InFiducial(tpcTrack.End(), 5, 5) && StoppingEnd(calos, tpcTrack.End(), tpcTrack.ID())){
        if(isNu) nNuRemovedStop++;
        if(isLep) nLepRemovedStop++;
        if(isCos) nCosmicRemovedStop++;
      }
      if(!removed && InFiducial(tpcTrack.Start(), 5, 5) && !InFiducial(tpcTrack.End(), 5, 5)){
        if(StoppingEnd(calos, tpcTrack.Start(), tpcTrack.ID())){
          removed = true;
          if(fVerbose) std::cout<<"REMOVED! - stopping muon\n";
          if(isNu) nNuRemoved++;
          if(isLep) nLepRemoved++;
          if(isCos) nCosmicRemoved++;
        }
      }
      if(InFiducial(tpcTrack.Start(), 5, 5) && !InFiducial(tpcTrack.End(), 5, 5) && StoppingEnd(calos, tpcTrack.Start(), tpcTrack.ID())){
        if(isNu) nNuRemovedStop++;
        if(isLep) nLepRemovedStop++;
        if(isCos) nCosmicRemovedStop++;
      }

      // TPC CUT
      // Remove any tracks that are detected in one TPC and reconstructed in another
      // Check if track is stitched
      // Loop over hits, get tpc, make sure it's the same for all
      // If it is check the start/end points are in same TPC
      // If not get the time from stitching, no time == in time with beam
      int tpc = DetectedInTPC(hits);
      double startX = tpcTrack.Start().X();
      double endX = tpcTrack.End().X();
      if(tpc < 0){
        if(fVerbose) std::cout<<"TRACK STITCHED\n";
      }
      else if(tpc == 0){
        if(!removed && (startX>0 || endX>0)){ 
          removed = true;
          if(fVerbose) std::cout<<"REMOVED! - diff tpc\n";
          if(isNu) nNuRemoved++;
          if(isLep) nLepRemoved++;
          // If primary cosmic flag
          if(isCos) nCosmicRemoved++;
        }
        if(startX>0 || endX>0){ 
          if(isNu) nNuRemovedDiffTPC++;
          if(isLep) nLepRemovedDiffTPC++;
          if(isCos) nCosmicRemovedDiffTPC++;
        }
      }
      else if(tpc == 1){
        if(!removed && (startX<0 || endX<0)){
          removed = true;
          if(fVerbose) std::cout<<"REMOVED! - diff tpc\n";
          if(isNu) nNuRemoved++;
          if(isLep) nLepRemoved++;
          // If primary cosmic flag
          if(isCos) nCosmicRemoved++;
        }
        if(startX<0 || endX<0){
          if(isNu) nNuRemovedDiffTPC++;
          if(isLep) nLepRemovedDiffTPC++;
          if(isCos) nCosmicRemovedDiffTPC++;
        }
      }
      else std::cout<<"SOMETHING WRONG\n";

      // TIME CUTS
      // Match CPA crossers, remove any outside of beam window
      double stitchTime = -99999;
      if(tpc == 0){
        stitchTime = T0FromCpaStitching(tpcTrack, tpcTracksTPC1, 50., 15.);
      }
      else if(tpc == 1){
        stitchTime = T0FromCpaStitching(tpcTrack, tpcTracksTPC0, 50., 15.);
      }
      if(stitchTime != -99999 && std::abs(stitchTime)>timeLim){
        if(isNu) nNuRemovedStitch++;
        if(isLep) nLepRemovedStitch++;
        if(isCos) nCosmicRemovedStitch++;
      }

      // Do CRT time matching
      // Try to get T0 from CRTTracks
      double trackTime = CRTAnaUtils::T0FromCRTTracks(tpcTrack, crtTracks, tpc, 0.2, 80.);
      if(trackTime != -99999 && std::abs(trackTime)>timeLim){
        if(isNu) nNuRemovedTrack++;
        if(isLep) nLepRemovedTrack++;
        if(isCos) nCosmicRemovedTrack++;
      }

      // Try to get T0 from CRTHits
      double hitTime = CRTAnaUtils::T0FromCRTHits(tpcTrack, crtHits, tpc, 20., 0.5, 15.);
      if(hitTime != -99999 && std::abs(hitTime)>timeLim){
        if(isNu) nNuRemovedHit++;
        if(isLep) nLepRemovedHit++;
        if(isCos) nCosmicRemovedHit++;
      }

      if(fVerbose) std::cout<<"True time = "<<trueTime<<" ticks, track time = "<<trackTime
                            <<" ticks, hit time = "<<hitTime<<" stitch time = "<<stitchTime<<"\n";

      double crossTime = T0FromApaCross(tpcTrack, crtT0, tpc);
      if(crossTime != -99999 && (!InFiducial(tpcTrack.Vertex(), 10., 10.) || !InFiducial(tpcTrack.End(), 10., 10.))) std::cout<<"Removed - cross & exit\n";
      if(crossTime != -99999 && (StoppingEnd(calos, tpcTrack.Start(), tpcTrack.ID()) || StoppingEnd(calos, tpcTrack.End(), tpcTrack.ID()))) std::cout<<"Removed - cross & stop\n";

      double matchedTime = -99999;
      if(trackTime != -99999) matchedTime = trackTime;
      else if(stitchTime != -99999) matchedTime = stitchTime;
      //else if(crossTime != -99999) matchedTime = crossTime;
      else if(hitTime != -99999) matchedTime = hitTime;

      if(fVerbose) std::cout<<"Matched time = "<<matchedTime<<"\n";

      if(!removed && matchedTime != -99999 && (std::abs(matchedTime) > timeLim)){
        removed = true;
        if(fVerbose) std::cout<<"REMOVED! - outside beam\n";
        if(isNu) nNuRemoved++;
        if(isLep) nLepRemoved++;
        // If primary cosmic flag
        if(isCos) nCosmicRemoved++;
      }

      if(fVerbose) std::cout<<"\n";

      for(size_t i = 0; i < 10; i++){
        double fidCut = fidStart + i*fidStep;
        double sa = saStart + i*saStep;
        double sd = sdStart + i*sdStep;
        double ta = taStart + i*taStep;
        double taRad = ta * TMath::Pi()/180.;
        double td = tdStart + i*tdStep;
        double hd = hdStart + i*hdStep;
        double hl = hlStart + i*hlStep;
        if(!InFiducial(tpcTrack.Start(), fidCut, fidCut) && !InFiducial(tpcTrack.End(), fidCut, fidCut)){
          if(isLep) hFidExitCutMu->Fill(fidCut);
          if(isCos) hFidExitCutCos->Fill(fidCut);
        }


        double time = -99999;
        if(tpc == 0){
          time = T0FromCpaStitching(tpcTrack, tpcTracksTPC1, 50., sa);
        }
        else if(tpc == 1){
          time = T0FromCpaStitching(tpcTrack, tpcTracksTPC0, 50., sa);
        }
        if(time != -99999 && std::abs(time)>timeLim){
          if(isLep) hStitchAngleCutMu->Fill(sa);
          if(isCos) hStitchAngleCutCos->Fill(sa);
        }

        time = -99999;
        if(tpc == 0){
          time = T0FromCpaStitching(tpcTrack, tpcTracksTPC1, sd, 15.);
        }
        else if(tpc == 1){
          time = T0FromCpaStitching(tpcTrack, tpcTracksTPC0, sd, 15.);
        }
        if(time != -99999 && std::abs(time) > timeLim){
          if(isLep) hStitchDistCutMu->Fill(sd);
          if(isCos) hStitchDistCutCos->Fill(sd);
        }

        time = CRTAnaUtils::T0FromCRTTracks(tpcTrack, crtTracks, tpc, taRad, 80.);
        if(time != -99999 && std::abs(time) > timeLim){
          if(isLep) hTrackAngleCutMu->Fill(ta);
          if(isCos) hTrackAngleCutCos->Fill(ta);
        }

        time = CRTAnaUtils::T0FromCRTTracks(tpcTrack, crtTracks, tpc, 0.2, td);
        if(time != -99999 && std::abs(time) > timeLim){
          if(isLep) hTrackDistCutMu->Fill(td);
          if(isCos) hTrackDistCutCos->Fill(td);
        }

        time = CRTAnaUtils::T0FromCRTHits(tpcTrack, crtHits, tpc, 20., 0.5, hd);
        if(time != -99999 && std::abs(time) > timeLim){
          if(isLep) hHitDistCutMu->Fill(hd);
          if(isCos) hHitDistCutCos->Fill(hd);
        }

        time = CRTAnaUtils::T0FromCRTHits(tpcTrack, crtHits, tpc, hl, 0.5, 30.);
        if(time != -99999 && std::abs(time) > timeLim){
          if(isLep) hHitLengthCutMu->Fill(hl);
          if(isCos) hHitLengthCutCos->Fill(hl);
        }
      }

      trackMatch.isNu = isNu;
      trackMatch.isLep = isLep;
      trackMatch.isRemoved = removed;
      matchingMap[tpcTrack.ID()] = trackMatch;

    }

    RecoTruth truthMatch;
    truthMatch.particles = parts;
    truthMatch.tpcTracks = tpcTracks;
    truthMatch.matchingMap = matchingMap;

    DrawTrueTracks(truthMatch, true, true, fTrackID);


  } // TPCCosmicRemoval::analyze()


  void TPCCosmicRemoval::endJob(){

    std::cout<<"Total nu tracks           = "<<nNuTracks<<"\n"
             <<"Tracks in fiducial volume = "<<nNuTracksFid<<"\n"
             <<"Removed in exit cut       = "<<nNuRemovedExit<<"\n"
             <<"Removed in stopping cut   = "<<nNuRemovedStop<<"\n"
             <<"Removed in diff TPC cut   = "<<nNuRemovedDiffTPC<<"\n"
             <<"Removed in stitch cut     = "<<nNuRemovedStitch<<"\n"
             <<"Removed in CRT track cut  = "<<nNuRemovedTrack<<"\n"
             <<"Removed in CRT hit cut    = "<<nNuRemovedHit<<"\n"
             <<"Total removed by cosmicID = "<<nNuRemoved<<"\n"
             <<"Percentage removed        = "<<(double)nNuRemoved/nNuTracksFid<<"\n"
             <<"-----------------------------------------------------------\n"
             <<"Total lepton tracks       = "<<nLepTracks<<"\n"
             <<"Tracks in fiducial volume = "<<nLepTracksFid<<"\n"
             <<"Removed in exit cut       = "<<nLepRemovedExit<<"\n"
             <<"Removed in stopping cut   = "<<nLepRemovedStop<<"\n"
             <<"Removed in diff TPC cut   = "<<nLepRemovedDiffTPC<<"\n"
             <<"Removed in stitch cut     = "<<nLepRemovedStitch<<"\n"
             <<"Removed in CRT track cut  = "<<nLepRemovedTrack<<"\n"
             <<"Removed in CRT hit cut    = "<<nLepRemovedHit<<"\n"
             <<"Total removed by cosmicID = "<<nLepRemoved<<"\n"
             <<"Percentage removed        = "<<(double)nLepRemoved/nLepTracksFid<<"\n"
             <<"-----------------------------------------------------------\n"
             <<"Total cosmic tracks       = "<<nCosmicTracks<<"\n"
             <<"Tracks in fiducial volume = "<<nCosmicTracksFid<<"\n"
             <<"Removed in exit cut       = "<<nCosmicRemovedExit<<"\n"
             <<"Removed in stopping cut   = "<<nCosmicRemovedStop<<"\n"
             <<"Removed in diff TPC cut   = "<<nCosmicRemovedDiffTPC<<"\n"
             <<"Removed in stitch cut     = "<<nCosmicRemovedStitch<<"\n"
             <<"Removed in CRT track cut  = "<<nCosmicRemovedTrack<<"\n"
             <<"Removed in CRT hit cut    = "<<nCosmicRemovedHit<<"\n"
             <<"Total removed by cosmicID = "<<nCosmicRemoved<<"\n"
             <<"Percentage removed        = "<<(double)nCosmicRemoved/nCosmicTracksFid<<"\n";

    std::ofstream myfile;
    myfile.open("results.txt");
  myfile<<nNuTracks<<","<<nNuTracksFid<<","<<nNuRemovedExit<<","<<nNuRemovedStop<<","<<nNuRemovedDiffTPC<<","<<nNuRemovedStitch<<","<<nNuRemovedTrack<<","<<nNuRemovedHit<<","<<nNuRemoved<<","<<(double)nNuRemoved/nNuTracksFid<<","<<nLepTracks<<","<<nLepTracksFid<<","<<nLepRemovedExit<<","<<nLepRemovedStop<<","<<nLepRemovedDiffTPC<<","<<nLepRemovedStitch<<","<<nLepRemovedTrack<<","<<nLepRemovedHit<<","<<nLepRemoved<<","<<(double)nLepRemoved/nLepTracksFid<<","<<nCosmicTracks<<","<<nCosmicTracksFid<<","<<nCosmicRemovedExit<<","<<nCosmicRemovedStop<<","<<nCosmicRemovedDiffTPC<<","<<nCosmicRemovedStitch<<","<<nCosmicRemovedTrack<<","<<nCosmicRemovedHit<<","<<nCosmicRemoved<<","<<(double)nCosmicRemoved/nCosmicTracksFid<<"\n";
    myfile.close();

    fFidExitMu->BayesDivide(hFidExitCutMu, hFidTotalMu);      
    fFidExitMu->Draw("ap"); 
    fStitchAngleMu->BayesDivide(hStitchAngleCutMu, hStitchAngleTotalMu);    
    fStitchAngleMu->Draw("ap");                                    
    fStitchDistMu->BayesDivide(hStitchDistCutMu, hStitchDistTotalMu);     
    fStitchDistMu->Draw("ap");                                     
    fTrackAngleMu->BayesDivide(hTrackAngleCutMu, hTrackAngleTotalMu);     
    fTrackAngleMu->Draw("ap");                                     
    fTrackDistMu->BayesDivide(hTrackDistCutMu, hTrackDistTotalMu);      
    fTrackDistMu->Draw("ap");                                      
    fHitDistMu->BayesDivide(hHitDistCutMu, hHitDistTotalMu);        
    fHitDistMu->Draw("ap");                                        
    fHitLengthMu->BayesDivide(hHitLengthCutMu, hHitLengthTotalMu);      
    fHitLengthMu->Draw("ap");                                      

    fFidExitCos->BayesDivide(hFidExitCutCos, hFidTotalCos);       
    fFidExitCos->Draw("ap");                                       
    fStitchAngleCos->BayesDivide(hStitchAngleCutCos, hStitchAngleTotalCos);   
    fStitchAngleCos->Draw("ap");                                   
    fStitchDistCos->BayesDivide(hStitchDistCutCos, hStitchDistTotalCos);    
    fStitchDistCos->Draw("ap");                                    
    fTrackAngleCos->BayesDivide(hTrackAngleCutCos, hTrackAngleTotalCos);    
    fTrackAngleCos->Draw("ap");                                    
    fTrackDistCos->BayesDivide(hTrackDistCutCos, hTrackDistTotalCos);     
    fTrackDistCos->Draw("ap");                                     
    fHitDistCos->BayesDivide(hHitDistCutCos, hHitDistTotalCos);       
    fHitDistCos->Draw("ap");                                       
    fHitLengthCos->BayesDivide(hHitLengthCutCos, hHitLengthTotalCos);
    fHitLengthCos->Draw("ap");

  } // TPCCosmicRemoval::endJob()

  bool TPCCosmicRemoval::InFiducial(geo::Point_t point, double fiducial, double fiducialTop){
    //
    double xmin = -200 + fiducial; //-2.0 * fGeometryService->DetHalfWidth() + fiducial;
    double xmax = 200 - fiducial; //2.0 * fGeometryService->DetHalfWidth() - fiducial;
    double ymin = -fGeometryService->DetHalfHeight() + fiducial;
    double ymax = fGeometryService->DetHalfHeight() - fiducialTop;
    double zmin = 0. + fiducial;
    double zmax = fGeometryService->DetLength() - fiducial;

    double x = point.X();
    double y = point.Y();
    double z = point.Z();
    if(x>xmin && x<xmax && y>ymin && y<ymax && z>zmin && z<zmax) return true;

    return false;
  }

  int TPCCosmicRemoval::DetectedInTPC(std::vector<art::Ptr<recob::Hit>> hits){
    //
    int tpc = hits[0]->WireID().TPC;
    for(size_t i = 0; i < hits.size(); i++){
      if((int)hits[i]->WireID().TPC != tpc) return -tpc;
    }
    return tpc;
  }

  // Function to draw true and reco tracks
  void TPCCosmicRemoval::DrawTrueTracks(RecoTruth rt, bool truth, bool tpcTracks, int id){

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
    int nparts = 0;
    size_t lim = rt.particles.size();

    if(truth){
      for(size_t i = 0; i < lim; i++){
        int trueID = rt.particles[i].TrackId();
        bool plot = false;
        bool isNu = false;
        for(auto& matching : (rt.matchingMap)){
          if(matching.second.trueID == trueID){ 
            plot = true;
            isNu = matching.second.isNu;
          }
        }
        if(plot){
          int nTraj = rt.particles[i].NumberTrajectoryPoints();
          trajectories[nparts] = new TPolyLine3D(nTraj);
          int ipt = 0;
          for(int j = 0; j < nTraj; j++){
            double px = rt.particles[i].Vx(j);
            double py = rt.particles[i].Vy(j);
            double pz = rt.particles[i].Vz(j);
            if(abs(px) < 250 && py < 250 && py > -250 && pz < 550 && pz > -50){
              trajectories[nparts]->SetPoint(ipt, px, py, pz);
              ipt++;
            }
          }
          trajectories[nparts]->SetLineColor(4);
          if(isNu) trajectories[nparts]->SetLineColor(41);
          trajectories[nparts]->SetLineWidth(2);
          if(id==-99999||rt.matchingMap[id].trueID==trueID){ 
            trajectories[nparts]->Draw();
            nparts++;
          }
        }
      }
    }

    if(tpcTracks){
      // Plot the tracks
      for(size_t i = 0; i < rt.tpcTracks.size(); i++){
        // Get the start and end points
        recob::Track tr = rt.tpcTracks[i];
        size_t npts = tr.NumberTrajectoryPoints();
        tpctrack[i] = new TPolyLine3D(npts);
        for(size_t j = 0; j < npts; j++){
          auto& pos = tr.LocationAtPoint(j);
          tpctrack[i]->SetPoint(j, pos.X(), pos.Y(), pos.Z());
        }
        // Draw a line between them
        tpctrack[i]->SetLineColor(3);
        if(rt.matchingMap[tr.ID()].isRemoved) tpctrack[i]->SetLineColor(31);
        if(rt.matchingMap[tr.ID()].isNu){ 
          tpctrack[i]->SetLineColor(2);
          if(rt.matchingMap[tr.ID()].isRemoved) tpctrack[i]->SetLineColor(6);
        }
        tpctrack[i]->SetLineWidth(2);
        if(id == -99999 || tr.ID() == id) tpctrack[i]->Draw();
      }
    }
    c1->SaveAs("removalPlot.root");

  } // TPCCosmicRemoval::DrawTrueTracks()

  double TPCCosmicRemoval::T0FromCpaStitching(recob::Track t1, std::vector<recob::Track> tracks, double stitchDist, double stitchAngle){
    //
    std::vector<std::pair<double, double>> matchCandidates;
    double matchedTime = -99999;

    TVector3 trk1Front = t1.Vertex<TVector3>();
    TVector3 trk1Back = t1.End<TVector3>();
    double closestX1 = std::min(std::abs(trk1Front.X()), std::abs(trk1Back.X()));

    for(auto & track : tracks){
      bool print = false;
      //if(t1.ID()==1 && track.ID()==5) print = true;

      TVector3 trk2Front = track.Vertex<TVector3>();
      TVector3 trk2Back = track.End<TVector3>();
      double closestX2 = std::min(std::abs(trk2Front.X()), std::abs(trk2Back.X()));

      if(print) std::cout<<"Closest X1 = "<<closestX1<<" closest X2 = "<<closestX2<<"\n";

      if(std::abs(closestX1-closestX2) < 10){
        TVector3 t1Pos = trk1Front;
        TVector3 t1Dir = t1.VertexDirection<TVector3>();
        if(std::abs(trk1Back.X()) == closestX1){ 
          t1Pos = trk1Back;
          t1Dir = t1.EndDirection<TVector3>();
        }

        TVector3 t2Pos = trk2Front;
        TVector3 t2Dir = track.VertexDirection<TVector3>();
        if(std::abs(trk2Back.X()) == closestX2){ 
          t2Pos = trk2Back;
          t2Dir = track.EndDirection<TVector3>();
        }

        double trkCos = std::abs(t1Dir.Dot(t2Dir));
        t1Pos[0] = 0.;
        t2Pos[0] = 0.;
        double dist = (t1Pos-t2Pos).Mag();

        if(print) std::cout<<"t1 pos = ("<<t1Pos.X()<<", "<<t1Pos.Y()<<", "<<t1Pos.Z()<<"), t2 pos = ("<<t2Pos.X()<<", "<<t2Pos.Y()<<", "<<t2Pos.Z()<<")\n";
        if(print) std::cout<<"Cos lim = "<<cos(TMath::Pi() * 15 / 180.)<<"Cos = "<<trkCos<<" dist = "<<dist<<" stitch dist = "<<stitchDist<<"\n";

        if(dist < stitchDist && trkCos > cos(TMath::Pi() * stitchAngle / 180.)){ 
          if(print) std::cout<<"Hello\n";
          matchCandidates.push_back(std::make_pair(trkCos, closestX1));
        }
      }
    }

    if(matchCandidates.size() > 0){
      std::sort(matchCandidates.begin(), matchCandidates.end(), [](auto& left, auto& right){
                return left.first < right.first;});
      double shiftX = matchCandidates[0].second;
      matchedTime = -(shiftX/fDetectorProperties->DriftVelocity()-17.); 
    }

    return matchedTime;
  }

  double TPCCosmicRemoval::T0FromApaCross(recob::Track track, std::vector<double> t0s, int tpc){

    bool print = true;
    //if(track.ID() == 1) print = true;

    double crossTime = -99999;
    //double xmin = -2.0 * fGeometryService->DetHalfWidth();
    double xmax = 2.0 * fGeometryService->DetHalfWidth();

    if(print) std::cout<<"t0 size = "<<t0s.size()<<std::endl;

    double minDist = 99999;
    double startX = track.Vertex().X();
    double endX = track.End().X();
    geo::Point_t point = track.Vertex();
    // Don't try to shift tracks near the Apa
    if(std::abs(startX) > xmax-10. || std::abs(endX) > xmax-10.) return crossTime;
    // If in tpc 0 use start/end with lowest X
    if(tpc == 0 && endX < startX) point = track.End();
    // If in tpc 1 use start/end with highest X
    if(tpc == 1 && endX > startX) point = track.End();
    //Shift track by all t0's
    for(auto const& t0 : t0s){
      // If particle crosses the APA before t = 0 the crossing point won't be reconstructed
      if(t0 < 0) continue;
      double shiftedX = point.X();
      double shift = t0 * fDetectorProperties->DriftVelocity();
      if(print) std::cout<<"t0 = "<<t0<<" shift = "<<shift<<" tpc = "<<tpc<<std::endl;
      if(print) std::cout<<"point x = "<<shiftedX<<std::endl;
      if(tpc == 0) shiftedX = point.X() - shift;
      if(tpc == 1) shiftedX = point.X() + shift;
      if(print) std::cout<<"shifted x = "<<shiftedX<<std::endl;
      //Check track still in TPC
      if(std::abs(shiftedX) > 201.){ 
        if(print) std::cout<<"shifted outside tpc\n\n"<<std::endl;
        continue;
      }
      //Calculate distance between start/end and APA
      double dist = std::abs(std::abs(shiftedX)-199.6);
      if(print) std::cout<<"dist = "<<dist<<"\n\n";
      if(dist < minDist) {
        minDist = dist;
        crossTime = t0;
      }
    }
    if(print) std::cout<<"minDist = "<<minDist<<" crossTime = "<<crossTime<<std::endl;
    //If distance < limit take CRT time
    if(minDist < 5.){
      return crossTime;
    }
    return -99999;
      
  }

  bool TPCCosmicRemoval::IsStoppingMuon(art::Ptr<anab::ParticleID> pid, art::Ptr<anab::Calorimetry> calo){
    /*TFile *file = TFile::Open("dEdxrestemplates.root");
    TProfile *dedx_range_mu = (TProfile*)file->Get("dedx_range_mu");
    size_t nhits = calo->dEdx().size();
    int npts = 0;
    double chi2 = 0;
    //Loop over hits in track trajectory
    for(size_t i = 0; i < nhits; i++){
      // If hit is on plane 2
      double dEdx = (calo->dEdx())[i];
      double resrg = (calo->ResidualRange())[i];
      // Make sure res range is right way round
      // If residual range is < 25
      // If dedx < 30
      if(resrg > 0.3 && resrg < 20){
        int bin = dedx_range_mu->FindBin(resrg);
        double err = 0.04231 + 0.0001783*dEdx*dEdx*dEdx;
        chi2 += pow((dEdx - dedx_range_mu->GetBinContent(bin))/std::sqrt(pow(dedx_range_mu->GetBinError(bin),2)+pow(err,2)),2);
        npts++;
      }
    }
    // Calculate chi2
    chi2 = chi2/npts;*/
    size_t nhits = calo->dEdx().size();
    std::vector<double> vdedx;
    std::vector<double> vresrg;
    for(size_t i = 0; i < nhits; i++){
      // If hit is on plane 2
      double dEdx = (calo->dEdx())[i];
      double resrg = (calo->ResidualRange())[i];
      if(resrg>1 && resrg<20 && dEdx < 20){
        vdedx.push_back(dEdx);
        vresrg.push_back(resrg);
      }
    }
    if(vdedx.size()==0) return false;
    
    TGraph *gresrgdedx = new TGraph(vdedx.size(), &vresrg[0], &vdedx[0]);
    gresrgdedx->Fit("pol0", "Q");
    TF1* polfit = gresrgdedx->GetFunction("pol0");
    double chi2pol = polfit->GetChisquare();
    gresrgdedx->Fit("expo", "Q");
    TF1* expfit = gresrgdedx->GetFunction("expo");
    double chi2exp = expfit->GetChisquare();
    std::cout<<"Pol chi2 = "<<chi2pol<<" exp chi2 = "<<chi2exp<<"\n";
    if(pid->PlaneID().Plane != 2) std::cout<<"NOT PLANE 2!\n";
    double chi2 = pid->MinChi2();
    std::cout<<"Muon Chi2 = "<<pid->Chi2Muon()<<" min chi2 = "<<pid->MinChi2()<<" PDG = "<<pid->Pdg()<<"\n";
    if((chi2exp < chi2pol && chi2exp > 0 && chi2exp < 10 && std::abs(chi2pol-chi2exp) > 1) || (chi2 > 0 && chi2 < 12)){ 
      //std::cout<<"STOPPING\n";
      return true;
    }
  
    return false;
  }

  double TPCCosmicRemoval::StoppingEnd(std::vector<art::Ptr<anab::Calorimetry>> calos, geo::Point_t end, int id){
    //Loop over residual range and dedx
    size_t nhits = 0;
    if(calos.size()==0) return false;
    art::Ptr<anab::Calorimetry> calo = calos[0];
    for( size_t i = calos.size(); i > 0; i--){
      if(calos[i-1]->dEdx().size()!=0 && nhits==0){
        nhits = calos[i-1]->dEdx().size();
        calo = calos[i-1];
      }
    }
    double distStart = (calo->XYZ()[0] - end).Mag2();
    double distEnd = (calo->XYZ()[nhits-1] - end).Mag2();
    std::vector<double> v_resrg;
    std::vector<double> v_dedx;
    for(size_t i = 0; i < nhits; i++){
      double dedx = calo->dEdx()[i];
      double resrg = calo->ResidualRange()[i];
      if(distStart < distEnd && calo->ResidualRange()[0] > calo->ResidualRange()[nhits-1]) resrg = calo->ResidualRange()[0] - calo->ResidualRange()[i];
      if(distStart > distEnd && calo->ResidualRange()[0] < calo->ResidualRange()[nhits-1]) resrg = calo->ResidualRange()[nhits-1] - calo->ResidualRange()[i];
      if(resrg > 1 && resrg < 20 && dedx < 30){
        v_resrg.push_back(resrg);
        v_dedx.push_back(dedx);
      }
    }
    if(v_dedx.size() < 10) return false;
    TGraph *gdedx = new TGraph(v_dedx.size(), &v_resrg[0], &v_dedx[0]);
    try{ gdedx->Fit("pol1", "Q"); } catch(...){ return false; }
    TF1* polfit = gdedx->GetFunction("pol1");
    double slope = polfit->GetParameter(1);
    double chi2 = polfit->GetChisquare()/polfit->GetNDF();
    if(id == fTrackID){
      TCanvas *c3 = new TCanvas("c3","",700,700);
      gdedx->SetMarkerStyle(3);
      gdedx->Draw("ap");
      c3->SaveAs("gdedx.root");
    }
    if(fVerbose) std::cout<<"resrg size = "<<v_resrg.size()<<" dedx size = "<<v_dedx.size()<<" slope = "<<slope<<" chi2 = "<<chi2<<"\n";
    if(slope < -0.06 && chi2 < 5) return true;
    return false;
    //return slope;
  }

  DEFINE_ART_MODULE(TPCCosmicRemoval)
} // namespace sbnd

