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
#include "TH2.h"
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

      fhicl::Atom<bool> Verbose {
        Name("Verbose"),
        Comment("Print information about what's going on")
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

    double T0FromApaCross(recob::Track track, std::vector<crt::CRTHit> crtHits, int tpc);
    
    // Function to draw true and reco tracks
    void DrawTrueTracks(RecoTruth truthMatch, bool truth, bool tpctracks, int id);

  private:

    // fcl file parameters
    art::InputTag fSimModuleLabel;      ///< name of detsim producer
    art::InputTag fCrtHitModuleLabel; ///< name of CRT track producer
    art::InputTag fTpcTrackModuleLabel; ///< name of TPC track producer
    bool          fVerbose;             ///< print information about what's going on
    
    // histograms
    // True time of neutrino particles to get beam window
    TH1D* hTrueNuMuTime;
    // Contained length of primary muon within fiducial volume
    TH1D* hTrueNuMuLength;

    TGraphAsymmErrors* fFidNoEnterMu;
    TH1D* hFidTotalMu;
    TH1D* hFidNoEnterCutMu;

    TGraphAsymmErrors* fFidNoEnterCos;
    TH1D* hFidTotalCos;
    TH1D* hFidNoEnterCutCos;

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

    TGraphAsymmErrors* fHitDirMu;
    TH1D* hHitDirTotalMu;
    TH1D* hHitDirCutMu;

    TGraphAsymmErrors* fHitDirCos;
    TH1D* hHitDirTotalCos;
    TH1D* hHitDirCutCos;

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
    int nCosmicRemovedNoEnter = 0;
    int nCosmicRemovedExit = 0;
    int nCosmicRemovedDiffTPC = 0;
    int nCosmicRemovedStitch = 0;
    int nCosmicRemovedTrack = 0;
    int nCosmicRemovedHit = 0;
    // Tracks from neutrino interactions
    int nNuTracks = 0;
    int nNuRemoved = 0;
    int nNuRemovedNoEnter = 0;
    int nNuRemovedExit = 0;
    int nNuRemovedDiffTPC = 0;
    int nNuRemovedStitch = 0;
    int nNuRemovedTrack = 0;
    int nNuRemovedHit = 0;

    int nLepTracks = 0;
    int nLepRemoved = 0;
    int nLepRemovedNoEnter = 0;
    int nLepRemovedExit = 0;
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
    double hdirStart = 0.025;
    double hdirStep = 0.05;
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
    , fVerbose              (config().Verbose())
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

    fFidNoEnterMu    = tfs->makeAndRegister<TGraphAsymmErrors>("FidNoEnterMu",   ";Fiducial volume cut (cm);Percentage removed");
    fFidExitMu       = tfs->makeAndRegister<TGraphAsymmErrors>("FidExitMu",      ";Fiducial volume cut (cm);Percentage removed");
    fStitchAngleMu   = tfs->makeAndRegister<TGraphAsymmErrors>("StitchAngleMu",  ";Stitch angle (deg);Percentage removed"      );
    fStitchDistMu    = tfs->makeAndRegister<TGraphAsymmErrors>("StitchDistMu",   ";Stitch dist (cm);Percentage removed"        );
    fTrackAngleMu    = tfs->makeAndRegister<TGraphAsymmErrors>("TrackAngleMu",   ";Max angle (deg);Percentage removed"         );
    fTrackDistMu     = tfs->makeAndRegister<TGraphAsymmErrors>("TrackDistMu",    ";Max dist (cm);Percentage removed"           );
    fHitDirMu        = tfs->makeAndRegister<TGraphAsymmErrors>("HitDirMu",       ";Direction average (%);Percentage removed"   );
    fHitDistMu       = tfs->makeAndRegister<TGraphAsymmErrors>("HitDistMu",      ";Max dist (cm);Percentage removed"           );
    fHitLengthMu     = tfs->makeAndRegister<TGraphAsymmErrors>("HitLengthMu",    ";Min length (cm);Percentage removed"         );

    fFidNoEnterCos    = tfs->makeAndRegister<TGraphAsymmErrors>("FidNoEnterCos",   ";Fiducial volume cut (cm);Percentage removed");
    fFidExitCos       = tfs->makeAndRegister<TGraphAsymmErrors>("FidExitCos",      ";Fiducial volume cut (cm);Percentage removed");
    fStitchAngleCos   = tfs->makeAndRegister<TGraphAsymmErrors>("StitchAngleCos",  ";Stitch angle (deg);Percentage removed"      );
    fStitchDistCos    = tfs->makeAndRegister<TGraphAsymmErrors>("StitchDistCos",   ";Stitch dist (cm);Percentage removed"        );
    fTrackAngleCos    = tfs->makeAndRegister<TGraphAsymmErrors>("TrackAngleCos",   ";Max angle (deg);Percentage removed"         );
    fTrackDistCos     = tfs->makeAndRegister<TGraphAsymmErrors>("TrackDistCos",    ";Max dist (cm);Percentage removed"           );
    fHitDirCos        = tfs->makeAndRegister<TGraphAsymmErrors>("HitDirCos",       ";Direction average (%);Percentage removed"   );
    fHitDistCos       = tfs->makeAndRegister<TGraphAsymmErrors>("HitDistCos",      ";Max dist (cm);Percentage removed"           );
    fHitLengthCos     = tfs->makeAndRegister<TGraphAsymmErrors>("HitLengthCos",    ";Min length (cm);Percentage removed"         );

    hFidTotalMu         = tfs->make<TH1D>("FidTotalMu",         "", 10, 0, 50);
    hStitchAngleTotalMu = tfs->make<TH1D>("StitchAngleTotalMu", "", 10, 0, 40);
    hStitchDistTotalMu  = tfs->make<TH1D>("StitchDistTotalMu",  "", 10, 0, 100);
    hTrackAngleTotalMu  = tfs->make<TH1D>("TrackAngleTotalMu",  "", 10, 0, 40);
    hTrackDistTotalMu   = tfs->make<TH1D>("TrackDistTotalMu",   "", 10, 0, 160);
    hHitDirTotalMu      = tfs->make<TH1D>("HitDirTotalMu",      "", 10, 0, 0.5);
    hHitDistTotalMu     = tfs->make<TH1D>("HitDistTotalMu",     "", 10, 0, 100);
    hHitLengthTotalMu   = tfs->make<TH1D>("HitLengthTotalMu",   "", 10, 0, 100);

    hFidTotalCos         = tfs->make<TH1D>("FidTotalCos",         "", 10, 0, 50);
    hStitchAngleTotalCos = tfs->make<TH1D>("StitchAngleTotalCos", "", 10, 0, 40);
    hStitchDistTotalCos  = tfs->make<TH1D>("StitchDistTotalCos",  "", 10, 0, 100);
    hTrackAngleTotalCos  = tfs->make<TH1D>("TrackAngleTotalCos",  "", 10, 0, 40);
    hTrackDistTotalCos   = tfs->make<TH1D>("TrackDistTotalCos",   "", 10, 0, 160);
    hHitDirTotalCos      = tfs->make<TH1D>("HitDirTotalCos",      "", 10, 0, 0.5);
    hHitDistTotalCos     = tfs->make<TH1D>("HitDistTotalCos",     "", 10, 0, 100);
    hHitLengthTotalCos   = tfs->make<TH1D>("HitLengthTotalCos",   "", 10, 0, 100);

    hFidNoEnterCutMu  = tfs->make<TH1D>("FidNoEnterCutMu",  "", 10, 0, 50);
    hFidExitCutMu     = tfs->make<TH1D>("FidExitCutMu",     "", 10, 0, 50);
    hStitchAngleCutMu = tfs->make<TH1D>("StitchAngleCutMu", "", 10, 0, 40);
    hStitchDistCutMu  = tfs->make<TH1D>("StitchDistCutMu",  "", 10, 0, 100);
    hTrackAngleCutMu  = tfs->make<TH1D>("TrackAngleCutMu",  "", 10, 0, 40);
    hTrackDistCutMu   = tfs->make<TH1D>("TrackDistCutMu",   "", 10, 0, 160);
    hHitDirCutMu      = tfs->make<TH1D>("HitDirCutMu",      "", 10, 0, 0.5);
    hHitDistCutMu     = tfs->make<TH1D>("HitDistCutMu",     "", 10, 0, 100);
    hHitLengthCutMu   = tfs->make<TH1D>("HitLengthCutMu",   "", 10, 0, 100);

    hFidNoEnterCutCos  = tfs->make<TH1D>("FidNoEnterCutCos",  "", 10, 0, 50);
    hFidExitCutCos     = tfs->make<TH1D>("FidExitCutCos",     "", 10, 0, 50);
    hStitchAngleCutCos = tfs->make<TH1D>("StitchAngleCutCos", "", 10, 0, 40);
    hStitchDistCutCos  = tfs->make<TH1D>("StitchDistCutCos",  "", 10, 0, 100);
    hTrackAngleCutCos  = tfs->make<TH1D>("TrackAngleCutCos",  "", 10, 0, 40);
    hTrackDistCutCos   = tfs->make<TH1D>("TrackDistCutCos",   "", 10, 0, 160);
    hHitDirCutCos      = tfs->make<TH1D>("HitDirCutCos",      "", 10, 0, 0.5);
    hHitDistCutCos     = tfs->make<TH1D>("HitDistCutCos",     "", 10, 0, 100);
    hHitLengthCutCos   = tfs->make<TH1D>("HitLengthCutCos",   "", 10, 0, 100);

    // Initial output
    if(fVerbose) std::cout<<"----------------- TPC Cosmic Removal Ana Module -------------------"<<std::endl;

  } // TPCCosmicRemoval::beginJob()


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

    if(fVerbose) std::cout<<"Number of CRT hits = "<<crtHitHandle->size()<<std::endl
                          <<"Number of TPC tracks = "<<tpcTrackHandle->size()<<std::endl;

    // Do track reconstruction from the hits
    std::vector<crt::CRTHit> crtHits;
    for(auto const& crtHit : (*crtHitHandle)){
      crtHits.push_back(crtHit);
    }
    std::vector<crt::CRTTrack> crtTracks = CRTAnaUtils::CreateCRTTracks(crtHitsPtr, 0.2, 30., true, 25.); 
    if(fVerbose) std::cout<<"Number of CRTTracks = "<<crtTracks.size()<<std::endl;

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

      if(fVerbose) std::cout<<"trueID = "<<trueId<<" pdg = "<<particles[trueId].PdgCode()<<" time = "<<trueTime<<" reco length = "<<tpcTrack.Length()
                            <<" reco start = "<<tpcTrack.Start()<<" end = "<<tpcTrack.End()<<"\n";

      // Is track from a neutrino interaction
      bool isNu = false;
      bool isLep = false;
      bool isCos = false;
      bool removed = false;
      if(std::find(nuParticleIds.begin(), nuParticleIds.end(), trueId) != nuParticleIds.end()){
        // Set neutrino flag
        isNu = true;
        nNuTracks++;
        if(fVerbose) std::cout<<"From neutrino!\n";
      }
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
          hHitDirTotalMu->Fill(hdirStart + i*hdirStep);     
          hHitDistTotalMu->Fill(hdStart + i*hdStep);    
          hHitLengthTotalMu->Fill(hlStart + i*hlStep);  
        }
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
          hHitDirTotalCos->Fill(hdirStart + i*hdirStep);     
          hHitDistTotalCos->Fill(hdStart + i*hdStep);    
          hHitLengthTotalCos->Fill(hlStart + i*hlStep);  
        }
      }


      // FIDUCIAL VOLUME CUT

      // Remove any tracks that don't enter the fiducial volume
      if(!removed && !EntersFiducial(tpcTrack, 10, 10)){ 
        removed = true;
        if(fVerbose) std::cout<<"REMOVED! - not in fiducial\n";
        if(isNu) nNuRemoved++;
        if(isLep) nLepRemoved++;
        // If primary cosmic flag
        if(isCos) nCosmicRemoved++;
      }
      if(!EntersFiducial(tpcTrack, 10, 10)){ 
        if(isNu) nNuRemovedNoEnter++;
        if(isLep) nLepRemovedNoEnter++;
        if(isCos) nCosmicRemovedNoEnter++;
      }

      // Remove any tracks that enter and exit the fiducial volume
      if(!removed && !InFiducial(tpcTrack.Start(), 10, 10) && !InFiducial(tpcTrack.End(), 10, 10)){ 
        removed = true;
        if(fVerbose) std::cout<<"REMOVED! - exits fiducial\n";
        if(isNu) nNuRemoved++;
        if(isLep) nLepRemoved++;
        // If primary cosmic flag
        if(isCos) nCosmicRemoved++;
      }
      if(!InFiducial(tpcTrack.Start(), 10, 10) && !InFiducial(tpcTrack.End(), 10, 10)){ 
        if(isNu) nNuRemovedExit++;
        if(isLep) nLepRemovedExit++;
        if(isCos) nCosmicRemovedExit++;
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
        stitchTime = T0FromCpaStitching(tpcTrack, tpcTracksTPC1, 15., 50.);
      }
      else if(tpc == 1){
        stitchTime = T0FromCpaStitching(tpcTrack, tpcTracksTPC0, 15., 50.);
      }
      if(stitchTime != -99999 && std::abs(stitchTime)>timeLim){
        if(isNu) nNuRemovedStitch++;
        if(isLep) nLepRemovedStitch++;
        if(isCos) nCosmicRemovedStitch++;
      }

      //double crossTime = T0FromApaCross(tpcTrack, crtHits, tpc);
      /*std::cout<<"------>Track "<<tpcTrack.ID()<<":\n";
      std::cout<<"trueID = "<<trueId<<" pdg = "<<particles[trueId].PdgCode()<<" time = "<<trueTime<<" reco length = "<<tpcTrack.Length()
               <<" reco start = "<<tpcTrack.Start()<<" end = "<<tpcTrack.End()<<"\n";
      std::cout<<"Cross time = "<<crossTime<<" Stitch time = "<<stitchTime<<"\n\n";*/

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
                            <<" ticks, hit time = "<<hitTime<<"\n";

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
        double hdir = hdirStart + i*hdirStep;
        double hd = hdStart + i*hdStep;
        double hl = hlStart + i*hlStep;
        if(!EntersFiducial(tpcTrack, fidCut, fidCut)){
          if(isLep) hFidNoEnterCutMu->Fill(fidCut);
          if(isCos) hFidNoEnterCutCos->Fill(fidCut);
        }
        if(!InFiducial(tpcTrack.Start(), fidCut, fidCut) && !InFiducial(tpcTrack.End(), fidCut, fidCut)){
          if(isLep) hFidExitCutMu->Fill(fidCut);
          if(isCos) hFidExitCutCos->Fill(fidCut);
        }

        double time = -99999;
        if(tpc == 0){
          time = T0FromCpaStitching(tpcTrack, tpcTracksTPC1, sa, 50.);
        }
        else if(tpc == 1){
          time = T0FromCpaStitching(tpcTrack, tpcTracksTPC0, sa, 50.);
        }
        if(time != -99999 && std::abs(time)>timeLim){
          if(isLep) hStitchAngleCutMu->Fill(sa);
          if(isCos) hStitchAngleCutCos->Fill(sa);
        }

        time = -99999;
        if(tpc == 0){
          time = T0FromCpaStitching(tpcTrack, tpcTracksTPC1, 15., sd);
        }
        else if(tpc == 1){
          time = T0FromCpaStitching(tpcTrack, tpcTracksTPC0, 15., sd);
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

        time = CRTAnaUtils::T0FromCRTHits(tpcTrack, crtHits, tpc, 20., hdir, 30.);
        if(time != -99999 && std::abs(time) > timeLim){
          if(isLep) hHitDirCutMu->Fill(hdir);
          if(isCos) hHitDirCutCos->Fill(hdir);
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

    DrawTrueTracks(truthMatch, true, true, -99999);


  } // TPCCosmicRemoval::analyze()


  void TPCCosmicRemoval::endJob(){

    std::cout<<"Total nu tracks           = "<<nNuTracks<<"\n"
             <<"Removed in no enter cut   = "<<nNuRemovedNoEnter<<"\n"
             <<"Removed in exit cut       = "<<nNuRemovedExit<<"\n"
             <<"Removed in diff TPC cut   = "<<nNuRemovedDiffTPC<<"\n"
             <<"Removed in stitch cut     = "<<nNuRemovedStitch<<"\n"
             <<"Removed in CRT track cut  = "<<nNuRemovedTrack<<"\n"
             <<"Removed in CRT hit cut    = "<<nNuRemovedHit<<"\n"
             <<"Total removed             = "<<nNuRemoved<<"\n"
             <<"Percentage removed        = "<<(double)nNuRemoved/nNuTracks<<"\n"
             <<"-----------------------------------------------------------\n"
             <<"Total lepton tracks       = "<<nLepTracks<<"\n"
             <<"Removed in no enter cut   = "<<nLepRemovedNoEnter<<"\n"
             <<"Removed in exit cut       = "<<nLepRemovedExit<<"\n"
             <<"Removed in diff TPC cut   = "<<nLepRemovedDiffTPC<<"\n"
             <<"Removed in stitch cut     = "<<nLepRemovedStitch<<"\n"
             <<"Removed in CRT track cut  = "<<nLepRemovedTrack<<"\n"
             <<"Removed in CRT hit cut    = "<<nLepRemovedHit<<"\n"
             <<"Total removed             = "<<nLepRemoved<<"\n"
             <<"Percentage removed        = "<<(double)nLepRemoved/nLepTracks<<"\n"
             <<"-----------------------------------------------------------\n"
             <<"Total cosmic tracks       = "<<nCosmicTracks<<"\n"
             <<"Removed in no enter cut   = "<<nCosmicRemovedNoEnter<<"\n"
             <<"Removed in exit cut       = "<<nCosmicRemovedExit<<"\n"
             <<"Removed in diff TPC cut   = "<<nCosmicRemovedDiffTPC<<"\n"
             <<"Removed in stitch cut     = "<<nCosmicRemovedStitch<<"\n"
             <<"Removed in CRT track cut  = "<<nCosmicRemovedTrack<<"\n"
             <<"Removed in CRT hit cut    = "<<nCosmicRemovedHit<<"\n"
             <<"Total removed             = "<<nCosmicRemoved<<"\n"
             <<"Percentage removed        = "<<(double)nCosmicRemoved/nCosmicTracks<<"\n";

    std::ofstream myfile;
    myfile.open("results.txt");
  myfile<<nNuTracks<<","<<nNuRemovedNoEnter<<","<<nNuRemovedExit<<","<<nNuRemovedDiffTPC<<","<<nNuRemovedStitch<<","<<nNuRemovedTrack<<","<<nNuRemovedHit<<","<<nNuRemoved<<","<<(double)nNuRemoved/nNuTracks<<","<<nLepTracks<<","<<nLepRemovedNoEnter<<","<<nLepRemovedExit<<","<<nLepRemovedDiffTPC<<","<<nLepRemovedStitch<<","<<nLepRemovedTrack<<","<<nLepRemovedHit<<","<<nLepRemoved<<","<<(double)nLepRemoved/nLepTracks<<","<<nCosmicTracks<<","<<nCosmicRemovedNoEnter<<","<<nCosmicRemovedExit<<","<<nCosmicRemovedDiffTPC<<","<<nCosmicRemovedStitch<<","<<nCosmicRemovedTrack<<","<<nCosmicRemovedHit<<","<<nCosmicRemoved<<","<<(double)nCosmicRemoved/nCosmicTracks<<"\n";
    myfile.close();

    fFidNoEnterMu->BayesDivide(hFidNoEnterCutMu, hFidTotalMu);
    fFidNoEnterMu->Draw("ap");
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
    fHitDirMu->BayesDivide(hHitDirCutMu, hHitDirTotalMu);         
    fHitDirMu->Draw("ap");                                         
    fHitDistMu->BayesDivide(hHitDistCutMu, hHitDistTotalMu);        
    fHitDistMu->Draw("ap");                                        
    fHitLengthMu->BayesDivide(hHitLengthCutMu, hHitLengthTotalMu);      
    fHitLengthMu->Draw("ap");                                      

    fFidNoEnterCos->BayesDivide(hFidNoEnterCutCos, hFidTotalCos);    
    fFidNoEnterCos->Draw("ap");                                    
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
    fHitDirCos->BayesDivide(hHitDirCutCos, hHitDirTotalCos);        
    fHitDirCos->Draw("ap");                                        
    fHitDistCos->BayesDivide(hHitDistCutCos, hHitDistTotalCos);       
    fHitDistCos->Draw("ap");                                       
    fHitLengthCos->BayesDivide(hHitLengthCutCos, hHitLengthTotalCos);
    fHitLengthCos->Draw("ap");

  } // TPCCosmicRemoval::endJob()

  bool TPCCosmicRemoval::EntersFiducial(recob::Track track, double fiducial, double fiducialTop){
    //
    size_t npts = track.NumberTrajectoryPoints();
    for(size_t i = 0; i < npts; i++){
      geo::Point_t point = track.LocationAtPoint(i);
      if(InFiducial(point, fiducial, fiducialTop)) return true;
    }

    return false;
  }

  bool TPCCosmicRemoval::InFiducial(geo::Point_t point, double fiducial, double fiducialTop){
    //
    double xmin = -2.0 * fGeometryService->DetHalfWidth() + fiducial;
    double xmax = 2.0 * fGeometryService->DetHalfWidth() - fiducial;
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
    TCanvas *c1 = new TCanvas("c1","",700,700);
    
    double xmin = -2.0 * fGeometryService->DetHalfWidth();
    double xmax = 2.0 * fGeometryService->DetHalfWidth();
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
      //if(t1.ID()==2 && track.ID()==10) print = true;

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
        if(print) std::cout<<"Cos lim = "<<cos(TMath::Pi() * 15 / 180.)<<"Cos = "<<trkCos<<" dist = "<<dist<<"\n";

        if(dist < stitchDist && trkCos > cos(TMath::Pi() * stitchAngle / 180.)){ 
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

  double TPCCosmicRemoval::T0FromApaCross(recob::Track track, std::vector<crt::CRTHit> crtHits, int tpc){

    bool print = false;
    //if(track.ID() == 15) print = true;

    double crossTime = -99999;
    double xmin = -2.0 * fGeometryService->DetHalfWidth();
    double xmax = 2.0 * fGeometryService->DetHalfWidth();

    std::vector<double> crtT0;
    //Get all unique t0's from CRT hits
    for(auto const& crtHit : crtHits){
      crtT0.push_back((double)(int)crtHit.ts1_ns*1e-4); //FIXME
    }
    crtT0.erase(std::unique(crtT0.begin(), crtT0.end()), crtT0.end());
    if(print) std::cout<<"t0 size = "<<crtT0.size()<<std::endl;

    double minDist = 99999;
    //Shift track by all t0's
    for(auto const& t0 : crtT0){
      geo::Point_t start = track.Vertex();
      geo::Point_t end = track.End();
      double shift = t0 * fDetectorProperties->DriftVelocity();
      if(print && t0>750 && t0 < 800) std::cout<<"t0 = "<<t0<<" shift = "<<shift<<" tpc = "<<tpc<<std::endl;
      if(print && t0>750 && t0 < 800) std::cout<<"start = "<<start<<std::endl;
      if(tpc == 0){
        start.SetX(start.X() - shift);
        end.SetX(end.X() - shift);
      }
      if(tpc == 1){
        start.SetX(start.X() + shift);
        end.SetX(end.X() + shift);
      }
      if(print && t0>750 && t0 < 800) std::cout<<"start = "<<start<<std::endl;
      //Check track still in TPC
      if(!InFiducial(start, -0.5, -0.5) || !InFiducial(end, -0.5, -0.5)) continue;
      if(start.X()*end.X() < 0) continue;
      //Calculate distance between start/end and APA
      double dist = std::abs(std::max(start.X(), end.X())-xmax);
      if(tpc == 0) dist = std::abs(std::min(start.X(), end.X())-xmin);
      if(dist < minDist) {
        minDist = dist;
        crossTime = t0;
      }
    }
    if(print) std::cout<<"minDist = "<<minDist<<std::endl;
    //If distance < limit take CRT time
    if(minDist < 5.){
      std::cout<<"Cross time = "<<crossTime<<" min dist = "<<minDist<<"\n";
      return crossTime;
    }
    return -99999;
      
  }

  DEFINE_ART_MODULE(TPCCosmicRemoval)
} // namespace sbnd

