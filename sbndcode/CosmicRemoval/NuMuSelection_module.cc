////////////////////////////////////////////////////////////////////////
// Class:       NuMuSelection
// Module Type: analyzer
// File:        NuMuSelection_module.cc
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbndcode/CosmicRemoval/CosmicRemovalUtils/CosmicRemovalUtils.h"
#include "sbndcode/CRT/CRTUtils/CRTAnaUtils.h"
#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"

// LArSoft includes
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
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

  class NuMuSelection : public art::EDAnalyzer {
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
    explicit NuMuSelection(Parameters const& config);
 
    // Called once, at start of the job
    virtual void beginJob() override;
 
    // Called once per event
    virtual void analyze(const art::Event& event) override;

    // Called once, at end of the job
    virtual void endJob() override;

  private:

    // fcl file parameters
    art::InputTag fSimModuleLabel;      ///< name of detsim producer
    art::InputTag fCrtHitModuleLabel; ///< name of CRT track producer
    art::InputTag fTpcTrackModuleLabel; ///< name of TPC track producer
    art::InputTag fCaloModuleLabel; ///< name of calorimetry producer
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

    // histograms

    // Performance Counters
    // Tracks from primary cosmic particles
    int nCosmicTracks = 0;
    int nCosmicTracksNotCut = 0;
    int nCosmicTracksFid = 0;
    int nCosSelected = 0;
    // Tracks from neutrino interactions
    int nNuTracks = 0;
    int nNuTracksNotCut = 0;
    int nNuTracksFid = 0;
    int nNuSelected = 0;
    int nLepTracks = 0;
    int nLepTracksNotCut = 0;
    int nLepTracksFid = 0;
    int nLepSelected = 0;
    int nOtherSelected = 0;
    int nNothingSelected = 0;

  }; // class NuMuSelection


  // Constructor
  NuMuSelection::NuMuSelection(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel       (config().SimModuleLabel())
    , fCrtHitModuleLabel    (config().CrtHitModuleLabel())
    , fTpcTrackModuleLabel  (config().TpcTrackModuleLabel())
    , fCaloModuleLabel      (config().CaloModuleLabel())
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

  } // NuMuSelection()


  void NuMuSelection::beginJob()
  {
    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;

    // Initial output
    if(fVerbose) std::cout<<"----------------- TPC Cosmic Removal Ana Module -------------------"<<std::endl;

  }// NuMuSelection::beginJob()


  void NuMuSelection::analyze(const art::Event& event)
  {
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

    // Fetch basic event info
    if(fVerbose){
      std::cout<<"============================================"<<std::endl
               <<"Run = "<<event.run()<<", SubRun = "<<event.subRun()<<", Event = "<<event.id().event()<<std::endl
               <<"============================================"<<std::endl;
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
    for(auto const& crtHit : (*crtHitHandle)){
      crtHits.push_back(crtHit);
    }
    std::vector<crt::CRTTrack> crtTracks = CRTAnaUtils::CreateCRTTracks(crtHitsPtr, 0.2, 30., true, 25.);
    if(fVerbose) std::cout<<"Number of CRTTracks = "<<crtTracks.size()<<std::endl;

    std::vector<double> stopT0 = CRTAnaUtils::ApaT0sFromCRTHits(crtHitsPtr, 2.);
    std::vector<double> throughT0 = CRTAnaUtils::ApaT0sFromCRTTracks(crtTracks);
    if(fVerbose) std::cout<<"Stopping CRT t0 size = "<<stopT0.size()<<" Through going CRT t0 size = "<<throughT0.size()<<"\n";

    //----------------------------------------------------------------------------------------------------------
    //                                          TRUTH MATCHING
    //----------------------------------------------------------------------------------------------------------
    // Retrieve all the truth info in the events
    auto particleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fSimModuleLabel);

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
        if(!CosmicRemovalUtils::InFiducial(vtx, 0, 0)) continue;
        if(std::abs(particle.PdgCode())==13 && particle.Mother()==0) lepParticleIds.push_back(partId);
        nuParticleIds.push_back(partId);

      }

    }

    if(fVerbose) std::cout<<"Number of true particles = "<<particles.size()<<"\n\n";

    //----------------------------------------------------------------------------------------------------------
    //                                          RECONSTRUCTED TRACKS
    //----------------------------------------------------------------------------------------------------------
    // Retrieve the TPC tracks
    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTpcTrackModuleLabel);

    // Get track to hit associations
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTpcTrackModuleLabel);
    art::FindManyP<anab::Calorimetry> findManyCalo(tpcTrackHandle, event, fCaloModuleLabel);

    std::vector<recob::Track> tpcTracksTPC0;
    std::vector<recob::Track> tpcTracksTPC1;

    // Loop over the tpc tracks
    for(auto const& tpcTrack : (*tpcTrackHandle)){

      std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
      int tpc = CosmicRemovalUtils::DetectedInTPC(hits);
      double startX = tpcTrack.Start().X();
      double endX = tpcTrack.End().X();

      if(tpc == 0 && !(startX>0 || endX>0)) tpcTracksTPC0.push_back(tpcTrack);
      else if(tpc == 1 && !(startX<0 || endX<0)) tpcTracksTPC1.push_back(tpcTrack);

    }

    std::vector<std::pair<recob::Track, int>> selectedTracks;

    for(auto const& tpcTrack : (*tpcTrackHandle)){

      if(fVerbose) std::cout<<"------>Track "<<tpcTrack.ID()<<":\n";

      // Match to the true particle
      std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
      int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
      if (particles.find(trueId) == particles.end()){ 
        if (fVerbose) std::cout<<"No valid true track!\n\n"; 
        continue; 
      }
      // Get the true T0
      double trueTime = particles[trueId].T()*1e-3;

      std::vector<art::Ptr<anab::Calorimetry>> calos = findManyCalo.at(tpcTrack.ID());

      if(fVerbose) std::cout<<"trueID = "<<trueId<<" pdg = "<<particles[trueId].PdgCode()<<" time = "<<trueTime<<" reco length = "<<tpcTrack.Length()
                            <<" reco start = "<<tpcTrack.Start()<<" end = "<<tpcTrack.End()<<"\n";

      //----------------------------------------------------------------------------------------------------------
      //                                          MATCH TRUE AND RECO
      //----------------------------------------------------------------------------------------------------------
      // Is track from a neutrino interaction
      bool isNu = false;
      bool isLep = false;
      bool isCos = false;
      if(std::find(lepParticleIds.begin(), lepParticleIds.end(), trueId) != lepParticleIds.end()){
        // Set neutrino flag
        isLep = true;
        nLepTracks++;
        if(fVerbose) std::cout<<"Primary lepton!\n";
      }
      if(std::find(nuParticleIds.begin(), nuParticleIds.end(), trueId) != nuParticleIds.end()){
        // Set neutrino flag
        isNu = true;
        nNuTracks++;
        if(fVerbose) std::cout<<"From neutrino!\n";
      }
      // Else is it from a primary cosmic (true length > 500)
      else if(particles[trueId].Trajectory().TotalLength() > 500. && std::abs(particles[trueId].PdgCode())==13){
        // Set primary cosmic flag
        isCos = true;
        nCosmicTracks++;
        if(fVerbose) std::cout<<"Primary cosmic!\n";
      }

      //----------------------------------------------------------------------------------------------------------
      //                                               CUTS
      //----------------------------------------------------------------------------------------------------------

      //---------------------------------  FIDUCIAL VOLUME CUT --------------------------------------
      // Remove any tracks that enter and exit the fiducial volume
      bool startInFiducial = CosmicRemovalUtils::InFiducial(tpcTrack.Vertex(), fFiducial, fFiducialTop);
      bool endInFiducial = CosmicRemovalUtils::InFiducial(tpcTrack.End(), fFiducial, fFiducialTop);
      if(!startInFiducial && !endInFiducial){ 
        continue;
      }
      if(isNu) nNuTracksFid++;
      if(isLep) nLepTracksFid++;
      if(isCos) nCosmicTracksFid++;

      //---------------------------------  STOPPING PARTICLE CUT --------------------------------------
      bool startInFiducialStop = CosmicRemovalUtils::InFiducial(tpcTrack.Vertex(), fFiducialStop, fFiducialStop);
      bool endInFiducialStop = CosmicRemovalUtils::InFiducial(tpcTrack.End(), fFiducialStop, fFiducialStop);
      bool startStops = CosmicRemovalUtils::StoppingEnd(calos, tpcTrack.Vertex(), fResRgMin, fResRgMax, fDEdxMax, fStoppingChi2Limit);
      bool endStops = CosmicRemovalUtils::StoppingEnd(calos, tpcTrack.End(), fResRgMin, fResRgMax, fDEdxMax, fStoppingChi2Limit);
      if(!startInFiducialStop && endInFiducialStop && endStops){
        continue;
      }
      if(startInFiducialStop && !endInFiducialStop && startStops){
        continue;
      }


      //---------------------------------  DIFFERENT TPC CUT --------------------------------------
      // Remove any tracks that are detected in one TPC and reconstructed in another
      int tpc = CosmicRemovalUtils::DetectedInTPC(hits);
      double startX = tpcTrack.Start().X();
      double endX = tpcTrack.End().X();
      // Check if track is stitched
      if(tpc < 0){
        // Try to get t0 from tracking algorithms, check against time
        if(fVerbose) std::cout<<"TRACK STITCHED\n";
      }
      // If it is check the start/end points are in same TPC
      else if(tpc == 0 && (startX>0 || endX>0)){ 
        continue;
      }
      else if(tpc == 1 && (startX<0 || endX<0)){
        continue;
      }

      //--------------------------------- CPA STITCHING CUT --------------------------------------
      double stitchTime = -99999;
      bool stitchExit = false;
      // Try to match tracks from CPA crossers
      if(tpc == 0){
        std::pair<double, bool> stitchResults = CosmicRemovalUtils::T0FromCpaStitching(tpcTrack, tpcTracksTPC1, fCpaStitchDistance, 
                                                                                       fCpaStitchAngle, fCpaXDifference, fFiducial, fFiducialTop);
        stitchTime = stitchResults.first;
        stitchExit = stitchResults.second;
      }
      else if(tpc == 1){
        std::pair<double, bool> stitchResults = CosmicRemovalUtils::T0FromCpaStitching(tpcTrack, tpcTracksTPC0, fCpaStitchDistance, 
                                                                                       fCpaStitchAngle, fCpaXDifference, fFiducial, fFiducialTop);
        stitchTime = stitchResults.first;
        stitchExit = stitchResults.second;
      }
      // If tracks are stitched, get time and remove any outside of beam window
      if(stitchTime != -99999 && std::abs(stitchTime) > fBeamTimeLimit){
        continue;
      }
      // If time within beam window, check if tracks exit the TPC
      if(stitchExit){
        continue;
      }

      //--------------------------------- CRTTRACK MATCHING --------------------------------------
      // Try to get T0 from CRTTracks, if there's a match then remove the track
      double trackTime = CRTAnaUtils::T0FromCRTTracks(tpcTrack, crtTracks, tpc, fMaxAngleDiff, fMaxDistance);
      if(trackTime != -99999){
        continue;
      }

      //---------------------------------- CRTHIT MATCHING ---------------------------------------
      // Try to get T0 from CRTHits, remove any matched outside the beam
      double hitTime = CRTAnaUtils::T0FromCRTHits(tpcTrack, crtHits, tpc, fMinTrackLength, fTrackDirectionFrac, fDistanceLimit);
      if(hitTime != -99999 && std::abs(hitTime) > fBeamTimeLimit){
        continue;
      }


      //----------------- APA CROSS MATCHING FOR THROUGH GOING PARTICLES -------------------------
      // Match APA crossers with times from CRT tracks that cross the APA
      double crossTimeThrough = CosmicRemovalUtils::T0FromApaCross(tpcTrack, throughT0, tpc, fFiducial, 2.);
      if(crossTimeThrough != -99999 && std::abs(crossTimeThrough) > fBeamTimeLimit){ 
        // Check that the end that doesn't cross the APA exits the TPC
        if(tpc == 0){
          if(tpcTrack.Vertex().X() < tpcTrack.End().X() && !endInFiducial) continue;
          else if(!startInFiducial) continue;
        }
        else if(tpc == 1){
          if(tpcTrack.Vertex().X() > tpcTrack.End().X() && !endInFiducial) continue;
          else if(!startInFiducial) continue;
        }
      }

      //-------------------- APA CROSS MATCHING FOR STOPPING PARTICLES ----------------------------
      // Match APA crossers with times from CRT hits that could result in a stopping particle that crosses the APA
      double crossTimeStop = CosmicRemovalUtils::T0FromApaCross(tpcTrack, stopT0, tpc, fFiducial, 2.);
      if(crossTimeStop != -99999 && std::abs(crossTimeStop) > fBeamTimeLimit){
        // Check that the end that doesn't cross the APA stops in the TPC
        if(tpc == 0){
          if(tpcTrack.Vertex().X() < tpcTrack.End().X() && endStops) continue;
          else if(startStops) continue;
        }
        else if(tpc == 1){
          if(tpcTrack.Vertex().X() > tpcTrack.End().X() && endStops) continue;
          else if(startStops) continue;
        }
      }

      if(fVerbose) std::cout<<"Not cut!\n";

      if(isNu) nNuTracksNotCut++;
      if(isLep) nLepTracksNotCut++;
      if(isCos) nCosmicTracksNotCut++;
        
      int type = 0;
      if(isLep) type = 1;
      else if(isNu) type = 2;
      else if(isCos) type = 3;
      selectedTracks.push_back(std::make_pair(tpcTrack, type));

    }

    double maxLength = 0;
    double maxType = -1;
    for(auto const& track : selectedTracks){
      if(track.first.Length() > maxLength){
        maxLength = track.first.Length();
        maxType = track.second;
      }
    }

    if(maxType == -1) nNothingSelected++;
    if(maxType == 0) nOtherSelected++;
    if(maxType == 1) nLepSelected++;
    if(maxType == 2) nNuSelected++;
    if(maxType == 3) nCosSelected++;

  } // NuMuSelection::analyze()


  void NuMuSelection::endJob(){

    std::cout<<"Total nu tracks             = "<<nNuTracks<<"\n"
             <<"Tracks in fiducial volume   = "<<nNuTracksFid<<"\n"
             <<"Not removed                 = "<<nNuTracksNotCut<<"\n"
             <<"Percentage removed by vol   = "<<(double)(nNuTracks-nNuTracksFid)/nNuTracks<<"\n"
             <<"Percentage removed by CosID = "<<(double)(nNuTracksFid-nNuTracksNotCut)/nNuTracksFid<<"\n"
             <<"-----------------------------------------------------------\n"
             <<"Total lepton tracks         = "<<nLepTracks<<"\n"
             <<"Tracks in fiducial volume   = "<<nLepTracksFid<<"\n"
             <<"Not removed                 = "<<nLepTracksNotCut<<"\n"
             <<"Percentage removed by vol   = "<<(double)(nLepTracks-nLepTracksFid)/nLepTracks<<"\n"
             <<"Percentage removed by CosID = "<<(double)(nLepTracksFid-nLepTracksNotCut)/nLepTracksFid<<"\n"
             <<"-----------------------------------------------------------\n"
             <<"Total cosmic tracks         = "<<nCosmicTracks<<"\n"
             <<"Tracks in fiducial volume   = "<<nCosmicTracksFid<<"\n"
             <<"Not removed                 = "<<nCosmicTracksNotCut<<"\n"
             <<"Percentage removed by vol   = "<<(double)(nCosmicTracks-nCosmicTracksFid)/nCosmicTracks<<"\n"
             <<"Percentage removed by CosID = "<<(double)(nCosmicTracksFid-nCosmicTracksNotCut)/nCosmicTracksFid<<"\n"
             <<"-----------------------------------------------------------\n"
             <<"Selected:\n"
             <<"Nu Muon      = "<<nLepSelected<<"\n"
             <<"Nu Particle  = "<<nNuSelected<<"\n"
             <<"Cosmic       = "<<nCosSelected<<"\n"
             <<"Other        = "<<nOtherSelected<<"\n"
             <<"Nothing      = "<<nNothingSelected<<"\n";

  } // NuMuSelection::endJob()


  DEFINE_ART_MODULE(NuMuSelection)
} // namespace sbnd

