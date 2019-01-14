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

    double T0FromCpaStitching(recob::Track track, std::vector<recob::Track> tracks);
    
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

    // Other variables shared between different methods.
    geo::GeometryCore const* fGeometryService;                 ///< pointer to Geometry provider
    detinfo::DetectorProperties const* fDetectorProperties;    ///< pointer to detector properties provider
    detinfo::DetectorClocks const* fDetectorClocks;            ///< pointer to detector clocks provider

    CRTTruthRecoAlg truthAlg;

    // Performance Counters
    // Tracks from primary cosmic particles
    int nCosmicTracks = 0;
    int nCosmicRemoved = 0;
    // Tracks from neutrino interactions
    int nNuTracks = 0;
    int nNuRemoved = 0;
    int nLepTracks = 0;
    int nLepRemoved = 0;

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

    // Loop over the tpc tracks
    for(auto const& tpcTrack : (*tpcTrackHandle)){
      tpcTracks.push_back(tpcTrack);
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
      }
      // Else is it from a primary cosmic (true length > 500)
      else if(particles[trueId].Trajectory().TotalLength() > 500.){
        // Set primary cosmic flag
        isCos = true;
        nCosmicTracks++;
        if(fVerbose) std::cout<<"Primary cosmic!\n";
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

      // Remove any tracks that enter and exit the fiducial volume
      if(!removed && !InFiducial(tpcTrack.Start(), 10, 10) && !InFiducial(tpcTrack.End(), 10, 10)){ 
        removed = true;
        if(fVerbose) std::cout<<"REMOVED! - exits fiducial\n";
        if(isNu) nNuRemoved++;
        if(isLep) nLepRemoved++;
        // If primary cosmic flag
        if(isCos) nCosmicRemoved++;
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
      }
      else std::cout<<"SOMETHING WRONG\n";

      // TIME CUTS
      // Match CPA crossers, remove any outside of beam window
      // Do CRT time matching
 
      // Try to get T0 from CRTTracks
      double trackTime = CRTAnaUtils::T0FromCRTTracks(tpcTrack, crtTracks, tpc, 0.2, 80.);

      // Try to get T0 from CRTHits
      double hitTime = CRTAnaUtils::T0FromCRTHits(tpcTrack, crtHits, tpc, 20., 0.5, 15.);

      if(fVerbose) std::cout<<"True time = "<<trueTime<<" ticks, track time = "<<trackTime
                            <<" ticks, hit time = "<<hitTime<<"\n";

      double matchedTime = -99999;
      if(trackTime != -99999) matchedTime = trackTime;
      else if(hitTime != -99999) matchedTime = hitTime;

      if(fVerbose) std::cout<<"Matched time = "<<matchedTime<<"\n";

      double timeLim = 20.;
      if(!removed && matchedTime != -99999 && (std::abs(matchedTime) > timeLim)){
        removed = true;
        if(fVerbose) std::cout<<"REMOVED! - outside beam\n";
        if(isNu) nNuRemoved++;
        if(isLep) nLepRemoved++;
        // If primary cosmic flag
        if(isCos) nCosmicRemoved++;
      }

      if(fVerbose) std::cout<<"\n";

      trackMatch.isNu = isNu;
      trackMatch.isLep = isLep;
      trackMatch.isRemoved = removed;
      matchingMap[tpcTrack.ID()] = trackMatch;

    }

    RecoTruth truthMatch;
    truthMatch.particles = parts;
    truthMatch.tpcTracks = tpcTracks;
    truthMatch.matchingMap = matchingMap;

    //DrawTrueTracks(truthMatch, true, true, -99999);


  } // TPCCosmicRemoval::analyze()


  void TPCCosmicRemoval::endJob(){

    std::cout<<"Total nu tracks = "<<nNuTracks<<"\n"
             <<"Percentage removed = "<<(double)nNuRemoved/nNuTracks<<"\n"
             <<"Total lepton tracks = "<<nLepTracks<<"\n"
             <<"Percentage removed = "<<(double)nLepRemoved/nLepTracks<<"\n"
             <<"Total cosmic tracks = "<<nCosmicTracks<<"\n"
             <<"Percentage removed = "<<(double)nCosmicRemoved/nCosmicTracks<<"\n";

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

  } // CRTTrackMatchingAna::DrawTrueTracks()

  double T0FromCpaStitching(recob::Track t1, std::vector<recob::Track> tracks, int tpc){
    //
    geo::Point_t trk1Front = t1.Start();
    geo::Point_t trk1Back = t1.End();
    double closestX1 = std::min(std::abs(trk1Front.X(), trk1Back.X()));
    for(auto & track : tracks){
      geo::Point_t trk2Front = track.Start();
      geo::Point_t trk2Back = track.End();
      double closestX2 = std::min(std::abs(trk2Front.X(), trk2Back.X());
      if(std::abs(closestX1-closestX2) < 
    

  }

  DEFINE_ART_MODULE(TPCCosmicRemoval)
} // namespace sbnd

