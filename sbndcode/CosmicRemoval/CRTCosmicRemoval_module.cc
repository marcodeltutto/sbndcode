////////////////////////////////////////////////////////////////////////
// Class:       CRTCosmicRemoval
// Module Type: analyzer
// File:        CRTCosmicRemovalAna_module.cc
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbndcode/CRT/CRTUtils/CRTAnaUtils.h"
#include "sbndcode/CRT/CRTUtils/CRTT0MatchAlg.h"
#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/CRTTruthRecoAlg.h"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"

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
#include "art_root_io/TFileService.h"
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

  class CRTCosmicRemoval : public art::EDAnalyzer {
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
    explicit CRTCosmicRemoval(Parameters const& config);
 
    // Called once, at start of the job
    virtual void beginJob() override;
 
    // Called once per event
    virtual void analyze(const art::Event& event) override;

    // Called once, at end of the job
    virtual void endJob() override;

    double DCA(TVector3 pos, TVector3 start, TVector3 end);

    bool isVisible(simb::MCParticle particle);

    bool entersTPC(simb::MCParticle particle);

    int numTaggerCross(simb::MCParticle particle);

  private:

    // fcl file parameters
    art::InputTag fSimModuleLabel;      ///< name of detsim producer
    art::InputTag fCrtHitModuleLabel; ///< name of CRT track producer
    art::InputTag fTpcTrackModuleLabel; ///< name of TPC track producer
    bool          fVerbose;             ///< print information about what's going on
    
    // histograms
    // Percentage of matched tracks with the correct time
    TGraphAsymmErrors* fPurity;
    TGraphAsymmErrors* fPurityHits;
    TGraphAsymmErrors* fPurityTracks;
    // Percentage of through going tracks with correct match
    TGraphAsymmErrors* fEfficiency;
    TGraphAsymmErrors* fEfficiencyHits;
    TGraphAsymmErrors* fEfficiencyTracks;
    TH1D* hCorrectMatch;
    TH1D* hCorrectMatchHits;
    TH1D* hCorrectMatchTracks;
    TH1D* hTotalMatch;
    TH1D* hTotalMatchHits;
    TH1D* hMatchHits;
    TH1D* hTotalMatchTracks;
    TH1D* hMatchTracks;
    TH1D* hTotalTracks;
    TH1D* hTotalTracksHits;
    TH1D* hTotalTracksTracks;

    int nL = 20;
    double startL = 0;
    double stepL = 10;
    // Percentage of cosmic tracks removed (maybe just primary?)
    TGraphAsymmErrors* fCosmicRemoved;
    // Number of removed cosmic tracks as a function of track length
    TH1D* hCosmicRemoved;
    // Total number of cosmic tracks as a function of track length
    TH1D* hCosmicTotal;
    TGraphAsymmErrors* fCosmicRemovedHits;
    TH1D* hCosmicRemovedHits;
    TGraphAsymmErrors* fCosmicRemovedTracks;
    TH1D* hCosmicRemovedTracks;
    TGraphAsymmErrors* fCosmicRemovedMinL;
    TH1D* hCosmicRemovedMinL;
    TH1D* hCosmicTotalMinL;

    // Percentage of neutrino tracks removed
    TGraphAsymmErrors* fNuRemoved;
    // Number of removed neutrino tracks as a function of track length
    TH1D* hNuRemoved;
    // Total number of neutrino tracks as a function of track length
    TH1D* hNuTotal;
    TGraphAsymmErrors* fNuRemovedHits;
    TH1D* hNuRemovedHits;
    TGraphAsymmErrors* fNuRemovedTracks;
    TH1D* hNuRemovedTracks;
    TGraphAsymmErrors* fNuRemovedMinL;
    TH1D* hNuRemovedMinL;
    TH1D* hNuTotalMinL;

    int nR = 20;
    double startR = 2.5;
    double stepR = 5;
    // Percentage of neutrino charge removed as a function of radius around CRT track
    TGraphAsymmErrors* fNuCharge;
    TH1D* hTotalNuCharge;
    TH1D* hRemovedNuCharge;
    // Percentage of cosmic charge removed as a function of radius around CRT track
    TGraphAsymmErrors* fCosmicCharge;
    TH1D* hTotalCosmicCharge;
    TH1D* hRemovedCosmicCharge;

    TGraphAsymmErrors* fNuChargeComp;
    TGraphAsymmErrors* fCosmicChargeComp;
    TH1D* hRemovedNuChargeComp;
    TH1D* hRemovedCosmicChargeComp;

    // Other variables shared between different methods.
    detinfo::DetectorProperties const* fDetectorProperties;    ///< pointer to detector properties provider
    detinfo::DetectorClocks const* fDetectorClocks;            ///< pointer to detector clocks provider
    TPCGeoAlg fTpcGeo;

    CRTTruthRecoAlg truthAlg;

    // Performance Counters
    int nTracks = 0;
    int nValTracks = 0;
    int nCorrectCross = 0;
    int nIncorrectCross = 0;
    int nNoMatchCross = 0;

    // Tracks from primary cosmic particles
    int nCosmicTracks = 0;
    int nCosmicRemoved = 0;
    // Tracks from neutrino interactions
    int nNuTracks = 0;
    int nNuRemoved = 0;

  }; // class CRTCosmicRemoval


  // Constructor
  CRTCosmicRemoval::CRTCosmicRemoval(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel       (config().SimModuleLabel())
    , fCrtHitModuleLabel    (config().CrtHitModuleLabel())
    , fTpcTrackModuleLabel  (config().TpcTrackModuleLabel())
    , fVerbose              (config().Verbose())
  {

    // Get a pointer to the geometry service provider
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 
    fDetectorClocks     = lar::providerFrom<detinfo::DetectorClocksService>(); 

  } // CRTCosmicRemoval()


  void CRTCosmicRemoval::beginJob()
  {
    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;
    fPurity     = tfs->makeAndRegister<TGraphAsymmErrors>("purity",     ";Reco track length (cm);Purity"    );
    fPurityHits = tfs->makeAndRegister<TGraphAsymmErrors>("purityhits",     ";Reco track length (cm);Purity"    );
    fPurityTracks = tfs->makeAndRegister<TGraphAsymmErrors>("puritytracks",     ";Reco track length (cm);Purity"    );
    fEfficiency = tfs->makeAndRegister<TGraphAsymmErrors>("efficiency", ";Reco track length (cm);Efficiency");
    fEfficiencyHits = tfs->makeAndRegister<TGraphAsymmErrors>("efficiencyhits", ";Reco track length (cm);Efficiency");
    fEfficiencyTracks = tfs->makeAndRegister<TGraphAsymmErrors>("efficiencytracks", ";Reco track length (cm);Efficiency");
    fCosmicRemoved = tfs->makeAndRegister<TGraphAsymmErrors>("cosmicremoved", ";Reco track length (cm);% tracks removed");
    fCosmicRemovedHits = tfs->makeAndRegister<TGraphAsymmErrors>("cosmicremovedhits", ";Reco track length (cm);% tracks removed");
    fCosmicRemovedTracks = tfs->makeAndRegister<TGraphAsymmErrors>("cosmicremovedtracks", ";Reco track length (cm);% tracks removed");
    fCosmicRemovedMinL = tfs->makeAndRegister<TGraphAsymmErrors>("cosmicremovedminl", ";Min reco track length (cm);% tracks removed");
    fNuRemoved = tfs->makeAndRegister<TGraphAsymmErrors>("nuremoved", ";Reco track length (cm);% tracks removed");
    fNuRemovedHits = tfs->makeAndRegister<TGraphAsymmErrors>("nuremovedhits", ";Reco track length (cm);% tracks removed");
    fNuRemovedTracks = tfs->makeAndRegister<TGraphAsymmErrors>("nuremovedtracks", ";Reco track length (cm);% tracks removed");
    fNuRemovedMinL = tfs->makeAndRegister<TGraphAsymmErrors>("nuremovedminl", ";Min reco track length (cm);% tracks removed");
    fCosmicCharge = tfs->makeAndRegister<TGraphAsymmErrors>("cosmiccharge", ";Radius around track (cm);% charge removed");
    fNuCharge = tfs->makeAndRegister<TGraphAsymmErrors>("nucharge", ";Radius around track (cm);% charge removed");
    fCosmicChargeComp = tfs->makeAndRegister<TGraphAsymmErrors>("cosmicchargecomp", ";Radius around track (cm);% charge removed");
    fNuChargeComp = tfs->makeAndRegister<TGraphAsymmErrors>("nuchargecomp", ";Radius around track (cm);% charge removed");

    hCorrectMatch = tfs->make<TH1D>("hCorrectMatch", "", 20, 0, 200);
    hCorrectMatchHits = tfs->make<TH1D>("hCorrectMatchHits", "", 20, 0, 200);
    hCorrectMatchTracks = tfs->make<TH1D>("hCorrectMatchTracks", "", 20, 0, 200);
    hTotalMatch = tfs->make<TH1D>("hTotalMatch", "", 20, 0, 200);
    hTotalMatchHits = tfs->make<TH1D>("hTotalMatchHits", "", 20, 0, 200);
    hMatchHits = tfs->make<TH1D>("hMatchHits", "", 20, 0, 200);
    hTotalMatchTracks = tfs->make<TH1D>("hTotalMatchTracks", "", 20, 0, 200);
    hMatchTracks = tfs->make<TH1D>("hMatchTracks", "", 20, 0, 200);
    hTotalTracks = tfs->make<TH1D>("hTotalTracks", "", 20, 0, 200);
    hTotalTracksHits = tfs->make<TH1D>("hTotalTracksHits", "", 20, 0, 200);
    hTotalTracksTracks = tfs->make<TH1D>("hTotalTracksTracks", "", 20, 0, 200);
    hCosmicRemoved = tfs->make<TH1D>("hCosmicRemoved", "", 20, 0 ,200);
    hCosmicRemovedHits = tfs->make<TH1D>("hCosmicRemovedHits", "", 20, 0 ,200);
    hCosmicRemovedTracks = tfs->make<TH1D>("hCosmicRemovedTracks", "", 20, 0 ,200);
    hCosmicRemovedMinL = tfs->make<TH1D>("hCosmicRemovedMinL", "", nL, startL-(stepL/2), startL+(nL-1/2)*stepL);
    hCosmicTotal = tfs->make<TH1D>("hCosmicTotal", "", 20, 0, 200);
    hCosmicTotalMinL = tfs->make<TH1D>("hCosmicTotalMinL", "", nL, startL-(stepL/2), startL+(nL-1/2)*stepL);
    hNuRemoved = tfs->make<TH1D>("hNuRemoved", "", 20, 0 ,200);
    hNuRemovedHits = tfs->make<TH1D>("hNuRemovedHits", "", 20, 0 ,200);
    hNuRemovedTracks = tfs->make<TH1D>("hNuRemovedTracks", "", 20, 0 ,200);
    hNuRemovedMinL = tfs->make<TH1D>("hNuRemovedMinL", "", nL, startL-(stepL/2), startL+(nL-1/2)*stepL);
    hNuTotal = tfs->make<TH1D>("hNuTotal", "", 20, 0, 200);
    hNuTotalMinL = tfs->make<TH1D>("hNuTotalMinL", "", nL, startL-(stepL/2), startL+(nL-1/2)*stepL);
    hTotalNuCharge = tfs->make<TH1D>("hTotalNuCharge", "", nR, startR-(stepR/2), startR+(nR-1/2)*stepR);
    hRemovedNuCharge = tfs->make<TH1D>("hRemovedNuCharge", "", nR, startR-(stepR/2), startR+(nR-1/2)*stepR);
    hTotalCosmicCharge = tfs->make<TH1D>("hTotalCosmicCharge", "", nR, startR-(stepR/2), startR+(nR-1/2)*stepR);
    hRemovedCosmicCharge = tfs->make<TH1D>("hRemovedCosmicCharge", "", nR, startR-(stepR/2), startR+(nR-1/2)*stepR);
    hRemovedNuChargeComp = tfs->make<TH1D>("hRemovedNuChargeComp", "", nR, startR-(stepR/2), startR+(nR-1/2)*stepR);
    hRemovedCosmicChargeComp = tfs->make<TH1D>("hRemovedCosmicChargeComp", "", nR, startR-(stepR/2), startR+(nR-1/2)*stepR);

    // Initial output
    if(fVerbose) std::cout<<"----------------- CRT Cosmic Removal Ana Module -------------------"<<std::endl;

  } // CRTCosmicRemoval::beginJob()


  void CRTCosmicRemoval::analyze(const art::Event& event)
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

    double readoutWindowMuS  = fDetectorClocks->TPCTick2Time((double)fDetectorProperties->ReadOutWindowSize()); // [us]
    double driftTimeMuS = fTpcGeo.MaxX()/fDetectorProperties->DriftVelocity(); // [us]

    // Loop over true particles in readout window
    std::map<int, simb::MCParticle> particles;
    std::vector<int> nuParticleIds;
    std::vector<std::pair<TVector3, double>> nuEnergy;
    std::vector<std::pair<TVector3, double>> cosmicEnergy;

    for (auto const& particle: (*particleHandle)){

      int partId = particle.TrackId();
      particles[partId] = particle;
      art::Ptr<simb::MCTruth> truth = pi_serv->TrackIdToMCTruth_P(partId);
      if(truth->Origin() == simb::kBeamNeutrino){
        nuParticleIds.push_back(partId);
      }

      double time = particle.T() * 1e-3; // [us]
      /*int pdg = std::abs(particle.PdgCode());
      double energy = particle.E();
      TVector3 start(particle.Vx(), particle.Vy(), particle.Vz());
      TVector3 end(particle.EndX(), particle.EndY(), particle.EndZ());*/

      if(time < -driftTimeMuS || time > readoutWindowMuS) continue;
      // If particle is visible
      if(!isVisible(particle)) continue;
      if(!entersTPC(particle)) continue;

      // Add up charge deposited in the detector
      int nTraj = particle.NumberTrajectoryPoints();
      double totEDep = 0;
      std::vector<std::pair<TVector3, double>> posEDep;
      for(int i = 0; i < nTraj-1; i++){
        TVector3 pos(particle.Vx(i), particle.Vy(i), particle.Vz(i));
        if(RecoUtils::IsInsideTPC(pos, 0)){
          double eDep = particle.E(i) - particle.E(i+1);
          totEDep += eDep;
          posEDep.push_back(std::make_pair(pos, eDep));
        }
      }

      // If from a neutrino interaction
      if(truth->Origin() == simb::kBeamNeutrino){
        // Add up total neutrino charge
        for(int j = 0; j < nR; j++){
          hTotalNuCharge->Fill(startR+j*stepR, totEDep);
        }
        nuEnergy.insert(nuEnergy.end(), posEDep.begin(), posEDep.end());
        //std::cout<<"Nu: ID = "<<partId<<" pdg = "<<pdg<<" time = "<<time<<" energy = "<<energy<<" pos = ("<<start.X()<<", "<<start.Y()<<", "<<start.Z()<<")->("<<end.X()<<", "<<end.Y()<<", "<<end.Z()<<") total EDep = "<<totEDep<<"\n";
      }

      // Else
      else{
        // Add up total cosmic charge
        for(int j = 0; j < nR; j++){
          hTotalCosmicCharge->Fill(startR+j*stepR, totEDep);
        }
        cosmicEnergy.insert(cosmicEnergy.end(), posEDep.begin(), posEDep.end());
        //std::cout<<"Cosmic: ID = "<<partId<<" pdg = "<<pdg<<" time = "<<time<<" energy = "<<energy<<" pos = ("<<start.X()<<", "<<start.Y()<<", "<<start.Z()<<")->("<<end.X()<<", "<<end.Y()<<", "<<end.Z()<<") total EDep = "<<totEDep<<"\n";
      }
    }

    if(fVerbose) std::cout<<"Number of true particles = "<<particles.size()<<std::endl;

    // Loop over neutrino energy
    for(auto const& nuE : nuEnergy){
      std::vector<bool> removed(nR, false);
      std::vector<bool> removedComp(nR, false);
      //std::cout<<"Nu EDep pos = ("<<nuE.first.X()<<", "<<nuE.first.Y()<<", "<<nuE.first.Z()<<") edep = "<<nuE.second<<"\n";

      // Loop over CRT tracks
      double minDca = 99999;
      for(auto const& crtTrack : crtTracks){
        // For each track create a cylinder of radius R around track (use incomplete?)
        TVector3 start(crtTrack.x1_pos, crtTrack.y1_pos, crtTrack.z1_pos);
        TVector3 end(crtTrack.x2_pos, crtTrack.y2_pos, crtTrack.z2_pos);
        // If trajectory point inside R then removed
        double dca = DCA(nuE.first, start, end);
        if(dca < minDca) minDca = dca;
        for(int i = 0; i < nR; i++){
          if(dca < startR + (i*stepR)) removed[i] = true;
          if(dca < startR + (i*stepR) && crtTrack.complete) removedComp[i] = true;
        }
      }

      //std::cout<<"Min dist = "<<minDca<<"\n";
      // Add up removed neutrino charge
      for(int i = 0; i < nR; i++){
        if(removed[i]) hRemovedNuCharge->Fill(startR+i*stepR, nuE.second);
        if(removedComp[i]) hRemovedNuChargeComp->Fill(startR+i*stepR, nuE.second);
      }

    }

    // Loop over cosmic energy
    for(auto const& cosE : cosmicEnergy){
      std::vector<bool> removed(nR, false);
      std::vector<bool> removedComp(nR, false);
      //std::cout<<"Cosmic EDep pos = ("<<cosE.first.X()<<", "<<cosE.first.Y()<<", "<<cosE.first.Z()<<") edep = "<<cosE.second<<"\n";

      // Loop over CRT tracks
      double minDca = 99999;
      for(auto const& crtTrack : crtTracks){
        // For each track create a cylinder of radius R around track (use incomplete?)
        //if(!crtTrack.complete) continue;

        TVector3 start(crtTrack.x1_pos, crtTrack.y1_pos, crtTrack.z1_pos);
        TVector3 end(crtTrack.x2_pos, crtTrack.y2_pos, crtTrack.z2_pos);
        // If trajectory point inside R then removed
        double dca = DCA(cosE.first, start, end);
        if(dca < minDca) minDca = dca;
        for(int i = 0; i < nR; i++){
          if(dca < startR + (i*stepR)) removed[i] = true;
          if(dca < startR + (i*stepR) && crtTrack.complete) removedComp[i] = true;
        }
      }

      //std::cout<<"Min dist = "<<minDca<<"\n";
      // Add up removed cosmic charge
      for(int i = 0; i < nR; i++){
        if(removed[i]) hRemovedCosmicCharge->Fill(startR+i*stepR, cosE.second);
        if(removedComp[i]) hRemovedCosmicChargeComp->Fill(startR+i*stepR, cosE.second);
      }

    }
    
    // Loop over the tpc tracks
    for(auto const& tpcTrack : (*tpcTrackHandle)){

      nTracks++;

      // Match to the true particle
      std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
      int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits);
      if (particles.find(trueId) == particles.end()){ 
        if (fVerbose) std::cout<<"No valid true track!\n"; 
        continue; 
      }
      // Get the true T0
      double trueTime = particles[trueId].T()*1e-3;

      if(fVerbose) std::cout<<"Track "<<tpcTrack.ID()<<" trueID = "<<trueId<<" pdg = "<<particles[trueId].PdgCode()<<" time = "<<trueTime<<" reco length = "<<tpcTrack.Length()<<"\n";

      // Is track from a neutrino interaction
      bool isNu = false;
      bool isCos = false;
      if(std::find(nuParticleIds.begin(), nuParticleIds.end(), trueId) != nuParticleIds.end()){
        // Fill total neutrino tracks with length
        hNuTotal->Fill(tpcTrack.Length());
        for(int i = 0; i < nL; i++){
          if(tpcTrack.Length()>(startL+i*stepL)) hNuTotalMinL->Fill(startL+i*stepL);
        }
        // Set neutrino flag
        isNu = true;
        if(tpcTrack.Length()>40) nNuTracks++;
        if(fVerbose) std::cout<<"From neutrino!\n";
      }
      // Else is it from a primary cosmic (true length > 500)
      else if(particles[trueId].Trajectory().TotalLength() > 500.){
        // Fill total cosmic tracks with length
        hCosmicTotal->Fill(tpcTrack.Length());
        for(int i = 0; i < nL; i++){
          if(tpcTrack.Length()>(startL+i*stepL)) hCosmicTotalMinL->Fill(startL+i*stepL);
        }
        // Set primary cosmic flag
        isCos = true;
        if(tpcTrack.Length()>40) nCosmicTracks++;
        if(fVerbose) std::cout<<"Primary cosmic!\n";
      }

      // Check if track is stitched
      int tpc = hits[0]->WireID().TPC;
 
      nValTracks++;
      hTotalTracks->Fill(tpcTrack.Length());
      int nCross = numTaggerCross(particles[trueId]);
      if(nCross >= 1) hTotalTracksHits->Fill(tpcTrack.Length());
      if(nCross >= 2) hTotalTracksTracks->Fill(tpcTrack.Length());

      // Try to get T0 from CRTTracks
      double trackTime = CRTAnaUtils::T0FromCRTTracks(tpcTrack, crtTracks, tpc, 0.2, 80.);
      if(trackTime != -99999){
        hTotalMatchTracks->Fill(tpcTrack.Length());
        if(nCross >= 2) hMatchTracks->Fill(tpcTrack.Length());
        if(std::abs(trueTime-trackTime) < 5.) hCorrectMatchTracks->Fill(tpcTrack.Length());
      }

      // Try to get T0 from CRTHits
      double hitTime = CRTAnaUtils::T0FromCRTHits(tpcTrack, crtHits, tpc, 0., 0.5, 30.);
      if(hitTime != -99999){
        hTotalMatchHits->Fill(tpcTrack.Length());
        if(nCross >= 1) hMatchHits->Fill(tpcTrack.Length());
        if(std::abs(trueTime-hitTime) < 5.) hCorrectMatchHits->Fill(tpcTrack.Length());
      }

      if(fVerbose) std::cout<<"True time = "<<trueTime<<" ticks, track time = "<<trackTime
                            <<" ticks, hit time = "<<hitTime<<"\n";

      double matchedTime = -99999;
      if(trackTime != -99999) matchedTime = trackTime;
      else if(hitTime != -99999) matchedTime = hitTime;

      if(fVerbose) std::cout<<"Matched time = "<<matchedTime<<"\n";

      double timeLim = 20.;
      // If neutrino flag
      if(isNu){
        if(matchedTime != -99999 && (std::abs(matchedTime) > timeLim)){
          if(fVerbose) std::cout<<"REMOVED!\n";
          // Fill removed neutrino tracks with length
          hNuRemoved->Fill(tpcTrack.Length());
          for(int i = 0; i < nL; i++){
            if(tpcTrack.Length()>(startL+i*stepL)) hNuRemovedMinL->Fill(startL+i*stepL);
          }
          if(tpcTrack.Length()>40) nNuRemoved++;
        }
        if(hitTime != -99999 && (std::abs(hitTime) > timeLim)){
          hNuRemovedHits->Fill(tpcTrack.Length());
        }
        if(trackTime != -99999 && (std::abs(hitTime) > timeLim)){
          hNuRemovedTracks->Fill(tpcTrack.Length());
        }
      }
      // If primary cosmic flag
      else if(isCos){
        if(matchedTime != -99999 && (std::abs(matchedTime) > timeLim)){
          if(fVerbose) std::cout<<"REMOVED!\n";
          // Fill removed cosmic tracks with length
          hCosmicRemoved->Fill(tpcTrack.Length());
          for(int i = 0; i < nL; i++){
            if(tpcTrack.Length()>(startL+i*stepL)) hCosmicRemovedMinL->Fill(startL+i*stepL);
          }
          if(tpcTrack.Length()>40) nCosmicRemoved++;
        }
        if(hitTime != -99999 && (std::abs(hitTime) > timeLim)){
          hCosmicRemovedHits->Fill(tpcTrack.Length());
        }
        if(trackTime != -99999 && (std::abs(hitTime) > timeLim)){
          hCosmicRemovedTracks->Fill(tpcTrack.Length());
        }
      }

      if(std::abs(trueTime-matchedTime) < 5.){
        nCorrectCross++;
        hCorrectMatch->Fill(tpcTrack.Length());
        hTotalMatch->Fill(tpcTrack.Length());
      }
      else if(matchedTime == -99999){
        nNoMatchCross++;
      }
      else{
        nIncorrectCross++;
        hTotalMatch->Fill(tpcTrack.Length());
      }

    }


  } // CRTCosmicRemoval::analyze()


  void CRTCosmicRemoval::endJob(){

    std::cout<<"Total tracks                   = "<<nTracks<<"\n"
             <<"Total truth-matched tracks     = "<<nValTracks<<"\n"
             <<"CRT Crossing tracks:\n"
             <<"Total tracks with correct T0   = "<<nCorrectCross<<"\n"
             <<"Total tracks with incorrect T0 = "<<nIncorrectCross<<"\n"
             <<"Total tracks with no T0 match  = "<<nNoMatchCross<<"\n"
             <<"Efficiency = "<<(double)(nCorrectCross+nIncorrectCross)/nValTracks<<"\n"
             <<"Purity     = "<<(double)nCorrectCross/(nCorrectCross+nIncorrectCross)<<"\n"
             <<"Total nu tracks = "<<nNuTracks<<"\n"
             <<"Percentage removed = "<<(double)nNuRemoved/nNuTracks<<"\n"
             <<"Total cosmic tracks = "<<nCosmicTracks<<"\n"
             <<"Percentage removed = "<<(double)nCosmicRemoved/nCosmicTracks<<"\n";

    //Purity = tracks with correct match / all matched tracks
    fPurity->SetMarkerStyle(8);
    fPurity->SetMarkerColor(1);
    fPurity->SetLineColor(1);
    fPurity->SetLineWidth(3);
    fPurity->BayesDivide(hCorrectMatch, hTotalMatch);
    fPurity->Draw("ap");

    fPurityHits->SetMarkerStyle(8);
    fPurityHits->SetMarkerColor(1);
    fPurityHits->SetLineColor(1);
    fPurityHits->SetLineWidth(3);
    fPurityHits->BayesDivide(hCorrectMatchHits, hTotalMatchHits);
    fPurityHits->Draw("ap");

    fPurityTracks->SetMarkerStyle(8);
    fPurityTracks->SetMarkerColor(1);
    fPurityTracks->SetLineColor(1);
    fPurityTracks->SetLineWidth(3);
    fPurityTracks->BayesDivide(hCorrectMatchTracks, hTotalMatchTracks);
    fPurityTracks->Draw("ap");

    //Efficiency = tracks with match / all through-going tracks
    fEfficiency->SetMarkerStyle(8);
    fEfficiency->SetMarkerColor(1);
    fEfficiency->SetLineColor(1);
    fEfficiency->SetLineWidth(3);
    fEfficiency->BayesDivide(hCorrectMatch, hTotalTracks);
    fEfficiency->Draw("ap");

    fEfficiencyHits->SetMarkerStyle(8);
    fEfficiencyHits->SetMarkerColor(1);
    fEfficiencyHits->SetLineColor(1);
    fEfficiencyHits->SetLineWidth(3);
    fEfficiencyHits->BayesDivide(hMatchHits, hTotalTracksHits);
    fEfficiencyHits->Draw("ap");

    fEfficiencyTracks->SetMarkerStyle(8);
    fEfficiencyTracks->SetMarkerColor(1);
    fEfficiencyTracks->SetLineColor(1);
    fEfficiencyTracks->SetLineWidth(3);
    fEfficiencyTracks->BayesDivide(hMatchTracks, hTotalTracksTracks);
    fEfficiencyTracks->Draw("ap");

    //CosmicRemoved = % of primary muon tracks removed as function of length
    fCosmicRemoved->SetMarkerStyle(8);
    fCosmicRemoved->SetMarkerColor(1);
    fCosmicRemoved->SetLineColor(1);
    fCosmicRemoved->SetLineWidth(3);
    fCosmicRemoved->BayesDivide(hCosmicRemoved, hCosmicTotal);
    fCosmicRemoved->Draw("ap");

    fCosmicRemovedHits->SetMarkerStyle(8);
    fCosmicRemovedHits->SetMarkerColor(1);
    fCosmicRemovedHits->SetLineColor(1);
    fCosmicRemovedHits->SetLineWidth(3);
    fCosmicRemovedHits->BayesDivide(hCosmicRemovedHits, hCosmicTotal);
    fCosmicRemovedHits->Draw("ap");

    fCosmicRemovedTracks->SetMarkerStyle(8);
    fCosmicRemovedTracks->SetMarkerColor(1);
    fCosmicRemovedTracks->SetLineColor(1);
    fCosmicRemovedTracks->SetLineWidth(3);
    fCosmicRemovedTracks->BayesDivide(hCosmicRemovedTracks, hCosmicTotal);
    fCosmicRemovedTracks->Draw("ap");

    fCosmicRemovedMinL->SetMarkerStyle(8);
    fCosmicRemovedMinL->SetMarkerColor(1);
    fCosmicRemovedMinL->SetLineColor(1);
    fCosmicRemovedMinL->SetLineWidth(3);
    fCosmicRemovedMinL->BayesDivide(hCosmicRemovedMinL, hCosmicTotalMinL);
    fCosmicRemovedMinL->Draw("ap");

    //NuRemoved = % of neutrino induced tracks removed as function of length
    fNuRemoved->SetMarkerStyle(8);
    fNuRemoved->SetMarkerColor(1);
    fNuRemoved->SetLineColor(1);
    fNuRemoved->SetLineWidth(3);
    fNuRemoved->BayesDivide(hNuRemoved, hNuTotal);
    fNuRemoved->Draw("ap");

    fNuRemovedHits->SetMarkerStyle(8);
    fNuRemovedHits->SetMarkerColor(1);
    fNuRemovedHits->SetLineColor(1);
    fNuRemovedHits->SetLineWidth(3);
    fNuRemovedHits->BayesDivide(hNuRemovedHits, hNuTotal);
    fNuRemovedHits->Draw("ap");

    fNuRemovedTracks->SetMarkerStyle(8);
    fNuRemovedTracks->SetMarkerColor(1);
    fNuRemovedTracks->SetLineColor(1);
    fNuRemovedTracks->SetLineWidth(3);
    fNuRemovedTracks->BayesDivide(hNuRemovedTracks, hNuTotal);
    fNuRemovedTracks->Draw("ap");

    fNuRemovedMinL->SetMarkerStyle(8);
    fNuRemovedMinL->SetMarkerColor(1);
    fNuRemovedMinL->SetLineColor(1);
    fNuRemovedMinL->SetLineWidth(3);
    fNuRemovedMinL->BayesDivide(hNuRemovedMinL, hNuTotalMinL);
    fNuRemovedMinL->Draw("ap");

    fCosmicCharge->SetMarkerStyle(8);
    fCosmicCharge->SetMarkerColor(1);
    fCosmicCharge->SetLineColor(1);
    fCosmicCharge->SetLineWidth(3);
    fCosmicCharge->BayesDivide(hRemovedCosmicCharge, hTotalCosmicCharge);
    fCosmicCharge->Draw("ap");

    fNuCharge->SetMarkerStyle(8);
    fNuCharge->SetMarkerColor(1);
    fNuCharge->SetLineColor(1);
    fNuCharge->SetLineWidth(3);
    fNuCharge->BayesDivide(hRemovedNuCharge, hTotalNuCharge);
    fNuCharge->Draw("ap");

    fCosmicChargeComp->SetMarkerStyle(8);
    fCosmicChargeComp->SetMarkerColor(1);
    fCosmicChargeComp->SetLineColor(1);
    fCosmicChargeComp->SetLineWidth(3);
    fCosmicChargeComp->BayesDivide(hRemovedCosmicChargeComp, hTotalCosmicCharge);
    fCosmicChargeComp->Draw("ap");

    fNuChargeComp->SetMarkerStyle(8);
    fNuChargeComp->SetMarkerColor(1);
    fNuChargeComp->SetLineColor(1);
    fNuChargeComp->SetLineWidth(3);
    fNuChargeComp->BayesDivide(hRemovedNuChargeComp, hTotalNuCharge);
    fNuChargeComp->Draw("ap");


  } // CRTCosmicRemoval::endJob()

  double CRTCosmicRemoval::DCA(TVector3 pos, TVector3 start, TVector3 end){
    double denominator = (end - start).Mag();
    double numerator = ((pos - start).Cross(pos - end)).Mag();
    return numerator/denominator;
  }

  bool CRTCosmicRemoval::isVisible(simb::MCParticle particle){
    int pdg = std::abs(particle.PdgCode());
    double energy = particle.E();
    if(energy < 0.01) return false;
    if(!(pdg==11||pdg==13||pdg==211||pdg==321||pdg==2212)) return false;
    return true;
  }

  bool CRTCosmicRemoval::entersTPC(simb::MCParticle particle){
    bool enters = false;
    int nTraj = particle.NumberTrajectoryPoints();
    for(int i = 0; i < nTraj-1; i++){
      TVector3 pos(particle.Vx(i), particle.Vy(i), particle.Vz(i));
      if(RecoUtils::IsInsideTPC(pos, 0)){
        enters = true;
      }
    }
    return enters;
  }

  int CRTCosmicRemoval::numTaggerCross(simb::MCParticle particle){

    int nCross = 0;
    for(int i = 0; i < 7; i++){
      if(truthAlg.CrossesTagger(particle, i)) nCross++;
    }
    return nCross;

  }

  DEFINE_ART_MODULE(CRTCosmicRemoval)
} // namespace sbnd

