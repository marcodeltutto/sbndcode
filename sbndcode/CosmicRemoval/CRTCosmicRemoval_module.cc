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
#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"

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
#include "larsim/Simulation/LArG4Parameters.h"
#include "larcorealg/CoreUtils/NumericUtils.h" // util::absDiff()

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

  private:

    // fcl file parameters
    art::InputTag fSimModuleLabel;      ///< name of detsim producer
    art::InputTag fCrtHitModuleLabel; ///< name of CRT track producer
    art::InputTag fTpcTrackModuleLabel; ///< name of TPC track producer
    bool          fVerbose;             ///< print information about what's going on
    
    // histograms
    TGraphAsymmErrors* fPurity;
    TGraphAsymmErrors* fEfficiency;
    TH1D* hCorrectMatch = new TH1D("hCorrectMatch", "", 20, 0, 200);
    TH1D* hTotalMatch = new TH1D("hTotalMatch", "", 20, 0, 200);
    TH1D* hTotalTracks = new TH1D("hTotalTracks", "", 20, 0, 200);

    // Other variables shared between different methods.
    geo::GeometryCore const* fGeometryService;                 ///< pointer to Geometry provider
    detinfo::DetectorProperties const* fDetectorProperties;    ///< pointer to detector properties provider
    detinfo::DetectorClocks const* fDetectorClocks;            ///< pointer to detector clocks provider

    // Performance Counters
    int nTracks = 0;
    int nLongTracks = 0;
    int nValTracks = 0;
    int nCorrectCross = 0;
    int nIncorrectCross = 0;
    int nNoMatchCross = 0;
    int nCorrectCont = 0;
    int nIncorrectCont = 0;

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
    fGeometryService    = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 
    fDetectorClocks     = lar::providerFrom<detinfo::DetectorClocksService>(); 

  } // CRTCosmicRemoval()


  void CRTCosmicRemoval::beginJob()
  {
    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;
    fPurity     = tfs->makeAndRegister<TGraphAsymmErrors>("purity",     ";Reco track length (cm);Purity"    );
    fEfficiency = tfs->makeAndRegister<TGraphAsymmErrors>("efficiency", ";Reco track length (cm);Efficiency");

    // Initial output
    if(fVerbose) std::cout<<"----------------- CRT Cosmic Removal Ana Module -------------------"<<std::endl;

  } // CRTCosmicRemoval::beginJob()


  void CRTCosmicRemoval::analyze(const art::Event& event)
  {

    // Fetch basic event info
    if(fVerbose){
      std::cout<<"============================================"<<std::endl
               <<"Run = "<<event.run()<<", SubRun = "<<event.subRun()<<", Event = "<<event.id().event()<<std::endl
               <<"============================================"<<std::endl;
    }

    // Retrieve all the truth info in the events
    auto particleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fSimModuleLabel);
    // Fill a map of true particles
    std::map<int, simb::MCParticle> particles;
    std::vector<simb::MCParticle> parts;
    for (auto const& particle: (*particleHandle)){
      int partId = particle.TrackId();
      particles[partId] = particle;
      parts.push_back(particle);
    }

    if(fVerbose) std::cout<<"Number of true particles = "<<particles.size()<<std::endl;

    // Retrieve list of CRT tracks
    auto crtHitHandle = event.getValidHandle<std::vector<crt::CRTHit>>(fCrtHitModuleLabel);

    // Retrieve the TPC tracks
    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTpcTrackModuleLabel);

    // Get track to hit associations
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTpcTrackModuleLabel);

    art::FindManyP<anab::T0> findManyT0(tpcTrackHandle, event, fTpcTrackModuleLabel);

    if(fVerbose) std::cout<<"Number of CRT hits = "<<crtHitHandle->size()<<std::endl
                          <<"Number of TPC tracks = "<<tpcTrackHandle->size()<<std::endl;

    // Do track reconstruction from the hits
    std::vector<crt::CRTHit> crtHits;
    for(auto const& crtHit : (*crtHitHandle)){
      crtHits.push_back(crtHit);
    }
    std::vector<crt::CRTTrack> crtTracks = CRTAnaUtils::CreateCRTTracks(crtHits, 0.2, 30., true, 25.); 
    if(fVerbose) std::cout<<"Number of CRTTracks = "<<crtTracks.size()<<std::endl;
    
    // Loop over the tpc tracks
    for(auto const& tpcTrack : (*tpcTrackHandle)){

      nTracks++;

      // Only look at tracks longer than some limit
      //if(tpcTrack.Length() < 20.) continue;
      nLongTracks++;

      // Match to the true particle
      std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
      int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits);
      if (particles.find(trueId) == particles.end()){ 
        if (fVerbose) std::cout<<"No valid true track!\n"; 
        continue; 
      }
      // Get the true T0
      double trueTime = fDetectorClocks->TPCG4Time2Tick(particles[trueId].T());

      //bool crtCrosser = false;
      // Loop over true particle trajectory
      if(particles[trueId].Trajectory().TotalLength() < 500) continue;
/*
      // Find any existing T0 associations
      if(findManyT0.at(tpcTrack.ID()).size() > 0){
        if(fVerbose) std::cout<<"Track has T0 associated! true time = "<<trueTime
                              <<" ticks, T0 size = "<<findManyT0.at(tpcTrack.ID()).size()<<"\n";
        continue;
      }*/

      // Check if track is stitched
      int tpc = hits[0]->WireID().TPC;
     /* if (tpc != (int)hits[hits.size()-1]->WireID().TPC){ 
        if(fVerbose) std::cout<<"Track stitched across CPA! true time = "<<trueTime<<" ticks\n";
        continue;
      }
*/
      nValTracks++;
      hTotalTracks->Fill(tpcTrack.Length());

      // Try to get T0 from CRTTracks
      double trackTime = CRTAnaUtils::T0FromCRTTracks(tpcTrack, crtTracks, tpc, 0.15, 60.);
      // Try to get T0 from CRTHits
      //double hitTime = CRTAnaUtils::T0FromCRTHits(tpcTrack, crtHits, tpc, 20., 0.5, 40.);
      double hitTime = CRTAnaUtils::T0FromCRTHits(tpcTrack, crtHits, tpc, 0., 0.5, 40.);
      if(fVerbose) std::cout<<"True time = "<<trueTime<<" ticks, track time = "<<trackTime
                            <<" ticks, hit time = "<<hitTime<<"\n";

      //if(crtCrosser){
        if(std::abs(trueTime-trackTime) < 5. || (trackTime==-99999 && std::abs(trueTime-hitTime) < 5.)){
          nCorrectCross++;
          hCorrectMatch->Fill(tpcTrack.Length());
          hTotalMatch->Fill(tpcTrack.Length());
        }
        else if(trackTime == -99999 && hitTime == -99999){
          nNoMatchCross++;
        }
        else{
          nIncorrectCross++;
          hTotalMatch->Fill(tpcTrack.Length());
        }
      /*}
      else{
        if(trackTime == -99999 && hitTime == -9999){
          nCorrectCont++;
        }
        else{
          nIncorrectCont++;
        }
      }*/
    }


  } // CRTCosmicRemoval::analyze()


  void CRTCosmicRemoval::endJob(){

    std::cout<<"Total tracks                   = "<<nTracks<<"\n"
             <<"Total tracks > 20 cm           = "<<nLongTracks<<"\n"
             <<"Total tracks with no T0        = "<<nValTracks<<"\n"
             <<"CRT Crossing tracks:\n"
             <<"Total tracks with correct T0   = "<<nCorrectCross<<"\n"
             <<"Total tracks with incorrect T0 = "<<nIncorrectCross<<"\n"
             <<"Total tracks with no T0 match  = "<<nNoMatchCross<<"\n"
             <<"Contained tracks:\n"
             <<"Total tracks with no T0 match  = "<<nCorrectCont<<"\n"
             <<"Total tracks with T0 match     = "<<nIncorrectCont<<"\n"
             <<"Efficiency = "<<(double)(nCorrectCross+nCorrectCont)/nValTracks<<"\n"
             <<"Purity = "<<(double)nCorrectCross/(nCorrectCross+nIncorrectCross+nIncorrectCont)<<"\n";

    //Purity = tracks with good match / all matched tracks
    fPurity->SetMarkerStyle(8);
    fPurity->SetMarkerColor(1);
    fPurity->SetLineColor(1);
    fPurity->SetLineWidth(3);
    fPurity->BayesDivide(hCorrectMatch, hTotalMatch);
    fPurity->Draw("ap");

    //Efficiency = tracks with good match / all matched tracks
    fEfficiency->SetMarkerStyle(8);
    fEfficiency->SetMarkerColor(1);
    fEfficiency->SetLineColor(1);
    fEfficiency->SetLineWidth(3);
    fEfficiency->BayesDivide(hCorrectMatch, hTotalTracks);
    fEfficiency->Draw("ap");

  } // CRTCosmicRemoval::endJob()

  DEFINE_ART_MODULE(CRTCosmicRemoval)
} // namespace sbnd

