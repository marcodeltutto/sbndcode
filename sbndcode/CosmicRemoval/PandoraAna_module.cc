////////////////////////////////////////////////////////////////////////
// Class:       PandoraAna
// Module Type: analyzer
// File:        PandoraAna_module.cc
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbndcode/CosmicRemoval/CosmicRemovalUtils/CosmicRemovalUtils.h"
#include "sbndcode/CRT/CRTUtils/CRTAnaUtils.h"
#include "sbndcode/CRT/CRTUtils/CRTTruthRecoAlg.h"
#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"

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

  class PandoraAna : public art::EDAnalyzer {
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

    }; // Config
 
    using Parameters = art::EDAnalyzer::Table<Config>;
 
    // Constructor: configures module
    explicit PandoraAna(Parameters const& config);
 
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

    // histograms
    TH1D *hTrueMuLength;
    TH1D *hTrueMuE;
    TH1D *hTrueMuLengthReco;
    TH1D *hTrueMuEReco;
    TH1D *hTrueCrLength;
    TH1D *hRecoCrLength;
    TH1D *hMuLengthDiff;
    TH1D *hTrueMuTheta;
    TH1D *hTrueMuPhi;
    TH1D *hRecoLongLength;
    TH1D *hRecoLongTheta;
    TH1D *hRecoLongPhi;
    TH1D *hRecoMuLength;
    TH1D *hRecoMuTheta;
    TH1D *hRecoMuPhi;
    TH1D *hTrueCrTheta;
    TH1D *hTrueCrPhi;
    TH1D *hRecoCrTheta;
    TH1D *hRecoCrPhi;
    TH1D *hTrueSecTrackAngle;
    TH1D *hRecoSecTrackAngle;
    TH1D *hTrueSecTrackLength;
    TH1D *hRecoSecTrackLength;
    TH1D *hRecoCrSecTrackAngle;
    TH1D *hRecoCrSecTrackLength;
    TH1D *hMuVertexDiff;

    TH2D *hTrueMuThetaPhi;
    TH2D *hRecoMuThetaPhi;
    TH2D *hTrueCrThetaPhi;
    TH2D *hRecoCrThetaPhi;
    TH2D *hTrueSecTrackLengthAngle;
    TH2D *hRecoSecTrackLengthAngle;
    TH2D *hRecoCrSecTrackLengthAngle;

    CRTTruthRecoAlg truthAlg;

    int nNuMuCC = 0;
    int nMu = 0;
    int nNuIsMu = 0;
    int nCrIsMu = 0;

    int nPandCr = 0;
    int nPandNu = 0;
    int nNuPandCut = 0;
    int nNuMuPandCut = 0;
    int nCrPandCut = 0;
    int nNuPandKeep = 0;
    int nNuMuPandKeep = 0;
    int nCrPandKeep = 0;

    void GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap);

    void PrintOutScores(const art::Event &evt, const PFParticleHandle &pfParticleHandle) const;

    void GetFinalStatePFParticleVectors(const PFParticleIdMap &pfParticleMap, PFParticleVector &crParticles, PFParticleVector &nuParticles);

    void CollectTracksAndShowers(const PFParticleVector &particles, const PFParticleHandle &pfParticleHandle, const art::Event &evt, TrackVector &tracks, ShowerVector &showers);

  }; // class PandoraAna


  // Constructor
  PandoraAna::PandoraAna(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel       (config().SimModuleLabel())
    , fCrtHitModuleLabel    (config().CrtHitModuleLabel())
    , fTpcTrackModuleLabel  (config().TpcTrackModuleLabel())
    , fCaloModuleLabel      (config().CaloModuleLabel())
    , fShowerModuleLabel    (config().ShowerModuleLabel())
    , fPandoraLabel         (config().PandoraLabel())
    , fVerbose              (config().Verbose())
  {

  } // PandoraAna()


  void PandoraAna::beginJob()
  {
    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;

    hTrueMuLength         = tfs->make<TH1D>("hTrueMuLength",          "", 20, 0,    500);
    hTrueMuE              = tfs->make<TH1D>("hTrueMuE",               "", 20, 0,    2);
    hTrueMuLengthReco     = tfs->make<TH1D>("hTrueMuLengthReco",      "", 20, 0,    500);
    hTrueMuEReco          = tfs->make<TH1D>("hTrueMuEReco",           "", 20, 0,    2);
    hTrueCrLength         = tfs->make<TH1D>("hTrueCrLength",          "", 20, 0,    500);
    hRecoCrLength         = tfs->make<TH1D>("hRecoCrLength",          "", 20, 0,    500);
    hMuLengthDiff         = tfs->make<TH1D>("hMuLengthDiff",          "", 50, -50, 50);
    hMuVertexDiff         = tfs->make<TH1D>("hMuVertexDiff",          "", 30, 0, 30);
    hTrueMuTheta          = tfs->make<TH1D>("hTrueMuTheta",           "", 20, 0,    3.2);
    hTrueMuPhi            = tfs->make<TH1D>("hTrueMuPhi",             "", 20, -3.2, 3.2);
    hRecoLongLength       = tfs->make<TH1D>("hRecoLongLength",        "", 20, 0,    500);
    hRecoLongTheta        = tfs->make<TH1D>("hRecoLongTheta",         "", 20, 0,    3.2);
    hRecoLongPhi          = tfs->make<TH1D>("hRecoLongPhi",           "", 20, -3.2, 3.2);
    hRecoMuLength         = tfs->make<TH1D>("hRecoMuLength",          "", 20, 0,    500);
    hRecoMuTheta          = tfs->make<TH1D>("hRecoMuTheta",           "", 20, 0,    3.2);
    hRecoMuPhi            = tfs->make<TH1D>("hRecoMuPhi",             "", 20, -3.2, 3.2);
    hTrueCrTheta          = tfs->make<TH1D>("hTrueCrTheta",           "", 20, 0,    3.2);
    hTrueCrPhi            = tfs->make<TH1D>("hTrueCrPhi",             "", 20, -3.2, 3.2);
    hRecoCrTheta          = tfs->make<TH1D>("hRecoCrTheta",           "", 20, 0,    3.2);
    hRecoCrPhi            = tfs->make<TH1D>("hRecoCrPhi",             "", 20, -3.2, 3.2);
    hTrueSecTrackAngle    = tfs->make<TH1D>("hTrueSecTrackAngle",     "", 20, 0,    3.2);
    hRecoSecTrackAngle    = tfs->make<TH1D>("hRecoSecTrackAngle",     "", 20, 0,    3.2);
    hTrueSecTrackLength   = tfs->make<TH1D>("hTrueSecTrackLength",    "", 20, 0,    250);
    hRecoSecTrackLength   = tfs->make<TH1D>("hRecoSecTrackLength",    "", 20, 0,    250);
    hRecoCrSecTrackAngle  = tfs->make<TH1D>("hRecoCrSecTrackAngle",   "", 20, 0,    3.2);
    hRecoCrSecTrackLength = tfs->make<TH1D>("hRecoCrSecTrackLength",  "", 20, 0,    250);
    
    hTrueMuThetaPhi             = tfs->make<TH2D>("hTrueMuThetaPhi",            "", 10, 0, 3.2, 10, -3.2, 3.2);
    hRecoMuThetaPhi             = tfs->make<TH2D>("hRecoMuThetaPhi",            "", 10, 0, 3.2, 10, -3.2, 3.2);
    hTrueCrThetaPhi             = tfs->make<TH2D>("hTrueCrThetaPhi",            "", 10, 0, 3.2, 10, -3.2, 3.2);
    hRecoCrThetaPhi             = tfs->make<TH2D>("hRecoCrThetaPhi",            "", 10, 0, 3.2, 10, -3.2, 3.2);
    hTrueSecTrackLengthAngle    = tfs->make<TH2D>("hTrueSecTrackLengthAngle",   "", 10, 0, 250, 10, 0,    3.2);
    hRecoSecTrackLengthAngle    = tfs->make<TH2D>("hRecoSecTrackLengthAngle",   "", 10, 0, 250, 10, 0,    3.2);
    hRecoCrSecTrackLengthAngle  = tfs->make<TH2D>("hRecoCrSecTrackLengthAngle", "", 10, 0, 250, 10, 0,    3.2);

    // Initial output
    if(fVerbose) std::cout<<"----------------- Pandora Removal Ana Module -------------------"<<std::endl;

  }// PandoraAna::beginJob()


  void PandoraAna::analyze(const art::Event& event)
  {

    // Fetch basic event info
    //if(fVerbose){
      std::cout<<"============================================"<<std::endl
               <<"Run = "<<event.run()<<", SubRun = "<<event.subRun()<<", Event = "<<event.id().event()<<std::endl
               <<"============================================"<<std::endl;
    //}

    PFParticleHandle pfParticleHandle;
    event.getByLabel(fPandoraLabel, pfParticleHandle);

    if( !pfParticleHandle.isValid() ){
      std::cout<<"Failed to find the PFParticles."<<std::endl;
      return;
    }

    PFParticleIdMap pfParticleMap;
    this->GetPFParticleIdMap(pfParticleHandle, pfParticleMap);

    std::vector<art::Ptr<recob::PFParticle>> crParticles;
    std::vector<art::Ptr<recob::PFParticle>> nuParticles;

    this->GetFinalStatePFParticleVectors(pfParticleMap, crParticles, nuParticles);

    std::vector<art::Ptr<recob::Track>> tracks;
    std::vector<art::Ptr<recob::Shower>> showers;
    this->CollectTracksAndShowers(nuParticles, pfParticleHandle, event, tracks, showers);

    std::vector<art::Ptr<recob::Track>> crTracks;
    std::vector<art::Ptr<recob::Shower>> crShowers;
    this->CollectTracksAndShowers(crParticles, pfParticleHandle, event, crTracks, crShowers);

    art::FindManyP< recob::Track > pfPartToTrackAssoc(pfParticleHandle, event, fTpcTrackModuleLabel);

    //----------------------------------------------------------------------------------------------------------
    //                                          TRUTH MATCHING
    //----------------------------------------------------------------------------------------------------------

    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    // Retrieve all the truth info in the events
    auto particleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fSimModuleLabel);

    // Loop over true particles in readout window
    std::map<int, simb::MCParticle> particles;
    std::vector<int> lepParticleIds;
    std::vector<int> nuParticleIds;

    std::map<double, std::vector<simb::MCParticle>> nuTruth;
    std::vector<simb::MCParticle> crTruth;

    std::vector<double> usedNuVtx;

    for (auto const& particle: (*particleHandle)){

      int partId = particle.TrackId();
      particles[partId] = particle;
      art::Ptr<simb::MCTruth> truth = pi_serv->TrackIdToMCTruth_P(partId);
      simb::MCNeutrino mcNu = truth->GetNeutrino();
      bool isNuMuCC = false;
      if(std::abs(mcNu.Lepton().PdgCode())==13 && mcNu.CCNC() == simb::kCC) isNuMuCC = true;
      if(truth->Origin() == simb::kBeamNeutrino){
        geo::Point_t vtx;
        vtx.SetX(mcNu.Nu().Vx()); vtx.SetY(mcNu.Nu().Vy()); vtx.SetZ(mcNu.Nu().Vz());


        double time = particle.T() * 1e-3; // [us]
        if(fVerbose && particle.Mother()==0 && isNuMuCC) 
          std::cout<<"Nu VTX = "<<vtx<<" ID = "<<partId<<" pdg = "<<particle.PdgCode()<<" time = "<<time
                   <<" length = "<<particle.Trajectory().TotalLength()
                   <<" start = ("<<particle.Vx()<<", "<<particle.Vy()<<", "<<particle.Vz()<<") end = ("
                   <<particle.EndX()<<", "<<particle.EndY()<<", "<<particle.EndZ()<<")\n";
        if(!CosmicRemovalUtils::InFiducial(vtx, 0, 0)) continue;
        nuParticleIds.push_back(partId);

        if(particle.Mother() == 0){
          nuTruth[vtx.X()].push_back(particle);
        }

        if(isNuMuCC && std::find(usedNuVtx.begin(), usedNuVtx.end(), vtx.X())==usedNuVtx.end()){ 
          std::cout<<"\n----> NuMuCC! id = "<<mcNu.Nu().TrackId()<<" VTX = "<<vtx<<" time = "<<time<<"\n";
          nNuMuCC++;
          usedNuVtx.push_back(vtx.X());
        }

        if(std::abs(particle.PdgCode())==13 && particle.Mother()==0){ 
          std::cout<<"--> Lepton ID = "<<partId<<" pdg = "<<particle.PdgCode()<<" time = "<<time
                   <<" length = "<<particle.Trajectory().TotalLength()
                   <<" start = ("<<particle.Vx()<<", "<<particle.Vy()<<", "<<particle.Vz()<<") end = ("
                   <<particle.EndX()<<", "<<particle.EndY()<<", "<<particle.EndZ()<<")\n";
          lepParticleIds.push_back(partId);
          nMu++;
        }

      }
      //Check enters TPC and is muon
      if(truth->Origin() == simb::kCosmicRay){
        if(particle.Mother() == 0 && std::abs(particle.PdgCode()) == 13 && particle.Trajectory().TotalLength() > 500){
          crTruth.push_back(particle);
        }
      }
    }

    auto nutruthHandle = event.getValidHandle<std::vector<simb::MCTruth>>("generator");
    for (auto const& nutruth : (*nutruthHandle)){
      if(nutruth.Origin() != simb::kBeamNeutrino) continue;
      simb::MCNeutrino nu = nutruth.GetNeutrino();
      geo::Point_t vtx;
      vtx.SetX(nu.Nu().Vx()); vtx.SetY(nu.Nu().Vy()); vtx.SetZ(nu.Nu().Vz());
      TVector3 vert(vtx.X(), vtx.Y(), vtx.Z());
      std::vector<simb::MCParticle> parts = nuTruth[vtx.X()];
      for(size_t i = 0; i < parts.size(); i++){
        simb::MCParticle part = parts[i];
        if(std::abs(part.PdgCode())==13 && part.Mother()==0){
          //
          std::pair<TVector3, TVector3> se = truthAlg.TpcCrossPoints(part);
          TVector3 end(part.EndX(), part.EndY(), part.EndZ());
          TVector3 diff = end - vert;
          hTrueMuLength->Fill(truthAlg.TpcLength(part));
          hTrueMuE->Fill(part.E());
          hTrueMuTheta->Fill(diff.Theta());
          hTrueMuPhi->Fill(diff.Phi());
          hTrueMuThetaPhi->Fill(diff.Theta(), diff.Phi());
          double secTrackLength = -99999;
          double secTrackAngle = -99999;
          for(size_t j = 0; j < parts.size(); j++){
            simb::MCParticle part2 = parts[j];
            int pdg = std::abs(part2.PdgCode());
            if(j!=i && (pdg==13||pdg==211||pdg==2212) && part2.Mother()==0){
              std::pair<TVector3, TVector3> se2 = truthAlg.TpcCrossPoints(part2);
              double len = truthAlg.TpcLength(part2);
              if(len > secTrackLength){
                secTrackLength = len;
                TVector3 end2(part2.EndX(), part2.EndY(), part2.EndZ());
                TVector3 diff2 = end2 - vert;
                secTrackAngle = diff.Angle(diff2);
              }
            }
          }
          if(secTrackLength != -99999){
            hTrueSecTrackLength->Fill(secTrackLength);
            hTrueSecTrackAngle->Fill(secTrackAngle);
            hTrueSecTrackLengthAngle->Fill(secTrackLength, secTrackAngle);
          }
        }
      }
    }

    for(size_t i = 0; i < crTruth.size(); i++){
      simb::MCParticle part = crTruth[i];
      if(std::abs(part.PdgCode())==13 && part.Mother()==0){
        std::pair<TVector3, TVector3> se = truthAlg.TpcCrossPoints(part);
        TVector3 start = se.first;
        TVector3 end = se.second;
        if(start.X()==end.X() && start.Y()==end.Y() && start.Z()==start.Z()) continue;
        TVector3 diff = end - start;
        if(end.Y() > start.Y()) diff = start - end;
        hTrueCrLength->Fill(truthAlg.TpcLength(part));
        hTrueCrTheta->Fill(diff.Theta());
        hTrueCrPhi->Fill(diff.Phi());
        hTrueCrThetaPhi->Fill(diff.Theta(), diff.Phi());
      }
    }

    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTpcTrackModuleLabel);
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTpcTrackModuleLabel);

    auto tpcShowerHandle = event.getValidHandle<std::vector<recob::Shower>>(fShowerModuleLabel);
    art::FindManyP<recob::Hit> findManyHitsShower(tpcShowerHandle, event, fShowerModuleLabel);

    art::FindManyP<anab::Calorimetry> findManyCalo(tpcTrackHandle, event, fCaloModuleLabel);

    std::vector<int> usedLepIds;
    std::vector<int> nuRecoTrackIds;
    std::vector<int> nuRecoShowerIds;
    std::vector<int> crRecoTrackIds;
    std::vector<int> crRecoShowerIds;

    for(auto const& track : tracks){
      std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(track->ID());
      int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
      nuRecoTrackIds.push_back(trueId);
      if(std::find(lepParticleIds.begin(), lepParticleIds.end(), trueId) != lepParticleIds.end()){
        std::cout<<"-> Track ID = "<<track->ID()<<", true ID = "<<trueId<<" length = "<<track->Length()<<" start = "<<track->Vertex()<<" end = "<<track->End()<<"\n";
        if(std::find(usedLepIds.begin(), usedLepIds.end(), trueId) == usedLepIds.end()){
          nNuIsMu++;
          usedLepIds.push_back(trueId);
        }
      }
      else{
        std::cout<<"Not Mu: Track ID = "<<track->ID()<<", true ID = "<<trueId<<" length = "<<track->Length()<<" start = "<<track->Vertex()<<" end = "<<track->End()<<"\n";
      } 
    }

    for(auto const& shower : showers){
      std::vector<art::Ptr<recob::Hit>> hits = findManyHitsShower.at(shower->ID());
      std::cout<<"Hits size = "<<hits.size()<<"\n";
      int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
      nuRecoShowerIds.push_back(trueId);
    }
/*
    for(auto const& shower : crShowers){
      std::vector<art::Ptr<recob::Hit>> hits = findManyHitsShower.at(shower->ID());
      int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
      crRecoShowerIds.push_back(trueId);
    }
*/
    for(auto const& crTrack : crTracks){
      std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(crTrack->ID());
      int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
      crRecoTrackIds.push_back(trueId);
      if(std::find(lepParticleIds.begin(), lepParticleIds.end(), trueId) != lepParticleIds.end()){
        std::cout<<"-> Track ID = "<<crTrack->ID()<<", true ID = "<<trueId<<" length = "<<crTrack->Length()<<" start = "<<crTrack->Vertex()<<" end = "<<crTrack->End()<<"\n";
        nCrIsMu++;
      }
    }

    for(size_t i = 0; i < lepParticleIds.size(); i++){
      int lepId = lepParticleIds[i];
      bool found = false;
      std::cout<<"True particle "<<lepId<<":\n";
      if(std::find(nuRecoTrackIds.begin(), nuRecoTrackIds.end(), lepId) != nuRecoTrackIds.end()){
        found = true;
        std::cout<<"ID AS NU TRACK\n";
      }
      /*if(std::find(nuRecoShowerIds.begin(), nuRecoShowerIds.end(), lepId) != nuRecoShowerIds.end()){
        found = true;
        std::cout<<"ID AS NU SHOWER\n";
      }*/
      if(std::find(crRecoTrackIds.begin(), crRecoTrackIds.end(), lepId) != crRecoTrackIds.end()){
        found = true;
        std::cout<<"ID AS CR TRACK\n";
      }
      /*if(std::find(crRecoShowerIds.begin(), nuRecoShowerIds.end(), lepId) != nuRecoShowerIds.end()){
        found = true;
        std::cout<<"ID AS CR SHOWER\n";
      }*/
      if(!found){
        std::cout<<"MISSED!\n";
      }
    }

    for (PFParticleIdMap::const_iterator it = pfParticleMap.begin(); it != pfParticleMap.end(); ++it){

      const art::Ptr<recob::PFParticle> pParticle(it->second);
      // Only look for primary particles
      if (!pParticle->IsPrimary()) continue;

      // Check if this particle is identified as the neutrino
      const int pdg(pParticle->PdgCode());
      const bool isNeutrino(std::abs(pdg) == pandora::NU_E || std::abs(pdg) == pandora::NU_MU || std::abs(pdg) == pandora::NU_TAU);
      //Find neutrino pfparticle
      //if(!isNeutrino) continue;

      //Get the tracks associated with thedaughters
      bool isTrueNuMu = false;
      bool isTrueNu = false;
      std::vector<recob::Track> nuTracks;
      for (const size_t daughterId : pParticle->Daughters()){

        art::Ptr<recob::PFParticle> pParticle = pfParticleMap.at(daughterId);

        const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pParticle.key()));
        if(associatedTracks.size() != 1) continue;

        //If cosmic ID then  veto pfparticle
        recob::Track tpcTrack = *associatedTracks.front();
        std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());

        int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
        if(std::find(lepParticleIds.begin(), lepParticleIds.end(), trueId) != lepParticleIds.end()) isTrueNuMu = true;
        if(std::find(nuParticleIds.begin(), nuParticleIds.end(), trueId) != nuParticleIds.end()) isTrueNu = true;
        nuTracks.push_back(tpcTrack);
      }

      if(!isNeutrino){ 
        nPandCr++;
        if(isTrueNu) nNuPandCut++;
        else nCrPandCut++;
        if(isTrueNuMu) nNuMuPandCut++;
        continue;
      }
      nPandNu++;
      if(isTrueNu) nNuPandKeep++;
      else nCrPandKeep++;
      if(isTrueNuMu) nNuMuPandKeep++;

      std::sort(nuTracks.begin(), nuTracks.end(), [](auto& left, auto& right){
                return left.Length() > right.Length();});
      if(nuTracks.size() == 0){
        continue;
      }

      TVector3 pandoraVtx = nuTracks[0].Vertex<TVector3>();

      if(nuTracks.size() > 1){
        // Get vertices of two longest tracks and see if they match
        size_t secTrack = 99999;
        TVector3 start = nuTracks[0].Vertex<TVector3>();
        TVector3 end = nuTracks[0].End<TVector3>();
        TVector3 diff1;
        TVector3 diff2;
        for(size_t i = 1; i < nuTracks.size(); i++){
          TVector3 start2 = nuTracks[i].Vertex<TVector3>();
          TVector3 end2 = nuTracks[i].End<TVector3>();
          double minDist = 5;
          if((start - start2).Mag() < minDist){
            secTrack = i;
            diff1 = end - start;
            diff2 = end2 - start2;
          }
          if(secTrack != 99999) continue;
        }
        // 
        if(isTrueNuMu && secTrack!=99999){
          
          bool nuFound = false;
          for(size_t i = 0; i < nuTracks.size(); i++){
            recob::Track nuTrack = nuTracks[i];
            std::cout<<"----------------------------------->LENGTH = "<<nuTrack.Length();
            std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(nuTrack.ID());
            int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
            if(std::find(lepParticleIds.begin(), lepParticleIds.end(), trueId) != lepParticleIds.end() && !nuFound){
              simb::MCParticle part = particles[trueId];
              nuFound = true;
              std::cout<<" IS NU MU!";
              double muLen = truthAlg.TpcLength(part);
              hMuLengthDiff->Fill(nuTrack.Length()-muLen);
              art::Ptr<simb::MCTruth> truth = pi_serv->TrackIdToMCTruth_P(trueId);
              simb::MCParticle nu = truth->GetNeutrino().Nu();
              TVector3 trueVtx(nu.Vx(), nu.Vy(), nu.Vz());
              hMuVertexDiff->Fill((pandoraVtx-trueVtx).Mag());
              hTrueMuLengthReco->Fill(muLen);
              hTrueMuEReco->Fill(part.E());
              hRecoLongLength->Fill(nuTrack.Length());
              hRecoLongTheta->Fill((nuTrack.End<TVector3>()-nuTrack.Vertex<TVector3>()).Theta());
              hRecoLongPhi->Fill((nuTrack.End<TVector3>()-nuTrack.Vertex<TVector3>()).Phi());
            }
            std::cout<<"\n";
          }

          hRecoLongLength->Fill(nuTracks[0].Length());
          hRecoLongTheta->Fill(diff1.Theta());
          hRecoLongPhi->Fill(diff1.Phi());
          hRecoSecTrackAngle->Fill(diff1.Angle(diff2));
          hRecoSecTrackLength->Fill(nuTracks[secTrack].Length());
          hRecoMuThetaPhi->Fill(diff1.Theta(), diff1.Phi());
          hRecoSecTrackLengthAngle->Fill(nuTracks[secTrack].Length(), diff1.Angle(diff2));
          continue;
        }
        else if(!isTrueNu && secTrack != 99999){
          //
          hRecoCrLength->Fill(nuTracks[0].Length());
          hRecoCrTheta->Fill(diff1.Theta());
          hRecoCrPhi->Fill(diff1.Phi());
          hRecoCrSecTrackAngle->Fill(diff1.Angle(diff2));
          hRecoCrSecTrackLength->Fill(nuTracks[secTrack].Length());
          hRecoCrThetaPhi->Fill(diff1.Theta(), diff1.Phi());
          hRecoCrSecTrackLengthAngle->Fill(nuTracks[secTrack].Length(), diff1.Angle(diff2));
          continue;
        }
      }
      // Treat as if only one track from vertex
      TVector3 diff = nuTracks[0].End<TVector3>() - nuTracks[0].Vertex<TVector3>();
      if(isTrueNuMu){
        std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(nuTracks[0].ID());
        int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
        if(std::find(lepParticleIds.begin(), lepParticleIds.end(), trueId) != lepParticleIds.end()){
          simb::MCParticle part = particles[trueId];
          double muLen = truthAlg.TpcLength(part);
          hMuLengthDiff->Fill(nuTracks[0].Length() - muLen);
          art::Ptr<simb::MCTruth> truth = pi_serv->TrackIdToMCTruth_P(trueId);
          simb::MCParticle nu = truth->GetNeutrino().Nu();
          TVector3 trueVtx(nu.Vx(), nu.Vy(), nu.Vz());
          hMuVertexDiff->Fill((pandoraVtx-trueVtx).Mag());
          hTrueMuLengthReco->Fill(muLen);
          hTrueMuEReco->Fill(part.E());
        }
        hRecoMuLength->Fill(nuTracks[0].Length());
        hRecoMuTheta->Fill(diff.Theta());
        hRecoMuPhi->Fill(diff.Phi());
        hRecoMuThetaPhi->Fill(diff.Theta(), diff.Phi());
      }
      else if(!isTrueNu){
        //
        hRecoCrLength->Fill(nuTracks[0].Length());
        hRecoCrTheta->Fill(diff.Theta());
        hRecoCrPhi->Fill(diff.Phi());
        hRecoCrThetaPhi->Fill(diff.Theta(), diff.Phi());
      }
    }
    
    //----------------------------------------------------------------------------------------------------------
    //                                          NUMU SELECTION
    //----------------------------------------------------------------------------------------------------------

  } // PandoraAna::analyze()


  void PandoraAna::endJob(){

    std::cout<<"Number of NuMuCC interactions in TPC = "<<nNuMuCC<<"\n"
             <<"Total true mu particles     = "<<nMu<<"\n"
             <<"Total pandora nu slices     = "<<nPandNu<<"\n"
             <<"Number of true cr           = "<<nCrPandKeep<<"\n"
             <<"Number of true nu           = "<<nNuPandKeep<<"\n"
             <<"Number of true numu         = "<<nNuMuPandKeep<<"\n"
             <<"Mu in neutrino collection   = "<<nNuIsMu<<"\n"
             <<"-----------------------------------------------------\n"
             <<"Total pandora cr slices     = "<<nPandCr<<"\n"
             <<"Number of true cr           = "<<nCrPandCut<<"\n"
             <<"Number of true nu           = "<<nNuPandCut<<"\n"
             <<"Number of true numu         = "<<nNuMuPandCut<<"\n"
             <<"Mu in cosmic collection     = "<<nCrIsMu<<"\n";
    
    std::ofstream myfile;
    myfile.open("results.txt");
    myfile<<nNuMuCC<<","<<nMu<<","<<nPandNu<<","<<nCrPandKeep<<","<<nNuPandKeep<<","<<nNuMuPandKeep<<","<<nNuIsMu<<","<<nPandCr<<","<<nCrPandCut<<","<<nNuPandCut<<","<<nNuMuPandCut<<","<<nCrIsMu;
    myfile.close();

  } // PandoraAna::endJob()

  void PandoraAna::GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap){
      for (unsigned int i = 0; i < pfParticleHandle->size(); ++i){
          const art::Ptr<recob::PFParticle> pParticle(pfParticleHandle, i);
          if (!pfParticleMap.insert(PFParticleIdMap::value_type(pParticle->Self(), pParticle)).second){
              std::cout << "  Unable to get PFParticle ID map, the input PFParticle collection has repeat IDs!" <<"\n";
          }
      }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void PandoraAna::PrintOutScores(const art::Event &evt, const PFParticleHandle &pfParticleHandle) const{

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
      
  void PandoraAna::GetFinalStatePFParticleVectors(const PFParticleIdMap &pfParticleMap, PFParticleVector &crParticles, PFParticleVector &nuParticles){

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

          // Add the daughters of the neutrino PFParticle to the nuPFParticles vector
          for (const size_t daughterId : pParticle->Daughters()){
              if (pfParticleMap.find(daughterId) == pfParticleMap.end())
                  std::cout << "  Invalid PFParticle collection!" <<"\n";

              nuParticles.push_back(pfParticleMap.at(daughterId));
          }
      }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------
      
  void PandoraAna::CollectTracksAndShowers(const PFParticleVector &particles, const PFParticleHandle &pfParticleHandle, const art::Event &evt, TrackVector &tracks, ShowerVector &showers)
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


  DEFINE_ART_MODULE(PandoraAna)
} // namespace sbnd

