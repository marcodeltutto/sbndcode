////////////////////////////////////////////////////////////////////////
// Class:       AnalyseEvents
// Plugin Type: analyzer (art v3_03_01)
// File:        AnalyseEvents_module.cc
//
// Generated at Mon Nov 25 02:54:01 2019 by Ala Zglam using cetskelgen
// from cetlib version v3_08_00.
////////////////////////////////////////////////////////////////////////

//Larsoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Track.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardataobj/Simulation/SimChannel.h"

// Additional framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

// Root includes
#include <TTree.h>
#include <TLorentzVector.h>
#include "TProfile.h"

// C++ includes
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <fstream>

namespace MCtest{
  class AnalyseEvents;
}


class MCtest::AnalyseEvents : public art::EDAnalyzer {
  public:
    explicit AnalyseEvents(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    AnalyseEvents(AnalyseEvents const&) = delete;
    AnalyseEvents(AnalyseEvents&&) = delete;
    AnalyseEvents& operator=(AnalyseEvents const&) = delete;
    AnalyseEvents& operator=(AnalyseEvents&&) = delete;

    // Required functions.
    void analyze(art::Event const& e) override;

    // Selected optional functions.
    void beginJob() override;
    void endJob() override;

  private:

    // Declare member data here.
    // Output tree declaration
    TTree *fTree;
    //Output TProfile declaration
    std::vector<TProfile*> fPurityPH;
    std::vector<TProfile*> fEffPH;

    // Variables to access and analyse
    // ----------integral variables-------------
    int     fWire;

    // ----------Fractional variables-----------
    float   fCharge;
    double  PlaneUEfficiency;
    double  PlaneVEfficiency;
    double  PlaneYEfficiency;
    double  PlaneUPurity;
    double  PlaneVPurity;
    double  PlaneYPurity;
    double  feffxpurityu;
    double  feffxpurityv;
    double  feffxpurityy;
    //double  efficiencyU = 0;
    //double  purityU = 0;
    double  effxpurityU = 0;
    double  effxpurityV = 0;
    double  effxpurityY = 0;

    //-----------unsigned integral variables----
    unsigned int  fEventID;
    unsigned int  fNHits;
    unsigned int  fnoevents = 0;
    //unsigned int  fplane;
    //-----------std library variables---------
    std::string   fHitsModuleLabel;
    std::string   fLArGeantModuleLabel;
    std::string   fHitsTrackLable;
    std::string   fTrackModuleLabel;
    //std::vector<int> wirefficiency;
    //std::vector<double> ChannelEnergy;

    // Variables to fill the output tree with
    // art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;
    art::ServiceHandle<geo::Geometry const> geom;
    art::ServiceHandle<cheat::BackTrackerService const> backtracker;
};


MCtest::AnalyseEvents::AnalyseEvents(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  // Get the producer name for the Hit from the configuration fhicl file
  //       analysisConfig.fcl
  fHitsModuleLabel     = p.get<std::string>("HitsModuleLabel");
  fLArGeantModuleLabel = p.get<std::string>("LArGeantModuleLabel");
  fHitsTrackLable      = p.get<std::string> ("HitsTrackLable");
  fTrackModuleLabel    = p.get<std::string>("TrackModuleLabel");
}

void MCtest::AnalyseEvents::analyze(art::Event const& e)
{
  //declaring variables for each event.
  // ----------integral variables-------------
    //-------- unsigned variables------
    //unsigned int fnoevents;
    // ----------double variables------
  //double  efficiencyU;
  //double  purityU;
  //double  effxpurityU;

  //int NuSimChan = 0;  //Number of the channel ides that have energy deposited
  //-----------std library variables---------
  std::vector<double> eventefficiency;
  std::vector<double> PlanePurity;
  std::vector<unsigned int> simchannelids;
  std::vector<art::Ptr<recob::Hit> > initialtrackhitslargest;
  // Implementation of required member function here.
  //Define our event ID variable
  fEventID = e.id().event();

  // Initialise our counters for this event
  fNHits     = 0;
  size_t initialtrackhitsize = 0;
  
  // Access the SimChannel to get SimWire Information
  //Get the simchannel to find the IDEs deposited on them.
  art::Handle<std::vector<sim::SimChannel> > simChannelHandle;
  std::vector<art::Ptr<sim::SimChannel> > simchannels;
  if(e.getByLabel(fLArGeantModuleLabel,simChannelHandle)) // make sure the handle is valid
  {art::fill_ptr_vector(simchannels, simChannelHandle);}

  // Access The Hits from Gaus
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector< art::Ptr<recob::Hit> > hits;
  if(e.getByLabel(fHitsModuleLabel, hitListHandle)) // Make sure the handle is valid
    art::fill_ptr_vector(hits, hitListHandle); // Fill the vector with art::Ptr hits
  if(!hits.size()) return; // If there are no reconstructed particles, skip the event
  fNHits = hits.size();
  //recob Track for the hits
  art::Handle<std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracks;
  if(e.getByLabel(fTrackModuleLabel, trackListHandle))
  {art::fill_ptr_vector(tracks, trackListHandle);}
  // Looking for the hists that associated with initial track
  art::FindManyP<recob::Hit> fmith(trackListHandle, e, fHitsTrackLable);
  //Git the Hits
  for(auto const& initialTrack: tracks){
    std::vector<art::Ptr<recob::Hit> > initialtrackhits = fmith.at(initialTrack.key());
    //Looking for the largest track
    if (initialtrackhits.size() > initialtrackhitsize){
      initialtrackhitslargest = fmith.at(initialTrack.key());
      initialtrackhitsize = initialtrackhits.size();
    }
  }
  if (initialtrackhitsize == 0) {
    //std::cout << " This Event has no hits which should not happen " << std::endl;
    return;
  }

  //Maps that created
  //Create a map to store the plane and the wireids that got hits
  std::map<geo::PlaneID, std::vector<art::Ptr<recob::Hit> > > planewhmap;
  // define the map of simChannel ides and energy
  //std::map<unsigned int, double> chan_eg_id_map;
  //**********************************************************************************
  ///////////////////////////////////////////////////////////////////////////////////
  //                          Looping over Hits                                   //
  /////////////////////////////////////////////////////////////////////////////////

  for(size_t nuHit = 0; nuHit < hitListHandle->size(); ++nuHit)
  //for(auto const& hit: initialtrackhitslargest)
  {
    // *** Finding Channel associated with the hit ***
    art::Ptr<recob::Hit> hit(hitListHandle, nuHit);
    //---------------
    fCharge = hit->Integral();
    fWire   = hit->WireID().Wire;
    //-----------------------
    geo::WireID const wire = hit->WireID();
    geo::PlaneID planeid = wire.asPlaneID();//Looking for all information from gwire "C:0 T:1 P:2"
    //-----------------------
    //Fill the map
    planewhmap[planeid].push_back(hit);
    //planewhmap[(hit->WireID().Plane)].push_back(hit);
    //std::cout << "WireID = " << alwireid << std::endl;
    //-------------------
  }//close the loop over the hits
  // Looping inside the map
  // --------------------
  for(auto mapiter : planewhmap) //loop on plane
  {
    unsigned int uptwc = 0;  // total number of wires between first and the last hit
    unsigned int uphwc = 0;  // number of wires that got hit between first and the last hit
    unsigned int hitc = 0; // the total count of hit from recob
    unsigned int thitc   = 0; // Number of hit that have true energy
    //unsigned int maxwire_id = 999999; // The last wire have been hit
    //unsigned int minwire_id = 0;     // The first wire have been hit
    // unsigned int thitsimc = 0; // the number of the hit that have true energy from sim channel
    auto MapPlaneId = mapiter.first; //The Plane number from the map "planewhmap"
    unsigned int tpcid = MapPlaneId.TPC;
    unsigned int MapPlane = MapPlaneId.Plane; //
    //std::vector<float> EffPeakAmp; //A vector of PeakAmplitude()
    //std::vector<unsigned int> wireID;
    std::vector<unsigned int> HitChannelID;
    //ChannelEnergy.clear();
    //wirefficiency.clear();
    //
    //-------------------loop inside the vecotr on hits--------------------//
    for (auto const& hit: mapiter.second)
    {
      // Check if a hit is a noise hit or a truth hit
      ++hitc;
      float PeakAmp = hit->PeakAmplitude();
      //EffPeakAmp.push_back(PeakAmp);
      unsigned int purity;
      {
        std::vector<sim::TrackIDE> trackides = backtracker->HitToTrackIDEs(hit);
        if(trackides.size() != 0){
          ++thitc;
          unsigned int hchanid = hit->Channel();
          HitChannelID.push_back(hchanid);
          purity = 1;

        }//Backtracker checker
        else {purity = 0;}
      } // Check the Purity
      fPurityPH[MapPlane]->Fill(PeakAmp, purity, 1.);
    }// Close the looping over the hits
    //std::cout << "P= " << MapPlane <<" NumHitsWithTrueEnergy = " << uphwc << " NumHitsWithOutTru    eEnergy = " << uptwc << std::endl;

    //////////////////////////////////////////////////////
    //                Loop over the wires              //
    ////////////////////////////////////////////////////
    for(geo::WireID const& wire_id: geom->IterateWireIDs()){
      // Declare Wire ID
      //unsigned int WirePlane = wire_id.Plane;
      unsigned int max_num_electrons = 0;
      unsigned int wireff;
      geo::PlaneID WirePlane = wire_id.asPlaneID();
      //------------------------------------------
      if(MapPlaneId == WirePlane /*&& WirePlane == 2*/){
        unsigned int WireChannelID = geom->PlaneWireToChannel(wire_id.Plane, wire_id.Wire, wire_id.TPC, wire_id.Cryostat);
        for(const art::Ptr<sim::SimChannel> &simchannel : simchannels){
          if((simchannel->Channel()) == WireChannelID){
            auto tdc_ide_map = simchannel->TDCIDEMap();
            for(auto const& tdc_pair : tdc_ide_map){
              auto const& ide_ne = tdc_pair.second;
              for(auto const& num_e : ide_ne){
                unsigned int num_electrons = num_e.numElectrons;
                if(num_electrons > max_num_electrons){
                  max_num_electrons = num_electrons;}
              }
            }
            if (max_num_electrons < 400) continue; 
            //std::cout << "P= " << t_num_electrons << std::endl;
            if(std::find(HitChannelID.begin(), HitChannelID.end(), WireChannelID)!=HitChannelID.end()){
              wireff = 1;
              //Count the number of Channels that have 1 efficiency
              ++uphwc;}  //uptwc = uptwc +1
            else {
              wireff = 0;
              //Count the number of Channels that have 0 efficiency
              ++uptwc; // uphwc = uphwc + 1
            }
            fEffPH[MapPlane]->Fill(max_num_electrons, wireff, 1.);
          }
        }//Close the if condition.
      }//Close the if statement that check if we are at the same plane
    }//close the looping over the wires
    //Fill a vector of the average of the efficiency per plane per event
    double ptwc = uphwc + uptwc;
    double aveEff = (double)uphwc/ptwc;
    if(tpcid==1){
    eventefficiency.push_back(aveEff);
    //Fill a vector of the average of the efficiency per plane per event
    double avepurity = (double)thitc/hitc;
    PlanePurity.push_back(avepurity);
    }
    //std::cout << " aveEff =  " << aveEff << " avepurity = " << avepurity << std::endl;
    //------------------------
  }//Close the looping on the map for plane
  //------------------------------
  if(eventefficiency.size() == 3){
    ++fnoevents;
    PlaneUEfficiency = eventefficiency[0];
    PlaneVEfficiency = eventefficiency[1];
    PlaneYEfficiency = eventefficiency[2];
    PlaneUPurity = PlanePurity[0];
    PlaneVPurity = PlanePurity[1];
    PlaneYPurity = PlanePurity[2];
    feffxpurityu = (PlanePurity[0] * eventefficiency[0]);
    feffxpurityv = (PlanePurity[1] * eventefficiency[1]);
    feffxpurityy = (PlanePurity[2] * eventefficiency[2]);
    // Printing the average Efficiency, Purity, and EffxPurity
    //efficiencyU += PlaneUEfficiency;
    //purityU += PlaneUPurity;
    effxpurityU += feffxpurityu;
    effxpurityV += feffxpurityv;
    effxpurityY += feffxpurityy;

  }
  //std::cout << " U = " << PlaneUEfficiency << " V = " << PlaneVEfficiency << " Y = " << PlaneYEf    ficiency  << std::endl;
  //--------------------------------------
  // Printing the average Efficiency, Purity, and EffxPurity
  if(fEventID == 100){
   // double avefficiencyU = (efficiencyU/fnoevents);
   // double avefficiencyV = (efficiencyV/fnoevents);
   // double avefficiencyY = (efficiencyY/fnoevents);
   // double avepurityU = (purityU/fnoevents);
    double aveffxpurityU = (effxpurityU/fnoevents);
    //double aveffxpurityV = (effxpurityV/fnoevents);
    //double aveffxpurityY = (effxpurityY/fnoevents);
    //double aveffxpurity = ((aveffxpurityU + aveffxpurityV + aveffxpurityY)/3.0);
    std::cout << " ave EffxPurity= " << aveffxpurityU << std::endl;
    std::cout << " ***************************************************************" << std::endl;
    {
      // save the average efficiency in txt file
      //std::ofstream MeanEffU;
      //MeanEffU.open ("MeanEffU.txt");
      //MeanEffU << avefficiency;
      //MeanEffU.close();
      // save the average Purity in txt file
      //std::ofstream MeanPurityU;
      //MeanPurityU.open ("MeanPurityU.txt");
      //MeanPurityU << avepurityU;
      //MeanPurityU.close();
      // save the average effxPurity in txt file
      std::ofstream MeanEffxPurityU;
      MeanEffxPurityU.open ("MeanEffxPurityU.txt");
      MeanEffxPurityU << aveffxpurityU;
      MeanEffxPurityU.close();
    }
  }
  //--------------------------------------
  //Fill the output TTree with all relevant variables 
  fTree->Fill();
 
  } // close the main function

  void MCtest::AnalyseEvents::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    {
      // Creating TProfile for Purity and Pulse Height
      for( unsigned int MapPlane=0; MapPlane<3; ++MapPlane){
    fPurityPH.push_back(tfs->make<TProfile> (Form("fPurityPH_%d",MapPlane),
            Form("Purity-Electron Sample-uBooNE Noise (Plane %d); Pulse Hight ADC; Purity",MapPlane),
            100,0,100,0,1));
    fEffPH.push_back(tfs->make<TProfile> (Form("fEffPH_%d",MapPlane),
         Form("Efficiency-Electron Sample-uBooNE Noise (Plane %d); Pulse Hight (N electrons);Efficiecny",MapPlane),
          100,0,20000,0,1));

      }
    }

    // Implementation of optional member function here.
    // The TFileService is used to define the TTree and writing it to the output f
   //art::ServiceHandle<art::TFileService> tfs;
   fTree = tfs->make<TTree>("tree", "Analyser Output Tree");

    // Add Branches to our TTree
    fTree->Branch("eventID", &fEventID, "eventID/I");
    fTree->Branch("nHits", &fNHits, "nHits/I");
    fTree->Branch("Hitcharge", &fCharge, "Hitcharge/f");
    fTree->Branch("Wire", &fWire, "Wire/I");
    fTree->Branch("EfficiencyU", &PlaneUEfficiency, "EfficiencyU/D");
    fTree->Branch("EfficiencyV", &PlaneVEfficiency, "EfficiencyV/D");
    fTree->Branch("EfficiencyY", &PlaneYEfficiency, "EfficiencyY/D");
    fTree->Branch("Purity_U", &PlaneUPurity,"Purity_U/D");
    fTree->Branch("Purity_V", &PlaneVPurity,"Purity_V/D");
    fTree->Branch("Purity_Y", &PlaneYPurity,"Purity_Y/D");
    fTree->Branch("EffxPurity_U", &feffxpurityu,"EffxPurity_U/D");
    fTree->Branch("EffxPurity_V", &feffxpurityv,"EffxPurity_V/D");
    fTree->Branch("EffxPurity_Y", &feffxpurityy,"EffxPurity_Y/D");
  }

  void MCtest::AnalyseEvents::endJob()
  {
    // Implementation of optional member function here.
  }

  DEFINE_ART_MODULE(MCtest::AnalyseEvents)
