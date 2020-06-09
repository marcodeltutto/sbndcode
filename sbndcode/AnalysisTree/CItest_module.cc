////////////////////////////////////////////////////////////////////////
// Class:       CItest
// Module Type: analyzer
// File:        CItest_module.cc
//
////////////////////////////////////////////////////////////////////////
#ifndef CItest_Module
#define CItest_Module

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "sbndcode/CRT/CRTProducts/CRTData.hh"
#include "sbndcode/CRT/CRTProducts/CRTHit.hh"
#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/CRTHitRecoAlg.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

// ROOT includes
#include "TTree.h"
#include "TTimeStamp.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGeoManager.h"

// C++ Includes
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>

const int kMaxHits       = 30000; ///< maximum number of hits

class CItest : public art::EDAnalyzer {
public:
  explicit CItest(fhicl::ParameterSet const & p);
  virtual ~CItest();

  // This method is called once, at the start of the job. In this
  // example, it will define the histograms and n-tuples we'll write.
  void beginJob();

  // This method is called once, at the start of each run. It's a
  // good place to read databases or files that may have
  // run-dependent information.
  // void beginRun(const art::Run& run);

  // This method reads in any parameters from the .fcl files. This
  // method is called 'reconfigure' because it might be called in the
  // middle of a job; e.g., if the user changes parameter values in an
  // interactive event display.
  void reconfigure(fhicl::ParameterSet const& pset);

  // The analysis routine, called once per event. 
  void analyze (const art::Event& evt); 

private:

  int _run;                            ///< The run number
  int _subrun;                         ///< The subrun number
  int _event;                          ///< The event number
  double _evttime;                     ///< The event time number
  int _t0;                             ///< The t0

  int    _nhits;                       ///< Number of reco hits in the event


  std::string fHitsModuleLabel;     ///< Label for Hit dataproduct (to be set via fcl)


  geo::GeometryCore const* fGeometryService;
  // detinfo::DetectorClocks const* fDetectorClocks;
  // detinfo::DetectorProperties const* fDetectorProperties;
  // detinfo::ElecClock fTrigClock;
  art::ServiceHandle<geo::AuxDetGeometry> fAuxDetGeoService;
  const geo::AuxDetGeometry* fAuxDetGeo;
  const geo::AuxDetGeometryCore* fAuxDetGeoCore;

  TH1D *hColl,*hIndFirst,*hIndMiddle;
  TH1D *hCollAmp,*hInd0Amp,*hInd1Amp;
  TH1D *hCollArea,*hInd0Area,*hInd1Area;
  TH1D *hNoiseCollArea,*hNoiseInd0Area,*hNoiseInd1Area;
  TH1D *hNoiseCollAmp,*hNoiseInd0Amp,*hNoiseInd1Amp;

};


CItest::CItest(fhicl::ParameterSet const& pset)
: EDAnalyzer(pset)
{

  fGeometryService = lar::providerFrom<geo::Geometry>();
  // fDetectorClocks = lar::providerFrom<detinfo::DetectorClocksService>();
  // fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
  // fTrigClock = fDetectorClocks->TriggerClock();
  fAuxDetGeo = &(*fAuxDetGeoService);
  fAuxDetGeoCore = fAuxDetGeo->GetProviderPtr();

  // Read in the parameters from the .fcl file.
  this->reconfigure(pset);
}

void CItest::reconfigure(fhicl::ParameterSet const& p)
{  


  fHitsModuleLabel     = p.get<std::string>("HitsModuleLabel","gaushit");

}


CItest::~CItest()
{
  // Clean up dynamic memory and other resources here.
}


void CItest::analyze(const art::Event& evt)
{

  // _run = evt.run();
  // _subrun = evt.subRun();
  // _event = evt.id().event();

  // _t0 = 0.;
  // t0 = detprop->TriggerOffset();  // units of TPC ticks

  //
  // Hits
  //
  art::Handle<std::vector<recob::Hit>> hitListHandle;
  std::vector<art::Ptr<recob::Hit>> hitlist;
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle)) {
    art::fill_ptr_vector(hitlist, hitListHandle);
    _nhits = hitlist.size();
  }
  else {
    std::cout << "Failed to get recob::Hit data product." << std::endl;
    _nhits = 0;
  }

  int timecutlow = 2000; 
  int timecuthigh = 2700;
  int hitcount[2000][3] = {0} ;
  int wirelow[3] = {301, 501, 101};
  int wirehigh[3] = {1400, 1700, 1500};

  // first count hits per wire in window
  for (int i = 0; i<_nhits; ++i){
    geo::WireID wireid = hitlist[i]->WireID();
    int this_tpc       = wireid.TPC;
    int this_plane       = wireid.Plane;
    int this_wire        = wireid.Wire;
    int this_hittime = hitlist[i]->PeakTime();
    if (this_tpc==1 && this_hittime<=timecuthigh && this_hittime>=timecutlow) hitcount[this_wire][this_plane]++;
  }
      
  // Now fill histograms
  for (int i = 0; i<_nhits; ++i){
    geo::WireID wireid = hitlist[i]->WireID();
    int this_tpc       = wireid.TPC;
    int this_plane       = wireid.Plane;
    int this_wire        = wireid.Wire;
    int this_hittime = hitlist[i]->PeakTime();
    float this_hitcharge      = hitlist[i]->Integral();
    float this_hitph      = hitlist[i]->PeakAmplitude();
    if (this_tpc==0) {
      if (this_plane==0) {hNoiseInd0Area->Fill(this_hitcharge); hNoiseInd0Amp->Fill(this_hitph);}
      else if (this_plane==1){ hNoiseInd1Area->Fill(this_hitcharge); hNoiseInd1Amp->Fill(this_hitph);}
      else  { hNoiseCollArea->Fill(this_hitcharge);hNoiseCollAmp->Fill(this_hitph);}
    }      
    else { 
	if (this_hittime<=timecuthigh && this_hittime>=timecutlow && this_wire<=wirehigh[this_plane] && this_wire>=wirelow[this_plane]  ) {
	if (this_plane==0) hIndFirst->Fill((double)this_wire); 
	else if (this_plane==1) hIndMiddle->Fill((double)this_wire);
	else  hColl->Fill((double)this_wire);
	if ( hitcount[this_wire][this_plane]==1) { 
	  if (this_plane==0) { hInd0Area->Fill(this_hitcharge); hInd0Amp->Fill(this_hitph);}
	  else if (this_plane==1) { hInd1Area->Fill(this_hitcharge);hInd1Amp->Fill(this_hitph);}
	  else  { hCollArea->Fill(this_hitcharge);hCollAmp->Fill(this_hitph);}
	}
      }
    }
    // _hit_ph[i]          = hitlist[i]->PeakAmplitude();
    // _hit_width[i]       = hitlist[i]->RMS();
  }//  loop over hits

}

 void CItest::beginJob()
 {
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  //  fTree = tfs->make<TTree>("hitdumpertree","analysis tree");
  hColl = tfs->make<TH1D>("WireColl","WireColl",1400,100.5,1500.5);
  hColl->GetXaxis()->SetTitle("wire number plane 2") ;
  hIndFirst = tfs->make<TH1D>("WireIndFirst","WireIndFirst",1100,300.5,1400.5);
  hIndFirst->GetXaxis()->SetTitle("wire number plane 0") ;
  hIndMiddle = tfs->make<TH1D>("WireIndMiddle","WireIndMiddle",1200,500.5,1700.5);
  hIndMiddle->GetXaxis()->SetTitle("wire number plane 1") ;
  hCollAmp = tfs->make<TH1D>("CollAmp","CollAmp",500,0.,200.);
  hCollAmp->GetXaxis()->SetTitle("muon hit amplitude plane 2") ;
  hInd0Amp = tfs->make<TH1D>("Ind0Amp","Ind0Amp",500,0.,200.);
  hInd0Amp->GetXaxis()->SetTitle("muon hit amplitude plane 0") ;
  hInd1Amp = tfs->make<TH1D>("Ind1Amp","Ind1Amp",500,0.,200.);
  hInd1Amp->GetXaxis()->SetTitle("muon hit amplitude plane 1") ;
  hNoiseCollAmp = tfs->make<TH1D>("NoiseCollAmp","NoiseCollAmp",100,0.,100.);
  hNoiseCollAmp->GetXaxis()->SetTitle("noise hit amplitude plane 2") ;
  hNoiseInd1Amp = tfs->make<TH1D>("NoiseInd1Amp","NoiseInd1Amp",100,0.,100.);
  hNoiseInd1Amp->GetXaxis()->SetTitle("noise hit amplitude plane 1") ;
  hNoiseInd0Amp = tfs->make<TH1D>("NoiseInd0Amp","NoiseInd0Amp",100,0.,100.);
  hNoiseInd0Amp->GetXaxis()->SetTitle("noise hit amplitude plane 0") ;
  hCollArea = tfs->make<TH1D>("CollArea","CollArea",500,0.,1000.);
  hCollArea->GetXaxis()->SetTitle("muon hit area plane 2") ;
  hInd0Area = tfs->make<TH1D>("Ind0Area","Ind0Area",500,0.,1000.);
  hInd0Area->GetXaxis()->SetTitle("muon hit area plane 0") ;
  hInd1Area = tfs->make<TH1D>("Ind1Area","Ind1Area",500,0.,1000.);
  hInd1Area->GetXaxis()->SetTitle("muon hit area plane 1") ;
  hNoiseCollArea = tfs->make<TH1D>("NoiseCollArea","NoiseCollArea",100,0.,200.);
  hNoiseCollArea->GetXaxis()->SetTitle("noise hit area plane 2") ;
  hNoiseInd1Area = tfs->make<TH1D>("NoiseInd1Area","NoiseInd1Area",100,0.,200.);
  hNoiseInd1Area->GetXaxis()->SetTitle("noise hit area plane 1") ;
  hNoiseInd0Area = tfs->make<TH1D>("NoiseInd0Area","NoiseInd0Area",100,0.,200.);
  hNoiseInd0Area->GetXaxis()->SetTitle("noise hit area plane 0") ;
}


DEFINE_ART_MODULE(CItest)

#endif // CItest_Module
