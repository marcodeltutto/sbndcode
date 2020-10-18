////////////////////////////////////////////////////////////////////////
// Class:       Playground
// Plugin Type: analyzer (art v3_05_01)
// File:        Playground_module.cc
//
// Generated at Wed Jul 15 15:08:53 2020 by Marco Del Tutto using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"


#include "TTree.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "larcore/Geometry/Geometry.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/Hit.h"

class Playground;


class Playground : public art::EDAnalyzer {
public:
  explicit Playground(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Playground(Playground const&) = delete;
  Playground(Playground&&) = delete;
  Playground& operator=(Playground const&) = delete;
  Playground& operator=(Playground&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:
  
  int _subtract_pedestal = false;

  detinfo::DetectorProperties const* _det_prop;
  geo::GeometryCore const* _geo_service;
  detinfo::LArProperties const* _lar_prop;

  // Declare member data here.
  TTree* tree;
  std::vector<double> _ides_x, _ides_y, _ides_z, _ides_num_electrons, _ides_energy;
  std::vector<float> _raw_digits;

  std::vector<float> _hit_time, _hit_wire, _hit_plane, _hit_integral, _hit_start_time, _hit_end_time;
};


Playground::Playground(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  _det_prop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  _geo_service = lar::providerFrom<geo::Geometry>();
  _lar_prop = lar::providerFrom<detinfo::LArPropertiesService>();

  art::ServiceHandle<art::TFileService> fs;
  tree = fs->make<TTree>("Tree","");

  tree->Branch("ides_x", "std::vector<double>", &_ides_x);
  tree->Branch("ides_y", "std::vector<double>", &_ides_y);
  tree->Branch("ides_z", "std::vector<double>", &_ides_z);
  tree->Branch("ides_num_electrons", "std::vector<double>", &_ides_num_electrons);
  tree->Branch("ides_energy", "std::vector<double>", &_ides_energy);
  tree->Branch("raw_digits", "std::vector<float>", &_raw_digits);
  tree->Branch("hit_time", "std::vector<float>", &_hit_time);
  tree->Branch("hit_wire", "std::vector<float>", &_hit_wire);
  tree->Branch("hit_plane", "std::vector<float>", &_hit_plane);
  tree->Branch("hit_integral", "std::vector<float>", &_hit_integral);
  tree->Branch("hit_start_time", "std::vector<float>", &_hit_start_time);
  tree->Branch("hit_end_time", "std::vector<float>", &_hit_end_time);

}

void Playground::analyze(art::Event const& e)
{
  // Implementation of required member function here.

  art::Handle< std::vector<sim::SimChannel> > simchannels_h;
  std::vector<art::Ptr<sim::SimChannel> > simchannels;
  if (e.getByLabel("largeant", simchannels_h)) {
    art::fill_ptr_vector(simchannels, simchannels_h);
  }

  _ides_x.clear();
  _ides_y.clear();
  _ides_z.clear();
  _ides_num_electrons.clear();
  _ides_energy.clear();

  /*
  for (auto & simch : *simchannels_h) {
    // get all the IDEs associated with this simch
    auto const& all_ide = simch.TrackIDsAndEnergies(0,19600);
    for (auto const& ide : all_ide){
      // std::cout << "3dpoint : [" << ide.x << ", " << ide.y << ", " << ide.z << "]" << std::endl;
      _ides_x.push_back(ide.x);
      _ides_y.push_back(ide.y);
      _ides_z.push_back(ide.z);
      _ides_num_electrons.push_back(ide.numElectrons);
      _ides_energy.push_back(ide.energy);
    }
  }
  */


  art::Handle< std::vector<raw::RawDigit> > raw_digit_h;
  std::vector<art::Ptr<raw::RawDigit> > raw_digits;
  if (e.getByLabel("daq", raw_digit_h)) {
    art::fill_ptr_vector(raw_digits, raw_digit_h);
  }

  size_t n_ticks = raw_digits.front()->ADCs().size();
  // _raw_digits.resize(n_ticks * _geo_service->Nwires(2, 0, 0), 0.);
  _raw_digits.resize(n_ticks * 500, 0.);

  for (auto const &rawdigit : raw_digits) {
      unsigned int ch  = rawdigit->Channel();
      float        ped = rawdigit->GetPedestal();

      std::vector<geo::WireID> widVec = _geo_service->ChannelToWire(ch);
      unsigned int wire = widVec[0].Wire;
      unsigned int plane = widVec[0].Plane;
      unsigned int tpc = widVec[0].TPC;
      unsigned int cryo = widVec[0].Cryostat;

      if (cryo != 0 || tpc != 0 || plane != 2) {
        continue;
      }

      if (wire < 100 || wire >= 600) continue;

      int offset = (wire - 100) * n_ticks;
      std::vector<float>::iterator startItr = _raw_digits.begin() + offset;
      
      float pedestal = 0;
      if (_subtract_pedestal) {
        pedestal = ped;
      }

      for(const auto& adcVal : rawdigit->ADCs()) {
        *startItr++ = adcVal - pedestal;
      }
  }



  std::cout << "number of ticks: " << n_ticks << std::endl;
  std::cout << "number of wires: " << _geo_service->Nwires(2, 0, 0) << std::endl;



  art::Handle<std::vector<recob::Hit>> hit_h;
  std::vector<art::Ptr<recob::Hit>> hits;
  if (e.getByLabel("gaushit", hit_h)) {
    art::fill_ptr_vector(hits, hit_h);
  }

  _hit_time.clear();
  _hit_wire.clear();
  _hit_plane.clear();
  _hit_integral.clear();
  _hit_start_time.clear();
  _hit_end_time.clear();

  for (auto const &hit : hits) {
    unsigned int plane = hit->WireID().Plane;
    unsigned int tpc = hit->WireID().TPC;
    unsigned int cryo = hit->WireID().Cryostat;
    
    if (cryo != 0 || tpc != 0 || plane != 2) {
      continue;
    }
    _hit_time.push_back(hit->PeakTime());
    _hit_wire.push_back(hit->WireID().Wire);
    _hit_plane.push_back(plane);
    _hit_integral.push_back(hit->Integral());
    _hit_start_time.push_back(hit->StartTick());
    _hit_end_time.push_back(hit->EndTick());
  }


  tree->Fill();



  // for(size_t sc = 0; sc < simchannels.size(); ++sc) {
  //   for (auto const& tdcide: simchannels[sc]->TDCIDEMap()) {
  //     std::vector<sim::IDE> const& idevec = tdcide.second;
  //     for(auto const& ide: idevec) {
  //       ides_x.push_back(ide.x);
  //       ides_y.push_back(ide.y);
  //       ides_z.push_back(ide.z);
  //       ides_num_electrons.push_back(ide.numElectrons);
  //       ides_energy.push_back(ide.energy);
  //       std::cout << "ide.x = " << ide.x << std::endl;
  //       tree->Fill();
  //     }
  //   }
  // }
}








DEFINE_ART_MODULE(Playground)
