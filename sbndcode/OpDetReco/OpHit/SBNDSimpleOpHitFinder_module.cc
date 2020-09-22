////////////////////////////////////////////////////////////////////////
// Class:       SBNDSimpleOpHitFinder
// Plugin Type: producer (art v3_05_01)
// File:        SBNDSimpleOpHitFinder_module.cc
//
// Generated at Wed Sep  9 23:12:35 2020 by Marco Del Tutto using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"

#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

#include <memory>
#include <numeric>
#include <fstream>

class SBNDSimpleOpHitFinder;


class SBNDSimpleOpHitFinder : public art::EDProducer {
public:
  explicit SBNDSimpleOpHitFinder(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SBNDSimpleOpHitFinder(SBNDSimpleOpHitFinder const&) = delete;
  SBNDSimpleOpHitFinder(SBNDSimpleOpHitFinder&&) = delete;
  SBNDSimpleOpHitFinder& operator=(SBNDSimpleOpHitFinder const&) = delete;
  SBNDSimpleOpHitFinder& operator=(SBNDSimpleOpHitFinder&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  void RunHitFinder(std::vector<raw::OpDetWaveform> const& waveform_v,
                    std::vector<recob::OpHit>& ophits,
                    geo::GeometryCore const& geo_service,
                    detinfo::DetectorClocksData const& clocks_data);
  bool FindPeak(std::vector<float>& wf, size_t& time_bin, float& area, float baseline_stdev);
  std::vector<int> PDNamesToList(std::vector<std::string>);

  std::string _input_module;
  std::vector< std::string > _input_labels;

  std::vector<std::string> _pd_to_use; ///< PDS to use (ex: "pmt", "barepmt")
  std::vector<int> _opch_to_use; ///< List of of opch (will be infered from _pd_to_use)

  int _baseline_sample; ///< Number of ticks used to estimate the baseline
  int _integration_window; ///< Number of ticks the pulse will be integrated for
  int _num_presample; ///< Number of ticks before maximum
  float _spe_area; ///< Single PhotoElectron area
  float _threshold; ///< Waveform threshold to claim a pulse (in ADC)
  int _dead_time; ///< Dead time in ticks after a pulse is claimed

  float _sampling_frequency; ///< Optical frequency of sampling

  opdet::sbndPDMapAlg _pds_map; ///< map for photon detector types

  int _n_wf_to_csvfile; ///< If greater than zero saves firsts waveforms with pedestal to csv file
  std::ofstream _csvfile; ///< To save waveforms
};


SBNDSimpleOpHitFinder::SBNDSimpleOpHitFinder(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
    _input_module = p.get<std::string>("InputModule");
    _input_labels = p.get<std::vector<std::string>>("InputLabels");

    _pd_to_use = p.get<std::vector<std::string>>("PD", _pd_to_use);
    _opch_to_use = this->PDNamesToList(_pd_to_use);

    _baseline_sample = p.get<int>("BaselineSample"); //in ticks
    _integration_window = p.get<int>("IntegrationWindow"); //in ticks
    _num_presample = p.get<int>("NumPreSample"); //in ticks
    _spe_area = p.get<float>("SPEArea");
    _threshold = p.get<float>("Threshold");
    _dead_time = p.get<int>("DeadTime"); //in ticks

    _n_wf_to_csvfile = p.get<int>("NWaveformsToFile", 0);
    if (_n_wf_to_csvfile > 0) {
      _csvfile.open ("wf_simpleophitfinder.csv", std::ofstream::out | std::ofstream::trunc);
      _csvfile << "n,time,wf,wf_ped_mean,wf_ped_rms" << std::endl;
    }

    produces<std::vector<recob::OpHit>>();
}

void SBNDSimpleOpHitFinder::produce(art::Event& e)
{
  std::unique_ptr<std::vector<recob::OpHit>> ophits (new std::vector<recob::OpHit>);

  auto const& geo_service (*lar::providerFrom< geo::Geometry >());
  auto const clocks_data = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  _sampling_frequency = clocks_data.OpticalClock().Frequency(); // MHz
  std::cout << "_sampling_frequency " << _sampling_frequency << std::endl;

  // Reserve a large enough array
  int totalsize = 0;
  for (auto label : _input_labels)
  {
    art::Handle<std::vector<raw::OpDetWaveform>> waveform_h;
    e.getByLabel(_input_module, label, waveform_h);
    if (!waveform_h.isValid()) continue; // Skip non-existent collections
    totalsize += waveform_h->size();
  }

  std::vector<raw::OpDetWaveform> waveform_v;
  waveform_v.reserve(totalsize);

  for (auto label : _input_labels)
  {
    art::Handle<std::vector<raw::OpDetWaveform>> waveform_h;
    e.getByLabel(_input_module, label, waveform_h);
    if (!waveform_h.isValid()) continue; // Skip non-existent collections

    for(auto const& wf : *waveform_h)
    {
      // If this channel is in the channel mask, ingore it
      // if ( fChannelMasks.find(wf.ChannelNumber()) != fChannelMasks.end() ) continue;
      // If this PDS in not in the list of PDS to use, ingore it
      if (std::find(_opch_to_use.begin(), _opch_to_use.end(), wf.ChannelNumber())
        == _opch_to_use.end()) continue;

      // if (wf.ChannelNumber() != 249) continue;

      waveform_v.push_back(wf);
    }
  }

  RunHitFinder(waveform_v, *ophits, geo_service, clocks_data);
  std::cout << "Done!" << std::endl;

  e.put(std::move(ophits));
}

void SBNDSimpleOpHitFinder::RunHitFinder(std::vector<raw::OpDetWaveform> const& waveform_v,
                                         std::vector<recob::OpHit>& ophits,
                                         geo::GeometryCore const& geo_service,
                                         detinfo::DetectorClocksData const& clocks_data) {

  float baseline_mean, baseline_stdev;

  int wf_saved = 0;

  for (auto const & w : waveform_v) {

    std::vector<float> wf(w.begin(), w.end());

    auto is_start = wf.begin();
    auto it_end = wf.begin() + _baseline_sample;

    baseline_mean = std::accumulate(is_start, it_end, 0.0);
    baseline_mean /= (float) _baseline_sample;


    for (auto & e : wf) {
      e = (e - baseline_mean) * -1;
    }

    float mean = 0;
    baseline_stdev = std::inner_product(is_start, it_end, is_start, 0.0);
    baseline_stdev = std::sqrt(baseline_stdev / (float) _baseline_sample - mean * mean);

    std::cout << "  ***** Baseline mean " << baseline_mean << ", std " << baseline_stdev << std::endl;

    if (wf_saved + 1 <= _n_wf_to_csvfile) {
      wf_saved ++;
      for (size_t i = 0; i < wf.size(); i++) {
        _csvfile << wf_saved-1 << "," << i << "," << wf[i] << "," << mean << "," << baseline_stdev << std::endl;
      }
    }

    size_t time_bin = 0;
    float area = 0.;

    std::vector<float> used_times;

    while (FindPeak(wf, time_bin, area, baseline_stdev)) {

      float time = w.TimeStamp() + (float) time_bin / _sampling_frequency;

      auto it = std::find(used_times.begin(), used_times.end(), time);
      if (it != used_times.end()) {
        mf::LogWarning("SBNDSimpleOpHitFinder")
          << "Got stuck trying to reconstruct OpHit for channel "
          << w.ChannelNumber() << std::endl;
        break;
      }

      recob::OpHit ophit(w.ChannelNumber(),   // channel
                         time,                // peaktime
                         time,                // peaktimeabs
                         1,                   // frame
                         1,                   // fwhm
                         area,                // area
                         0,                   // amplitude
                         area / _spe_area,    // PEs
                         3./4.);              // fast-to-total

      used_times.emplace_back(time);

      std::cout << "Found ophit for ch " << w.ChannelNumber()
                << " time stamp " << w.TimeStamp()
                << " time_bin " << time_bin
                << " at time " << time
                << " with PE " << area / _spe_area  << std::endl;

      ophits.emplace_back(ophit);
    }
  }
}


bool SBNDSimpleOpHitFinder::FindPeak(std::vector<float>& wf,
                                     size_t& time_bin,
                                     float& area,
                                     float baseline_stdev) {


  std::vector<float>::iterator max_element_it = std::max_element(wf.begin(), wf.end());

  auto amplitude = *max_element_it;

  if(amplitude < _threshold) return false; // stop if there's no more peaks

  time_bin = std::distance(wf.begin(), max_element_it);


  // it_end contains the iterator to the last element in the peak
  // where waveform is above threshold
  // auto it_end = std::find_if(max_element_it, wf.end(),
  //                             [&baseline_stdev](const double& x)->bool
  //                             {return x < baseline_stdev;} );

  // it_start contains the iterator to the first element in the peak
  // where waveform is above threshold
  auto it_start = std::find_if(std::make_reverse_iterator(max_element_it),
                               std::make_reverse_iterator(wf.begin()),
                               [&baseline_stdev](const double& x)->bool
                               {return x < baseline_stdev;} ).base();

  int presample = _num_presample;
  if (std::distance(wf.begin(), it_start) <= presample) {
    presample = std::distance(wf.begin(), it_start);
  }

  int integration_window = _integration_window;
  if (std::distance(it_start, wf.end()) <= integration_window) {
    integration_window = std::distance(it_start, wf.end());
  }

  int dead_time = _dead_time;
  if (std::distance(it_start, wf.end()) <= dead_time) {
    dead_time = std::distance(it_start, wf.end());
  }

  // Integrate for a predefined range only
  area = std::accumulate(it_start - presample, it_start + integration_window, 0.0);
  area = area / (_sampling_frequency / 1000.);

  // zero the peak so it won't be picked up again
  std::fill(it_start - presample, it_start + dead_time, 0.0);

  return true;

}

std::vector<int> SBNDSimpleOpHitFinder::PDNamesToList(std::vector<std::string> pd_names) {

  std::vector<int> out_ch_v;

  for (auto name : pd_names) {
    auto ch_v = _pds_map.getChannelsOfType(name);
    out_ch_v.insert(out_ch_v.end(), ch_v.begin(), ch_v.end());
  }

  return out_ch_v;

}


DEFINE_ART_MODULE(SBNDSimpleOpHitFinder)
