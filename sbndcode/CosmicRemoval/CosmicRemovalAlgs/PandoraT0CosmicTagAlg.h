#ifndef PANDORAT0COSMICTAGALG_H_SEEN
#define PANDORAT0COSMICTAGALG_H_SEEN


///////////////////////////////////////////////
// PandoraT0CosmicTagAlg.h
//
// Functions for pandora t0 cosmic tagger
// T Brooks (tbrooks@fnal.gov), November 2018
///////////////////////////////////////////////

#include "sbndcode/CosmicRemoval/CosmicRemovalUtils/CosmicRemovalUtils.h"

// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art_root_io/TFileService.h" 
#include "art_root_io/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

// c++
#include <vector>
#include <utility>


namespace sbnd{

  class PandoraT0CosmicTagAlg {
  public:

    struct BeamTime {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<double> BeamTimeMin {
        Name("BeamTimeMin"),
        Comment("")
      };

      fhicl::Atom<double> BeamTimeMax {
        Name("BeamTimeMax"),
        Comment("")
      };

    };

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> PandoraLabel {
        Name("PandoraLabel"),
        Comment("")
      };

      fhicl::Atom<art::InputTag> TpcTrackModuleLabel {
        Name("TpcTrackModuleLabel"),
        Comment("")
      };

      fhicl::Table<BeamTime> BeamTimeLimits {
        Name("BeamTimeLimits"),
        Comment("")
      };

    };

    PandoraT0CosmicTagAlg(const Config& config);

    PandoraT0CosmicTagAlg(const fhicl::ParameterSet& pset) :
      PandoraT0CosmicTagAlg(fhicl::Table<Config>(pset, {})()) {}

    PandoraT0CosmicTagAlg();

    ~PandoraT0CosmicTagAlg();

    void reconfigure(const Config& config);

    // Finds any t0s associated with track by pandora, tags if outside beam
    bool PandoraT0CosmicTag(recob::Track track, const art::Event& event);

    // Finds any t0s associated with pfparticle by pandora, tags if outside beam
    bool PandoraT0CosmicTag(recob::PFParticle pfparticle, std::map< size_t, art::Ptr<recob::PFParticle> > pfParticleMap, const art::Event& event);

  private:

    art::InputTag fPandoraLabel;
    art::InputTag fTpcTrackModuleLabel;
    double fBeamTimeMin;
    double fBeamTimeMax;

  };

}

#endif