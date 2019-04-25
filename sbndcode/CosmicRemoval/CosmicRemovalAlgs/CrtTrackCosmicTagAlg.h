#ifndef CRTTRACKCOSMICTAGALG_H_SEEN
#define CRTTRACKCOSMICTAGALG_H_SEEN


///////////////////////////////////////////////
// CrtTrackCosmicTagAlg.h
//
// Functions for fiducial volume cosmic tagger
// T Brooks (tbrooks@fnal.gov), November 2018
///////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/CRTTrackMatchAlg.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

// c++
#include <vector>
#include <utility>


namespace sbnd{

  class CrtTrackCosmicTagAlg {
  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Table<CRTTrackMatchAlg::Config> TrackMatchAlg {
        Name("TrackMatchAlg"),
        Comment("")
      };

      fhicl::Atom<double> BeamTimeLimit {
        Name("BeamTimeLimit"),
        Comment("")
      };

    };

    CrtTrackCosmicTagAlg(const Config& config);

    CrtTrackCosmicTagAlg(const fhicl::ParameterSet& pset) :
      CrtTrackCosmicTagAlg(fhicl::Table<Config>(pset, {})()) {}

    CrtTrackCosmicTagAlg();

    ~CrtTrackCosmicTagAlg();

    void reconfigure(const Config& config);

    bool CrtTrackCosmicTag(recob::Track track, std::vector<crt::CRTTrack> crtTracks, int tpc);

  private:

    CRTTrackMatchAlg trackMatchAlg;
    double fBeamTimeLimit;

  };

}

#endif
