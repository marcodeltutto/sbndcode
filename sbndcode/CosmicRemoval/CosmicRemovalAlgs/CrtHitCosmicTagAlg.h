#ifndef CRTHITCOSMICTAGALG_H_SEEN
#define CRTHITCOSMICTAGALG_H_SEEN


///////////////////////////////////////////////
// CrtHitCosmicTagAlg.h
//
// Functions for fiducial volume cosmic tagger
// T Brooks (tbrooks@fnal.gov), November 2018
///////////////////////////////////////////////

#include "sbndcode/CosmicRemoval/CosmicRemovalUtils/CosmicRemovalUtils.h"
#include "sbndcode/CRT/CRTUtils/CRTAnaUtils.h"

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
#include "sbndcode/CRT/CRTProducts/CRTHit.hh"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

// c++
#include <vector>
#include <utility>


namespace sbnd{

  class CrtHitCosmicTagAlg {
  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Table<CRTT0MatchAlg::Config> T0Alg {
        Name("T0Alg"),
        Comment("")
      };

      fhicl::Atom<double> BeamTimeLimit {
        Name("BeamTimeLimit"),
        Comment("")
      };

    };

    CrtHitCosmicTagAlg(const Config& config);

    CrtHitCosmicTagAlg(const fhicl::ParameterSet& pset) :
      CrtHitCosmicTagAlg(fhicl::Table<Config>(pset, {})()) {}

    CrtHitCosmicTagAlg();

    ~CrtHitCosmicTagAlg();

    void reconfigure(const Config& config);

    bool CrtHitCosmicTag(recob::Track track, std::vector<crt::CRTHit> crtHits, int tpc);

  private:

    CRTT0MatchAlg t0Alg;
    double fBeamTimeLimit;

  };

}

#endif
