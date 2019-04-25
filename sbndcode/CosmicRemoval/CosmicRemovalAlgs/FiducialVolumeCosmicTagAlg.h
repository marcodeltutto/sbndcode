#ifndef FIDUCIALVOLUMECOSMICTAGALG_H_SEEN
#define FIDUCIALVOLUMECOSMICTAGALG_H_SEEN


///////////////////////////////////////////////
// FiducialVolumeCosmicTagAlg.h
//
// Functions for fiducial volume cosmic tagger
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
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft
#include "lardataobj/RecoBase/Track.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

// c++
#include <vector>


namespace sbnd{

  class FiducialVolumeCosmicTagAlg {
  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<double> Fiducial {
        Name("Fiducial"),
        Comment("Fiducial volume cut (cm) from all sides and bottom")
      };

      fhicl::Atom<double> FiducialTop {
        Name("FiducialTop"),
        Comment("Fiducial volume cut from top plane")
      };

    };

    FiducialVolumeCosmicTagAlg(const Config& config);

    FiducialVolumeCosmicTagAlg(const fhicl::ParameterSet& pset) :
      FiducialVolumeCosmicTagAlg(fhicl::Table<Config>(pset, {})()) {}

    FiducialVolumeCosmicTagAlg();

    ~FiducialVolumeCosmicTagAlg();

    void reconfigure(const Config& config);

    bool InFiducial(geo::Point_t point);

    bool FiducialVolumeCosmicTag(recob::Track track);

  private:

    double fFiducial;
    double fFiducialTop;

  };

}

#endif
