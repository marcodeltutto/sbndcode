#ifndef CPACROSSCOSMICTAGALG_H_SEEN
#define CPACROSSCOSMICTAGALG_H_SEEN


///////////////////////////////////////////////
// CpaCrossCosmicTagAlg.h
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
#include "lardataobj/RecoBase/Hit.h"
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

  class CpaCrossCosmicTagAlg {
  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<double> CpaStitchDistance {
        Name("CpaStitchDistance"),
        Comment("")
      };

      fhicl::Atom<double> CpaStitchAngle {
        Name("CpaStitchAngle"),
        Comment("")
      };

      fhicl::Atom<double> CpaXDifference {
        Name("CpaXDifference"),
        Comment("")
      };

      fhicl::Atom<double> Fiducial {
        Name("Fiducial"),
        Comment("")
      };

      fhicl::Atom<double> FiducialTop {
        Name("FiducialTop"),
        Comment("")
      };

      fhicl::Atom<double> BeamTimeLimit {
        Name("BeamTimeLimit"),
        Comment("")
      };

    };

    CpaCrossCosmicTagAlg(const Config& config);

    CpaCrossCosmicTagAlg(const fhicl::ParameterSet& pset) :
      CpaCrossCosmicTagAlg(fhicl::Table<Config>(pset, {})()) {}

    CpaCrossCosmicTagAlg();

    ~CpaCrossCosmicTagAlg();

    void reconfigure(const Config& config);

    std::pair<double, bool> T0FromCpaStitching(recob::Track t1, std::vector<recob::Track> tracks);

    bool CpaCrossCosmicTag(recob::Track track, std::vector<recob::Track> tracks, art::FindManyP<recob::Hit> hitAssoc);

  private:

    double fCpaStitchDistance;
    double fCpaStitchAngle;
    double fCpaXDifference;
    double fFiducial;
    double fFiducialTop;
    double fBeamTimeLimit;

    detinfo::DetectorProperties const* fDetectorProperties;

  };

}

#endif
