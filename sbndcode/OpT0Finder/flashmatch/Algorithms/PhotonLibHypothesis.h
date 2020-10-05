/**
 * \file PhotonLibHypothesis.h
 *
 * \ingroup Algorithms
 *
 * \brief Class def header for a class PhotonLibHypothesis
 *
 * @author yuntse
 */

/** \addtogroup Algorithms

    @{*/

#ifndef PHOTONLIBHYPOTHESIS_H
#define PHOTONLIBHYPOTHESIS_H

#ifndef USING_LARSOFT
#define USING_LARSOFT 1
#endif

#if USING_LARSOFT == 0
#include "flashmatch/Base/OpT0FinderException.h"
#include "flashmatch/Base/FMWKInterface.h"
#include "flashmatch/Base/BaseFlashFilter.h"
#include "flashmatch/Base/FlashHypothesisFactory.h"
#else
#include "sbndcode/OpT0Finder/flashmatch/Base/OpT0FinderException.h"
#include "sbndcode/OpT0Finder/flashmatch/Base/FMWKInterface.h"
#include "sbndcode/OpT0Finder/flashmatch/Base/BaseFlashFilter.h"
#include "sbndcode/OpT0Finder/flashmatch/Base/FlashHypothesisFactory.h"
#endif

#include <iostream>
#include <cassert>
#include <chrono>

namespace flashmatch {
  /**
     \class PhotonLibHypothesis
     User custom analysis class made by SHELL_USER_NAME
   */
  class PhotonLibHypothesis : public BaseFlashHypothesis {

  public:

    /// Default constructor
    PhotonLibHypothesis(const std::string name="PhotonLibHypothesis");

    /// Default destructor
    virtual ~PhotonLibHypothesis(){}

    void FillEstimate(const QCluster_t&, Flash_t&) const;

  protected:

    void _Configure_(const Config_t &pset);

    double _global_qe;             ///< Global QE
    double _sigma_qe;              ///< Sigma for Gaussian centered on Global QE
    std::vector<double> _qe_v;     ///< PMT-wise relative QE
  };

  /**
     \class flashmatch::PhotonLibHypothesisFactory
  */
  class PhotonLibHypothesisFactory : public FlashHypothesisFactoryBase {
  public:
    /// ctor
    PhotonLibHypothesisFactory() { FlashHypothesisFactory::get().add_factory("PhotonLibHypothesis",this); }
    /// dtor
    ~PhotonLibHypothesisFactory() {}
    /// creation method
    BaseFlashHypothesis* create(const std::string instance_name) { return new PhotonLibHypothesis(instance_name); }
  };
}
#endif

/** @} */ // end of doxygen group