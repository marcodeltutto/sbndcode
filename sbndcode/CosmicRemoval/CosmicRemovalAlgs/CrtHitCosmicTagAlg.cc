#include "CrtHitCosmicTagAlg.h"

namespace sbnd{

CrtHitCosmicTagAlg::CrtHitCosmicTagAlg(const Config& config){

  this->reconfigure(config);

}


CrtHitCosmicTagAlg::CrtHitCosmicTagAlg(){

}


CrtHitCosmicTagAlg::~CrtHitCosmicTagAlg(){

}


void CrtHitCosmicTagAlg::reconfigure(const Config& config){

  t0Alg = config.T0Alg();
  fBeamTimeLimit = config.BeamTimeLimit();

  return;
}


bool CrtHitCosmicTagAlg::CrtHitCosmicTag(recob::Track track, std::vector<crt::CRTHit> crtHits, int tpc){

  double crtHitTime = t0Alg.T0FromCRTHits(track, crtHits, tpc);

  if(crtHitTime != -99999 && (crtHitTime < 0 || crtHitTime > fBeamTimeLimit)) return true;

  return false;

}
 
}
