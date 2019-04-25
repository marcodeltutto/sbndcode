#include "CrtTrackCosmicTagAlg.h"

namespace sbnd{

CrtTrackCosmicTagAlg::CrtTrackCosmicTagAlg(const Config& config){

  this->reconfigure(config);

}


CrtTrackCosmicTagAlg::CrtTrackCosmicTagAlg(){

}


CrtTrackCosmicTagAlg::~CrtTrackCosmicTagAlg(){

}


void CrtTrackCosmicTagAlg::reconfigure(const Config& config){

  trackMatchAlg = config.TrackMatchAlg();
  fBeamTimeLimit = config.BeamTimeLimit();

  return;
}


bool CrtTrackCosmicTagAlg::CrtTrackCosmicTag(recob::Track track, std::vector<crt::CRTTrack> crtTracks, int tpc){

  double crtTrackTime = trackMatchAlg.T0FromCRTTracks(track, crtTracks, tpc);

  if(crtTrackTime != -99999) return true;

  return false;

}
 
}
