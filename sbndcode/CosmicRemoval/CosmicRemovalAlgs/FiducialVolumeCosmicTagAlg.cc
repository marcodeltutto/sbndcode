#include "FiducialVolumeCosmicTagAlg.h"

namespace sbnd{

FiducialVolumeCosmicTagAlg::FiducialVolumeCosmicTagAlg(const Config& config){

  this->reconfigure(config);

}


FiducialVolumeCosmicTagAlg::FiducialVolumeCosmicTagAlg(){

}


FiducialVolumeCosmicTagAlg::~FiducialVolumeCosmicTagAlg(){

}


void FiducialVolumeCosmicTagAlg::reconfigure(const Config& config){

  fFiducial = config.Fiducial(); 
  fFiducialTop = config.FiducialTop();

  return;
}

bool FiducialVolumeCosmicTagAlg::InFiducial(geo::Point_t point){
  return CosmicRemovalUtils::InFiducial(point, fFiducial, fFiducialTop);
}

bool FiducialVolumeCosmicTagAlg::FiducialVolumeCosmicTag(recob::Track track){
  
  bool startInFiducial = CosmicRemovalUtils::InFiducial(track.Vertex(), fFiducial, fFiducialTop);

  bool endInFiducial = CosmicRemovalUtils::InFiducial(track.End(), fFiducial, fFiducialTop);

  if(!startInFiducial && !endInFiducial)  return true;
  
  return false;
  
}


}
