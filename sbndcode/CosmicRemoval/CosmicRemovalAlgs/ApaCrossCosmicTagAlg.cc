#include "ApaCrossCosmicTagAlg.h"

namespace sbnd{

ApaCrossCosmicTagAlg::ApaCrossCosmicTagAlg(const Config& config){

  this->reconfigure(config);

  fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
  fGeometryService = lar::providerFrom<geo::Geometry>();

}


ApaCrossCosmicTagAlg::ApaCrossCosmicTagAlg(){

  fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
  fGeometryService = lar::providerFrom<geo::Geometry>();

}


ApaCrossCosmicTagAlg::~ApaCrossCosmicTagAlg(){

}


void ApaCrossCosmicTagAlg::reconfigure(const Config& config){

  fApaDistance = config.ApaDistance(); 
  fFiducial = config.Fiducial();
  fBeamTimeLimit = config.BeamTimeLimit();

  return;
}

double ApaCrossCosmicTagAlg::T0FromApaCross(recob::Track track, std::vector<double> t0List, int tpc){

  double crossTime = -99999;
  double xmax = 2.0 * fGeometryService->DetHalfWidth();

  double minDist = 99999;
  double startX = track.Vertex().X();
  double endX = track.End().X();
  geo::Point_t point = track.Vertex();

  // Don't try to shift tracks near the Apa
  if(std::abs(startX) > xmax-fFiducial || std::abs(endX) > xmax-fFiducial) return crossTime;

  // If in tpc 0 use start/end with lowest X
  if(tpc == 0 && endX < startX) point = track.End();

  // If in tpc 1 use start/end with highest X
  if(tpc == 1 && endX > startX) point = track.End();

  //Shift track by all t0's
  for(auto const& t0 : t0List){
    // If particle crosses the APA before t = 0 the crossing point won't be reconstructed
    if(t0 < 0) continue;
    double shiftedX = point.X();
    double shift = t0 * fDetectorProperties->DriftVelocity();
    if(tpc == 0) shiftedX = point.X() - shift;
    if(tpc == 1) shiftedX = point.X() + shift;

    //Check track still in TPC
    if(std::abs(shiftedX) > 201.) continue; //FIXME
    //Calculate distance between start/end and APA
    double dist = std::abs(std::abs(shiftedX)-199.6);//FIXME
    if(dist < minDist) {
      minDist = dist;
      crossTime = t0;
    }
  }

  if(minDist < fApaDistance) return crossTime;

  return -99999;

}

bool ApaCrossCosmicTagAlg::ApaCrossCosmicTag(recob::Track track, std::vector<art::Ptr<recob::Hit>> hits, std::vector<double> t0Tpc0, std::vector<double> t0Tpc1){

  int tpc = CosmicRemovalUtils::DetectedInTPC(hits);

  if(tpc == 0){
    double crossTimeThrough = T0FromApaCross(track, t0Tpc0, tpc);
    if(crossTimeThrough != -99999 && (crossTimeThrough < 0 || crossTimeThrough > fBeamTimeLimit)) return true;
  }
  if(tpc == 1){
    double crossTimeThrough = T0FromApaCross(track, t0Tpc1, tpc);
    if(crossTimeThrough != -99999 && (crossTimeThrough < 0 || crossTimeThrough > fBeamTimeLimit)) return true;
  }

  return false;

}

}
