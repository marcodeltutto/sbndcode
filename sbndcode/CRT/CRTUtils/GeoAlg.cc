#include "GeoAlg.h"

namespace sbnd{


GeoAlg::GeoAlg(){

  fMinX = 99999;
  fMinY = 99999;
  fMinZ = 99999;
  fMaxX = -99999;
  fMaxY = -99999;
  fMaxZ = -99999;
  fCpaWidth = 0;

  fGeometryService = lar::providerFrom<geo::Geometry>();

  for(size_t cryo_i = 0; cryo_i < fGeometryService->Ncryostats(); cryo_i++){
    const geo::CryostatGeo& cryostat = fGeometryService->Cryostat(cryo_i);

    for (size_t tpc_i = 0; tpc_i < cryostat.NTPC(); tpc_i++)
    {
      const geo::TPCGeo& tpcg = cryostat.TPC(tpc_i);
      if (tpcg.MinX() < fMinX) fMinX = tpcg.MinX();
      if (tpcg.MaxX() > fMaxX) fMaxX = tpcg.MaxX();
      if (tpcg.MinY() < fMinY) fMinY = tpcg.MinY();
      if (tpcg.MaxY() > fMaxY) fMaxY = tpcg.MaxY();
      if (tpcg.MinZ() < fMinZ) fMinZ = tpcg.MinZ();
      if (tpcg.MaxZ() > fMaxZ) fMaxZ = tpcg.MaxZ();
      fCpaWidth = std::min(std::abs(tpcg.MinX()), std::abs(tpcg.MaxX()));
    }
  }
}


GeoAlg::~GeoAlg(){

}

double GeoAlg::MinX() const{
  return fMinX;
}

double GeoAlg::MinY() const{
  return fMinY;
}

double GeoAlg::MinZ() const{
  return fMinZ;
}

double GeoAlg::MaxX() const{
  return fMaxX;
}

double GeoAlg::MaxY() const{
  return fMaxY;
}

double GeoAlg::MaxZ() const{
  return fMaxZ;
}

double GeoAlg::CpaWidth() const{
  return fCpaWidth;
}

bool GeoAlg::InFiducial(geo::Point_t point, double fiducial){
  return InFiducial(point, fiducial, fiducial);
}

bool GeoAlg::InFiducial(geo::Point_t point, double fiducial, double fiducialTop){
  return InFiducial(point, fiducial, fiducial, fiducial, fiducial, fiducialTop, fiducial);
}

bool GeoAlg::InFiducial(geo::Point_t point, double minXCut, double minYCut, double minZCut, double maxXCut, double maxYCut, double maxZCut){
  
  double xmin = fMinX + minXCut;
  double xmax = fMaxX - maxXCut;
  double ymin = fMinY + minYCut;
  double ymax = fMaxY - maxYCut;
  double zmin = fMinZ + minZCut;
  double zmax = fMaxZ - maxZCut;

  double x = point.X();
  double y = point.Y();
  double z = point.Z();
  if(x>xmin && x<xmax && y>ymin && y<ymax && z>zmin && z<zmax) return true;

  return false;
}

int GeoAlg::DetectedInTPC(std::vector<art::Ptr<recob::Hit>> hits){
  // Return tpc of hit collection or -1 if in multiple
  int tpc = hits[0]->WireID().TPC;
  for(size_t i = 0; i < hits.size(); i++){
    if((int)hits[i]->WireID().TPC != tpc) return -1;
  }
  return tpc;
}


bool GeoAlg::EntersVolume(simb::MCParticle particle){
  for(size_t i = 0; i < particle.NumberTrajectoryPoints(); i++){
    double x = particle.Vx(i); 
    double y = particle.Vy(i);
    double z = particle.Vz(i);
    if(x > fMinX && y > fMinY && z > fMinZ && x < fMaxX && y < fMaxY && z < fMaxZ){
      return true;
    }
  }
  return false;
}

bool GeoAlg::CrossesVolume(simb::MCParticle particle){
  bool enters = false;
  bool startOutside = false;
  bool endOutside = false;
  for(size_t i = 0; i < particle.NumberTrajectoryPoints(); i++){
    double x = particle.Vx(i); 
    double y = particle.Vy(i);
    double z = particle.Vz(i);
    if(x > fMinX && y > fMinY && z > fMinZ && x < fMaxX && y < fMaxY && z < fMaxZ){
      enters = true;
    }
    else if(i == 0) startOutside = true;
    else if(i == particle.NumberTrajectoryPoints()-1) endOutside = true;
  }
  if(startOutside && enters && endOutside) return true;
  return false;
}

}
