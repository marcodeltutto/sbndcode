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

}
