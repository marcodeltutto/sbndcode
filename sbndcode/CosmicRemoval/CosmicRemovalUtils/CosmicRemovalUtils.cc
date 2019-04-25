#include "CosmicRemovalUtils.h"

namespace sbnd{

// =============================== UTILITY FUNCTIONS ==============================

  bool CosmicRemovalUtils::InFiducial(geo::Point_t point, double fiducial, double fiducialTop){
    
    geo::GeometryCore const* fGeometryService = lar::providerFrom<geo::Geometry>();
    double xmin = -200 + fiducial; //-2.0 * fGeometryService->DetHalfWidth() + fiducial;
    double xmax = 200 - fiducial; //2.0 * fGeometryService->DetHalfWidth() - fiducial;
    double ymin = -fGeometryService->DetHalfHeight() + fiducial;
    double ymax = fGeometryService->DetHalfHeight() - fiducialTop;
    double zmin = 0. + fiducial;
    double zmax = fGeometryService->DetLength() - fiducial;

    double x = point.X();
    double y = point.Y();
    double z = point.Z();
    if(x>xmin && x<xmax && y>ymin && y<ymax && z>zmin && z<zmax) return true;

    return false;
  }

  int CosmicRemovalUtils::DetectedInTPC(std::vector<art::Ptr<recob::Hit>> hits){
    //
    int tpc = hits[0]->WireID().TPC;
    for(size_t i = 0; i < hits.size(); i++){
      if((int)hits[i]->WireID().TPC != tpc) return -1;
    }
    return tpc;
  }

  std::pair<double,bool> CosmicRemovalUtils::T0FromCpaStitching(recob::Track t1, std::vector<recob::Track> tracks, double stitchDist, double stitchAngle, double xDiff, double fiducial, double fiducialTop){

    detinfo::DetectorProperties const* fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    std::vector<std::pair<double, std::pair<double, bool>>> matchCandidates;
    double matchedTime = -99999;
    std::pair<double, bool> returnVal = std::make_pair(matchedTime, false);

    TVector3 trk1Front = t1.Vertex<TVector3>();
    TVector3 trk1Back = t1.End<TVector3>();
    double closestX1 = std::min(std::abs(trk1Front.X()), std::abs(trk1Back.X()));

    for(auto & track : tracks){

      TVector3 trk2Front = track.Vertex<TVector3>();
      TVector3 trk2Back = track.End<TVector3>();
      double closestX2 = std::min(std::abs(trk2Front.X()), std::abs(trk2Back.X()));

      if(std::abs(closestX1-closestX2) < xDiff){
        TVector3 t1Pos = trk1Front;
        TVector3 t1PosEnd = trk1Back;
        TVector3 t1Dir = t1.VertexDirection<TVector3>();
        if(std::abs(trk1Back.X()) == closestX1){ 
          t1Pos = trk1Back;
          t1PosEnd = trk1Front;
          t1Dir = t1.EndDirection<TVector3>();
        }

        TVector3 t2Pos = trk2Front;
        TVector3 t2PosEnd = trk2Back;
        TVector3 t2Dir = track.VertexDirection<TVector3>();
        if(std::abs(trk2Back.X()) == closestX2){ 
          t2Pos = trk2Back;
          t2PosEnd = trk2Front;
          t2Dir = track.EndDirection<TVector3>();
        }

        double trkCos = std::abs(t1Dir.Dot(t2Dir));
        t1Pos[0] = 0.;
        t2Pos[0] = 0.;
        double dist = (t1Pos-t2Pos).Mag();

        geo::Point_t mergeStart {t1PosEnd.X(), t1PosEnd.Y(), t1PosEnd.Z()};
        geo::Point_t mergeEnd {t2PosEnd.X(), t2PosEnd.Y(), t2PosEnd.Z()};
        bool exits = false;
        if(!CosmicRemovalUtils::InFiducial(mergeStart, fiducial, fiducialTop) && !CosmicRemovalUtils::InFiducial(mergeEnd, fiducial, fiducialTop)) exits = true;

        if(dist < stitchDist && trkCos > cos(TMath::Pi() * stitchAngle / 180.)){ 
          matchCandidates.push_back(std::make_pair(trkCos, std::make_pair(closestX1, exits)));
        }
      }
    }

    if(matchCandidates.size() > 0){
      std::sort(matchCandidates.begin(), matchCandidates.end(), [](auto& left, auto& right){
                return left.first < right.first;});
      double shiftX = matchCandidates[0].second.first;
      matchedTime = -(shiftX/fDetectorProperties->DriftVelocity()-17.); //FIXME
      returnVal = std::make_pair(matchedTime, matchCandidates[0].second.second);
    }

    return returnVal;
  }

  double CosmicRemovalUtils::T0FromApaCross(recob::Track track, std::vector<double> t0s, int tpc, double fiducial, double distLimit){

    geo::GeometryCore const* fGeometryService = lar::providerFrom<geo::Geometry>();
    detinfo::DetectorProperties const* fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    double crossTime = -99999;
    double xmax = 2.0 * fGeometryService->DetHalfWidth();

    double minDist = 99999;
    double startX = track.Vertex().X();
    double endX = track.End().X();
    geo::Point_t point = track.Vertex();
    // Don't try to shift tracks near the Apa
    if(std::abs(startX) > xmax-fiducial || std::abs(endX) > xmax-fiducial) return crossTime;
    // If in tpc 0 use start/end with lowest X
    if(tpc == 0 && endX < startX) point = track.End();
    // If in tpc 1 use start/end with highest X
    if(tpc == 1 && endX > startX) point = track.End();

    //Shift track by all t0's
    for(auto const& t0 : t0s){
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
    if(minDist < distLimit) return crossTime;

    return -99999;
      
  }

  double CosmicRemovalUtils::StoppingEnd(std::vector<art::Ptr<anab::Calorimetry>> calos, geo::Point_t end, double rangeMin, double rangeMax, double dedxMax, double chi2Lim){

    //Loop over residual range and dedx
    if(calos.size()==0) return false;
    size_t nhits = 0;
    art::Ptr<anab::Calorimetry> calo = calos[0];
    for( size_t i = calos.size(); i > 0; i--){
      if(calos[i-1]->dEdx().size() > nhits*1.5){
        nhits = calos[i-1]->dEdx().size();
        calo = calos[i-1];
      }
    }

    double distStart = (calo->XYZ()[0] - end).Mag2();
    double distEnd = (calo->XYZ()[nhits-1] - end).Mag2();

    double maxDedx = 0;
    double resrgStart = 0;
    std::vector<double> v_resrg;
    std::vector<double> v_dedx;

    for(size_t i = 0; i < nhits; i++){
      double dedx = calo->dEdx()[i];
      double resrg = calo->ResidualRange()[i];

      if(distStart < distEnd && calo->ResidualRange()[0] > calo->ResidualRange()[nhits-1]) resrg = calo->ResidualRange()[0] - calo->ResidualRange()[i];
      if(distStart > distEnd && calo->ResidualRange()[0] < calo->ResidualRange()[nhits-1]) resrg = calo->ResidualRange()[nhits-1] - calo->ResidualRange()[i];

      //if(resrg < rangeMin && dedx > maxDedx && dedx < dedxMax){
      if(resrg < 10 && dedx > maxDedx && dedx < 30){
        maxDedx = dedx;
        resrgStart = resrg;
      }

    }

    for(size_t i = 0; i < nhits; i++){
      double dedx = calo->dEdx()[i];
      double resrg = calo->ResidualRange()[i];

      if(distStart < distEnd && calo->ResidualRange()[0] > calo->ResidualRange()[nhits-1]) resrg = calo->ResidualRange()[0] - calo->ResidualRange()[i];
      if(distStart > distEnd && calo->ResidualRange()[0] < calo->ResidualRange()[nhits-1]) resrg = calo->ResidualRange()[nhits-1] - calo->ResidualRange()[i];

      //if(resrg > resrgStart && resrg < resrgStart+rangeMax && dedx < dedxMax){
      if(resrg > resrgStart && resrg < resrgStart+20 && dedx < 30){
        v_resrg.push_back(resrg);
        v_dedx.push_back(dedx);
      }
    }

    if(v_dedx.size() < 10) return false;

    TGraph *gdedx = new TGraph(v_dedx.size(), &v_resrg[0], &v_dedx[0]);
    try{ gdedx->Fit("pol0", "Q"); } catch(...){ return false; }
    TF1* polfit = gdedx->GetFunction("pol0");
    double polchi2 = polfit->GetChisquare();

    try{ gdedx->Fit("expo", "Q"); } catch(...){ return false; }
    TF1* expfit = gdedx->GetFunction("expo");
    double expchi2 = expfit->GetChisquare();

    //if(polchi2/expchi2 > chi2Lim) return true;
    if(polchi2/expchi2 > 1.2) return true;

    return false;

  }

  std::pair<std::vector<double>, std::vector<double>> CosmicRemovalUtils::FakeTpcFlashes(std::vector<simb::MCParticle> particles){
    //
    geo::GeometryCore const* fGeometryService    = lar::providerFrom<geo::Geometry>();
    detinfo::DetectorProperties const* fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 
    detinfo::DetectorClocks const* fDetectorClocks = lar::providerFrom<detinfo::DetectorClocksService>(); 

    // Create fake flashes in each tpc
    std::vector<double> fakeTpc0Flashes;
    std::vector<double> fakeTpc1Flashes;

    double readoutWindowMuS  = fDetectorClocks->TPCTick2Time((double)fDetectorProperties->ReadOutWindowSize()); // [us]
    double driftTimeMuS = (2.*fGeometryService->DetHalfWidth()+3.)/fDetectorProperties->DriftVelocity(); // [us]

    // Loop over all true particles
    for (auto const particle: particles){

      // Get particle info
      int pdg = std::abs(particle.PdgCode());
      double time = particle.T() * 1e-3;

      //Check if time is in reconstructible window
      if(time < -driftTimeMuS || time > readoutWindowMuS) continue; 
      //Check if particle is visible, electron, muon, proton, pion, kaon, photon
      if(!(pdg==13||pdg==11||pdg==22||pdg==2212||pdg==211||pdg==321||pdg==111)) continue;

      //Loop over the trajectory
      int npts = particle.NumberTrajectoryPoints();
      double TPC0Energy = 0;
      double TPC1Energy = 0;
      for(int i = 1; i < npts; i++){
        geo::Point_t pt;
        pt.SetX(particle.Vx(i)); pt.SetY(particle.Vy(i)); pt.SetZ(particle.Vz(i));
        if(!CosmicRemovalUtils::InFiducial(pt, 0, 0)) continue;
        // Add up the energy deposited in each tpc
        if(pt.X() < 0) TPC0Energy += particle.E(i-1) - particle.E(i);
        else TPC1Energy += particle.E(i-1) - particle.E(i);
      }
      // If the total energy deposited is > 50 MeV then create fake flash
      if(TPC0Energy > 0.05) fakeTpc0Flashes.push_back(time);
      else if(TPC1Energy > 0.05) fakeTpc1Flashes.push_back(time);
    }

    std::sort(fakeTpc0Flashes.begin(), fakeTpc0Flashes.end());
    double previousTime = -99999;
    // Loop over flashes in tpc 0
    for(size_t i = 0; i < fakeTpc0Flashes.size(); i++){
      double time = fakeTpc0Flashes[i];
      // Combine flashes within 0.1 us
      if(std::abs(time-previousTime) < 0.1){
        fakeTpc0Flashes.erase(fakeTpc0Flashes.begin()+i);
      }
      else previousTime = time;
    }

    std::sort(fakeTpc1Flashes.begin(), fakeTpc1Flashes.end());
    previousTime = -99999;
    // Loop over flashes in tpc 0
    for(size_t i = 0; i < fakeTpc1Flashes.size(); i++){
      double time = fakeTpc1Flashes[i];
      // Combine flashes within 0.1 us
      if(std::abs(time-previousTime) < 0.1){
        fakeTpc1Flashes.erase(fakeTpc1Flashes.begin()+i);
      }
      else previousTime = time;
    }

    return std::make_pair(fakeTpc0Flashes, fakeTpc1Flashes);
  }

  bool CosmicRemovalUtils::BeamFlash(std::vector<double> flashes, double beamTimeLimit){
    //
    bool beamFlash = false;
    std::sort(flashes.begin(), flashes.end());
    // Loop over flashes in tpc 0
    for(size_t i = 0; i < flashes.size(); i++){
      double time = flashes[i];
      if(time > 0 && time < beamTimeLimit) beamFlash = true;
    }

    return beamFlash;
  }

}
