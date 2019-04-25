#include "PandoraT0CosmicTagAlg.h"

namespace sbnd{

PandoraT0CosmicTagAlg::PandoraT0CosmicTagAlg(const Config& config){

  this->reconfigure(config);

}


PandoraT0CosmicTagAlg::PandoraT0CosmicTagAlg(){

}


PandoraT0CosmicTagAlg::~PandoraT0CosmicTagAlg(){

}


void PandoraT0CosmicTagAlg::reconfigure(const Config& config){

  fPandoraLabel = config.PandoraLabel();
  fTpcTrackModuleLabel = config.TpcTrackModuleLabel();
  fBeamTimeLimit = config.BeamTimeLimit();

  return;
}

bool PandoraT0CosmicTagAlg::PandoraT0CosmicTag(recob::Track track, const art::Event& event){

  art::Handle< std::vector<recob::PFParticle> > pfParticleHandle;
  event.getByLabel(fPandoraLabel, pfParticleHandle);
  art::FindManyP<recob::Track> pfPartToTrackAssoc(pfParticleHandle, event, fTpcTrackModuleLabel);
  art::FindManyP<anab::T0> findManyT0(pfParticleHandle, event, fPandoraLabel);
  for(auto const pfp : (*pfParticleHandle)){
    const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pfp.Self()));
    if(associatedTracks.size() != 1) continue;
    recob::Track trk = *associatedTracks.front();
    if(trk.ID() != track.ID()) continue;
    const std::vector< art::Ptr<anab::T0> > associatedT0s(findManyT0.at(pfp.Self()));
    for(size_t i = 0; i < associatedT0s.size(); i++){
      double pandoraTime = associatedT0s[i]->Time()*1e-3;
      if(pandoraTime < 0 || pandoraTime > fBeamTimeLimit) return true;
    }
  }

  return false;

}

bool PandoraT0CosmicTagAlg::PandoraT0CosmicTag(recob::PFParticle pfparticle, std::map< size_t, art::Ptr<recob::PFParticle> > pfParticleMap, const art::Event& event){

  art::Handle< std::vector<recob::PFParticle> > pfParticleHandle;
  event.getByLabel(fPandoraLabel, pfParticleHandle);
  art::FindManyP<anab::T0> findManyT0(pfParticleHandle, event, fPandoraLabel);
  for (const size_t daughterId : pfparticle.Daughters()){
    art::Ptr<recob::PFParticle> pParticle = pfParticleMap.at(daughterId);
    const std::vector< art::Ptr<anab::T0> > associatedT0s(findManyT0.at(pParticle.key()));
    for(size_t i = 0; i < associatedT0s.size(); i++){
      double pandoraTime = associatedT0s[i]->Time()*1e-3;
      if(pandoraTime < 0 || pandoraTime > fBeamTimeLimit) return true;
    }
  }

  return false;

}


}
