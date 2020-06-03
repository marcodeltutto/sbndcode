////////////////////////////////////////////////////////////////////////
// Class:       PFPSliceValidation
// Plugin Type: analyzer (art v3_02_06)
// File:        PFPSliceValidation_module.cc
//
// Generated at Wed Oct  2 03:27:09 2019 by Edward Tyley using cetskelgen
// from cetlib version v3_07_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"

//LArSoft Includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

//Root Includes
#include "TMath.h"
#include "TTree.h"
#include <vector>
#include <iostream>

namespace ana {
  class PFPSliceValidation;
}


class ana::PFPSliceValidation : public art::EDAnalyzer {
  public:
    explicit PFPSliceValidation(fhicl::ParameterSet const& pset);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    PFPSliceValidation(PFPSliceValidation const&) = delete;
    PFPSliceValidation(PFPSliceValidation&&) = delete;
    PFPSliceValidation& operator=(PFPSliceValidation const&) = delete;
    PFPSliceValidation& operator=(PFPSliceValidation&&) = delete;

    // Required functions.
    void analyze(art::Event const& evt) override;
    void beginJob();

    std::map<art::Ptr<simb::MCTruth>, int> GetTruthHitMap(
        const sim::ParticleList& trueParticlesMap,
        const std::map<int, art::Ptr<simb::MCTruth> >& particleTruthMap,
        const std::vector< art::Ptr< recob::Hit> >& allHits);

    art::Ptr<simb::MCTruth> GetSliceTruthMatchHits(
        const std::vector< art::Ptr< recob::Hit> >& sliceHits,
        const std::map<int, art::Ptr<simb::MCTruth> >& particleTruthMap,
        const std::map<art::Ptr<simb::MCTruth>, int>& truthHitMap,
        float& completeness, float& purity);

    void ClearTrueTree();
    void ClearEventTree();

    template <class T>
    void initTree(TTree* Tree, std::string branchName,
        std::map<std::string, T>& Metric,
        std::vector<std::string> fPFParticleLabels);

  private:

    int fVerbose;
    std::string fHitLabel, fGenieGenModuleLabel;
    std::vector<std::string> fPFParticleLabels;
    art::ServiceHandle<art::TFileService> tfs;
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;
    // Declare member data here.
    TTree* eventTree;
    TTree* trueTree;

    int eventTrueNeutrinos;
    std::map<std::string, int>  eventPFPSlices, eventPFPNeutrinos;
    std::map<std::string, std::vector<float> > eventCosmicScores, eventNeutrinoScores;

    std::map<std::string, bool> nuMatchNeutrino;
    std::map<std::string, int> nuSlices, nuNeutrinos;
    std::map<std::string, float> bestNuPurity, bestNuComp, bestNuScore;

    int intType, CCNC, neutrinoPDG, numProtons, numNeutrons, numPi, numPi0;
    double W, X, Y, QSqr, Pt, Theta, neutrinoE, leptonP;
    float trueVertexX, trueVertexY, trueVertexZ;

    std::map<std::string, float> pfpVertexX, pfpVertexY, pfpVertexZ;
    std::map<std::string, float> pfpVertexDistX, pfpVertexDistY, pfpVertexDistZ, pfpVertexDistMag;
};


ana::PFPSliceValidation::PFPSliceValidation(fhicl::ParameterSet const& pset)
  : EDAnalyzer{pset}
{
  fVerbose             = pset.get<int>("Verbose", 0);
  fHitLabel            = pset.get<std::string>("HitLabel");
  fGenieGenModuleLabel = pset.get<std::string>("GenieGenModuleLabel");
  fPFParticleLabels    = pset.get<std::vector<std::string> >("PFParticleLabels");
}

void ana::PFPSliceValidation::beginJob() {
  trueTree  = tfs->make<TTree>("trueTree", "Tree with true neutrino metrics");
  eventTree = tfs->make<TTree>("eventTree", "Tree with event wide metrics");

  eventTree->Branch("trueNeutrinos", &eventTrueNeutrinos);

  initTree(eventTree, "pfpNeutrinos", eventPFPNeutrinos, fPFParticleLabels);
  initTree(eventTree, "pfpSlices", eventPFPSlices, fPFParticleLabels);
  initTree(eventTree, "cosmicScores", eventCosmicScores, fPFParticleLabels);
  initTree(eventTree, "nuScores", eventNeutrinoScores, fPFParticleLabels);

  trueTree->Branch("intType",&intType);
  trueTree->Branch("CCNC",&CCNC);
  trueTree->Branch("neutrinoPDG",&neutrinoPDG);
  trueTree->Branch("numProtons",&numProtons);
  trueTree->Branch("numPi",&numPi);
  trueTree->Branch("numPi0",&numPi0);

  trueTree->Branch("W",&W);
  trueTree->Branch("X",&X);
  trueTree->Branch("Y",&Y);
  trueTree->Branch("QSqr",&QSqr);
  trueTree->Branch("Pt",&Pt);
  trueTree->Branch("Theta",&Theta);
  trueTree->Branch("neutrinoE",&neutrinoE);
  trueTree->Branch("leptonP",&leptonP);

  initTree(trueTree, "bestMatchNeutrino", nuMatchNeutrino, fPFParticleLabels);
  initTree(trueTree, "numSlices", nuSlices, fPFParticleLabels);
  initTree(trueTree, "numNeutrinos", nuNeutrinos, fPFParticleLabels);
  initTree(trueTree, "purity", bestNuPurity, fPFParticleLabels);
  initTree(trueTree, "comp", bestNuComp, fPFParticleLabels);
  initTree(trueTree, "score", bestNuScore, fPFParticleLabels);

  // Throw some vertex reco stuff into the tree
  trueTree->Branch("trueVertexX",&trueVertexX);
  trueTree->Branch("trueVertexY",&trueVertexY);
  trueTree->Branch("trueVertexZ",&trueVertexZ);

  initTree(trueTree, "pfpVertexX", pfpVertexX, fPFParticleLabels);
  initTree(trueTree, "pfpVertexY", pfpVertexY, fPFParticleLabels);
  initTree(trueTree, "pfpVertexZ", pfpVertexZ, fPFParticleLabels);
  initTree(trueTree, "pfpVertexDistX", pfpVertexDistX, fPFParticleLabels);
  initTree(trueTree, "pfpVertexDistY", pfpVertexDistY, fPFParticleLabels);
  initTree(trueTree, "pfpVertexDistZ", pfpVertexDistZ, fPFParticleLabels);
  initTree(trueTree, "pfpVertexDistMag", pfpVertexDistMag, fPFParticleLabels);
}

void ana::PFPSliceValidation::analyze(art::Event const& evt)
{
  ClearEventTree();

  // Get the truths in the event:
  const std::vector<art::Ptr<simb::MCTruth> > truthVec = particleInventory->MCTruthVector_Ps();
  // for (auto const& truth: truthVec){
  // }

  // Get a map of each true particle to the MC Truth
  std::map<int, art::Ptr<simb::MCTruth> > particleTruthMap;
  const sim::ParticleList& trueParticlesMap= particleInventory->ParticleList();
  for (auto const& [trackId, particle]: trueParticlesMap){
    particleTruthMap[trackId] = particleInventory->ParticleToMCTruth_P(particle);
  }
  eventTrueNeutrinos = truthVec.size();

  // Get reco
  // Initialse some stuff???
  std::vector<art::Handle<std::vector<recob::Hit> > > hitHandles;
  evt.getManyByType(hitHandles);
  std::vector<art::Handle<std::vector<recob::Vertex> > > vertexHandles;
  evt.getManyByType(vertexHandles);

  // Set the handles
  art::Handle<std::vector<recob::Hit> > hitHandle;
  art::Handle<std::vector<recob::PFParticle> > pfpHandle;
  art::Handle<std::vector<recob::Slice> > sliceHandle;
  art::Handle<std::vector<recob::Vertex > > vertexHandle;

  // Get all the hits
  std::vector<art::Ptr<recob::Hit> > allHits;
  if(evt.getByLabel(fHitLabel,hitHandle))
  {art::fill_ptr_vector(allHits, hitHandle);}

  // Get map of true primary particle to number of reco hits / energy in reco hits
  std::map<art::Ptr<simb::MCTruth>, int> truthHitMap = GetTruthHitMap(trueParticlesMap,
      particleTruthMap, allHits);

  std::map<std::string, std::map<art::Ptr<simb::MCTruth>, unsigned int> > pfpTruthNuCounterMap;
  std::map<std::string, std::map<art::Ptr<simb::MCTruth>, unsigned int> > pfpTruthSliceCounterMap;
  std::map<std::string, std::map<art::Ptr<simb::MCTruth>, float> > pfpTruthCompMap;
  std::map<std::string, std::map<art::Ptr<simb::MCTruth>, float> > pfpTruthPurityMap;
  std::map<std::string, std::map<art::Ptr<simb::MCTruth>, float> > pfpTruthScoreMap;
  std::map<std::string, std::map<art::Ptr<simb::MCTruth>, int> > pfpTruthNuMap;
  std::map<std::string, std::map<art::Ptr<simb::MCTruth>, double> > pfpTruthVtxMapX;
  std::map<std::string, std::map<art::Ptr<simb::MCTruth>, double> > pfpTruthVtxMapY;
  std::map<std::string, std::map<art::Ptr<simb::MCTruth>, double> > pfpTruthVtxMapZ;

  for (auto const fPFParticleLabel: fPFParticleLabels){
    for (auto const& truth: truthVec){
      pfpTruthNuCounterMap[fPFParticleLabel][truth]    = 0;
      pfpTruthSliceCounterMap[fPFParticleLabel][truth] = 0;
    }
  }

  for (auto const fPFParticleLabel: fPFParticleLabels){

    if (fVerbose){
      std::cout << "On PFParticleLabel: " << fPFParticleLabel <<std::endl;
    }

    // Get all the PFPs
    std::vector<art::Ptr<recob::Slice> > pfpSliceVec;
    if(evt.getByLabel(fPFParticleLabel, sliceHandle))
    {art::fill_ptr_vector(pfpSliceVec, sliceHandle);}

    std::vector<art::Ptr<recob::PFParticle> > pfps;
    if(evt.getByLabel(fPFParticleLabel, pfpHandle))
    {art::fill_ptr_vector(pfps, pfpHandle);}

    art::FindOneP<recob::Vertex> fopfv(pfpHandle, evt, fPFParticleLabel);
    if (fopfv.isValid() && fopfv.size()>0){
      evt.get(fopfv.at(0).id(),vertexHandle);
      if(!vertexHandle.isValid()) {
        std::cout<<"Vertex handle not valid"<<std::endl;
        return;
      }
    }

    art::FindManyP<recob::Hit> fmSliceHits(pfpSliceVec, evt, fPFParticleLabel);
    if (!fmSliceHits.isValid() || fmSliceHits.size()==0){
      std::cout<<"FindMany Slice Hits not valid"<<std::endl;
      return;
    }
    art::FindManyP<recob::PFParticle> fmSlicePFPs(pfpSliceVec, evt, fPFParticleLabel);
    if (!fmSlicePFPs.isValid() || fmSlicePFPs.size()==0){
      std::cout<<"FindMany Slice PFPs not valid"<<std::endl;
      return;
    }
    // Create a map between PFParticles and their IDs
    art::FindManyP<larpandoraobj::PFParticleMetadata> fmpfpmd(pfps, evt, fPFParticleLabel);
    if (!fmpfpmd.isValid() || fmpfpmd.size()==0){
      std::cout<<"PFP MetaData handle not valid"<<std::endl;
      return;
    }

    std::map<long unsigned int, art::Ptr<recob::PFParticle> > pfpMap;
    std::map<long unsigned int, float > pfpNuScoreMap;
    std::vector<art::Ptr<recob::PFParticle> > pfpNeutrinoVec;
    for (auto const& pfp: pfps){
      long unsigned int pfpID = pfp->Self();
      pfpMap[pfpID] = pfp;
      if ((pfp->PdgCode()==12) ||(pfp->PdgCode()==14)){
        pfpNeutrinoVec.push_back(pfp);

        const std::vector< art::Ptr<larpandoraobj::PFParticleMetadata> > pfpMetaVec = fmpfpmd.at(pfpID);
        for (auto const pfpMeta: pfpMetaVec)
        {
          larpandoraobj::PFParticleMetadata::PropertiesMap propertiesMap = pfpMeta->GetPropertiesMap();
          pfpNuScoreMap[pfpID] = propertiesMap.at("NuScore");
        }
      }
    }

    eventPFPSlices[fPFParticleLabel] = pfpSliceVec.size();
    eventPFPNeutrinos[fPFParticleLabel] = pfpNeutrinoVec.size();

    for (const auto& pfpSlice: pfpSliceVec){

      ++eventPFPSlices[fPFParticleLabel];

      std::vector<art::Ptr<recob::Hit> > sliceHits = fmSliceHits.at(pfpSlice.key());
      std::vector<art::Ptr<recob::PFParticle> > slicePFPs = fmSlicePFPs.at(pfpSlice.key());

      bool isNeutrinoSlice(false);
      float nuScore(-999);
      long unsigned int pfpNu(-999);
      for (auto const& pfp: slicePFPs){
        if ((pfp->PdgCode()==12) ||(pfp->PdgCode()==14)){
          pfpNu = pfp->Self();
          nuScore = pfpNuScoreMap[pfpNu];
          isNeutrinoSlice = true;
          ++eventPFPNeutrinos[fPFParticleLabel];
          break;
        }
      }

      float purity(-999), completeness(-999);
      art::Ptr<simb::MCTruth> trueMatch = GetSliceTruthMatchHits(sliceHits, particleTruthMap,
          truthHitMap, completeness, purity);

      if (trueMatch.isNull()){
        if (isNeutrinoSlice){
          eventCosmicScores[fPFParticleLabel].push_back(nuScore);
        }
        continue;
      }

      if (fVerbose){
        std::cout << "True Match: "  << trueMatch << " with completeness: " << completeness
          << " and purity: " << purity    << " and score: "         << nuScore
          << std::endl;;
      }

      ++pfpTruthSliceCounterMap[fPFParticleLabel][trueMatch];
      if (isNeutrinoSlice) {
        ++pfpTruthNuCounterMap[fPFParticleLabel][trueMatch];
        eventNeutrinoScores[fPFParticleLabel].push_back(nuScore);
      }
      if (completeness > pfpTruthCompMap[fPFParticleLabel][trueMatch]){
        pfpTruthCompMap[fPFParticleLabel][trueMatch]   = completeness;
        pfpTruthPurityMap[fPFParticleLabel][trueMatch] = purity;
        pfpTruthScoreMap[fPFParticleLabel][trueMatch]  = nuScore;
        pfpTruthNuMap[fPFParticleLabel][trueMatch]     = pfpNu;

        art::Ptr<recob::PFParticle> pfpNeutrino = pfpMap.at(pfpNu);;
        art::Ptr<recob::Vertex> pfpVertex = fopfv.at(pfpNeutrino.key());

        double pfpVtx[3];
        pfpVertex->XYZ(pfpVtx);
        pfpTruthVtxMapX[fPFParticleLabel][trueMatch] = pfpVtx[0];
        pfpTruthVtxMapY[fPFParticleLabel][trueMatch] = pfpVtx[1];
        pfpTruthVtxMapZ[fPFParticleLabel][trueMatch] = pfpVtx[2];
      }
    }
  }

  eventTree->Fill();

  // Get the true neutrinos
  for (auto const& truth: truthVec){

    ClearTrueTree();

    const simb::MCNeutrino neutrino = truth->GetNeutrino();
    const simb::MCParticle nu = neutrino.Nu();
    const simb::MCParticle lepton  = neutrino.Lepton();

    intType = neutrino.Mode();
    CCNC = neutrino.CCNC();
    neutrinoPDG = nu.PdgCode();

    W = neutrino.W();
    X = neutrino.X();
    Y = neutrino.Y();
    QSqr = neutrino.QSqr();
    Pt = neutrino.Pt();
    Theta = neutrino.Theta();
    neutrinoE = nu.E();
    leptonP = lepton.P();

    trueVertexX = nu.Vx();
    trueVertexY = nu.Vy();
    trueVertexZ = nu.Vz();

    for (auto const& [primary, truthIter]: particleTruthMap){
      if (truthIter!=truth) continue;
      const simb::MCParticle* particle = trueParticlesMap.at(primary);
      if (particle->PdgCode() == 2212 && (particle->E()-particle->Mass())>0.021) {
        ++numProtons;
      }  else if (TMath::Abs(particle->PdgCode()) == 211) {
        ++numPi;
      }  else if (TMath::Abs(particle->PdgCode()) == 111) {
        ++numPi0;
      }
    }

    for (auto const fPFParticleLabel: fPFParticleLabels){
      // if (pfpTruthSliceCounterMap[fPFParticleLabel].find(truth)
      //     != pfpTruthSliceCounterMap[fPFParticleLabel].end()){

        nuSlices[fPFParticleLabel] = pfpTruthSliceCounterMap[fPFParticleLabel][truth];
        nuNeutrinos[fPFParticleLabel] = pfpTruthNuCounterMap[fPFParticleLabel][truth];
        bestNuComp[fPFParticleLabel] = pfpTruthCompMap[fPFParticleLabel][truth];
        bestNuPurity[fPFParticleLabel] = pfpTruthPurityMap[fPFParticleLabel][truth];
        bestNuScore[fPFParticleLabel] = pfpTruthScoreMap[fPFParticleLabel][truth];
        nuMatchNeutrino[fPFParticleLabel] = (pfpTruthNuMap[fPFParticleLabel][truth]!=-999);

        if (!nuMatchNeutrino[fPFParticleLabel])
          continue;

        pfpVertexX[fPFParticleLabel] = pfpTruthVtxMapX[fPFParticleLabel][truth];
        pfpVertexY[fPFParticleLabel] = pfpTruthVtxMapY[fPFParticleLabel][truth];
        pfpVertexZ[fPFParticleLabel] = pfpTruthVtxMapZ[fPFParticleLabel][truth];

        pfpVertexDistX[fPFParticleLabel] = pfpVertexX[fPFParticleLabel] - nu.Vx();
        pfpVertexDistY[fPFParticleLabel] = pfpVertexY[fPFParticleLabel] - nu.Vy();
        pfpVertexDistZ[fPFParticleLabel] = pfpVertexZ[fPFParticleLabel] - nu.Vz();

        pfpVertexDistMag[fPFParticleLabel] = sqrt(
            pfpVertexDistX[fPFParticleLabel] * pfpVertexDistX[fPFParticleLabel] +
            pfpVertexDistY[fPFParticleLabel] * pfpVertexDistY[fPFParticleLabel] +
            pfpVertexDistZ[fPFParticleLabel] * pfpVertexDistZ[fPFParticleLabel]);
      }
    // }
    trueTree->Fill();
  }

  eventTree->Fill();
  std::cout<<"\n"<<std::endl;

}

std::map<art::Ptr<simb::MCTruth>, int> ana::PFPSliceValidation::GetTruthHitMap(
    const sim::ParticleList& trueParticlesMap,
    const std::map<int, art::Ptr<simb::MCTruth> >& particleTruthMap,
    const std::vector< art::Ptr< recob::Hit> >& allHits){

  std::map<int,int> trueParticleHits;
  for (const auto& hit: allHits){
    int trackID     = 0;
    float hitEnergy = 0;
    std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(hit);
    for (const auto& ide: trackIDEs) {
      if (ide.energy > hitEnergy){
        hitEnergy = ide.energy;
        trackID   = TMath::Abs(ide.trackID);
      }
    }
    ++trueParticleHits[trackID];
  }

  std::map<art::Ptr<simb::MCTruth>, int> truthHitMap;
  for (const auto& [trueParticle, truth]: particleTruthMap){
    truthHitMap[truth] += trueParticleHits[trueParticle];
  }
  return truthHitMap;
}

art::Ptr<simb::MCTruth> ana::PFPSliceValidation::GetSliceTruthMatchHits(
    const std::vector< art::Ptr< recob::Hit> >& sliceHits,
    const std::map<int, art::Ptr<simb::MCTruth> >& particleTruthMap,
    const std::map<art::Ptr<simb::MCTruth>, int>& truthHitMap,
    float& completeness, float& purity){

  std::map<int,int> particleHits;
  for (const auto& hit: sliceHits){
    int trackID     = 0;
    float hitEnergy = 0;
    std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(hit);
    for (const auto& ide: trackIDEs) {
      if (ide.energy > hitEnergy){
        hitEnergy = ide.energy;
        trackID   = TMath::Abs(ide.trackID);
      }
    }
    ++particleHits[trackID];
  }

  std::map<art::Ptr<simb::MCTruth>, int> sliceTruthHitMap;
  for (const auto& [particle, truth]: particleTruthMap){
    sliceTruthHitMap[truth] += particleHits[particle];
  }

  int maxHits = 0;
  art::Ptr<simb::MCTruth> bestTruthMatch;
  for (const auto& [truth, truthHits]: sliceTruthHitMap){
    if (truthHits > maxHits){
      maxHits        = truthHits;
      bestTruthMatch = truth;
    }
  }

  purity       = (float) maxHits / sliceHits.size();
  completeness = (float) maxHits / truthHitMap.at(bestTruthMatch);
  return bestTruthMatch;
}

void ana::PFPSliceValidation::ClearTrueTree(){

  intType     = -999;
  CCNC        = -999;
  neutrinoPDG = -999;
  numProtons  = -999;
  numNeutrons = -999;
  numPi       = -999;
  numPi0      = -999;

  W         = -999;
  X         = -999;
  Y         = -999;
  QSqr      = -999;
  Pt        = -999;
  Theta     = -999;
  neutrinoE = -999;
  leptonP   = -999;

  trueVertexX = -999;
  trueVertexY = -999;
  trueVertexZ = -999;

  for (auto const fPFParticleLabel: fPFParticleLabels){

    nuMatchNeutrino[fPFParticleLabel] = false;
    nuSlices[fPFParticleLabel]       = -999;
    nuNeutrinos[fPFParticleLabel]    = -999;
    bestNuPurity[fPFParticleLabel]   = -999;
    bestNuComp[fPFParticleLabel]     = -999;
    bestNuScore[fPFParticleLabel]    = -999;

    pfpVertexX[fPFParticleLabel] = -999;
    pfpVertexY[fPFParticleLabel] = -999;
    pfpVertexZ[fPFParticleLabel] = -999;

    pfpVertexDistX[fPFParticleLabel]   = -999;
    pfpVertexDistY[fPFParticleLabel]   = -999;
    pfpVertexDistZ[fPFParticleLabel]   = -999;
    pfpVertexDistMag[fPFParticleLabel] = -999;
  }
}

void ana::PFPSliceValidation::ClearEventTree(){
  eventTrueNeutrinos = -999;
  for (auto const fPFParticleLabel: fPFParticleLabels){
    eventPFPNeutrinos[fPFParticleLabel]   = -999;
    eventPFPSlices[fPFParticleLabel]      = -999;
    eventCosmicScores[fPFParticleLabel].clear();
    eventNeutrinoScores[fPFParticleLabel].clear();
  }
}

template <class T>
void ana::PFPSliceValidation::initTree(TTree* Tree, std::string branchName,
    std::map<std::string, T>& Metric, std::vector<std::string> fPFParticleLabels){

  for (auto const& fPFParticleLabel: fPFParticleLabels){
    std::string branchString = branchName + "_" + fPFParticleLabel;
    const char* branchChar   = branchString.c_str();
    Tree->Branch(branchChar, &Metric[fPFParticleLabel], 32000, 0);
  }
}

DEFINE_ART_MODULE(ana::PFPSliceValidation)
