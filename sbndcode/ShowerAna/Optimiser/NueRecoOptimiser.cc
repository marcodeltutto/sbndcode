#include "NueRecoOptimiser.hh"

//C++ Includes 
#include <iostream>
#include <typeinfo>

//Root Includes
#include "TFile.h"
#include "TMVA/MethodBDT.h"
#include "TMVA/CrossValidation.h"
#include "TMVA/Tools.h"

#include "TError.h"

optimiser::NueRecoOptimiser::NueRecoOptimiser(TTree* signaltree, TTree* backgroundtree, float bkpot, float sigpot){

  SignalTree.AddBranch<std::vector<float>>("number_of_showers_per_neutrino",signaltree);
  SignalTree.AddBranch<std::vector<float>>("vertex_recoK",signaltree);
  SignalTree.AddBranch<std::vector<float>>("vertex_trueK",signaltree);
  SignalTree.AddBranch<std::vector<float>>("vertex_reco",signaltree);
  SignalTree.AddBranch<std::vector<float>>("nu_reco_energy",signaltree);
  SignalTree.AddBranch<std::vector<float>>("nu_truth_energy",signaltree);
  SignalTree.AddBranch<std::vector<float>>("nu_interaction_type",signaltree);
  SignalTree.AddBranch<std::vector<float>>("nu_mode",signaltree);
  SignalTree.AddBranch<std::vector<float>>("nu_E",signaltree);
  SignalTree.AddBranch<std::vector<float>>("nu_E_numtrue",signaltree);
  SignalTree.AddBranch<std::vector<float>>("nu_distance",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("truepionE",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("trueprotonE",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("truekaonE",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("truetrackE",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_energy",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("truth_pid",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("true_energy",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_coversion_gap",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_residual_dist",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_length",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_density",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_length_perp",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_density_perp",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_density_3D",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_density_grad_perp",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_density_grad",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_density_ratio",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_density_grad_ratio",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_density_grad_sq",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_density_grad_perp_sq",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_density_grad_ratio_sq",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_dEdx",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("track_lengths",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("track_PIDA_branch",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("track_E",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("track_trueE",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("track_pdg",signaltree); 
  SignalTree.AddBranch<std::vector<std::vector<float> >>("track_resE",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_density_grad_new",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_open_angle",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_density_pw_new",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_density_pwgrad_new",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_tracklength",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_trackwidth",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_hits",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_len",signaltree);

  SignalTree.AddBranch<std::vector<float> >("in_FV",signaltree);
  SignalTree.AddBranch<std::vector<float> >("nu_osc_prob",signaltree);
  SignalTree.AddBranch<std::vector<float> >("nu_pdg",signaltree);
  SignalTree.AddBranch<std::vector<float> >("nu_truepdg",signaltree);
  SignalTree.AddBranch<std::vector<float> >("nu_cc",signaltree);
  SignalTree.AddBranch<std::vector<float> >("nu_leptonE",signaltree);
  SignalTree.AddBranch<std::vector<float> >("nu_X",signaltree);
  SignalTree.AddBranch<std::vector<float> >("nu_W",signaltree);  
  SignalTree.AddBranch<std::vector<float> >("nu_Y",signaltree);
  SignalTree.AddBranch<std::vector<float> >("nu_QSqr",signaltree);
  SignalTree.AddBranch<std::vector<float> >("nu_Pt",signaltree);
  SignalTree.AddBranch<std::vector<float> >("nu_Theta",signaltree);

  SignalTree.AddBranch<int>("trueShower_num",signaltree);
  SignalTree.AddBranch<int>("numtrueVtx_branch",signaltree);
  SignalTree.AddBranch<float>("POT",signaltree);

  BackgroundTree.AddBranch<std::vector<float>>("number_of_showers_per_neutrino",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float>>("vertex_recoK",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float>>("vertex_trueK",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float>>("vertex_reco",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float>>("nu_reco_energy",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float>>("nu_truth_energy",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float>>("nu_interaction_type",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float>>("nu_mode",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float>>("nu_E",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float>>("nu_E_numtrue",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float>>("nu_distance",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("truepionE",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("trueprotonE",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("truekaonE",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("truetrackE",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_energy",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("truth_pid",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("true_energy",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_coversion_gap",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_residual_dist",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_length",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_density",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_length_perp",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_density_perp",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_density_3D",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_density_grad_perp",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_density_grad",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_density_ratio",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_density_grad_ratio",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_density_grad_perp_sq",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_density_grad_sq",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_density_grad_ratio_sq",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_dEdx",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("track_lengths",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("track_PIDA_branch",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("track_E",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("track_trueE",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("track_pdg",backgroundtree); 
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("track_resE",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_density_grad_new",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_open_angle",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_density_pw_new",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_density_pwgrad_new",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_tracklength",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_trackwidth",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_hits",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_len",backgroundtree);

  BackgroundTree.AddBranch<std::vector<float> >("in_FV",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float> >("nu_osc_prob",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float> >("nu_pdg",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float> >("nu_truepdg",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float> >("nu_cc",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float> >("nu_leptonE",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float> >("nu_X",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float> >("nu_W",backgroundtree);  
  BackgroundTree.AddBranch<std::vector<float> >("nu_Y",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float> >("nu_QSqr",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float> >("nu_Pt",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float> >("nu_Theta",backgroundtree);

  BackgroundTree.AddBranch<int>("trueShower_num",backgroundtree);
  BackgroundTree.AddBranch<int>("numtrueVtx_branch",backgroundtree);
  BackgroundTree.AddBranch<float>("POT",backgroundtree);

  MVAFile = new TFile("MVAFile.root","RECREATE");
  MVAFile->cd();
  BackgroundTreeMVA = new TTree("BackgroundTreeMVA","BackgroundTreeMVA");
  SignalTreeMVA     = new TTree("SignalTreeMVA","SignalTreeMVA");

  BackgroundTreeMVA->Branch("weight",&BackgrndWeight,"weight/D");
  SignalTreeMVA->Branch("weight",&SignalWeight,"weight/D");


  TotalSigPOT = sigpot;
  TotalBKPOT  = bkpot;

}

void optimiser::NueRecoOptimiser::InitialiseMetrics(){
  
  //MVA Analsysis
  ApplyOscWht = true;
  ApplyPOTWht = true;
  fApplyFVCut = true;
  InitialiseMetric("dEdxMVA","shower_dEdx",50,0,10,"StandardShowerCutFinder",true,true,true,-1);
  InitialiseMetric("CoversionGapMVA","shower_coversion_gap",50,0,10,"StandardShowerCutFinder",true,true,true,-1);
  InitialiseMetric("MaxTrackLenghMVA","track_lengths",50,0,200,"StandardMaxDaughterCutFinder",true,true,true,-1); 
  //  InitialiseMetric("MaxTrackLenghMVANoScale","track_lengths",50,0,200,"StandardMaxDaughterCutFinder",true,true,true,-1); 
  //  InitialiseMetric("CoversionGapMVANoScale","shower_coversion_gap",50,0,10,"StandardShowerCutFinder",true,true,true,-1);

  InitialiseMetric("MaxTrackPIDAMVA","track_PIDA_branch",50,0,25,"StandardMaxDaughterCutFinder",true,true,true,-1); 
  InitialiseMetric("ShowerEnergyMVA","shower_energy",50,0,1500,"StandardShowerCutFinder",true,true,true,-100);
  InitialiseMetric("ShowerLengthEMVA","shower_length",50,0,150,"StandardShowerCutFinder",true,true,true,-1);
  InitialiseMetric("NumberOfNeutrinosMVA","nu_reco_energy",5,0,5,"OneShowerSizeNeutrinoAnalysisCutFinder",true,true,true,-1);
  InitialiseMetric("ShowerResidualNumShowersMVA","shower_residual_dist",6,0,6,"ShowerResidualNumShowersCutFinder",true,true,true,-1);
  InitialiseMetric("ShowerDensityPWMVA","shower_density_pw_new",50,1,2,"StandardShowerCutFinder",true,true,true,-1);
  InitialiseMetric("ShowerDensityGradNewMVA","shower_density_grad_new",50,0,1,"StandardShowerCutFinder",true,true,true,-1);
  InitialiseMetric("ShowerDensityOpeningAngleMVA","shower_open_angle",50,0,1.5,"StandardShowerCutFinder",true,true,true,-1);
  ////  InitialiseMetric("NuetrinoPdGMVA","nu_pdg",20,0,20,"NeutrinoComparision",true,true,true,-1);
  InitialiseMetric("NuetrinoPdGMVA","nu_pdg",2,12,14,"StandardNeutrinoCutFinder",true,true,true,-1);

  InitialiseMetric("NuetrinoPdGMVAPOther","nu_pdg",3,12,15,"StandardNeutrinoCutFinder",true,true,false,-1);

  InitialiseMetric("ShowerTrackWidthMVA","shower_trackwidth",50,0,5,"StandardShowerCutFinder",true,true,true,-1);
  InitialiseMetric("ShowerTrackLengthMVA","shower_tracklength",50,0,10,"StandardShowerCutFinder",true,true,true,-1);
  //  InitialiseMetric("ShowerRes","shower_residual_dist",50,0,100,"StandardShowerCutFinder",true,true,true,-1);

  InitialiseMetric("dEdxAboveE","shower_dEdx",50,0,10,"StandardAboveEDaughterCutFinder",true,true,false,-1);




  //  InitialiseMetric("NumberOfShowersIntegralMVA","shower_energy",50,0,1000,"number_of_showers_per_neutrino",10,0,10,"StandardNumShowerAnalysis",true,true,false,-1);

  //  InitialiseMetric("OneShowerResidualAnalysis","shower_energy",100,10,260,"shower_residual_dist",100,0,3.8,"OneShower2DEnergyAnalysis",true,true,false,-1);
  //   InitialiseMetric("OneShowerResidualAnalysisHits","shower_hits",100,0,800,"shower_residual_dist",100,0,5,"OneShower2DEnergyAnalysis");
  //  InitialiseMetric("ShowerResiudal","shower_residual_dist",100,0,10);
  //  InitialiseMetric("ShowerResiudalEcut","shower_residual_dist",100,0,10,"StandardAboveEDaughterCutFinder");
  //  InitialiseMetric("ShowerdEdx","shower_dEdx",50,0,10,"StandardAboveEDaughterCutFinder",true,true,true,-1);
  //  InitialiseMetric("TrackRes","track_resE",50,-2,2,"StandardMaxDaughterCutFinder");

  //  InitialiseMetric("NumberOfShowersIntegral","shower_energy",50,0,1000,"number_of_showers_per_neutrino",10,0,10,"StandardNumShowerAnalysis",false,false,false,-1);
  InitialiseMetric("NumberOfShowers","number_of_showers_per_neutrino",5,0,5,"StandardNeutrinoCutFinder",true,true,false,10);
  //InitialiseMetric("MaxTrackLengh2D","track_lengths",200,0,100,"track_PIDA_branch",200,0,25,"StandardMaxDaughterCutFinder",true,true,false,-1); 
  InitialiseMetric("zBDTGCrossVal","nu_reco_energy", 40,-1,1,"MVAAnalysisCutFinder",true,true,false,-999,"BDTGCrossVal");
  InitialiseMetric("zBDTGCrossValScaled","nu_reco_energy",40,-1,1,"MVAAnalysisCutFinder",true,true,false,-999,"BDTGCrossValScaled");

  //  InitialiseMetric("zBDTGCrossValNotNorm","nu_reco_energy", 200,-1,1,"MVAAnalysisCutFinder",false,false,false,-999,"BDTGCrossVal");
  // InitialiseMetric("zBDTGCrossValScaledNorNorm","nu_reco_energy",200,-1,1,"MVAAnalysisCutFinder",false,false,false,-999,"BDTGCrossValScaled");


  //  InitialiseMetric("zBDTGCrossValOsc","nu_reco_energy",1000,-1,1,"MVAAnalysisCutFinder",true,true,false,-999,"BDTGCrossValOsc");
  //  InitialiseMetric("zBDTGCrossValStandSet","nu_reco_energy",1000,-1,1,"MVAAnalysisCutFinder",true,true,false,-999,"BDTGCrossValStandSet");
  //  InitialiseMetric("zBDTG","nu_reco_energy",1000,-1,1,"MVAAnalysisCutFinder",true,true,false,-999,"BDTG");

  //  InitialiseMetric("VertexReco","vertex_trueK",50,0,300,"vertex_reco",10,0,10,"EfficiencyNeutrinoMaker",false,false,true,-1);

  //Big Boi Full selection 
  //  InitialiseMetric("ShowerAllsCutTrueE","nu_truth_energy",50,0,5000,"ShowerAllsCut",true,true,false,-1);
  //  InitialiseMetric("ShowerAllsCutTrueENoScale","nu_truth_energy",50,0,5000,"ShowerAllsCut",false,false,false,-1);

  //  InitialiseMetric("BDTAllsCutTrueE","nu_truth_energy",50000,-1.01,1,"BDTAllsCutFinder",false,false,false,-1.02);
  
  //InitialiseMetric("ShowerAllsCutRecoE","nu_reco_energy",50,0,5,"ShowerAllsCut");
  // InitialiseMetric("ShowerAllsCutLeptonE","nu_leptonE",50,0,5,"ShowerAllsCut",true,true,false,-1);
  // InitialiseMetric("ShowerAllsX","nu_X",50,0,5,"ShowerAllsCut",true,true,false,-1);
  // InitialiseMetric("ShowerAllsW","nu_W",50,0,5,"ShowerAllsCut",true,true,false,-1);
  // InitialiseMetric("ShowerAllsY","nu_Y",50,0,5,"ShowerAllsCut",true,true,false,-1);
  // InitialiseMetric("ShowerAllsQSqr","nu_QSqr",50,0,5,"ShowerAllsCut",true,true,false,-1);
  // InitialiseMetric("ShowerAllsPt","nu_Pt",50,0,0.2,"ShowerAllsCut",true,true,false,-1);
  // InitialiseMetric("ShowerAllsTheta","nu_Theta",50,0,5,"ShowerAllsCut",true,true,false,-1);
  //  InitialiseMetric("ShowerAllRecoE","nu_reco_energy",50,0,5000,"ShowerAllsCut",true,true,false,-1);
  //  InitialiseMetric("ShowerAllRecoENoScale","nu_reco_energy",50,0,5000,"ShowerAllsCut",false,false,false,-1);


  //InitialiseMetric("StandardNeutrinoEnergyResolutionProfileTrueE","nu_truth_energy",50,0,5,"nu_truth_energy",50,-2,1.1,"StandardNeutrinoEnergyResolutionProfile",false,false,false,-1);
  // InitialiseMetric("StandardNeutrinoEnergyResolutionProfileLeptonE","nu_leptonE",50,0,5,"nu_truth_energy",50,-2,1.1,"StandardNeutrinoEnergyResolutionProfile",false,false,false,-1);
  // InitialiseMetric("StandardNeutrinoEnergyResolutionProfileX","nu_X",50,0,5,"nu_truth_energy",50,-2,1.1,"StandardNeutrinoEnergyResolutionProfile",false,false,false,-1);
  // InitialiseMetric("StandardNeutrinoEnergyResolutionProfileW","nu_W",50,0,5,"nu_truth_energy",50,-2,1.1,"StandardNeutrinoEnergyResolutionProfile",false,false,false,-1);
  // InitialiseMetric("StandardNeutrinoEnergyResolutionProfileY","nu_Y",50,0,5,"nu_truth_energy",50,-2,1.1,"StandardNeutrinoEnergyResolutionProfile",false,false,false,-1);
  // InitialiseMetric("StandardNeutrinoEnergyResolutionProfileQSqr","nu_QSqr",50,0,5,"nu_truth_energy",50,-2,1.1,"StandardNeutrinoEnergyResolutionProfile",false,false,false,-1);
  // InitialiseMetric("StandardNeutrinoEnergyResolutionProfilePt","nu_Pt",50,0,0.05,"nu_truth_energy",50,-2,1.1,"StandardNeutrinoEnergyResolutionProfile",false,false,false,-1);
  // InitialiseMetric("StandardNeutrinoEnergyResolutionProfileTheta","nu_Theta",50,0,5,"nu_truth_energy",50,-2,1.1,"StandardNeutrinoEnergyResolutionProfile",false,false,false,-1);
  //  InitialiseMetric("StandardNeutrinoEnergyResolutionRecoEnergy","nu_reco_energy",25,0,5000,"nu_truth_energy",25,-2,1.1,"StandardNeutrinoEnergyResolutionProfile",false,false,false,-1);
  //  InitialiseMetric("StandardNeutrinoEnergyResolutionRecoEnergy","nu_reco_energy",25,0,5000,"nu_truth_energy",50,-2,2,"StandardNeutrinoEnergyResolutionProfile",false,false,false,-1);


  //  InitialiseMetric("zBDTGCrossVal","MVAAllsCut",true,true,false,-999,"zBDTGCrossVal");
  //  InitialiseMetric("zBDTGCrossValScaled","MVAAllsCut",true,true,false,-999,"zBDTGCrossValScaled");

  //InitialiseMetric("zBDTGAllsCut","MVAAllsCut",true,true,false,-999,"BDTG");
  // SetLogic("zBDT",false,false,false);  
  // SetLogic("zBDTG",false,false,false);  
  // SetLogic("zBDTGNotNorm",false,false,false);  
  SetLogic("zBDTGCrossVal",false,false,false);
  SetLogic("zBDTGCrossValScaled",false,false,false);
  SetLogic("BDTAllsCutTrueE",false,false,false);
  SetLogic("ShowerEnergyMVA",false,false,false);
  SetLogic("ShowerLengthEMVA",false,false,false);
  SetLogic("MaxTrackPIDAMVA",false,false,false);
  SetLogic("MaxTrackLengh2D",true,false,false);
  // SetLogic("zBDTGCrossValOsc",false,false,false);
  //  SetLogic("zBDTGCrossValStandSet",false,false,false);

  //  SetLogic("zBDTG",false,false,false);

}

void optimiser::NueRecoOptimiser::InitialiseMetric(std::string Name, std::string branchname,int numbins, float xmin, float xmax, std::string AnalysisName, bool ApplyPOT, 
						   bool ApplyOscProb, bool fFillMVA, float MVAErrorVal, TString MVAName){
  MetricMap[Name] = optimiser::MetricHolder(Name,branchname,numbins,xmin,xmax,TotalSigPOT,TotalBKPOT,ApplyPOT,ApplyOscProb,AnalysisName,fFillMVA, MVAName,MVAErrorVal);
  if(fFillMVA){
    std::string OverF = Name + "/F";
    BackgroundTreeMVA->Branch(Name.c_str(),&MetricMap[Name].GetMVAValBackground(),OverF.c_str());
    SignalTreeMVA->Branch(Name.c_str(),&MetricMap[Name].GetMVAValSignal(),OverF.c_str());
  }
  return;
}

void optimiser::NueRecoOptimiser::InitialiseMetric(std::string Name, std::string AnalysisName, bool ApplyPOT,bool ApplyOscProb, bool fFillMVA, float MVAErrorVal, TString MVAName){
  MetricMap[Name] = optimiser::MetricHolder(Name,TotalSigPOT,TotalBKPOT,ApplyPOT,ApplyOscProb,AnalysisName,fFillMVA,MVAName,MVAErrorVal);
  if(fFillMVA){
    std::string OverF = Name + "/F";
    BackgroundTreeMVA->Branch(Name.c_str(),&MetricMap[Name].GetMVAValBackground(),OverF.c_str());
    SignalTreeMVA->Branch(Name.c_str(),&MetricMap[Name].GetMVAValSignal(),OverF.c_str());
  }
  return;
}


void optimiser::NueRecoOptimiser::InitialiseMetric(std::string Name, std::string branchname,int numbins, float xmin, float xmax, std::string metricpartner,int numbins_partner, float xmin_partner, float xmax_partner, std::string AnalysisName, bool ApplyPOT, bool ApplyOscProb, bool fFillMVA, float MVAErrorVal, TString MVAName){
  MetricMap[Name] = optimiser::MetricHolder(Name,branchname,numbins,xmin,xmax,metricpartner,numbins_partner,xmin_partner,xmax_partner,TotalSigPOT,TotalBKPOT,ApplyPOT,ApplyOscProb,AnalysisName,fFillMVA,MVAName,MVAErrorVal);

  if(fFillMVA){
    std::string OverF = Name + "/F";
    BackgroundTreeMVA->Branch(Name.c_str(),&MetricMap[Name].GetMVAValBackground(),OverF.c_str());
    SignalTreeMVA->Branch(Name.c_str(),&MetricMap[Name].GetMVAValSignal(),OverF.c_str());
  }
  return;
}

void optimiser::NueRecoOptimiser::SetLogic(std::string Name, bool ls, bool pls, bool andor){
  MetricMap[Name].SetLogic(ls,pls,andor);
  return;
}

float optimiser::NueRecoOptimiser::GetOscProb(int& iter, std::string sample, bool applyoscprob){

  if(!applyoscprob){return 1;}

  std::string branchname = "nu_osc_prob";
  std::vector<float>* oscprob;  
  int err = GetBranchPointer(branchname,sample,oscprob);
  if(err){
    std::cerr << "Osclillation probability broken" << std::endl;
    return 0;
  }

  if(oscprob->size() <= iter){
    //    std::cerr << "Osclillation probability broken" << std::endl;
    return 0;
  }
  return oscprob->at(iter);
}

void optimiser::NueRecoOptimiser::FillData(TTree* tree_signal, TTree* tree_background){

  //Fill The Background
  std::cout << " tree_background->GetEntries(): " <<  tree_background->GetEntries() << std::endl;
  for (Long64_t evt = 0; evt < tree_background->GetEntries(); evt++) {
    
    tree_background->GetEntry(evt);

    std::string branchname = "nu_truepdg";
    std::vector<float>* branch_vals;
    int electron = 0;
    std::string bk = "background";
    int err = GetBranchPointer(branchname,bk,branch_vals);
    if(err){return;}
    for(auto const& neutrino: *branch_vals){
      if(neutrino == 12){
	++electron;
      }
    }
    if(electron !=0){continue;}



    for(auto& Metric: MetricMap){
      
      std::vector<float> vals;
      std::vector<std::pair<float,float> > TwoDvals;
      TString MVAName = Metric.second.GetMVAName();
      int err = PerformCuts(Metric.second,"background",vals,TwoDvals,MVAName);
      if(err){return;}
      for(int i=0; i<vals.size(); ++i){

	if(vals[i] == Metric.second.GetErrorVal()){continue;}

	double oscprob = GetOscProb(i,"background", Metric.second.ApplyOscProb());
	double potweight = 1;
        if(Metric.second.ApplyPOT()){potweight = Metric.second.ReturnPOT("background");}
	else{potweight = potweight/tree_background->GetEntries();}

	double oscpot = oscprob*potweight;  
	Metric.second.FillBackground(vals[i],oscpot);
      }
      for(int i=0; i<TwoDvals.size(); ++i){

	int iter = i;
	if(Metric.second.GetAnalysisName() == "StandardNumShowerAnalysis" || Metric.second.GetAnalysisName() == "ShowerResidualAnalysis"){
	  iter = i/(Metric.second.GetSignalHistPartner()->GetNbinsX());
	}
	else if(Metric.second.GetAnalysisName() == "OneShower2DEnergyAnalysis"){
	  iter = i/(Metric.second.GetSignalHistPartner()->GetNbinsX()*Metric.second.GetSignalHistPartner()->GetNbinsY());
	}

	//	if(TwoDvals[i].first == Metric.second.GetErrorVal() || TwoDvals[i].second == Metric.second.GetErrorVal()){continue;}
	double oscprob = GetOscProb(iter,"background", Metric.second.ApplyOscProb());
	double potweight = 1;
        if(Metric.second.ApplyPOT()){potweight = Metric.second.ReturnPOT("background");}
	//	else{potweight = potweight/tree_background->GetEntries();}
	double oscpot = oscprob*potweight;
	Metric.second.FillBackground(TwoDvals[i].first,TwoDvals[i].second,oscpot);
      }

      if(Metric.second.FillForMVA()){
	Metric.second.MVAValBackgroundVec(vals);
      }
    }

    int MVAVectorSize = -999;
    for(auto& Metric: MetricMap){
      
      if(!Metric.second.FillForMVA()){continue;}
      if(MVAVectorSize < Metric.second.GetMVAVectorBackgroundSize()){
	MVAVectorSize =  Metric.second.GetMVAVectorBackgroundSize();
      }
    }

    for(int neutrino=0; neutrino<MVAVectorSize; ++neutrino){

      if(fApplyFVCut){
	std::vector<float>* fv_vals; 
	std::string fv_branch_name = "in_FV";
	std::string sig= "background";
	  int err = GetBranchPointer(fv_branch_name,sig,fv_vals);
	  if(err){ std::cout << "Did not get the FV branch" << std::endl; return;}
	  if(fv_vals->at(neutrino) == 0){
	    continue;
	    
	  }
      }

      for(auto& Metric: MetricMap){
	if(!Metric.second.FillForMVA()){continue;}
	float val = Metric.second.GetMVAValBackground(neutrino);
	Metric.second.FillMVAValBackground(val);
      }
      
     if(ApplyOscWht){BackgrndWeight =  GetOscProb(neutrino,"background", true);}
     if(ApplyPOTWht){BackgrndWeight /= (TotalBKPOT/6.6e20);}
     else{BackgrndWeight /= tree_background->GetEntries();}

     for(auto& Metric: MetricMap){

	if(Metric.second.GetMVAName() != ""){
	  TString MVAName = Metric.second.GetMVAName();
	  float val = reader->EvaluateMVA(MVAName);
	  Metric.second.FillBackground(val,BackgrndWeight);
	}
      }
     BackgroundTreeMVA->Fill();
     BackgrndWeight = 1;

    }
  }

  //Fill The Signal
  std::cout << "tree_signal->GetEntries(): " << tree_signal->GetEntries() << std::endl;
  for (Long64_t evt = 0; evt < tree_signal->GetEntries(); evt++) {
    tree_signal->GetEntry(evt);
   
    std::string branchname = "nu_truepdg";
    std::vector<float>* branch_vals;
    std::string sig = "signal";
    int err = GetBranchPointer(branchname,sig,branch_vals);
    int muon = 0;
    std::vector<int> muons;
    if(err){return;}
    for(auto const& neutrino: *branch_vals){
      if(neutrino == 14){
	muons.push_back(muon);
      }
	++muon;
    }
    //    if(muon !=0){continue;}

    std::string branchnamecc = "nu_cc";
    std::vector<float>* branch_valscc;
    err = GetBranchPointer(branchnamecc,sig,branch_valscc);
    int nc = 0;
    int cc = 0;
    std::vector<int> ncs;
    if(err){return;}
    for(auto const& neutrino: *branch_valscc){
      if(neutrino == 1){
	ncs.push_back(nc);
      }
      ++nc;
    }
    //if(nc !=0 && cc==0){continue;}

   
    
    for(auto& Metric: MetricMap){

      std::vector<float> vals;
      std::vector<std::pair<float,float> > TwoDvals;
      TString MVAName = Metric.second.GetMVAName();
      int err = PerformCuts(Metric.second,"signal",vals,TwoDvals,MVAName);
      if(err){return;}
      for(int i=0; i<vals.size(); ++i){

	if(std::find(muons.begin(),muons.end(),i) != muons.end()){continue;}
	if(std::find(ncs.begin(),ncs.end(),i) != ncs.end()){continue;}

	if(vals[i] == Metric.second.GetErrorVal()){continue;}
	double oscprob = GetOscProb(i,"signal",Metric.second.ApplyOscProb());
	double potweight = 1;
	if(Metric.second.ApplyPOT()){potweight = Metric.second.ReturnPOT("signal");}
	else{potweight = potweight/tree_signal->GetEntries();}

	double oscpot = oscprob*potweight;  
	Metric.second.FillSignal(vals[i],oscpot);
      }
      for(int i=0; i<TwoDvals.size(); ++i){

	if(std::find(muons.begin(),muons.end(),i) != muons.end()){continue;}
	if(std::find(ncs.begin(),ncs.end(),i) != ncs.end()){continue;}

	int iter = i;
	if(Metric.second.GetAnalysisName() == "StandardNumShowerAnalysis" || Metric.second.GetAnalysisName() == "ShowerResidualAnalysis"){
	  iter = i/(Metric.second.GetSignalHistPartner()->GetNbinsX());
	}
	else if(Metric.second.GetAnalysisName() == "OneShower2DEnergyAnalysis"){
	  iter = i/(Metric.second.GetSignalHistPartner()->GetNbinsX()*Metric.second.GetSignalHistPartner()->GetNbinsY());
	}
	
	//	if(TwoDvals[i].first == Metric.second.GetErrorVal() || TwoDvals[i].second == Metric.second.GetErrorVal()){continue;}
	double oscprob = GetOscProb(iter,"signal",Metric.second.ApplyOscProb());
	double potweight = 1;  
	if(Metric.second.ApplyPOT()){potweight = Metric.second.ReturnPOT("signal");}
	//	else{potweight = potweight/tree_signal->GetEntries();}
	double oscpot = oscprob*potweight;
	Metric.second.FillSignal(TwoDvals[i].first,TwoDvals[i].second,oscpot);
      }
      
      if(Metric.second.FillForMVA()){
	Metric.second.FillMVAValSignalVec(vals);
      }
    }

    int MVAVectorSize = -999;
    for(auto& Metric: MetricMap){
      if(!Metric.second.FillForMVA()){continue;}
      if(MVAVectorSize < Metric.second.GetMVAVectorSignalSize()){
	MVAVectorSize = Metric.second.GetMVAVectorSignalSize();
      }
    }
  
    for(int neutrino=0; neutrino<MVAVectorSize; ++neutrino){
      
      if(std::find(muons.begin(),muons.end(),neutrino) != muons.end()){continue;}
      if(std::find(ncs.begin(),ncs.end(),neutrino) != ncs.end()){continue;}


      if(fApplyFVCut){
	std::vector<float>* fv_vals; 
	std::string fv_branch_name = "in_FV";
	std::string sig= "signal";
	int err = GetBranchPointer(fv_branch_name,sig,fv_vals);
	if(err){ std::cout << "Did not get the FV branch" << std::endl; return;}
	if(fv_vals->at(neutrino) == 0){
	  continue;
	}
      }

      for(auto& Metric: MetricMap){
	if(!Metric.second.FillForMVA()){continue;}
	float val = Metric.second.GetMVAValSignal(neutrino);
	Metric.second.FillMVAValSignal(val);
      }

      if(ApplyOscWht){SignalWeight =  GetOscProb(neutrino,"signal", true);}
      if(ApplyPOTWht){SignalWeight /= (TotalSigPOT/6.6e20);}
      else{SignalWeight /= tree_signal->GetEntries();}

      for(auto& Metric: MetricMap){
	if(Metric.second.GetMVAName() != ""){
	  TString MVAName = Metric.second.GetMVAName();
	  float val = reader->EvaluateMVA(MVAName);
	  Metric.second.FillSignal(val,SignalWeight);
	}
      }

      SignalTreeMVA->Fill();
      SignalWeight = 1;
    }
  }
  return;
} 

void optimiser::NueRecoOptimiser::AnalyseData(){

  TFile* file = new TFile("CutFile.root", "RECREATE");

  //Make new directory for plots 
  for(auto const& Metric: MetricMap){
    gDirectory->mkdir((Metric.first).c_str());
  }

  //Loop over metrics and make the graphs.
  for(auto& Metric: MetricMap){

    file->cd((Metric.first).c_str());

    std::cout << "##################" << std::endl;
    std::cout << "On Metric: " << Metric.first << std::endl;
    std::cout << "##################" << std::endl;

    if(Metric.second.GetAnalysisName().find("CutFinder") != std::string::npos){
      //Perform the 1D analysis
      
      if(Metric.second.TotalEntries() != 0){
	Metric.second.Perform1DCutFinder();
      }

      //Perform the 2D analysis if any.
      if(Metric.second.GetPartnerMetricName() != "" 
	 && Metric.second.Total2DEntries() != 0){
	Metric.second.Perform2DCutFinder();
      }
    }

    if(Metric.second.GetAnalysisName().find("Resolution") != std::string::npos){
      if(Metric.second.GetPartnerMetricName() != "" 
	 && Metric.second.Total2DEntries() != 0){
	Metric.second.Perform2DCutFinder();
      }
    }


    //Now for the special analyses
    if(Metric.second.GetAnalysisName() == "StandardNumShowerAnalysis"){

      int xbins = Metric.second.GetSignalHistPartner()->GetNbinsX();
      int ybins  = Metric.second.GetSignalHistPartner()->GetNbinsY();

      float totalsig = Metric.second.GetSignalHistPartner()->Integral(0,xbins+1,0,ybins+1)/(Metric.second.GetSignalHistPartner()->GetNbinsX());
      float totalbk  = Metric.second.GetBackgroundHistPartner()->Integral(0,xbins+1,0,ybins+1)/(Metric.second.GetSignalHistPartner()->GetNbinsX());
      Metric.second.MakeStandardNumShowerAnalysisGraphs(totalsig,totalbk);
    }
    if(Metric.second.GetAnalysisName() == "OneShower2DEnergyAnalysis"){
      
      int xbins = Metric.second.GetSignalHistPartner()->GetNbinsX();
      int ybins  = Metric.second.GetSignalHistPartner()->GetNbinsY();

      float totalsig = Metric.second.GetSignalHistPartner()->Integral(0,xbins+1,0,ybins+1)/(Metric.second.GetSignalHistPartner()->GetNbinsX()*Metric.second.GetSignalHistPartner()->GetNbinsY());
      float totalbk  = Metric.second.GetBackgroundHistPartner()->Integral(0,xbins+1,0,ybins+1)/(Metric.second.GetSignalHistPartner()->GetNbinsX()*Metric.second.GetSignalHistPartner()->GetNbinsY());

      std::cout << " Integral: " <<  Metric.second.GetSignalHistPartner()->Integral(0,xbins+1,0,ybins+1) << " Entries: " << Metric.second.GetSignalHistPartner()->GetEntries() << " totalsig: " << totalsig << std::endl;

      Metric.second.MakeStandardNumShowerAnalysisGraphs(totalsig,totalbk);
    }
    if(Metric.second.GetAnalysisName() == "ShowerResidualAnalysis"){

      int xbins = Metric.second.GetSignalHistPartner()->GetNbinsX();
      int ybins  = Metric.second.GetSignalHistPartner()->GetNbinsY();


      float totalsig = Metric.second.GetSignalHistPartner()->Integral(0,xbins+1,0,ybins+1)/(Metric.second.GetSignalHistPartner()->GetNbinsX());
      float totalbk  = Metric.second.GetBackgroundHistPartner()->Integral(0,xbins+1,0,ybins+1)/(Metric.second.GetSignalHistPartner()->GetNbinsX());
      Metric.second.MakeStandardNumShowerAnalysisGraphs(totalsig,totalbk);
    }
    if(Metric.second.GetAnalysisName().find("Comparision") != std::string::npos){
      
      if(Metric.second.TotalEntries() != 0){
	Metric.second.Plot1D();
      }

      if(Metric.second.GetPartnerMetricName() != "" 
	 && Metric.second.Total2DEntries() != 0){
	Metric.second.Plot2D();
      }
    }

    if(Metric.second.GetAnalysisName() == "EfficiencyNeutrinoMaker"){
      Metric.second.MakeEfficiencyPlots(Metric.second.GetHistogram("NeutrinoRecoEAfter","signal"),Metric.second.GetHistogram("NeutrinoRecoEBefore","signal"),Metric.second.GetHistogram("NeutrinoRecoEAfterExtra","signal"));
      
      Metric.second.MakeEfficiencyPlots(Metric.second.GetHistogram("NeutrinoRecoEAfter","background"),Metric.second.GetHistogram("NeutrinoRecoEBefore","background"),Metric.second.GetHistogram("NeutrinoRecoEAfterExtra","background"));
    }

    if(Metric.second.GetAnalysisName() == "ShowerAllsCut"){


      std::string NeutrinoTrueEBefore = Metric.second.GetMetricName() + "NeutrinoTrueEBefore";
      std::string NeutrinoRecoEBefore = Metric.second.GetMetricName() + "NeutrinoRecoEBefore";
      std::string NeutrinoRecoEBeforeExtra = Metric.second.GetMetricName() + "NeutrinoRecoEBeforeExtra";
      std::string NeutrinoRecoEAfter  = Metric.second.GetMetricName() + "NeutrinoRecoEAfter"; 
      std::string NeutrinoTrueEAfter = Metric.second.GetMetricName() + "NeutrinoTrueEAfter";
      std::string NeutrinoRecoEAfterExtra = Metric.second.GetMetricName() + "NeutrinoRecoEAfterExtra";
      std::string NeutrinoTrueEAfterExtra =  Metric.second.GetMetricName() + "NeutrinoTrueEAfterExtra";

      std::string NeutrinoRecoEBeforeNorm = Metric.second.GetMetricName() + "NeutrinoRecoEBeforeNorm";
      std::string NeutrinoRecoEAfterNorm  = Metric.second.GetMetricName() + "NeutrinoRecoEAfterNorm"; 


      // Metric.second.MakeEfficiencyPlots(Metric.second.GetHistogram("NeutrinoRecoEAfter","signal"),Metric.second.GetHistogram("NeutrinoRecoEBefore","signal"),Metric.second.GetHistogram("NeutrinoRecoEAfterExtra","signal"));

      // Metric.second.MakeEfficiencyPlots(Metric.second.GetHistogram("NeutrinoRecoEAfter","background"),Metric.second.GetHistogram("NeutrinoRecoEBefore","background"),Metric.second.GetHistogram("NeutrinoRecoEAfterExtra","background"));
      // Metric.second.MakeEfficiencyPlots(Metric.second.GetHistogram("NeutrinoTrueEAfter","signal"),Metric.second.GetHistogram("NeutrinoTrueEBefore","signal"),Metric.second.GetHistogram("NeutrinoTrueEAfterExtra","signal"));
      // Metric.second.MakeEfficiencyPlots(Metric.second.GetHistogram("NeutrinoTrueEAfter","background"),Metric.second.GetHistogram("NeutrinoTrueEBefore","background"),Metric.second.GetHistogram("NeutrinoTrueEAfterExtra","background"));
      
      Metric.second.MakeEfficiencyPlots(Metric.second.GetHistogram(NeutrinoRecoEAfter,"signal"),Metric.second.GetHistogram(NeutrinoRecoEBefore,"signal"),Metric.second.GetHistogram(NeutrinoRecoEAfterExtra,"signal"),Metric.second.GetHistogram(NeutrinoRecoEAfter,"background"),Metric.second.GetHistogram(NeutrinoRecoEBefore,"background"),Metric.second.GetHistogram(NeutrinoRecoEAfterExtra,"background"),Metric.second.GetHistogram(NeutrinoRecoEAfterNorm,"signal"),Metric.second.GetHistogram(NeutrinoRecoEBeforeNorm,"signal"),Metric.second.GetHistogram(NeutrinoRecoEAfterNorm,"background"),Metric.second.GetHistogram(NeutrinoRecoEBeforeNorm,"background"));
      
      //Metric.second.MakeEfficiencyPlots(Metric.second.GetHistogram(NeutrinoTrueEAfter,"signal"),Metric.second.GetHistogram(NeutrinoTrueEBefore,"signal"),Metric.second.GetHistogram(NeutrinoTrueEAfterExtra,"signal"), Metric.second.GetHistogram(NeutrinoTrueEAfter,"background"),Metric.second.GetHistogram(NeutrinoTrueEBefore,"background"),Metric.second.GetHistogram(NeutrinoTrueEAfterExtra,"background"));



    }

    if(Metric.second.GetAnalysisName() == "MVAAllsCut"){

      TString MVAMethod = Metric.second.GetMVAName();

      std::string NeutrinoTrueEBeforeName = "NeutrinoTrueEBefore" + (std::string) MVAMethod; 
      std::string NeutrinoRecoEBeforeName = "NeutrinoRecoEBeforeName" +(std::string) MVAMethod;
      std::string NeutrinoRecoEBeforeExtraName = "NeutrinoRecoEBeforeExtraName" +(std::string) MVAMethod;
      std::string NeutrinoRecoEAfterName = "NeutrinoRecoEAfterName" +(std::string) MVAMethod;
      std::string NeutrinoTrueEAfterName = "NeutrinoTrueEAfterName" +(std::string) MVAMethod;
      std::string NeutrinoRecoEAfterExtraName = "NeutrinoRecoEAfterExtraName" +(std::string) MVAMethod;
      std::string NeutrinoTrueEAfterExtraName = "NeutrinoTrueEAfterExtraName" +(std::string) MVAMethod;
      

      Metric.second.MakeEfficiencyPlots(Metric.second.GetHistogram(NeutrinoRecoEAfterName,"signal"),Metric.second.GetHistogram(NeutrinoRecoEBeforeName,"signal"),Metric.second.GetHistogram(NeutrinoRecoEAfterExtraName,"signal"));
      
      Metric.second.MakeEfficiencyPlots(Metric.second.GetHistogram(NeutrinoRecoEAfterName,"background"),Metric.second.GetHistogram(NeutrinoRecoEBeforeName,"background"),Metric.second.GetHistogram(NeutrinoRecoEAfterExtraName,"background"));
      Metric.second.MakeEfficiencyPlots(Metric.second.GetHistogram(NeutrinoTrueEAfterName,"signal"),Metric.second.GetHistogram(NeutrinoTrueEBeforeName,"signal"),Metric.second.GetHistogram(NeutrinoTrueEAfterExtraName,"signal"));
      Metric.second.MakeEfficiencyPlots(Metric.second.GetHistogram(NeutrinoTrueEAfterName,"background"),Metric.second.GetHistogram(NeutrinoTrueEBeforeName,"background"),Metric.second.GetHistogram(NeutrinoTrueEAfterExtraName,"background"));
    }


  }
  
  file->Close();
  return;
}



int optimiser::NueRecoOptimiser::PerformCuts(optimiser::MetricHolder& Metric, std::string sample,std::vector<float>& vals, std::vector<std::pair<float,float> >& TwoDvals, TString& MVAName){

  int err=0;

  //Make Efficiency Hists
  //Make Resolutions Hits

  if(Metric.GetAnalysisName() == "StandardShowerCutFinder"){err = StandardDaughterAnalysis(Metric,sample,vals,TwoDvals);}
  else if(Metric.GetAnalysisName() == "StandardNeutrinoCutFinder"){err = StandardNeutrinoAnalysis(Metric,sample,vals,TwoDvals);}
  else if(Metric.GetAnalysisName() == "StandardNeutrinoDaughterCutFinder"){err = StandardNeutrinoDaughterAnalysis(Metric,sample,vals,TwoDvals);}
  else if(Metric.GetAnalysisName() == "StandardNumShowerAnalysis"){err = StandardNumShowerAnalysis(Metric,sample,TwoDvals);}
  else if(Metric.GetAnalysisName() == "OneShower2DEnergyAnalysis"){err = OneShower2DEnergyAnalysis(Metric,sample,TwoDvals);}  
  else if(Metric.GetAnalysisName() == "ShowerResidualAnalysis"){ err = ShowerResidualAnalysis(Metric,sample,TwoDvals);}
  else if(Metric.GetAnalysisName() == "DaughterComparision"){err = StandardDaughterAnalysis(Metric,sample,vals,TwoDvals);}
  else if(Metric.GetAnalysisName() == "NeutrinoComparision"){err = StandardNeutrinoAnalysis(Metric,sample,vals,TwoDvals);}
  else if(Metric.GetAnalysisName() == "NeutrinoDaughterComparision"){err = StandardNeutrinoDaughterAnalysis(Metric,sample,vals,TwoDvals);}
  else if(Metric.GetAnalysisName() == "ShowerResidualCutFinder"){err = ShowerResidualCut(Metric,sample,vals);}
  else if(Metric.GetAnalysisName() == "ShowerAllsCut"){err = ShowerAllsCut(Metric,sample);}
  else if(Metric.GetAnalysisName() == "OneShowerSizeNeutrinoAnalysisCutFinder"){err =StandardSizeNeutrinoAnalysis(Metric,sample,vals,TwoDvals);}
  else if(Metric.GetAnalysisName() == "StandardOverEDaughterCutFinder"){err = StandardOverEDaughterAnalysis(Metric,sample,vals,TwoDvals);}
  else if(Metric.GetAnalysisName() == "StandardAboveEDaughterCutFinder"){err = StandardAboveEDaughterAnalysis(Metric,sample,vals,TwoDvals);}
  else if(Metric.GetAnalysisName() == "StandardAboveOverEDaughterCutFinder"){err = StandardAboveOverEDaughterAnalysis(Metric,sample,vals,TwoDvals);}
  else if(Metric.GetAnalysisName() == "StandardMaxDaughterCutFinder"){err = StandardMaxDaughterAnalysis(Metric,sample,vals,TwoDvals);}
  else if(Metric.GetAnalysisName() == "ShowerResidualNumShowersCutFinder"){err = ShowerResidualNumShowers(Metric,sample,vals);}
  //else if(Metric.GetAnalysisName() == "MVAAnalysisCutFinder"){err = MVAAnalysis(Metric,sample,vals,MVAName);}
  else if(Metric.GetAnalysisName() == "MVAAllsCut"){err = MVAAnalysisAllsCut(Metric,sample,MVAName);}
  else if(Metric.GetAnalysisName() == "EfficiencyNeutrinoMaker"){err = EfficiencyNeutrinoMaker(Metric,sample);}
  else if(Metric.GetAnalysisName() == "OneShower1DCutFinder"){err = OneShower1DCutAnalysis(Metric,sample,vals);}
  else if(Metric.GetAnalysisName() == "StandardNeutrinoEnergyResolutionProfile"){err  = StandardNeutrinoEnergyResolutionProfile(Metric,sample,TwoDvals);}
  else if(Metric.GetAnalysisName() == "BDTAllsCutFinder"){err = BDTAllsCut(Metric,sample,vals);}

  return err;
}

template <class T>
int optimiser::NueRecoOptimiser::GetBranchPointer(std::string& branchname, std::string& sample, T* &vals){

  optimiser::BranchTypeBase* branchbase;

  //Decide if it is a signal or background branch
  if(sample == "signal"){
    if(SignalTree.Tree.find(branchname) == SignalTree.Tree.end()){
      std::cerr << "Signal Tree branch " << branchname << " doesn't exist. You probably mistypes the name in the metric initialser" << std::endl;
      return 1;
    }
    branchbase = SignalTree.Tree[branchname];
  }
  else{
    if(BackgroundTree.Tree.find(branchname) == BackgroundTree.Tree.end()){
      std::cerr << "Background Tree branch " << branchname << " doesn't exist. You probably mistypes the name in the metric initialser" << std::endl;
      return 1;
    }
    branchbase = BackgroundTree.Tree[branchname];
  }

  //Get the branch.
  optimiser::BranchType<T>* branch  =  dynamic_cast<optimiser::BranchType<T>*>(branchbase);
  if(branch == NULL){
    std::cerr << "could not find branch: " << branchname << " with type given " << typeid(*vals).name() << std::endl;
    return 1;
  }
  //Actual values.
  vals = branch->GetPointer();
  return 0;
}

//#######################################
//######### Standard Daughter ###########
//#######################################


int optimiser::NueRecoOptimiser::StandardDaughterAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals,  std::vector<std::pair<float,float> >& TwoDvals){

  int err = StandardDaughter1DAnalysis(Metric,sample,vals);
  if(err){
    std::cerr << "Error in 1D daughter analysis and bailing" << std::endl;
    return 1;
  }

  if(Metric.GetPartnerMetricName() == ""){return 0;}

  err =  StandardDaughter2DAnalysis(Metric,sample,TwoDvals);
  if(err){
    std::cerr << "Error in 2D daughter analysis and bailing" << std::endl;
    return 1;
  }
  return 0;
}


int optimiser::NueRecoOptimiser::StandardDaughter1DAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals){

  std::string branchname = Metric.GetBranchName();
  std::vector<std::vector<float> >* branch_vals;  
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  //Loop over the neutrinos
  int i=0;
  for(auto const& neutrino: *branch_vals){
    //Add the largest shower in the event. The metric are ordered this way
    if(neutrino.size() > 0){
      
      if(fApplyFVCut){
	std::vector<float>* fv_vals; 
	std::string fv_branch_name = "in_FV";
	int err = GetBranchPointer(fv_branch_name,sample,fv_vals);
	if(err){ std::cout << "Did not get the FV branch" << std::endl; return 1;}
	if(fv_vals->size() < branch_vals->size()){std::cout << "branch val not the right size" << std::endl; return 1;}
	if(fv_vals->at(i) == 0){
	  vals.push_back(Metric.GetErrorVal());
	  ++i;
	  continue;
	}
	++i;
      }

      // int j=0;
      // for(auto const& val: neutrino){
      // 	if(j==0){++j; continue;}
      // 	vals.push_back(val);
      // 	++j;
      // }

      float maxval = neutrino[0];

      //      maxval = 1;
      //if(neutrino.size() !=1){maxval = 4;}
      
      
      //      float maxval = -999;
      //      for(auto const& val: neutrino){
      //      	if(maxval < val){
      //      	  maxval = val;
      //      	}
      //      }

      vals.push_back(maxval);
    }
    else{
      vals.push_back(Metric.GetErrorVal());
    }
  }
  return 0;
}

int optimiser::NueRecoOptimiser::StandardDaughter2DAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<std::pair<float,float> >& TwoDvals){

  //Get the metric branch 
  std::string branchname = Metric.GetBranchName();
  std::vector<std::vector<float> >* branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  //Get the partner metric branch
  std::string partnername = Metric.GetPartnerMetricName();
  std::vector<std::vector<float> >* partner_vals;
  err = GetBranchPointer(partnername,sample,partner_vals);
  if(err){return 1;}

  if(partner_vals->size() != branch_vals->size()){
    std::cerr << "expecting same vector length for daughter and daughter info. Please Check" << std::endl;
  }

  //Loop over the neutrinos
  for(int iter=0; iter<branch_vals->size(); ++iter){

      if(fApplyFVCut){
	std::vector<float>* fv_vals; 
	std::string fv_branch_name = "in_FV";
	int err = GetBranchPointer(fv_branch_name,sample,fv_vals);
	if(err){ std::cout << "Did not get the FV branch" << std::endl; return 1;}
	if(fv_vals->size() < branch_vals->size()){std::cout << "branch val not the right size" << std::endl; return 1;}
	if(fv_vals->at(iter) == 0){
	  std::pair<float,float> vals = std::make_pair(Metric.GetErrorVal(),Metric.GetErrorVal());
	  TwoDvals.push_back(vals);
	  continue;
	}
      }

    //Add the largest shower in the event. The metric are ordered this way
    if(branch_vals->at(iter).size() > 0 && partner_vals->at(iter).size() > 0){
      std::pair<float,float> vals = std::make_pair(branch_vals->at(iter).at(0),partner_vals->at(iter).at(0));
      TwoDvals.push_back(vals);
    }
    else{
      std::pair<float,float> vals = std::make_pair(Metric.GetErrorVal(),Metric.GetErrorVal());
      TwoDvals.push_back(vals);
    }

  }
  return 0;
}

//########################
//### Biggest Analysis ###
//########################


int optimiser::NueRecoOptimiser::StandardMaxDaughterAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals,  std::vector<std::pair<float,float> >& TwoDvals){

  int err = StandardMaxDaughter1DAnalysis(Metric,sample,vals);
  if(err){
    std::cerr << "Error in 1D daughter analysis and bailing" << std::endl;
    return 1;
  }

  if(Metric.GetPartnerMetricName() == ""){return 0;}

  err =  StandardMaxDaughter2DAnalysis(Metric,sample,TwoDvals);
  if(err){
    std::cerr << "Error in 2D daughter analysis and bailing" << std::endl;
    return 1;
  }
  return 0;
}


int optimiser::NueRecoOptimiser::StandardMaxDaughter1DAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals){

  std::string branchname = Metric.GetBranchName();
  std::vector<std::vector<float> >* branch_vals;  
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  //Loop over the neutrinos
  int i=0;
  for(auto const& neutrino: *branch_vals){

    if(fApplyFVCut){
      std::vector<float> * fv_vals; 
      std::string fv_branch_name = "in_FV";
      int err = GetBranchPointer(fv_branch_name,sample,fv_vals);
      if(err){return 1;}
      if(fv_vals->size()-1 < i){return 1;}
      if(fv_vals->at(i) == 0){
	vals.push_back(Metric.GetErrorVal());
	++i;
	continue;
      }
      ++i;
    }


    //Add the largest shower in the event. The metric are ordered this way
    if(neutrino.size() > 0){
      
      float maxval = -999;
      int max_val_iter = -999;
      int i = 0;
      for(auto const& val: neutrino){
      	if(maxval < val){
      	  maxval = val;
	  max_val_iter = i;
      	}
	++i;
      }
      vals.push_back(maxval);
    }
    else{
      //      vals.push_back(Metric.GetErrorVal());
      vals.push_back(0);
    }
  }
  return 0;
}

int optimiser::NueRecoOptimiser::StandardMaxDaughter2DAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<std::pair<float,float> >& TwoDvals){

  //Get the metric branch 
  std::string branchname = Metric.GetBranchName();
  std::vector<std::vector<float> >* branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  //Get the partner metric branch
  std::string partnername = Metric.GetPartnerMetricName();
  std::vector<std::vector<float> >* partner_vals;
  err = GetBranchPointer(partnername,sample,partner_vals);
  if(err){return 1;}

  if(partner_vals->size() != branch_vals->size()){
    std::cerr << "expecting same vector length for daughter and daughter info. Please Check" << std::endl;
  }

  //Loop over the neutrinos
  for(int iter=0; iter<branch_vals->size(); ++iter){
    
    if(fApplyFVCut){
      std::vector<float> * fv_vals; 
      std::string fv_branch_name = "in_FV";
      int err = GetBranchPointer(fv_branch_name,sample,fv_vals);
      if(err){return 1;}
      if(fv_vals->size()-1 < iter){return 1;}
      if(fv_vals->at(iter) == 0){
	//	std::pair<float,float> vals = std::make_pair(Metric.GetErrorVal(),Metric.GetErrorVal());
	//	TwoDvals.push_back(vals);
	continue;
      }
    }


    //Add the largest shower in the event. The metric are ordered this way
    if(branch_vals->at(iter).size() > 0 && partner_vals->at(iter).size() > 0){
      float maxval = -999;
      int max_iter = -999;
      int i = 0;
      for(auto const& val: branch_vals->at(iter)){
      	if(maxval < val){
      	  maxval = val;
	  max_iter = i;
      	}
	++i;
      }

      if(max_iter == -999 || max_iter > (branch_vals->at(iter).size()-1) || max_iter > (partner_vals->at(iter).size()-1)){
	//	std::pair<float,float> vals = std::make_pair(Metric.GetErrorVal(),Metric.GetErrorVal());
	std::pair<float,float> vals = std::make_pair(0,0);
	TwoDvals.push_back(vals);
      }
      else{
	std::pair<float,float> vals = std::make_pair(branch_vals->at(iter).at(max_iter),partner_vals->at(iter).at(max_iter));
	TwoDvals.push_back(vals);
      }
    }
    else{
      //      std::pair<float,float> vals = std::make_pair(Metric.GetErrorVal(),Metric.GetErrorVal());
      std::pair<float,float> vals = std::make_pair(0,0);  
      TwoDvals.push_back(vals);
    }

  }
  return 0;
}




//#############################################
//######### Standard OverE Daughter ###########
//#############################################


int optimiser::NueRecoOptimiser::StandardOverEDaughterAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals, std::vector<std::pair<float,float> >& TwoDvals){

  int err = StandardOverEDaughter1DAnalysis(Metric,sample,vals);
  if(err){
    std::cerr << "Error in 1D daughter analysis and bailing" << std::endl;
    return 1;
  }

  if(Metric.GetPartnerMetricName() == ""){return 0;}

  err =  StandardOverEDaughter2DAnalysis(Metric,sample,TwoDvals);
  if(err){
    std::cerr << "Error in 2D daughter analysis and bailing" << std::endl;
    return 1;
  }
  return 0;
}


int optimiser::NueRecoOptimiser::StandardOverEDaughter1DAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals){

  std::string branchname = Metric.GetBranchName();
  std::vector<std::vector<float> >* branch_vals;  
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  //Get the Energy
  std::string shower_energy = "shower_energy";
  std::vector<std::vector<float> >* shower_energy_vals;
  err = GetBranchPointer(shower_energy,sample,shower_energy_vals);
  if(err){return 1;}

  //Loop over the neutrinos
  for(int iter=0; iter<branch_vals->size(); ++iter){
    
    if(fApplyFVCut){
      std::vector<float> * fv_vals; 
      std::string fv_branch_name = "in_FV";
      int err = GetBranchPointer(fv_branch_name,sample,fv_vals);
      if(err){return 1;}
      if(fv_vals->size() < branch_vals->size()){return 1;}
      if(fv_vals->at(iter) == 0){
	vals.push_back(Metric.GetErrorVal());
	continue;
      }
    }


    //Add the largest shower in the event. The metric are ordered this way
    if(branch_vals->at(iter).size() > 0){
      vals.push_back(branch_vals->at(iter).at(0)/shower_energy_vals->at(iter).at(0));
    }
    else{
      vals.push_back(Metric.GetErrorVal());
    }
  }
  return 0;
}

int optimiser::NueRecoOptimiser::StandardOverEDaughter2DAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<std::pair<float,float> >& TwoDvals){

  //Get the metric branch 
  std::string branchname = Metric.GetBranchName();
  std::vector<std::vector<float> >* branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  //Get the partner metric branch
  std::string partnername = Metric.GetPartnerMetricName();
  std::vector<std::vector<float> >* partner_vals;
  err = GetBranchPointer(partnername,sample,partner_vals);
  if(err){return 1;}

  //Get the Energy
  std::string shower_energy = "shower_energy";
  std::vector<std::vector<float> >* shower_energy_vals;
  err = GetBranchPointer(shower_energy,sample,shower_energy_vals);
  if(err){return 1;}


  if(partner_vals->size() != branch_vals->size()){
    std::cerr << "expecting same vector length for daughter and daughter info. Please Check" << std::endl;
  }

  //Loop over the neutrinos
  for(int iter=0; iter<branch_vals->size(); ++iter){

    if(fApplyFVCut){
      std::vector<float> * fv_vals; 
      std::string fv_branch_name = "in_FV";
      int err = GetBranchPointer(fv_branch_name,sample,fv_vals);
      if(err){return 1;}
      if(fv_vals->size() < branch_vals->size()){return 1;}
      if(fv_vals->at(iter) == 0){
	std::pair<float,float> vals = std::make_pair(Metric.GetErrorVal(),Metric.GetErrorVal());
	TwoDvals.push_back(vals);
	continue;
      }
    }


    //Add the largest shower in the event. The metric are ordered this way
    if(branch_vals->at(iter).size() > 0 && partner_vals->at(iter).size() > 0){
      std::pair<float,float> vals = std::make_pair(branch_vals->at(iter).at(0)/shower_energy_vals->at(iter).at(0),partner_vals->at(iter).at(0));
      TwoDvals.push_back(vals);
    }
    else{
      std::pair<float,float> vals = std::make_pair(Metric.GetErrorVal(),Metric.GetErrorVal());

      TwoDvals.push_back(vals);
    }
  }
  return 0;
}

//###########################
//######## Above E ##########
//###########################


int optimiser::NueRecoOptimiser::StandardAboveEDaughterAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals, std::vector<std::pair<float,float> >& TwoDvals){

  int err = StandardAboveEDaughter1DAnalysis(Metric,sample,vals);
  if(err){
    std::cerr << "Error in 1D daughter analysis and bailing" << std::endl;
    return 1;
  }

  if(Metric.GetPartnerMetricName() == ""){return 0;}

  err =  StandardAboveEDaughter2DAnalysis(Metric,sample,TwoDvals);
  if(err){
    std::cerr << "Error in 2D daughter analysis and bailing" << std::endl;
    return 1;
  }
  return 0;
}


int optimiser::NueRecoOptimiser::StandardAboveEDaughter1DAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals){

  float fEnergyCut=225;

  std::string branchname = Metric.GetBranchName();
  std::vector<std::vector<float> >* branch_vals;  
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  //Get the Energy
  std::string shower_energy = "shower_energy";
  std::vector<std::vector<float> >* shower_energy_vals;
  err = GetBranchPointer(shower_energy,sample,shower_energy_vals);
  if(err){return 1;}

  //Loop over the neutrinos
  for(int iter=0; iter<branch_vals->size(); ++iter){
    //Add the largest shower in the event. The metric are ordered this way
    
    if(fApplyFVCut){
      std::vector<float> * fv_vals; 
      std::string fv_branch_name = "in_FV";
      int err = GetBranchPointer(fv_branch_name,sample,fv_vals);
      if(err){return 1;}
      if(fv_vals->size() < branch_vals->size()){return 1;}
      if(fv_vals->at(iter) == 0){
	vals.push_back(Metric.GetErrorVal());
	//	vals.push_back(10); 
	continue;
      }
    }


    if(branch_vals->at(iter).size() > 0){

      int numshowers = 0; 
      for(auto const& shower: shower_energy_vals->at(iter)){
	if(shower > fEnergyCut){
	  ++numshowers;
	}
      }
      //      if(numshowers == 0){
      //	vals.push_back(5);
      //      }
      //      else{
      //	vals.push_back(numshowers);
      //      }
      if(shower_energy_vals->at(iter).at(0) > fEnergyCut){
       	vals.push_back(branch_vals->at(iter).at(0));	
      }
      else{
       	vals.push_back(Metric.GetErrorVal());
      }
    }
    else{
      vals.push_back(Metric.GetErrorVal());
      //vals.push_back(9); 
    }
  }
  return 0;
}

int optimiser::NueRecoOptimiser::StandardAboveEDaughter2DAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<std::pair<float,float> >& TwoDvals){

  float fEnergyCut=10;//210

  //Get the metric branch 
  std::string branchname = Metric.GetBranchName();
  std::vector<std::vector<float> >* branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  //Get the partner metric branch
  std::string partnername = Metric.GetPartnerMetricName();
  std::vector<std::vector<float> >* partner_vals;
  err = GetBranchPointer(partnername,sample,partner_vals);
  if(err){return 1;}

  //Get the Energy
  std::string shower_energy = "shower_energy";
  std::vector<std::vector<float> >* shower_energy_vals;
  err = GetBranchPointer(shower_energy,sample,shower_energy_vals);
  if(err){return 1;}


  if(partner_vals->size() != branch_vals->size()){
    std::cerr << "expecting same vector length for daughter and daughter info. Please Check" << std::endl;
  }

  //Loop over the neutrinos
  for(int iter=0; iter<branch_vals->size(); ++iter){

    if(fApplyFVCut){
      std::vector<float> * fv_vals; 
      std::string fv_branch_name = "in_FV";
      int err = GetBranchPointer(fv_branch_name,sample,fv_vals);
      if(err){return 1;}
      if(fv_vals->size() < branch_vals->size()){return 1;}
      if(fv_vals->at(iter) == 0){
	std::pair<float,float> vals = std::make_pair(Metric.GetErrorVal(),Metric.GetErrorVal());
	TwoDvals.push_back(vals);
	continue;
      }
    }
    

    //Add the largest shower in the event. The metric are ordered this way
    if(branch_vals->at(iter).size() > 0 && partner_vals->at(iter).size() > 0){
      if(shower_energy_vals->at(iter).at(0) > fEnergyCut){
	std::pair<float,float> vals = std::make_pair(branch_vals->at(iter).at(0),partner_vals->at(iter).at(0));
	TwoDvals.push_back(vals);
      }
      else{
	std::pair<float,float> vals = std::make_pair(Metric.GetErrorVal(),Metric.GetErrorVal());
	TwoDvals.push_back(vals);
      }
    }
    else{
      std::pair<float,float> vals = std::make_pair(Metric.GetErrorVal(),Metric.GetErrorVal());

      TwoDvals.push_back(vals);
    }
  }
  return 0;
}


//###############################
//######## AboveOver E ##########
//###############################


int optimiser::NueRecoOptimiser::StandardAboveOverEDaughterAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals, std::vector<std::pair<float,float> >& TwoDvals){

  int err = StandardAboveOverEDaughter1DAnalysis(Metric,sample,vals);
  if(err){
    std::cerr << "Error in 1D daughter analysis and bailing" << std::endl;
    return 1;
  }

  if(Metric.GetPartnerMetricName() == ""){return 0;}

  err =  StandardAboveOverEDaughter2DAnalysis(Metric,sample,TwoDvals);
  if(err){
    std::cerr << "Error in 2D daughter analysis and bailing" << std::endl;
    return 1;
  }
  return 0;
}


int optimiser::NueRecoOptimiser::StandardAboveOverEDaughter1DAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals){

  float fEnergyCut = 0.1;

  std::string branchname = Metric.GetBranchName();
  std::vector<std::vector<float> >* branch_vals;  
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  //Get the Energy
  std::string shower_energy = "shower_energy";
  std::vector<std::vector<float> >* shower_energy_vals;
  err = GetBranchPointer(shower_energy,sample,shower_energy_vals);
  if(err){return 1;}

  //Loop over the neutrinos
  for(int iter=0; iter<branch_vals->size(); ++iter){

    if(fApplyFVCut){
      std::vector<float> * fv_vals; 
      std::string fv_branch_name = "in_FV";
      int err = GetBranchPointer(fv_branch_name,sample,fv_vals);
      if(err){return 1;}
      if(fv_vals->size() < branch_vals->size()){return 1;}
      if(fv_vals->at(iter) == 0){
	vals.push_back(Metric.GetErrorVal());
	continue;
      }
    }


    //Add the largest shower in the event. The metric are ordered this way
    if(branch_vals->at(iter).size() > 0){
      if(shower_energy_vals->at(iter).at(0) > fEnergyCut){
	vals.push_back(branch_vals->at(iter).at(0)/shower_energy_vals->at(iter).at(0));	
      }
      else{
	vals.push_back(Metric.GetErrorVal());
      }
    }
    else{
      vals.push_back(Metric.GetErrorVal());
    }
  }
  return 0;
}

int optimiser::NueRecoOptimiser::StandardAboveOverEDaughter2DAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<std::pair<float,float> >& TwoDvals){

  float fEnergyCut=210;

  //Get the metric branch 
  std::string branchname = Metric.GetBranchName();
  std::vector<std::vector<float> >* branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  //Get the partner metric branch
  std::string partnername = Metric.GetPartnerMetricName();
  std::vector<std::vector<float> >* partner_vals;
  err = GetBranchPointer(partnername,sample,partner_vals);
  if(err){return 1;}

  //Get the Energy
  std::string shower_energy = "shower_energy";
  std::vector<std::vector<float> >* shower_energy_vals;
  err = GetBranchPointer(shower_energy,sample,shower_energy_vals);
  if(err){return 1;}


  if(partner_vals->size() != branch_vals->size()){
    std::cerr << "expecting same vector length for daughter and daughter info. Please Check" << std::endl;
  }

  //Loop over the neutrinos
  for(int iter=0; iter<branch_vals->size(); ++iter){

    if(fApplyFVCut){
      std::vector<float> * fv_vals; 
      std::string fv_branch_name = "in_FV";
      int err = GetBranchPointer(fv_branch_name,sample,fv_vals);
      if(err){return 1;}
      if(fv_vals->size() < branch_vals->size()){return 1;}
      if(fv_vals->at(iter) == 0){
	std::pair<float,float> vals = std::make_pair(Metric.GetErrorVal(),Metric.GetErrorVal());
	TwoDvals.push_back(vals);
	continue;
      }
    }


    //Add the largest shower in the event. The metric are ordered this way
    if(branch_vals->at(iter).size() > 0 && partner_vals->at(iter).size() > 0){
      if(shower_energy_vals->at(iter).at(0) > fEnergyCut){
	std::pair<float,float> vals = std::make_pair(branch_vals->at(iter).at(0)/shower_energy_vals->at(iter).at(0),partner_vals->at(iter).at(0));
	TwoDvals.push_back(vals);
      }
      else{
	std::pair<float,float> vals = std::make_pair(Metric.GetErrorVal(),Metric.GetErrorVal());
	TwoDvals.push_back(vals);
      }
    }
    else{
      std::pair<float,float> vals = std::make_pair(Metric.GetErrorVal(),Metric.GetErrorVal());

      TwoDvals.push_back(vals);
    }
  }
  return 0;
}




//#######################################
//######### Standard Neutrino ###########
//#######################################

int optimiser::NueRecoOptimiser::StandardNeutrinoAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals, std::vector<std::pair<float,float> >& TwoDvals){

  int err = StandardNeutrino1DAnalysis(Metric,sample,vals);
  if(err){
    std::cerr << "Error or in 1D neutrino analysis and bailing" << std::endl;
    return 1;
  }

  if(Metric.GetPartnerMetricName() == ""){return 0;}

  err =  StandardNeutrino2DAnalysis(Metric,sample,TwoDvals);
  if(err){
    std::cerr << "Error or in 2D neutrino analysis and bailing" << std::endl;
    return 1;
  }
  return 0;
}


int optimiser::NueRecoOptimiser::StandardNeutrino1DAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals){

  //Get the values
  std::string branchname = Metric.GetBranchName();
  std::vector<float>* branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}
  
  //Loop over the neutrinos
  int i = 0;
  for(auto const& neutrino: *branch_vals){
    
    if(fApplyFVCut){
      std::vector<float> * fv_vals; 
      std::string fv_branch_name = "in_FV";
      int err = GetBranchPointer(fv_branch_name,sample,fv_vals);
      if(err){return 1;}
      if(fv_vals->size() < branch_vals->size()){return 1;}
      if(fv_vals->at(i) == 0){
	vals.push_back(Metric.GetErrorVal());
	++i;continue;
      }
    }
    ++i;
    
    //    if(neutrino != 1){vals.push_back(4); continue;}

    vals.push_back(neutrino);
  }
  return 0;
}

int optimiser::NueRecoOptimiser::StandardNeutrino2DAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<std::pair<float,float> >& TwoDvals){

  //Get the metric branch 
  std::string branchname = Metric.GetBranchName();
  std::vector<float>* branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  //Get the partner metric branch
  std::string partnername = Metric.GetPartnerMetricName();
  std::vector<float>* partner_vals;
  err = GetBranchPointer(partnername,sample,partner_vals);
  if(err){return 1;}

  if(partner_vals->size() != branch_vals->size()){
    std::cerr << "expecting same vector length for neutrino and neutrino info. Please Check" << std::endl;
  }


  //Loop over the neutrinos
  int i=0;
  for(int iter=0; iter<branch_vals->size(); ++iter){

    if(fApplyFVCut){
      std::vector<float> * fv_vals; 
      std::string fv_branch_name = "in_FV";
      int err = GetBranchPointer(fv_branch_name,sample,fv_vals);
      if(err){return 1;}
      if(fv_vals->size() < branch_vals->size()){return 1;}
      if(fv_vals->at(i) == 0){++i;continue;}
    }
    ++i;


    std::pair<float,float> vals = std::make_pair(branch_vals->at(iter),partner_vals->at(iter));
    TwoDvals.push_back(vals);
  }
  return 0;
}


int optimiser::NueRecoOptimiser::StandardNeutrinoEnergyResolutionProfile(optimiser::MetricHolder& Metric, std::string& sample, std::vector<std::pair<float,float> >& TwoDvals){

  //Get the metric branch 
  std::string branchname = Metric.GetBranchName();
  std::vector<float>* branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  //Get the partner metric branch
  std::string partnername = Metric.GetPartnerMetricName();
  std::vector<float>* partner_vals;
  err = GetBranchPointer(partnername,sample,partner_vals);
  if(err){return 1;}

  if(partner_vals->size() != branch_vals->size()){
    std::cerr << "expecting same vector length for neutrino and neutrino info. Please Check" << std::endl;
  }
  

  //Get the partner metric branch
  std::string partnername_reco = "nu_reco_energy";
  std::vector<float>* partner_vals_recoE;
  err = GetBranchPointer(partnername_reco,sample,partner_vals_recoE);
  if(err){return 1;}

  if(partner_vals->size() != branch_vals->size()){
    std::cerr << "expecting same vector length for neutrino and neutrino info. Please Check" << std::endl;
  }

    //Get the partner metric branch
  std::string partnername_int = "nu_interaction_type";
  std::vector<float>* partner_vals_inttype;
  err = GetBranchPointer(partnername_int,sample,partner_vals_inttype);
  if(err){return 1;}



  //Loop over the neutrinos
  int i=0;
  for(int iter=0; iter<branch_vals->size(); ++iter){

    if(fApplyFVCut){
      std::vector<float> * fv_vals; 
      std::string fv_branch_name = "in_FV";
      int err = GetBranchPointer(fv_branch_name,sample,fv_vals);
      if(err){return 1;}
      if(fv_vals->size() < branch_vals->size()){return 1;}
      if(fv_vals->at(i) == 0){++i;continue;}
    }
    ++i;
    
    //    if(partner_vals_recoE->at(iter) < 225){continue;}
    //if(partner_vals_inttype->at(iter) != 1001){continue;}

    float resolution = (partner_vals->at(iter)*1000 -  partner_vals_recoE->at(iter))/(partner_vals_recoE->at(iter));
    //    float resolution = partner_vals_recoE->at(iter);

    std::pair<float,float> vals = std::make_pair(branch_vals->at(iter),resolution);
    TwoDvals.push_back(vals);
  }
  return 0;
}


//############################################
//######### Standard Size Neutrino ###########
//############################################

int optimiser::NueRecoOptimiser::StandardSizeNeutrinoAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals, std::vector<std::pair<float,float> >& TwoDvals){

  int err = StandardSizeNeutrino1DAnalysis(Metric,sample,vals);
  if(err){
    std::cerr << "Error or in 1D neutrino analysis and bailing" << std::endl;
    return 1;
  }

  if(Metric.GetPartnerMetricName() == ""){return 0;}

  err =  StandardSizeNeutrino2DAnalysis(Metric,sample,TwoDvals);
  if(err){
    std::cerr << "Error or in 2D neutrino analysis and bailing" << std::endl;
    return 1;
  }
  return 0;
}




// int optimiser::NueRecoOptimiser::StandardShowerEnergyResolutionProfile(optimiser::MetricHolder& Metric, std::string& sample, std::vector<std::pair<float,float> >& TwoDvals){

//   //Get the metric branch 
//   std::string branchname = Metric.GetBranchName();
//   std::vector<float>* branch_vals;
//   int err = GetBranchPointer(branchname,sample,branch_vals);
//   if(err){return 1;}

//   //Get the partner metric branch
//   std::string partnername = Metric.GetPartnerMetricName();
//   std::vector<float>* partner_vals;
//   err = GetBranchPointer(partnername,sample,partner_vals);
//   if(err){return 1;}

//   if(partner_vals->size() != branch_vals->size()){
//     std::cerr << "expecting same vector length for neutrino and neutrino info. Please Check" << std::endl;
//   }
  

//   //Get the partner metric branch
//   std::string partnername_reco = "nu_reco_energy";
//   std::vector<float>* partner_vals_recoE;
//   err = GetBranchPointer(partnername_reco,sample,partner_vals_recoE);
//   if(err){return 1;}

//   if(partner_vals->size() != branch_vals->size()){
//     std::cerr << "expecting same vector length for neutrino and neutrino info. Please Check" << std::endl;
//   }

//     //Get the partner metric branch
//   std::string partnername_int = "nu_interaction_type";
//   std::vector<float>* partner_vals_inttype;
//   err = GetBranchPointer(partnername_int,sample,partner_vals_inttype);
//   if(err){return 1;}



//   //Loop over the neutrinos
//   int i=0;
//   for(int iter=0; iter<branch_vals->size(); ++iter){

//     if(fApplyFVCut){
//       std::vector<float> * fv_vals; 
//       std::string fv_branch_name = "in_FV";
//       int err = GetBranchPointer(fv_branch_name,sample,fv_vals);
//       if(err){return 1;}
//       if(fv_vals->size() < branch_vals->size()){return 1;}
//       if(fv_vals->at(i) == 0){++i;continue;}
//     }
//     ++i;
    
//     if(partner_vals_recoE->at(iter) < 225){continue;}
//     //    if(partner_vals_inttype->at(iter) != 1001){continue;}

//     float resolution = (partner_vals->at(iter)*1000 -  partner_vals_recoE->at(iter))/(1000*partner_vals->at(iter));

//     std::pair<float,float> vals = std::make_pair(branch_vals->at(iter),resolution);
//     TwoDvals.push_back(vals);
//   }
//   return 0;
// }


int optimiser::NueRecoOptimiser::StandardSizeNeutrino1DAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals){

  //Get the values
  std::string branchname = Metric.GetBranchName();
  std::vector<float>* branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  //Loop over the neutrinos
  int i=0;
  for(int iter=0; iter<branch_vals->size(); ++iter){

    if(fApplyFVCut){
      std::vector<float> * fv_vals; 
      std::string fv_branch_name = "in_FV";
      int err = GetBranchPointer(fv_branch_name,sample,fv_vals);
      if(err){return 1;}
      if(fv_vals->size() < branch_vals->size()){return 1;}
      if(fv_vals->at(iter) == 0){++i;continue;}
    }
    ++i;
  
    vals.push_back(branch_vals->size());
  }

  return 0;
}

int optimiser::NueRecoOptimiser::StandardSizeNeutrino2DAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<std::pair<float,float> >& TwoDvals){

  //Get the metric branch 
  std::string branchname = Metric.GetBranchName();
  std::vector<float>* branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  //Get the partner metric branch
  std::string partnername = Metric.GetPartnerMetricName();
  std::vector<float>* partner_vals;
  err = GetBranchPointer(partnername,sample,partner_vals);
  if(err){return 1;}

  //Loop over the neutrinos
  std::pair<float,float> vals = std::make_pair(branch_vals->size(),partner_vals->size());
  TwoDvals.push_back(vals);

  return 0;
}


//##################################################
//######### Standard Neutrino & Daughter ###########
//##################################################

int optimiser::NueRecoOptimiser::StandardNeutrinoDaughterAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals, std::vector<std::pair<float,float> >& TwoDvals){

  if(Metric.GetPartnerMetricName() == ""){return 1;}

  int err =  StandardNeutrinoDaughter2DAnalysis(Metric,sample,TwoDvals);
  if(err){
    std::cerr << "Error or in 2D daughter+neutrino analysis and bailing" << std::endl;
    return 1;
  }
  return 0;
}


int optimiser::NueRecoOptimiser::StandardNeutrinoDaughter2DAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<std::pair<float,float> >& TwoDvals){

  //Get the metric branch 
  std::string branchname = Metric.GetBranchName();
  std::vector<std::vector<float> >* branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  //Get the partner metric branch
  std::string partnername = Metric.GetPartnerMetricName();
  std::vector<float>* partner_vals;
  err = GetBranchPointer(partnername,sample,partner_vals);
  if(err){return 1;}

  if(partner_vals->size() != branch_vals->size()){
    std::cerr << "expecting same vector length for neutrino and neutrino daughter info. Please Check" << std::endl;
    return 1;
  }

  //Loop over the neutrinos
  for(int iter=0; iter<branch_vals->size(); ++iter){

    if(fApplyFVCut){
      std::vector<float> * fv_vals; 
      std::string fv_branch_name = "in_FV";
      int err = GetBranchPointer(fv_branch_name,sample,fv_vals);
      if(err){return 1;}
      if(fv_vals->size() < branch_vals->size()){return 1;}
      if(fv_vals->at(iter) == 0){      
	std::pair<float,float> vals = std::make_pair(Metric.GetErrorVal(),Metric.GetErrorVal());
	TwoDvals.push_back(vals);
	continue;
      }
    }

    //Add the largest shower in the event. The metric are ordered this way
    if(branch_vals->at(iter).size() > 0){
      std::pair<float,float> vals = std::make_pair(branch_vals->at(iter).at(0),partner_vals->at(iter));
      TwoDvals.push_back(vals);
    }
    else{
      std::pair<float,float> vals = std::make_pair(Metric.GetErrorVal(),Metric.GetErrorVal());

      TwoDvals.push_back(vals);
    }
  }
  return 0;
}

//######################
//### Num Shower Ana ###
//######################

int optimiser::NueRecoOptimiser::StandardNumShowerAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<std::pair<float,float> >& TwoDvals){

  //Get the metric branch 
  std::string branchname = Metric.GetBranchName();
  std::vector<std::vector<float> >* branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  int numshowerbins = Metric.GetSignalHist()->GetNbinsX();
  
  //Loop over the neutrinos
  for(int bin=0; bin<numshowerbins; ++bin){
  
    float binenergy = ((TAxis*)Metric.GetSignalHist()->GetXaxis())->GetBinCenter(bin);
  
    for(auto const& neutrino: *branch_vals){

      int  numshowers = 0;
      int iter =0;
      for(auto const& shower_energy: neutrino){

	if(iter == 0){
	  ++iter;
	  continue;
	}
	++iter;
	
	if(shower_energy > binenergy){
	++numshowers;
	}
	//	if(shower_energy < binenergy){
	  //break;
	  //	}
	++numshowers;
      }

      if(numshowers == 0){numshowers =8;}
      if(neutrino.size() == 0){numshowers =8;}

      std::pair<float,float> vals = std::make_pair(binenergy,numshowers);
      TwoDvals.push_back(vals);
    }
  }
  return 0;
}

//######################
//### Num Shower Ana ###
//######################

int optimiser::NueRecoOptimiser::OneShower1DCutAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals){

  float fEnergyCut = 210;

  //Get the shower energy 
  std::string shower_energy = "shower_energy";
  std::vector<std::vector<float> >* shower_energy_vals;
  int err = GetBranchPointer(shower_energy,sample,shower_energy_vals);
  if(err){return 1;}
  
  //Get the metric branch 
  std::string branchname = Metric.GetBranchName();
  std::vector<std::vector<float> >* branch_vals;
  err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  for(int iter=0; iter<branch_vals->size(); ++iter){

    std::vector<float> shower_energy_vals_neutrino = shower_energy_vals->at(iter);
    int  numshowers = 0;
    for(auto const& shower_energy: shower_energy_vals_neutrino){
      if(shower_energy > fEnergyCut){
        ++numshowers;
      }
    }
    
    if(numshowers != 1){
      vals.push_back(Metric.GetErrorVal());
      continue;
    }
    
    if(fApplyFVCut){
      std::vector<float> * fv_vals; 
      std::string fv_branch_name = "in_FV";
      int err = GetBranchPointer(fv_branch_name,sample,fv_vals);
      if(err){return 1;}
      if(fv_vals->size() < branch_vals->size()){return 1;}
      if(fv_vals->at(iter) == 0){
	vals.push_back(Metric.GetErrorVal());
	continue;
      }
    }

  
    //Add the largest shower in the event. The metric are ordered this way
    float maxmetric = -999;
    for(auto const& metric: branch_vals->at(iter)){
      if(metric > maxmetric){
	maxmetric = metric;
      }
    }
     
    if(maxmetric == -999){
      vals.push_back(Metric.GetErrorVal());
      continue;
    }
    vals.push_back(maxmetric);
    
    //    if(branch_vals->at(iter).size() > 0){
    //      vals.push_back(branch_vals->at(iter).at(0));
    //    }
  }
  return 0;



}

int optimiser::NueRecoOptimiser::OneShower2DEnergyAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<std::pair<float,float> >& TwoDvals){

  //Not optimised use sparingly

  //Get the metric branch 
  std::string branchname = Metric.GetBranchName();
  std::vector<std::vector<float> >* branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  //Get the partner metric branch
  std::string partnername = Metric.GetPartnerMetricName();
  std::vector<std::vector<float> >* partner_vals;
  err = GetBranchPointer(partnername,sample,partner_vals);
  if(err){return 1;}

  //Get the partner metric branch
  std::string partnername_len = "shower_len";
  std::vector<std::vector<float> >* partner_vals_len;
  err = GetBranchPointer(partnername_len,sample,partner_vals_len);
  if(err){return 1;}

  //Get the partner metric branch
  std::string partnername_energy = "shower_energy";
  std::vector<std::vector<float> >* partner_vals_energy;
  err = GetBranchPointer(partnername_energy,sample,partner_vals_energy);
  if(err){return 1;}


  //Get the partner metric branch
  std::string partnername_openingangle = "shower_open_angle";
  std::vector<std::vector<float> >* partner_vals_openingangle;
  err = GetBranchPointer(partnername_openingangle,sample,partner_vals_openingangle);
  if(err){return 1;}

  std::string branchnamecc = "nu_cc";
  std::vector<float>* branch_valscc;
  err = GetBranchPointer(branchnamecc,sample,branch_valscc);
  if(err){return 1;}

  if(partner_vals->size() != branch_vals->size()){
    std::cerr << "expecting same vector length for neutrino and neutrino daughter info. Please Check" << std::endl;
  }

  int numshowerbins = Metric.GetSignalHist()->GetNbinsX();
  int numshowerbins_partner = Metric.GetSignalHistPartner()->GetNbinsY();
  
  //Loop over the neutrinos
  for(int bin_p=0; bin_p<numshowerbins_partner; ++bin_p){
    float fCut = ((TAxis*)Metric.GetSignalHistPartner()->GetYaxis())->GetBinCenter(bin_p);
    for(int bin=0; bin<numshowerbins; ++bin){
  
      float binenergy = ((TAxis*)Metric.GetSignalHist()->GetXaxis())->GetBinCenter(bin);
 
      for(int neutrino_iter=0; neutrino_iter<branch_vals->size(); ++neutrino_iter){

	if(branch_valscc->at(neutrino_iter) != 0 && sample=="signal"){continue;}
	
	if(fApplyFVCut){
	  std::vector<float> * fv_vals; 
	  std::string fv_branch_name = "in_FV";
	  int err = GetBranchPointer(fv_branch_name,sample,fv_vals);
	  if(err){return 1;}
	  if(fv_vals->size() < branch_vals->size()){return 1;}
	  if(fv_vals->at(neutrino_iter) == 0){
	    //	    std::pair<float,float> vals = std::make_pair(Metric.GetErrorVal(),Metric.GetErrorVal());
	    //	    TwoDvals.push_back(vals);
	    continue;
	  }
	}
	
	//No shower no cut continue;
	if(branch_vals->at(neutrino_iter).size() == 0){
	  std::pair<float,float> vals = std::make_pair(Metric.GetErrorVal(),Metric.GetErrorVal());
	  TwoDvals.push_back(vals);
	  continue;
	}

	//Count the number of showers.
	int numshowers = 0;
	int shower_iter_to_ana = -999;
	for(int shower_iter=0; shower_iter<branch_vals->at(neutrino_iter).size(); ++shower_iter){

	 
	  float shower_energy =  partner_vals_energy->at(neutrino_iter).at(shower_iter); //branch_vals->at(neutrino_iter).at(shower_iter);
	  float showerdist = partner_vals_len->at(neutrino_iter).at(shower_iter) * TMath::Tan(0.5*partner_vals_openingangle->at(neutrino_iter).at(0));

	  if(shower_iter==0){continue;}
	  
	  if(shower_energy > binenergy && partner_vals->at(neutrino_iter).at(shower_iter) > TMath::Abs(fCut*showerdist)+2){
	  //  if(partner_vals->at(neutrino_iter).at(shower_iter) > TMath::Abs(fCut*showerdist)+binenergy){
	    ++numshowers;
	    shower_iter_to_ana = shower_iter;
	  }
	}

	//Only Keep one shower events
	if(numshowers > 0){
	  std::pair<float,float> vals = std::make_pair(Metric.GetErrorVal(),Metric.GetErrorVal());
	  TwoDvals.push_back(vals);
	  continue;
	}
	

	//	//Only analyse match showers
	//	if(shower_iter_to_ana == -999 || partner_vals->at(neutrino_iter).size() < shower_iter_to_ana){
	//	  std::pair<float,float> vals = std::make_pair(Metric.GetErrorVal(),Metric.GetErrorVal());
	//	  TwoDvals.push_back(vals);
	//	  continue;
	//	}

	//Cut is unfortately done here 
	//if(partner_vals->at(neutrino_iter).at(shower_iter_to_ana) > fCut){

	//std::pair<float,float> vals = std::make_pair(Metric.GetErrorVal(),Metric.GetErrorVal());
	//	  TwoDvals.push_back(vals);
	//	  continue;
	//	}
      
	//Add the highest energy shower and its counter part 
	std::pair<float,float> vals = std::make_pair(binenergy,fCut);
	TwoDvals.push_back(vals);
      }
    }
  }
  return 0;
}


//######################
//### Shower Res Ana ###
//######################

int optimiser::NueRecoOptimiser::ShowerResidualAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<std::pair<float,float> >& TwoDvals){

  //Get the metric branch 
  std::string branchname = Metric.GetBranchName();
  std::vector<std::vector<float> >* branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  int numshowerbins = Metric.GetSignalHist()->GetNbinsX();

  std::string shower_energy = "shower_energy";
  std::vector<std::vector<float> >* energy_vals;
  err = GetBranchPointer(shower_energy,sample,energy_vals);
  if(err){return 1;}

  if(energy_vals->size() != branch_vals->size()){
    std::cerr << "expecting same vector length for neutrino and neutrino daughter info. Please Check" << std::endl;
  }

  //Get the partner metric branch
  std::string partnername_len = "shower_len";
  std::vector<std::vector<float> >* partner_vals_len;
  err = GetBranchPointer(partnername_len,sample,partner_vals_len);
  if(err){return 1;}

  //Get the partner metric branch
  std::string partnername_openingangle = "shower_open_angle";
  std::vector<std::vector<float> >* partner_vals_openingangle;
  err = GetBranchPointer(partnername_openingangle,sample,partner_vals_openingangle);
  if(err){return 1;}


  
  //Loop over the neutrinos
  for(int bin=0; bin<numshowerbins; ++bin){
  
    float binresidual = ((TAxis*)Metric.GetSignalHist()->GetXaxis())->GetBinCenter(bin);
  
    for(int n_iter=0; n_iter<branch_vals->size(); ++n_iter){
      int  numshowers = 0;

      if(energy_vals->at(n_iter).size() != branch_vals->at(n_iter).size()){
	std::cerr << "energy_vals  size: " << energy_vals->at(n_iter).size() << " branch_vals->at(n_iter).size(): " << branch_vals->at(n_iter).size() << std::endl;
	std::pair<float,float> vals = std::make_pair(Metric.GetErrorVal(),Metric.GetErrorVal());
	TwoDvals.push_back(vals);
	continue;
      }

      
      for(int s_iter=0; s_iter<branch_vals->at(n_iter).size(); ++s_iter){
	//    for(auto const& neutrino: *branch_vals){
	//      for(auto const& shower_residual: neutrino){
	//if(shower_residual > binresidual){ //|| shower_residual<1e-5){


	if(branch_vals->at(n_iter).at(s_iter) > binresidual && energy_vals->at(n_iter).at(s_iter) > 50){
	  ++numshowers;
	}
      }
      
      //Only Keep one shower events
      // if(numshowers > 1){
      // 	std::pair<float,float> vals = std::make_pair(-999,-999);
      // 	TwoDvals.push_back(vals);
      // 	continue;
      // }
      
      std::pair<float,float> vals = std::make_pair(binresidual,numshowers);
      TwoDvals.push_back(vals);
    }
  }
  return 0;
}

int optimiser::NueRecoOptimiser::ShowerResidualCut(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals){

  //Get the metric branch 
  std::string branchname = Metric.GetBranchName();
  std::vector<std::vector<float> >* branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  int numshowerbins = Metric.GetSignalHist()->GetNbinsX();

  std::string shower_energy = "shower_energy";
  std::vector<std::vector<float> >* energy_vals;
  err = GetBranchPointer(shower_energy,sample,energy_vals);
  if(err){return 1;}

  std::string shower_residual_dist = "shower_residual_dist";
  std::vector<std::vector<float> >* shower_residual_vals;
  err = GetBranchPointer(shower_residual_dist,sample,shower_residual_vals);
  if(err){return 1;}


  if(energy_vals->size() != branch_vals->size()){
    std::cerr << "expecting same vector length for neutrino and neutrino daughter info. Please Check" << std::endl;
  }

  
  //Loop over the neutrinos
  int i=0;
  for(int n_iter=0; n_iter<branch_vals->size(); ++n_iter){
    int  numshowers = 0;
    
    if(fApplyFVCut){
      std::vector<float> * fv_vals; 
      std::string fv_branch_name = "in_FV";
      int err = GetBranchPointer(fv_branch_name,sample,fv_vals);
      if(err){return 1;}
      if(fv_vals->size()-1 < i){return 1;}
      if(fv_vals->at(n_iter) == 0){
	//	vals.push_back(Metric.GetErrorVal());
	++i;
	continue;
      }
      ++i;
    }

    
    if(energy_vals->at(n_iter).size() != branch_vals->at(n_iter).size()){
      std::cerr << "energy_vals  size: " << energy_vals->at(n_iter).size() << " branch_vals->at(n_iter).size(): " << branch_vals->at(n_iter).size() << std::endl;
      vals.push_back(Metric.GetErrorVal());
      continue;
    }
    
    if(branch_vals->at(n_iter).size() == 0){continue;}


    for(int s_iter=0; s_iter<branch_vals->at(n_iter).size(); ++s_iter){
      if(shower_residual_vals->at(n_iter).at(s_iter) > 1.25 && energy_vals->at(n_iter).at(s_iter) > 10){
	++numshowers;
      }
    }
    
    //Only Keep one shower events
    if(numshowers > 1){
      vals.push_back(Metric.GetErrorVal());
      continue;
    }
    
    vals.push_back(branch_vals->at(n_iter).at(0));
  }

return 0;
}

int optimiser::NueRecoOptimiser::ShowerResidualNumShowers(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals){

  //Get the metric branch
  std::string branchname = Metric.GetBranchName();
  std::vector<std::vector<float> >* branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  int numshowerbins = Metric.GetSignalHist()->GetNbinsX();

  std::string shower_energy = "shower_energy";
  std::vector<std::vector<float> >* energy_vals;
  err = GetBranchPointer(shower_energy,sample,energy_vals);
  if(err){return 1;}

  std::string shower_residual_dist = "shower_residual_dist";
  std::vector<std::vector<float> >* shower_residual_vals;
  err = GetBranchPointer(shower_residual_dist,sample,shower_residual_vals);
  if(err){return 1;}


  if(energy_vals->size() != branch_vals->size()){
    std::cerr << "expecting same vector length for neutrino and neutrino daughter info. Please Check" << std::endl;
  }

  //Get the partner metric branch
  std::string partnername_len = "shower_len";
  std::vector<std::vector<float> >* partner_vals_len;
  err = GetBranchPointer(partnername_len,sample,partner_vals_len);
  if(err){return 1;}

  //Get the partner metric branch
  std::string partnername_openingangle = "shower_open_angle";
  std::vector<std::vector<float> >* partner_vals_openingangle;
  err = GetBranchPointer(partnername_openingangle,sample,partner_vals_openingangle);
  if(err){return 1;}



  
  //Loop over the neutrinos
  int i=0;
  for(int n_iter=0; n_iter<branch_vals->size(); ++n_iter){
    int  numshowers = 0;

    if(fApplyFVCut){
      std::vector<float> * fv_vals; 
      std::string fv_branch_name = "in_FV";
      int err = GetBranchPointer(fv_branch_name,sample,fv_vals);
      if(err){return 1;}
      if(fv_vals->size()-1 < i){return 1;}
      if(fv_vals->at(n_iter) == 0){
	vals.push_back(Metric.GetErrorVal());
	++i;
	continue;
      }
      ++i;
    }

    
    if(energy_vals->at(n_iter).size() != branch_vals->at(n_iter).size()){
      std::cerr << "energy_vals  size: " << energy_vals->at(n_iter).size() << " branch_vals->at(n_iter).size(): " << branch_vals->at(n_iter).size() << std::endl;
      vals.push_back(Metric.GetErrorVal());
      continue;
    }
    
    if(branch_vals->at(n_iter).size() == 0){
      //vals.push_back(Metric.GetErrorVal());
       vals.push_back(8); 
      continue;
    }
    
    int showere =0;
    for(int s_iter=0; s_iter<branch_vals->at(n_iter).size(); ++s_iter){

      if(energy_vals->at(n_iter).at(s_iter) > 10){
	++showere;
      }


      if(s_iter == 0){continue;}
      float showerdist = partner_vals_len->at(n_iter).at(s_iter) * TMath::Tan(0.5*partner_vals_openingangle->at(n_iter).at(0));

      ///      if(energy_vals->at(n_iter).at(s_iter) > 10 && shower_residual_vals->at(n_iter).at(s_iter) > TMath::Abs(1*showerdist)+2){
      if(energy_vals->at(n_iter).at(s_iter) > 10 && shower_residual_vals->at(n_iter).at(s_iter) > TMath::Abs(0.25*showerdist)+2)
	//if(shower_residual_vals->at(n_iter).at(s_iter) > 1 ){
	++numshowers;
      }
      //      vals.push_back(shower_residual_vals->at(n_iter).at(s_iter));
    //      std::cout << "showerdist: " << showerdist << std::endl;
    //    vals.push_back(shower_residual_vals->at(n_iter).at(s_iter));
    
    if(showere == 0){
      vals.push_back(8);
      continue;
    }

    vals.push_back(numshowers);
  }

  return 0;
}


int optimiser::NueRecoOptimiser::BDTAllsCut(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals){

  //Cuts 
  const float fResidualCutVal          = 0.025;//3.135;//0.025;
  const std::string fMVAMethod         = "BDTGCrossVal";
  const float fResidualEnergyCut       = 10;//10;23.75
  const float fEnergyCut               = 0;

  //Get the shower energy 
  std::string shower_energy = "shower_energy";
  std::vector<std::vector<float> >* shower_energy_vals;
  int err = GetBranchPointer(shower_energy,sample,shower_energy_vals);
  if(err){return 1;}

  //Get the residual distance 
  std::string shower_residual_dist = "shower_residual_dist";
  std::vector<std::vector<float> >* shower_residual_vals;
  err = GetBranchPointer(shower_residual_dist,sample,shower_residual_vals);
  if(err){return 1;}

  //Get the dEdx
  std::string shower_dEdx = "shower_dEdx";
  std::vector<std::vector<float> >* shower_dEdx_vals;
  err = GetBranchPointer(shower_dEdx,sample,shower_dEdx_vals);
  if(err){return 1;}

  //Get the conversion distance
  std::string shower_coversion_gap = "shower_coversion_gap";
  std::vector<std::vector<float> >* shower_coversion_gap_vals;
  err = GetBranchPointer(shower_coversion_gap,sample,shower_coversion_gap_vals);
  if(err){return 1;}

  //Get the max track length 
  std::string track_lengths = "track_lengths";
  std::vector<std::vector<float> >* shower_track_length_vals;
  err = GetBranchPointer(track_lengths,sample,shower_track_length_vals);
  if(err){return 1;}

  //Get the nu reco energy  
  std::string nu_reco_energy = "nu_reco_energy";
  std::vector<float>* nu_reco_energy_vals;
  err = GetBranchPointer(nu_reco_energy,sample,nu_reco_energy_vals);
  if(err){return 1;}

  //Get the nu truth energy  
  std::string nu_truth_energy = "nu_truth_energy";
  std::vector<float> * nu_truth_energy_vals;
  err = GetBranchPointer(nu_truth_energy,sample,nu_truth_energy_vals);
  if(err){return 1;}

  //Get the nu truth energy  
  std::string track_PIDAs = "track_PIDA_branch";
  std::vector<std::vector<float> >* track_PIDA_vals;
  err = GetBranchPointer(track_PIDAs,sample,track_PIDA_vals);
  if(err){return 1;}

  std::string shower_density_pw_news = "shower_density_pw_new";
  std::vector<std::vector<float> >* shower_density_pw_new_vals;
  err = GetBranchPointer(shower_density_pw_news,sample,shower_density_pw_new_vals);
  if(err){return 1;}

  std::string shower_density_grad_news = "shower_density_grad_new";
  std::vector<std::vector<float> >* shower_density_grad_new_vals;
  err = GetBranchPointer(shower_density_grad_news,sample,shower_density_grad_new_vals);
  if(err){return 1;}

  std::string shower_trackwidths = "shower_trackwidth";
  std::vector<std::vector<float> >* shower_trackwidth_vals;
  err = GetBranchPointer(shower_trackwidths,sample,shower_trackwidth_vals);
  if(err){return 1;}

  std::string shower_tracklengths = "shower_tracklength";
  std::vector<std::vector<float> >* shower_tracklength_vals;
  err = GetBranchPointer(shower_tracklengths,sample,shower_tracklength_vals);
  if(err){return 1;}

  std::string shower_lengths = "shower_length";
  std::vector<std::vector<float> >* shower_length_vals;
  err = GetBranchPointer(shower_lengths,sample,shower_length_vals);
  if(err){return 1;}

  std::string shower_open_angles = "shower_open_angle";
  std::vector<std::vector<float> >* shower_open_angle_vals;
  err = GetBranchPointer(shower_open_angles,sample,shower_open_angle_vals);
  if(err){return 1;}



  //Get the nu truth energy  
  std::string nu_E = "nu_E";
  std::vector<float> * nu_E_vals;
  err = GetBranchPointer(nu_E,sample,nu_E_vals);
  if(err){return 1;}

  //Get the nu truth energy  
  std::string nu_pdg = "nu_pdg";
  std::vector<float> * nu_pdg_vals;
  err = GetBranchPointer(nu_pdg,sample,nu_pdg_vals);
  if(err){return 1;}


  //Get the metric branch
  std::string branchname = Metric.GetBranchName();
  std::vector<float> * branch_vals;
  err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  //Get the partner metric branch
  std::string partnername_len = "shower_len";
  std::vector<std::vector<float> >* partner_vals_len;
  err = GetBranchPointer(partnername_len,sample,partner_vals_len);
  if(err){return 1;}

  //Get the partner metric branch
  std::string partnername_openingangle = "shower_open_angle";
  std::vector<std::vector<float> >* partner_vals_openingangle;
  err = GetBranchPointer(partnername_openingangle,sample,partner_vals_openingangle);
  if(err){return 1;}



  std::string branchnamecc = "nu_cc";
  std::vector<float>* branch_valscc;
  err = GetBranchPointer(branchnamecc,sample,branch_valscc);
  



  // //Loop over the neutrinos and get their true energy
  // for(int n_iter=0; n_iter<nu_E_vals->size(); ++n_iter){

  //   double oscprob = GetOscProb(n_iter,sample,Metric.ApplyOscProb());
  //   double potweight = 1;  
  //   if(Metric.ApplyPOT()){potweight = Metric.ReturnPOT(sample);}
  //   double oscpot = oscprob*potweight;
   
    
  //   Metric.GetHistogram(NeutrinoTrueEBefore,sample)->Fill(nu_truth_energy_vals->at(n_iter)*1000,oscpot);
  // }

  //Loop over the reco neutrinos and idnetify the reco neutrino with the most energy
  std::map<float,std::vector<float> > NeutrinoE;
  std::map<float,std::vector<float> > NeutrinoMVA;
  std::map<float,float > MaxNeutrinoE;
  std::map<float, int >  MaxNeutrinoE_iter;

  // for(int n_iter=0; n_iter<nu_E_vals->size(); ++n_iter){
  //   if(branch_valscc->at(n_iter) == 1 && sample == "signal"){continue;}

  //   MaxNeutrinoE[nu_E_vals->at(n_iter)]      = -99999;
  //   MaxNeutrinoE_iter[nu_E_vals->at(n_iter)] = -99999;   
  // }

  // for(int n_iter=0; n_iter<nu_reco_energy_vals->size(); ++n_iter){
  //   if(branch_valscc->at(n_iter) == 1 && sample == "signal"){continue;}

  //   //check if the true neutrino is already be identified 
  //   if(MaxNeutrinoE[nu_truth_energy_vals->at(n_iter)] < nu_reco_energy_vals->at(n_iter)){
  //     MaxNeutrinoE[nu_truth_energy_vals->at(n_iter)] = nu_reco_energy_vals->at(n_iter);
  //     MaxNeutrinoE_iter[nu_truth_energy_vals->at(n_iter)] = n_iter;
  //   }
  // }

  //Re initialise the maps
  // for(int n_iter=0; n_iter<nu_E_vals->size(); ++n_iter){
  //   if(branch_valscc->at(n_iter) == 1 && sample == "signal"){continue;}

  //   MaxNeutrinoE[nu_E_vals->at(n_iter)]      = -99999;
  //   MaxNeutrinoE_iter[nu_E_vals->at(n_iter)] = -99999;   
  //   NeutrinoMVA[nu_E_vals->at(n_iter)] = -99999; 
  // }


  //Loop over the neutrinos
  for(int n_iter=0; n_iter<nu_reco_energy_vals->size(); ++n_iter){
    
    if(branch_valscc->at(n_iter) == 1 && sample == "signal"){continue;}

    if(fApplyFVCut){
      std::vector<float> * fv_vals; 
      std::string fv_branch_name = "in_FV";
      int err = GetBranchPointer(fv_branch_name,sample,fv_vals);
      if(err){return 1;}
      if(fv_vals->size() < nu_reco_energy_vals->size()){return 1;}
      if(fv_vals->at(n_iter) == 0){
	continue;
      }
    }

    if(shower_energy_vals->at(n_iter).size() == 0){
      continue;
    }

    int  numshowers = 0;

    for(int s_iter=0; s_iter<shower_residual_vals->at(n_iter).size(); ++s_iter){
      float showerdist = partner_vals_len->at(n_iter).at(s_iter) * TMath::Tan(0.5*partner_vals_openingangle->at(n_iter).at(0));
      if(shower_energy_vals->at(n_iter).at(s_iter) > fResidualEnergyCut && shower_residual_vals->at(n_iter).at(s_iter) > TMath::Abs(fResidualCutVal*showerdist)+2){
	//	if(shower_residual_vals->at(n_iter).at(s_iter) > fResidualCut && shower_energy_vals->at(n_iter).at(s_iter) > 25){
	++numshowers;
      }
    }
      	  
    //Get the biggest shower 
    float MaxEnergy = -99999;
    int MaxEnergy_iter = -999;
    for(int s_iter=0; s_iter<shower_energy_vals->at(n_iter).size(); ++s_iter){
      if(shower_energy_vals->at(n_iter).at(s_iter) > MaxEnergy){
	MaxEnergy      = shower_energy_vals->at(n_iter).at(s_iter);
	MaxEnergy_iter = s_iter;
      }
    }
    //If there is no shower continue;
    if(MaxEnergy_iter == -999){continue;}
 
    if(MaxEnergy < fEnergyCut){continue;}

    float max_track_length = -99999;
    float max_pida         = -99999;
    int track_iter = 0;
    for(auto const& shower_track_length_val: shower_track_length_vals->at(n_iter)){
      if(max_track_length < shower_track_length_val){
	max_track_length = shower_track_length_val;
	max_pida = track_PIDA_vals->at(n_iter).at(track_iter);
      }
      ++track_iter;
    }

    for(auto& Metric: MetricMap){
      if(!Metric.second.FillForMVA()){continue;}
      
      float numneutrinos = shower_length_vals->size();
      float numshowers_mva = numshowers;
      
      if(Metric.first.c_str() == "dEdxMVA"){Metric.second.SetMVAVal(shower_dEdx_vals->at(n_iter).at(MaxEnergy_iter));}
      if(Metric.first.c_str() == "CoversionGapMVA"){Metric.second.SetMVAVal(shower_coversion_gap_vals->at(n_iter).at(MaxEnergy_iter));}
      if(Metric.first.c_str() == "MaxTrackLenghMVA"){Metric.second.SetMVAVal(max_track_length);}
      if(Metric.first.c_str() == "MaxTrackPIDAMVA"){Metric.second.SetMVAVal(max_pida);}
      if(Metric.first.c_str() == "ShowerEnergyMVA"){Metric.second.SetMVAVal(shower_energy_vals->at(n_iter).at(MaxEnergy_iter));}
      if(Metric.first.c_str() == "ShowerLengthEMVA"){Metric.second.SetMVAVal(shower_length_vals->at(n_iter).at(MaxEnergy_iter));}
      if(Metric.first.c_str() == "NumberOfNeutrinosMVA"){Metric.second.SetMVAVal(numneutrinos);}
      if(Metric.first.c_str() == "ShowerResidualNumShowersMVA"){Metric.second.SetMVAVal(numshowers_mva);}
      if(Metric.first.c_str() == "ShowerDensityPWMVA"){Metric.second.SetMVAVal(shower_density_pw_new_vals->at(n_iter).at(MaxEnergy_iter));}
      if(Metric.first.c_str() == "ShowerDensityGradNewMVA"){Metric.second.SetMVAVal(shower_density_grad_new_vals->at(n_iter).at(MaxEnergy_iter));}
      if(Metric.first.c_str() == "ShowerDensityOpeningAngleMVA"){Metric.second.SetMVAVal(shower_open_angle_vals->at(n_iter).at(MaxEnergy_iter));}
      if(Metric.first.c_str() == "NuetrinoPdGMVA"){Metric.second.SetMVAVal(nu_pdg_vals->at(n_iter));}
      if(Metric.first.c_str() == "ShowerTrackWidthMVA"){Metric.second.SetMVAVal(shower_trackwidth_vals->at(n_iter).at(MaxEnergy_iter));}
      if(Metric.first.c_str() == "ShowerTrackLengthMVA"){Metric.second.SetMVAVal(shower_length_vals->at(n_iter).at(MaxEnergy_iter));}
    }

    float MVAvalue = reader->EvaluateMVA(fMVAMethod);
    

    //    if(nu_reco_energy_vals->at(n_iter)<225 && sample=="signal"){std::cout << "nu_reco_energy_vals: " << nu_reco_energy_vals->at(n_iter) << "shower energy: " << MaxEnergy << " fMinEnergyCutVal: " << fMinEnergyCut  << " sample: " << sample << std::endl;}
    //Save the Neutrino E
    NeutrinoE[nu_truth_energy_vals->at(n_iter)].push_back(nu_reco_energy_vals->at(n_iter));
    NeutrinoMVA[nu_truth_energy_vals->at(n_iter)].push_back(MVAvalue);

  }  


  //Identify if this neutrino is the largest
  for(int nt_iter=0; nt_iter<nu_truth_energy_vals->size(); ++nt_iter){
    
    int added = 0;
    for(auto const& trueNeutrino: NeutrinoE){

      if(nu_truth_energy_vals->at(nt_iter) != trueNeutrino.first){continue;}

      float MaxNeutrinoEval      = -99999;
      int   MaxNeutrinoEval_iter = -99999; 

      //Identify the remianing reco neutrino with the max energy
      for(int n_iter=0; n_iter<(trueNeutrino.second.size()); ++n_iter){
	if(MaxNeutrinoEval < trueNeutrino.second.at(n_iter)){
	  MaxNeutrinoEval = trueNeutrino.second.at(n_iter);
	  MaxNeutrinoEval_iter = n_iter;
	}
      }
      vals.push_back(NeutrinoMVA.at(trueNeutrino.first).at(MaxNeutrinoEval_iter));
      ++added;
    }

    if(added == 0){
      vals.push_back(-1);
    }
  }
  return 0;
}



int optimiser::NueRecoOptimiser::ShowerAllsCut(optimiser::MetricHolder& Metric, std::string& sample){

  //Cuts 
  const float fResidualCutVal          = 0.025;//3.135;//0.025;
  const float fdEdxCutVal              = 2.5; //2.59
  const float fShowerTrackLengthCutVal = 100; //37.25;//54;//7.75; //40 
  const float fConversionGapVal        = 3;//3.5;//2.1;
  const float fSingleShowerEnergyVal   = 345;//345;
  const float fMinEnergyCutVal         = 200;//225;
  const float fScaleFactor             = 1000;//1.1949;
  const float fTrackPIDAScore          = 9.5625;//9.5625;//9.94;
  const std::string fMVAMethod         = "BDTGCrossValScaled";
  const float fMVAcut                  = -0.985;//0;
  const float fResidualEnergyCut       = 10;//10;23.75

  //Bools to swithc cut on l
  const bool fResidualCut          = false;
  const bool fSingleShowerCut      = false;
  const bool fdEdxCut              = true;
  const bool fShowerTrackLengthCut = true;
  const bool fConversionGap        = true;
  const bool fMinEnergyCut         = true; 
  const bool fOneShowerCut         = true;
  const bool fShowerTrackPIDACut   = false;
  const bool fNeutrinoPdgCut       = false;
  const bool fDoMVACut             = false;

  std::string NeutrinoTrueEBefore = Metric.GetMetricName() + "NeutrinoTrueEBefore";
  std::string NeutrinoRecoEBefore = Metric.GetMetricName() + "NeutrinoRecoEBefore";
  std::string NeutrinoRecoEBeforeExtra = Metric.GetMetricName() + "NeutrinoRecoEBeforeExtra";
  std::string NeutrinoRecoEAfter  = Metric.GetMetricName() + "NeutrinoRecoEAfter"; 
  std::string NeutrinoTrueEAfter = Metric.GetMetricName() + "NeutrinoTrueEAfter";
  std::string NeutrinoRecoEAfterExtra = Metric.GetMetricName() + "NeutrinoRecoEAfterExtra";
  std::string NeutrinoTrueEAfterExtra =  Metric.GetMetricName() + "NeutrinoTrueEAfterExtra";

  std::string NeutrinoRecoEBeforeNorm = Metric.GetMetricName() + "NeutrinoRecoEBeforeNorm";
  std::string NeutrinoRecoEAfterNorm = Metric.GetMetricName() + "NeutrinoRecoEAfterNorm";

  if(!Metric.CheckHistogram(NeutrinoTrueEBefore,sample)){Metric.SetHistogram(NeutrinoTrueEBefore,30,0,3000,sample);}
  if(!Metric.CheckHistogram(NeutrinoRecoEBefore,sample)){Metric.SetHistogram(NeutrinoRecoEBefore,Metric.Numbins,Metric.Xmin,Metric.Xmax,sample);}
  if(!Metric.CheckHistogram(NeutrinoRecoEBeforeExtra,sample)){Metric.SetHistogram(NeutrinoRecoEBeforeExtra,Metric.Numbins,Metric.Xmin,Metric.Xmax,sample);}
  if(!Metric.CheckHistogram(NeutrinoRecoEAfter,sample)){Metric.SetHistogram(NeutrinoRecoEAfter,Metric.Numbins,Metric.Xmin,Metric.Xmax,sample);}
  if(!Metric.CheckHistogram(NeutrinoTrueEAfter,sample)){Metric.SetHistogram(NeutrinoTrueEAfter,30,0,3000,sample);}
  if(!Metric.CheckHistogram(NeutrinoRecoEAfterExtra,sample)){Metric.SetHistogram(NeutrinoRecoEAfterExtra,Metric.Numbins,Metric.Xmin,Metric.Xmax,sample);}
  if(!Metric.CheckHistogram(NeutrinoTrueEAfterExtra,sample)){Metric.SetHistogram(NeutrinoTrueEAfterExtra,30,0,3000,sample);}
  if(!Metric.CheckHistogram(NeutrinoRecoEBeforeNorm,sample)){Metric.SetHistogram(NeutrinoRecoEBeforeNorm,Metric.Numbins,Metric.Xmin,Metric.Xmax,sample);}
  if(!Metric.CheckHistogram(NeutrinoRecoEAfterNorm,sample)){Metric.SetHistogram(NeutrinoRecoEAfterNorm,Metric.Numbins,Metric.Xmin,Metric.Xmax,sample);}

  //Get the shower energy 
  std::string shower_energy = "shower_energy";
  std::vector<std::vector<float> >* shower_energy_vals;
  int err = GetBranchPointer(shower_energy,sample,shower_energy_vals);
  if(err){return 1;}

  //Get the residual distance 
  std::string shower_residual_dist = "shower_residual_dist";
  std::vector<std::vector<float> >* shower_residual_vals;
  err = GetBranchPointer(shower_residual_dist,sample,shower_residual_vals);
  if(err){return 1;}

  //Get the dEdx
  std::string shower_dEdx = "shower_dEdx";
  std::vector<std::vector<float> >* shower_dEdx_vals;
  err = GetBranchPointer(shower_dEdx,sample,shower_dEdx_vals);
  if(err){return 1;}

  //Get the conversion distance
  std::string shower_coversion_gap = "shower_coversion_gap";
  std::vector<std::vector<float> >* shower_coversion_gap_vals;
  err = GetBranchPointer(shower_coversion_gap,sample,shower_coversion_gap_vals);
  if(err){return 1;}

  //Get the max track length 
  std::string track_lengths = "track_lengths";
  std::vector<std::vector<float> >* shower_track_length_vals;
  err = GetBranchPointer(track_lengths,sample,shower_track_length_vals);
  if(err){return 1;}

  //Get the nu reco energy  
  std::string nu_reco_energy = "nu_reco_energy";
  std::vector<float>* nu_reco_energy_vals;
  err = GetBranchPointer(nu_reco_energy,sample,nu_reco_energy_vals);
  if(err){return 1;}

  //Get the nu truth energy  
  std::string nu_truth_energy = "nu_truth_energy";
  std::vector<float> * nu_truth_energy_vals;
  err = GetBranchPointer(nu_truth_energy,sample,nu_truth_energy_vals);
  if(err){return 1;}

  //Get the nu truth energy  
  std::string track_PIDAs = "track_PIDA_branch";
  std::vector<std::vector<float> >* track_PIDA_vals;
  err = GetBranchPointer(track_PIDAs,sample,track_PIDA_vals);
  if(err){return 1;}

  std::string shower_density_pw_news = "shower_density_pw_new";
  std::vector<std::vector<float> >* shower_density_pw_new_vals;
  err = GetBranchPointer(shower_density_pw_news,sample,shower_density_pw_new_vals);
  if(err){return 1;}

  std::string shower_density_grad_news = "shower_density_grad_new";
  std::vector<std::vector<float> >* shower_density_grad_new_vals;
  err = GetBranchPointer(shower_density_grad_news,sample,shower_density_grad_new_vals);
  if(err){return 1;}

  std::string shower_trackwidths = "shower_trackwidth";
  std::vector<std::vector<float> >* shower_trackwidth_vals;
  err = GetBranchPointer(shower_trackwidths,sample,shower_trackwidth_vals);
  if(err){return 1;}

  std::string shower_tracklengths = "shower_tracklength";
  std::vector<std::vector<float> >* shower_tracklength_vals;
  err = GetBranchPointer(shower_tracklengths,sample,shower_tracklength_vals);
  if(err){return 1;}

  std::string shower_lengths = "shower_length";
  std::vector<std::vector<float> >* shower_length_vals;
  err = GetBranchPointer(shower_lengths,sample,shower_length_vals);
  if(err){return 1;}

  std::string shower_open_angles = "shower_open_angle";
  std::vector<std::vector<float> >* shower_open_angle_vals;
  err = GetBranchPointer(shower_open_angles,sample,shower_open_angle_vals);
  if(err){return 1;}



  //Get the nu truth energy  
  std::string nu_E = "nu_E";
  std::vector<float> * nu_E_vals;
  err = GetBranchPointer(nu_E,sample,nu_E_vals);
  if(err){return 1;}

  //Get the nu truth energy  
  std::string nu_pdg = "nu_pdg";
  std::vector<float> * nu_pdg_vals;
  err = GetBranchPointer(nu_pdg,sample,nu_pdg_vals);
  if(err){return 1;}


  //Get the metric branch
  std::string branchname = Metric.GetBranchName();
  std::vector<float> * branch_vals;
  err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  //Get the partner metric branch
  std::string partnername_len = "shower_len";
  std::vector<std::vector<float> >* partner_vals_len;
  err = GetBranchPointer(partnername_len,sample,partner_vals_len);
  if(err){return 1;}

  //Get the partner metric branch
  std::string partnername_openingangle = "shower_open_angle";
  std::vector<std::vector<float> >* partner_vals_openingangle;
  err = GetBranchPointer(partnername_openingangle,sample,partner_vals_openingangle);
  if(err){return 1;}



  std::string branchnamecc = "nu_cc";
  std::vector<float>* branch_valscc;
  err = GetBranchPointer(branchnamecc,sample,branch_valscc);
  



  // //Loop over the neutrinos and get their true energy
  // for(int n_iter=0; n_iter<nu_E_vals->size(); ++n_iter){

  //   double oscprob = GetOscProb(n_iter,sample,Metric.ApplyOscProb());
  //   double potweight = 1;  
  //   if(Metric.ApplyPOT()){potweight = Metric.ReturnPOT(sample);}
  //   double oscpot = oscprob*potweight;
   
    
  //   Metric.GetHistogram(NeutrinoTrueEBefore,sample)->Fill(nu_truth_energy_vals->at(n_iter)*1000,oscpot);
  // }

  //Loop over the reco neutrinos and idnetify the reco neutrino with the most energy
  std::map<float,std::vector<float> > NeutrinoE;
  std::map<float,float > MaxNeutrinoE;
  std::map<float, int >  MaxNeutrinoE_iter;
  for(int n_iter=0; n_iter<nu_E_vals->size(); ++n_iter){
    if(branch_valscc->at(n_iter) == 1 && sample == "signal"){continue;}

    MaxNeutrinoE[nu_E_vals->at(n_iter)]      = -99999;
    MaxNeutrinoE_iter[nu_E_vals->at(n_iter)] = -99999;   
  }

  for(int n_iter=0; n_iter<nu_reco_energy_vals->size(); ++n_iter){
    if(branch_valscc->at(n_iter) == 1 && sample == "signal"){continue;}

    //check if the true neutrino is already be identified 
    if(MaxNeutrinoE[nu_truth_energy_vals->at(n_iter)] < nu_reco_energy_vals->at(n_iter)){
      MaxNeutrinoE[nu_truth_energy_vals->at(n_iter)] = nu_reco_energy_vals->at(n_iter);
      MaxNeutrinoE_iter[nu_truth_energy_vals->at(n_iter)] = n_iter;
    }
  }
  //Add the max reco values to the efficiency and the others to the additional graph
  for(int n_iter=0; n_iter<nu_reco_energy_vals->size(); ++n_iter){
    if(branch_valscc->at(n_iter) == 1 && sample == "signal"){continue;}

    double oscprob = GetOscProb(n_iter,sample,Metric.ApplyOscProb());
    double potweight = 1;  
    if(Metric.ApplyPOT()){potweight = Metric.ReturnPOT(sample);}
    double oscpot = oscprob*potweight;
    //    if(sample=="signal")std::cout << "presel: " << oscpot  << " oscprob: " << oscprob << " potweight: " << potweight << std::endl;

    if(MaxNeutrinoE[nu_truth_energy_vals->at(n_iter)] == nu_reco_energy_vals->at(n_iter)){
      Metric.GetHistogram(NeutrinoTrueEBefore,sample)->Fill(nu_truth_energy_vals->at(n_iter)*1000,oscpot);
      Metric.GetHistogram(NeutrinoRecoEBefore,sample)->Fill(branch_vals->at(n_iter)*fScaleFactor,oscpot);
      Metric.GetHistogram(NeutrinoRecoEBeforeNorm,sample)->Fill(branch_vals->at(n_iter)*fScaleFactor);

    }
    else{
      Metric.GetHistogram(NeutrinoRecoEBeforeExtra,sample)->Fill(branch_vals->at(n_iter)*fScaleFactor,oscpot); 
    }
  }

  //Re initialise the maps
  for(int n_iter=0; n_iter<nu_E_vals->size(); ++n_iter){
    if(branch_valscc->at(n_iter) == 1 && sample == "signal"){continue;}

    MaxNeutrinoE[nu_E_vals->at(n_iter)]      = -99999;
    MaxNeutrinoE_iter[nu_E_vals->at(n_iter)] = -99999;   
  }


  //Loop over the neutrinos
  for(int n_iter=0; n_iter<nu_reco_energy_vals->size(); ++n_iter){
    
    if(branch_valscc->at(n_iter) == 1 && sample == "signal"){continue;}

    if(fApplyFVCut){
      std::vector<float> * fv_vals; 
      std::string fv_branch_name = "in_FV";
      int err = GetBranchPointer(fv_branch_name,sample,fv_vals);
      if(err){return 1;}
      if(fv_vals->size() < nu_reco_energy_vals->size()){return 1;}
      if(fv_vals->at(n_iter) == 0){
	continue;
      }
    }

    if(nu_pdg_vals->at(n_iter) != 12 && fNeutrinoPdgCut){
      continue;
    }
    
    if(shower_energy_vals->at(n_iter).size() == 0){
      continue;
    }

    int  numshowers = 0;

    
    //Residual Cut
    if(fResidualCut){
      numshowers = 0;
      for(int s_iter=0; s_iter<shower_residual_vals->at(n_iter).size(); ++s_iter){
	float showerdist = partner_vals_len->at(n_iter).at(s_iter) * TMath::Tan(0.5*partner_vals_openingangle->at(n_iter).at(0));
	if(shower_energy_vals->at(n_iter).at(s_iter) > fResidualEnergyCut && shower_residual_vals->at(n_iter).at(s_iter) > TMath::Abs(fResidualCutVal*showerdist)+2){
	  //	if(shower_residual_vals->at(n_iter).at(s_iter) > fResidualCut && shower_energy_vals->at(n_iter).at(s_iter) > 25){
	  ++numshowers;
	}
      }
      //Only Keep one shower events
      if(numshowers > 1 && fResidualCut){
	continue;
      }
    }
    else if(fSingleShowerCut){
      numshowers = 0;
      for(int s_iter=0; s_iter<shower_residual_vals->at(n_iter).size(); ++s_iter){
	if(shower_energy_vals->at(n_iter).at(s_iter) > fSingleShowerEnergyVal){
	  ++numshowers;
	}
      }
      //Only Keep one shower events
      if(numshowers != 1 && fSingleShowerCut){
	continue;
      }
    }
    else if (fOneShowerCut){
      int showers = 0;
      for(int s_iter=0; s_iter<shower_energy_vals->at(n_iter).size(); ++s_iter){
	if(shower_energy_vals->at(n_iter).at(s_iter) > 10){
	  ++showers;
	}
      }
      if(showers != 1) {continue;}
    }
  	  
    //Get the biggest shower 
    float MaxEnergy = -99999;
    int MaxEnergy_iter = -999;
    for(int s_iter=0; s_iter<shower_energy_vals->at(n_iter).size(); ++s_iter){
      if(shower_energy_vals->at(n_iter).at(s_iter) > MaxEnergy){
	MaxEnergy      = shower_energy_vals->at(n_iter).at(s_iter);
	MaxEnergy_iter = s_iter;
      }
    }
    //If there is no shower continue;
    if(MaxEnergy_iter == -999){continue;}
 

    //Energy Cut
    if(shower_energy_vals->at(n_iter).at(MaxEnergy_iter) < fMinEnergyCutVal && fMinEnergyCut){
      continue;
    } 

    //dEdx Cut
    if(shower_dEdx_vals->at(n_iter).at(MaxEnergy_iter) > fdEdxCutVal && fdEdxCut){
      continue;
    }

    //Conversion Gap Cut
    if(shower_coversion_gap_vals->at(n_iter).at(MaxEnergy_iter) > fConversionGapVal && fConversionGap){
      continue;
    }

    //Track length cut
    float max_track_length = -99999;
    float max_pida         = -99999;
    int track_iter = 0;
    for(auto const& shower_track_length_val: shower_track_length_vals->at(n_iter)){
      if(max_track_length < shower_track_length_val){
	max_track_length = shower_track_length_val;
	max_pida = track_PIDA_vals->at(n_iter).at(track_iter);
      }
      ++track_iter;
    }

    if(max_track_length > fShowerTrackLengthCutVal && fShowerTrackLengthCut){continue;}
    if(max_track_length > fShowerTrackLengthCutVal && max_pida < fTrackPIDAScore && fShowerTrackPIDACut){continue;}

    if(fDoMVACut){
      for(auto& Metric: MetricMap){
	if(!Metric.second.FillForMVA()){continue;}
	
	float numneutrinos = shower_length_vals->size();
	float numshowers_mva = numshowers;

	if(Metric.first.c_str() == "dEdxMVA"){Metric.second.SetMVAVal(shower_dEdx_vals->at(n_iter).at(MaxEnergy_iter));}
	if(Metric.first.c_str() == "CoversionGapMVA"){Metric.second.SetMVAVal(shower_coversion_gap_vals->at(n_iter).at(MaxEnergy_iter));}
	if(Metric.first.c_str() == "MaxTrackLenghMVA"){Metric.second.SetMVAVal(max_track_length);}
	if(Metric.first.c_str() == "MaxTrackPIDAMVA"){Metric.second.SetMVAVal(max_pida);}
	if(Metric.first.c_str() == "ShowerEnergyMVA"){Metric.second.SetMVAVal(shower_energy_vals->at(n_iter).at(MaxEnergy_iter));}
	if(Metric.first.c_str() == "ShowerLengthEMVA"){Metric.second.SetMVAVal(shower_length_vals->at(n_iter).at(MaxEnergy_iter));}
	if(Metric.first.c_str() == "NumberOfNeutrinosMVA"){Metric.second.SetMVAVal(numneutrinos);}
	if(Metric.first.c_str() == "ShowerResidualNumShowersMVA"){Metric.second.SetMVAVal(numshowers_mva);}
	if(Metric.first.c_str() == "ShowerDensityPWMVA"){Metric.second.SetMVAVal(shower_density_pw_new_vals->at(n_iter).at(MaxEnergy_iter));}
	if(Metric.first.c_str() == "ShowerDensityGradNewMVA"){Metric.second.SetMVAVal(shower_density_grad_new_vals->at(n_iter).at(MaxEnergy_iter));}
	if(Metric.first.c_str() == "ShowerDensityOpeningAngleMVA"){Metric.second.SetMVAVal(shower_open_angle_vals->at(n_iter).at(MaxEnergy_iter));}
	if(Metric.first.c_str() == "NuetrinoPdGMVA"){Metric.second.SetMVAVal(nu_pdg_vals->at(n_iter));}
	if(Metric.first.c_str() == "ShowerTrackWidthMVA"){Metric.second.SetMVAVal(shower_trackwidth_vals->at(n_iter).at(MaxEnergy_iter));}
	if(Metric.first.c_str() == "ShowerTrackLengthMVA"){Metric.second.SetMVAVal(shower_length_vals->at(n_iter).at(MaxEnergy_iter));}
      }
      if(reader->EvaluateMVA(fMVAMethod) < fMVAcut){continue;}
    }

    //    if(nu_reco_energy_vals->at(n_iter)<225 && sample=="signal"){std::cout << "nu_reco_energy_vals: " << nu_reco_energy_vals->at(n_iter) << "shower energy: " << MaxEnergy << " fMinEnergyCutVal: " << fMinEnergyCut  << " sample: " << sample << std::endl;}

    //Save the Neutrino E
    NeutrinoE[nu_truth_energy_vals->at(n_iter)].push_back(nu_reco_energy_vals->at(n_iter));
  }  


  //Identify if this neutrino is the largest
  int neut = 0;
  for(auto const& trueNeutrino: NeutrinoE){
    float MaxNeutrinoEval      = -99999;
    int   MaxNeutrinoEval_iter = -99999; 

    //Identify the remianing reco neutrino with the max energy
    for(int n_iter=0; n_iter<(trueNeutrino.second.size()); ++n_iter){
      if(MaxNeutrinoEval < trueNeutrino.second.at(n_iter)){
	MaxNeutrinoEval = trueNeutrino.second.at(n_iter);
	MaxNeutrinoEval_iter = n_iter;

	}

      }

    double oscpot;
    for(int n_iter=0; n_iter<nu_truth_energy_vals->size(); ++n_iter){
      if(nu_truth_energy_vals->at(n_iter) == trueNeutrino.first){
	double oscprob = GetOscProb(n_iter,sample,Metric.ApplyOscProb());
	double potweight = 1;  
	if(Metric.ApplyPOT()){potweight = Metric.ReturnPOT(sample);}
	oscpot = oscprob*potweight;
	//	if(sample == "signal"){std::cout << "postsel: " << oscpot  << " oscprob: " << oscprob << " potweight: " << potweight << std::endl;}
	break;
      }
    }
    
    for(int n_iter=0; n_iter<(trueNeutrino.second.size()); ++n_iter){ 

      int   filling_iter         = -99999;
      for(int reco_iter=0; reco_iter<nu_reco_energy_vals->size(); ++reco_iter){
	if(nu_reco_energy_vals->at(reco_iter) == trueNeutrino.second.at(n_iter)){
	  filling_iter = reco_iter;
	}
      }


      //Add the biggest to account for the efficiency.
      if(n_iter == MaxNeutrinoEval_iter){
	//	if(sample=="signal")std::cout << "oscpot: " << oscpot << std::endl;
	Metric.GetHistogram(NeutrinoRecoEAfter,sample)->Fill(branch_vals->at(filling_iter)*fScaleFactor,oscpot);
	Metric.GetHistogram(NeutrinoRecoEAfterNorm,sample)->Fill(branch_vals->at(filling_iter)*fScaleFactor);
	Metric.GetHistogram(NeutrinoTrueEAfter,sample)->Fill(trueNeutrino.first*1000,oscpot);
      }
      //Add the extra reco
      else{
	Metric.GetHistogram(NeutrinoRecoEAfterExtra,sample)->Fill(branch_vals->at(n_iter)*fScaleFactor,oscpot);
	Metric.GetHistogram(NeutrinoTrueEAfterExtra,sample)->Fill(trueNeutrino.first*1000,oscpot);
      }
    }
  }
return 0;
}


// int optimiser::NueRecoOptimiser::ShowerAllsCut(optimiser::MetricHolder& Metric, std::string& sample){

//   //Cuts 
//   const float fResidualCutVal          = 0.025;//0.025;
//   const float fdEdxCutVal              = 2.5; //2.59
//   const float fShowerTrackLengthCutVal = 7.75;//54;//7.75; //40 
//   const float fConversionGapVal        = 2.1;//3.5;
//   const float fSingleShowerEnergyVal   = 225;//345;
//   const float fMinEnergyCutVal         = 225;//345;
//   const float fScaleFactor             = 1.1949;//1.1949;
//   const float fTrackPIDAScore          = 9.5625;//9.94;
//   const std::string fMVAMethod         = "BDTGCrossVal";
//   const float fMVAcut                  = 0;

//   //Bools to swithc cut on l
//   const bool fResidualCut          = false;
//   const bool fSingleShowerCut      = false;
//   const bool fdEdxCut              = false;
//   const bool fShowerTrackLengthCut = false;
//   const bool fConversionGap        = false;
//   const bool fMinEnergyCut         = false; 
//   const bool fOneShowerCut         = false;
//   const bool fShowerTrackPIDACut   = false;
//   const bool fNeutrinoPdgCut       = false;

//   std::string NeutrinoTrueEBefore = Metric.GetMetricName() + "NeutrinoTrueEBefore";
//   std::string NeutrinoRecoEBefore = Metric.GetMetricName() + "NeutrinoRecoEBefore";
//   std::string NeutrinoRecoEBeforeExtra = Metric.GetMetricName() + "NeutrinoRecoEBeforeExtra";
//   std::string NeutrinoRecoEAfter  = Metric.GetMetricName() + "NeutrinoRecoEAfter"; 
//   std::string NeutrinoTrueEAfter = Metric.GetMetricName() + "NeutrinoTrueEAfter";
//   std::string NeutrinoRecoEAfterExtra = Metric.GetMetricName() + "NeutrinoRecoEAfterExtra";
//   std::string NeutrinoTrueEAfterExtra =  Metric.GetMetricName() + "NeutrinoTrueEAfterExtra";

//   std::string NeutrinoRecoEBeforeNorm = Metric.GetMetricName() + "NeutrinoRecoEBeforeNorm";
//   std::string NeutrinoRecoEAfterNorm = Metric.GetMetricName() + "NeutrinoRecoEAfterNorm";

//   if(!Metric.CheckHistogram(NeutrinoTrueEBefore,sample)){Metric.SetHistogram(NeutrinoTrueEBefore,30,0,3000,sample);}
//   if(!Metric.CheckHistogram(NeutrinoRecoEBefore,sample)){Metric.SetHistogram(NeutrinoRecoEBefore,Metric.Numbins,Metric.Xmin,Metric.Xmax,sample);}
//   if(!Metric.CheckHistogram(NeutrinoRecoEBeforeExtra,sample)){Metric.SetHistogram(NeutrinoRecoEBeforeExtra,Metric.Numbins,Metric.Xmin,Metric.Xmax,sample);}
//   if(!Metric.CheckHistogram(NeutrinoRecoEAfter,sample)){Metric.SetHistogram(NeutrinoRecoEAfter,Metric.Numbins,Metric.Xmin,Metric.Xmax,sample);}
//   if(!Metric.CheckHistogram(NeutrinoTrueEAfter,sample)){Metric.SetHistogram(NeutrinoTrueEAfter,30,0,3000,sample);}
//   if(!Metric.CheckHistogram(NeutrinoRecoEAfterExtra,sample)){Metric.SetHistogram(NeutrinoRecoEAfterExtra,Metric.Numbins,Metric.Xmin,Metric.Xmax,sample);}
//   if(!Metric.CheckHistogram(NeutrinoTrueEAfterExtra,sample)){Metric.SetHistogram(NeutrinoTrueEAfterExtra,30,0,3000,sample);}
//   if(!Metric.CheckHistogram(NeutrinoRecoEBeforeNorm,sample)){Metric.SetHistogram(NeutrinoRecoEBeforeNorm,Metric.Numbins,Metric.Xmin,Metric.Xmax,sample);}
//   if(!Metric.CheckHistogram(NeutrinoRecoEAfterNorm,sample)){Metric.SetHistogram(NeutrinoRecoEAfterNorm,Metric.Numbins,Metric.Xmin,Metric.Xmax,sample);}

//   //Get the shower energy 
//   std::string shower_energy = "shower_energy";
//   std::vector<std::vector<float> >* shower_energy_vals;
//   int err = GetBranchPointer(shower_energy,sample,shower_energy_vals);
//   if(err){return 1;}

//   //Get the residual distance 
//   std::string shower_residual_dist = "shower_residual_dist";
//   std::vector<std::vector<float> >* shower_residual_vals;
//   err = GetBranchPointer(shower_residual_dist,sample,shower_residual_vals);
//   if(err){return 1;}

//   //Get the dEdx
//   std::string shower_dEdx = "shower_dEdx";
//   std::vector<std::vector<float> >* shower_dEdx_vals;
//   err = GetBranchPointer(shower_dEdx,sample,shower_dEdx_vals);
//   if(err){return 1;}

//   //Get the conversion distance
//   std::string shower_coversion_gap = "shower_coversion_gap";
//   std::vector<std::vector<float> >* shower_coversion_gap_vals;
//   err = GetBranchPointer(shower_coversion_gap,sample,shower_coversion_gap_vals);
//   if(err){return 1;}

//   //Get the max track length 
//   std::string track_lengths = "track_lengths";
//   std::vector<std::vector<float> >* shower_track_length_vals;
//   err = GetBranchPointer(track_lengths,sample,shower_track_length_vals);
//   if(err){return 1;}

//   //Get the nu reco energy  
//   std::string nu_reco_energy = "nu_reco_energy";
//   std::vector<float>* nu_reco_energy_vals;
//   err = GetBranchPointer(nu_reco_energy,sample,nu_reco_energy_vals);
//   if(err){return 1;}

//   //Get the nu truth energy  
//   std::string nu_truth_energy = "nu_truth_energy";
//   std::vector<float> * nu_truth_energy_vals;
//   err = GetBranchPointer(nu_truth_energy,sample,nu_truth_energy_vals);
//   if(err){return 1;}

//   //Get the nu truth energy  
//   std::string track_PIDAs = "track_PIDA_branch";
//   std::vector<std::vector<float> >* track_PIDA_vals;
//   err = GetBranchPointer(track_PIDAs,sample,track_PIDA_vals);
//   if(err){return 1;}


//   //Get the nu truth energy  
//   std::string nu_E = "nu_E";
//   std::vector<float> * nu_E_vals;
//   err = GetBranchPointer(nu_E,sample,nu_E_vals);
//   if(err){return 1;}

//   //Get the nu truth energy  
//   std::string nu_pdg = "nu_pdg";
//   std::vector<float> * nu_pdg_vals;
//   err = GetBranchPointer(nu_pdg,sample,nu_pdg_vals);
//   if(err){return 1;}


//   //Get the metric branch
//   std::string branchname = Metric.GetBranchName();
//   std::vector<float> * branch_vals;
//   err = GetBranchPointer(branchname,sample,branch_vals);
//   if(err){return 1;}

//   //Get the partner metric branch
//   std::string partnername_len = "shower_len";
//   std::vector<std::vector<float> >* partner_vals_len;
//   err = GetBranchPointer(partnername_len,sample,partner_vals_len);
//   if(err){return 1;}

//   //Get the partner metric branch
//   std::string partnername_openingangle = "shower_open_angle";
//   std::vector<std::vector<float> >* partner_vals_openingangle;
//   err = GetBranchPointer(partnername_openingangle,sample,partner_vals_openingangle);
//   if(err){return 1;}



//   std::string branchnamecc = "nu_cc";
//   std::vector<float>* branch_valscc;
//   err = GetBranchPointer(branchnamecc,sample,branch_valscc);
  



//   // //Loop over the neutrinos and get their true energy
//   // for(int n_iter=0; n_iter<nu_E_vals->size(); ++n_iter){

//   //   double oscprob = GetOscProb(n_iter,sample,Metric.ApplyOscProb());
//   //   double potweight = 1;  
//   //   if(Metric.ApplyPOT()){potweight = Metric.ReturnPOT(sample);}
//   //   double oscpot = oscprob*potweight;
   
    
//   //   Metric.GetHistogram(NeutrinoTrueEBefore,sample)->Fill(nu_truth_energy_vals->at(n_iter)*1000,oscpot);
//   // }

//   //Loop over the reco neutrinos and idnetify the reco neutrino with the most energy
//   std::map<float,std::vector<float> > NeutrinoE;
//   std::map<float,float > MaxNeutrinoE;
//   std::map<float, int >  MaxNeutrinoE_iter;
//   for(int n_iter=0; n_iter<nu_E_vals->size(); ++n_iter){
//     if(branch_valscc->at(n_iter) == 1 && sample == "signal"){continue;}

//     MaxNeutrinoE[nu_E_vals->at(n_iter)]      = -99999;
//     MaxNeutrinoE_iter[nu_E_vals->at(n_iter)] = -99999;   
//   }

//   for(int n_iter=0; n_iter<nu_reco_energy_vals->size(); ++n_iter){
//     if(branch_valscc->at(n_iter) == 1 && sample == "signal"){continue;}

//     //check if the true neutrino is already be identified 
//     if(MaxNeutrinoE[nu_truth_energy_vals->at(n_iter)] < nu_reco_energy_vals->at(n_iter)){
//       MaxNeutrinoE[nu_truth_energy_vals->at(n_iter)] = nu_reco_energy_vals->at(n_iter);
//       MaxNeutrinoE_iter[nu_truth_energy_vals->at(n_iter)] = n_iter;
//     }
//   }
//   //Add the max reco values to the efficiency and the others to the additional graph
//   for(int n_iter=0; n_iter<nu_reco_energy_vals->size(); ++n_iter){
//     if(branch_valscc->at(n_iter) == 1 && sample == "signal"){continue;}

//     double oscprob = GetOscProb(n_iter,sample,Metric.ApplyOscProb());
//     double potweight = 1;  
//     if(Metric.ApplyPOT()){potweight = Metric.ReturnPOT(sample);}
//     double oscpot = oscprob*potweight;
//     //    if(sample=="signal")std::cout << "presel: " << oscpot  << " oscprob: " << oscprob << " potweight: " << potweight << std::endl;

//     if(MaxNeutrinoE[nu_truth_energy_vals->at(n_iter)] == nu_reco_energy_vals->at(n_iter)){
//       Metric.GetHistogram(NeutrinoTrueEBefore,sample)->Fill(nu_truth_energy_vals->at(n_iter)*1000,oscpot);
//       Metric.GetHistogram(NeutrinoRecoEBefore,sample)->Fill(branch_vals->at(n_iter)*fScaleFactor,oscpot);
//       Metric.GetHistogram(NeutrinoRecoEBeforeNorm,sample)->Fill(branch_vals->at(n_iter)*fScaleFactor);

//     }
//     else{
//       Metric.GetHistogram(NeutrinoRecoEBeforeExtra,sample)->Fill(branch_vals->at(n_iter)*fScaleFactor,oscpot); 
//     }
//   }

//   //Re initialise the maps
//   for(int n_iter=0; n_iter<nu_E_vals->size(); ++n_iter){
//     if(branch_valscc->at(n_iter) == 1 && sample == "signal"){continue;}

//     MaxNeutrinoE[nu_E_vals->at(n_iter)]      = -99999;
//     MaxNeutrinoE_iter[nu_E_vals->at(n_iter)] = -99999;   
//   }


//   //Loop over the neutrinos
//   for(int n_iter=0; n_iter<nu_reco_energy_vals->size(); ++n_iter){
    
//     if(branch_valscc->at(n_iter) == 1 && sample == "signal"){continue;}

//     if(fApplyFVCut){
//       std::vector<float> * fv_vals; 
//       std::string fv_branch_name = "in_FV";
//       int err = GetBranchPointer(fv_branch_name,sample,fv_vals);
//       if(err){return 1;}
//       if(fv_vals->size() < nu_reco_energy_vals->size()){return 1;}
//       if(fv_vals->at(n_iter) == 0){
// 	continue;
//       }
//     }

//     if(nu_pdg_vals->at(n_iter) != 12 && fNeutrinoPdgCut){
//       continue;
//     }
    
//     if(shower_energy_vals->at(n_iter).size() == 0){
//       continue;
//     }
    
//     //Residual Cut
//     if(fResidualCut){
//       int  numshowers = 0;
//       for(int s_iter=0; s_iter<shower_residual_vals->at(n_iter).size(); ++s_iter){
// 	float showerdist = partner_vals_len->at(n_iter).at(s_iter) * TMath::Tan(0.5*partner_vals_openingangle->at(n_iter).at(0));
// 	if(shower_energy_vals->at(n_iter).at(s_iter) > 10 && shower_residual_vals->at(n_iter).at(s_iter) > TMath::Abs(fResidualCutVal*showerdist)+2){
// 	  //	if(shower_residual_vals->at(n_iter).at(s_iter) > fResidualCut && shower_energy_vals->at(n_iter).at(s_iter) > 25){
// 	  ++numshowers;
// 	}
//       }
//       //Only Keep one shower events
//       if(numshowers > 1 && fResidualCut){
// 	continue;
//       }
//     }
//     else if(fSingleShowerCut){
//       int  numshowers = 0;
//       for(int s_iter=0; s_iter<shower_residual_vals->at(n_iter).size(); ++s_iter){
// 	if(shower_energy_vals->at(n_iter).at(s_iter) > fSingleShowerEnergyVal){
// 	  ++numshowers;
// 	}
//       }
//       //Only Keep one shower events
//       if(numshowers != 1 && fSingleShowerCut){
// 	continue;
//       }
//     }
//     else if (fOneShowerCut){
//       int showers = 0;
//       for(int s_iter=0; s_iter<shower_energy_vals->at(n_iter).size(); ++s_iter){
// 	if(shower_energy_vals->at(n_iter).at(s_iter) > 10){
// 	  ++showers;
// 	}
//       }
//       if(showers != 1) {continue;}
//     }
  	  
//     //Get the biggest shower 
//     float MaxEnergy = -99999;
//     int MaxEnergy_iter = -999;
//     for(int s_iter=0; s_iter<shower_energy_vals->at(n_iter).size(); ++s_iter){
//       if(shower_energy_vals->at(n_iter).at(s_iter) > MaxEnergy){
// 	MaxEnergy      = shower_energy_vals->at(n_iter).at(s_iter);
// 	MaxEnergy_iter = s_iter;
//       }
//     }
//     //If there is no shower continue;
//     if(MaxEnergy_iter == -999){continue;}
 

//     //Energy Cut
//     if(shower_energy_vals->at(n_iter).at(MaxEnergy_iter) < fMinEnergyCutVal && fMinEnergyCut){
//       continue;
//     } 

//     //dEdx Cut
//     if(shower_dEdx_vals->at(n_iter).at(MaxEnergy_iter) > fdEdxCutVal && fdEdxCut){
//       continue;
//     }

//     //Conversion Gap Cut
//     if(shower_coversion_gap_vals->at(n_iter).at(MaxEnergy_iter) > fConversionGapVal && fConversionGap){
//       continue;
//     }

//     //Track length cut
//     float max_track_length = -99999;
//     float max_pida         = -99999;
//     int track_iter = 0;
//     for(auto const& shower_track_length_val: shower_track_length_vals->at(n_iter)){
//       if(max_track_length < shower_track_length_val){
// 	max_track_length = shower_track_length_val;
// 	max_pida = track_PIDA_vals->at(n_iter).at(track_iter);
//       }
//       ++track_iter;
//     }

//     if(reader->EvaluateMVA(fMVAMethod) < fMVAcut){continue;}

//     if(max_track_length > fShowerTrackLengthCutVal && fShowerTrackLengthCut){continue;}
//     if(max_track_length > fShowerTrackLengthCutVal && max_pida < fTrackPIDAScore && fShowerTrackPIDACut){continue;}

//     //    if(nu_reco_energy_vals->at(n_iter)<225 && sample=="signal"){std::cout << "nu_reco_energy_vals: " << nu_reco_energy_vals->at(n_iter) << "shower energy: " << MaxEnergy << " fMinEnergyCutVal: " << fMinEnergyCut  << " sample: " << sample << std::endl;}

//     //Save the Neutrino E
//     NeutrinoE[nu_truth_energy_vals->at(n_iter)].push_back(nu_reco_energy_vals->at(n_iter));
//   }  


//   //Identify if this neutrino is the largest
//   int neut = 0;
//   for(auto const& trueNeutrino: NeutrinoE){
//     float MaxNeutrinoEval      = -99999;
//     int   MaxNeutrinoEval_iter = -99999; 

//     //Identify the remianing reco neutrino with the max energy
//     for(int n_iter=0; n_iter<(trueNeutrino.second.size()); ++n_iter){
//       if(MaxNeutrinoEval < trueNeutrino.second.at(n_iter)){
// 	MaxNeutrinoEval = trueNeutrino.second.at(n_iter);
// 	MaxNeutrinoEval_iter = n_iter;

// 	}

//       }

//     double oscpot;
//     for(int n_iter=0; n_iter<nu_truth_energy_vals->size(); ++n_iter){
//       if(nu_truth_energy_vals->at(n_iter) == trueNeutrino.first){
// 	double oscprob = GetOscProb(n_iter,sample,Metric.ApplyOscProb());
// 	double potweight = 1;  
// 	if(Metric.ApplyPOT()){potweight = Metric.ReturnPOT(sample);}
// 	oscpot = oscprob*potweight;
// 	//	if(sample == "signal"){std::cout << "postsel: " << oscpot  << " oscprob: " << oscprob << " potweight: " << potweight << std::endl;}
// 	break;
//       }
//     }
    
//     for(int n_iter=0; n_iter<(trueNeutrino.second.size()); ++n_iter){ 

//       int   filling_iter         = -99999;
//       for(int reco_iter=0; reco_iter<nu_reco_energy_vals->size(); ++reco_iter){
// 	if(nu_reco_energy_vals->at(reco_iter) == trueNeutrino.second.at(n_iter)){
// 	  filling_iter = reco_iter;
// 	}
//       }


//       //Add the biggest to account for the efficiency.
//       if(n_iter == MaxNeutrinoEval_iter){
// 	//	if(sample=="signal")std::cout << "oscpot: " << oscpot << std::endl;
// 	Metric.GetHistogram(NeutrinoRecoEAfter,sample)->Fill(branch_vals->at(filling_iter)*fScaleFactor,oscpot);
// 	Metric.GetHistogram(NeutrinoRecoEAfterNorm,sample)->Fill(branch_vals->at(filling_iter)*fScaleFactor);
// 	Metric.GetHistogram(NeutrinoTrueEAfter,sample)->Fill(trueNeutrino.first*1000,oscpot);
//       }
//       //Add the extra reco
//       else{
// 	Metric.GetHistogram(NeutrinoRecoEAfterExtra,sample)->Fill(branch_vals->at(n_iter)*fScaleFactor,oscpot);
// 	Metric.GetHistogram(NeutrinoTrueEAfterExtra,sample)->Fill(trueNeutrino.first*1000,oscpot);
//       }
//     }
//   }
// return 0;
// }

int optimiser::NueRecoOptimiser::TrainMVA(){

  MVAFile->cd();
  SignalTreeMVA->Write();
  BackgroundTreeMVA->Write();

  std::cout << "sig entries: " << SignalTreeMVA->GetEntries() << " bk entries: " << BackgroundTreeMVA->GetEntries() << std::endl;

  //Load the libaries 
  TMVA::Tools::Instance();

    
  TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");

  //Define the variables.
  for(auto& Metric: MetricMap){
    if(!Metric.second.FillForMVA()){continue;}
    dataloader->AddVariable(Metric.first.c_str(),Metric.first.c_str(),"units",'F');
  }

  double signalWeight = 1;
  double backgroundWeight = 1;
  
  // You can add an arbitrary number of signal or background trees
  dataloader->AddSignalTree    (SignalTreeMVA,     signalWeight);
  dataloader->AddBackgroundTree(BackgroundTreeMVA, backgroundWeight);

  //event weights 
  dataloader->SetSignalWeightExpression("weight");
  dataloader->SetBackgroundWeightExpression("weight");


  TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
  TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

  //  dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
  //  					  "nTrain_Signal=14000:nTrain_Background=7000:SplitMode=Random:NormMode=None:!V" );

  dataloader->PrepareTrainingAndTestTree(mycuts, mycutb, "nTest_Signal=1:nTest_Background=1:NormMode=None:!V");

  TString outfileName( "TMVA.root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  //TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
  //					      "!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification:VerboseLevel=Debug" );


  //State the training methods. 

  //Standard Cut Based
  //  factory->BookMethod( dataloader, TMVA::Types::kCuts, "Cuts",
  //		       "!H:!V:FitMethod=MC:EffSel:SampleSize=10000:VarProp=FSmart" );
   //Some whack deep learning
  //  factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

  //Mr BDT
  //  factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
  //		       "!H:!V:NTrees=500:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );
  //  factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTSimp",
  //		       "!H:!V:NTrees=200:MaxDepth=2" );

  //  factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTGNotNorm",
  //		       "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2");


  float NTrees=100;
  float MinNodeSize=2.0;
  float Shrinkage=1.0;
  float BaggedSampleFraction=1.0;
  float MaxDepth=10.0;

  std::string options = "!H:NTrees=" + std::to_string(NTrees) + ":MinNodeSize=" + std::to_string(MinNodeSize) + "%:BoostType=Grad:Shrinkage=" + std::to_string(Shrinkage) + ":UseBaggedBoost:BaggedSampleFraction=" + std::to_string(BaggedSampleFraction) + ":nCuts=20:MaxDepth=" + std::to_string(MaxDepth) + ":SkipNormalization=False";
  
  //TMVA::MethodBase* methodbase = factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTGNotNorm",options.c_str());

  //  factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTGNotNorm",
  //			 "!H:!V:NTrees=500:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2:SkipNormalization=False");

  //  factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG",
  //		       "!H:!V:NTrees=500:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2:SkipNormalization=False");


  UInt_t numFolds = 10;
  TString analysisType = "Classification";
  TString splitType = "Random";
  TString splitExpr = "";

  TString cvOptions = Form(
			   ":!V"
			   ":!Silent"
			   ":ModelPersistence"                             
			   ":AnalysisType=%s"
			   ":SplitType=%s"
			   ":NumFolds=%i"
			   ":SplitExpr=%s",
			   analysisType.Data(), splitType.Data(), numFolds,
			   splitExpr.Data());
  
  cvOptions += ":OutputEnsembling=Avg";
  
  TMVA::CrossValidation cv{"TMVACrossValidation", dataloader,outputFile, cvOptions};
  
  cv.BookMethod(TMVA::Types::kBDT, "BDTGCrossVal", options.c_str());
  cv.Evaluate();
     
  const std::vector<TMVA::CrossValidationResult> results = cv.GetResults();
  for(auto const& result: results){ 
    result.Print();
  }
  
  
  //TMVA::Factory cv_factory = cv.GetFactory();
  //  cv.GetFactory().TestAllMethods();
  for(auto const& result: results){
    std::ofstream myfile;
    std::cout << "result.GetROCAverage(): " << result.GetROCAverage() << " ROCStdDeV: " << result.GetROCStandardDeviation() << std::endl;
    if(result.GetROCStandardDeviation() > 0.01){
          myfile.open ("histcomp.txt");
	  myfile << 0;
	  myfile.close();
    }
    else{
      myfile.open ("histcomp.txt");
      myfile << result.GetROCAverage();
      myfile.close();
      std::cout << "result.GetROCAverage(): " << result.GetROCAverage() << " ROCStdDeV: " << result.GetROCStandardDeviation() << std::endl;
      std::cout << "ROCAverage/ROCStdDeV: " << result.GetROCAverage() << std::endl;
      result.Print();
    }
   }


  //  cv.GetFactory().EvaluateAllMethods();

  // Now you can tell the factory to train, test, and evaluate the MVAs
  //
  //  Train MVAs using the set of training events
  //   factory->TrainAllMethods();
  //  Evaluate all MVAs using the set of test events
  //    factory->TestAllMethods();
  //Evaluate and compare performance of all configured MVAs
  //   factory->EvaluateAllMethods();
  
  //TMVA::MethodBDT* method = dynamic_cast<TMVA::MethodBDT*>(methodbase);
  

  //  double background_ks = method->GetKSTrainingVsTest('b',"XD");
  //  double signal_ks     = method->GetKSTrainingVsTest('s',"XD");

  //  std::cout << "signal_ks: " << signal_ks << " background_ks: " << background_ks << std::endl;


  // --------------------------------------------------------------
  // Save the output

  outputFile->Close();
  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;
  //delete factory;
  //  delete dataloader;
  
  return 0;
}

int optimiser::NueRecoOptimiser::SetupMVA(){

  reader = new TMVA::Reader( "!Color:!Silent" );

  //Define the variables.
  for(auto& Metric: MetricMap){
    if(!Metric.second.FillForMVA()){continue;}
    std::cout << "adding variable: " << Metric.first << std::endl;
    reader->AddVariable(Metric.first.c_str(),&Metric.second.GetMVAVal());
  }

  //  std::vector<TString> mvamethods = {"BDTG"};
  std::vector<TString> mvamethods = {};  
  for(auto const method: mvamethods){
    TString dir    = "dataset/weights/";
    TString prefix = "TMVAClassification";
    TString weightfile = dir + prefix + TString("_") + method + TString(".weights.xml");
    reader->BookMVA( method, weightfile);
  }
  
  std::vector<TString> cross_val_mvamethods = {"BDTGCrossVal","BDTGCrossValScaled"};
  for(auto const method: cross_val_mvamethods){
    TString dir    = "dataset/weights/";
    TString prefix = "TMVACrossValidation";
    TString weightfile = dir + prefix + TString("_") + method + TString(".weights.xml");
    reader->BookMVA( method, weightfile);
  }


  return 0; 
}

int optimiser::NueRecoOptimiser::MVAAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals, TString& MVAMethod){

  //Get the values
  std::string branchname = Metric.GetBranchName();
  std::vector<float>* branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  
  //Loop over the neutrinos
  vals.push_back(reader->EvaluateMVA(MVAMethod));

  return 0;
}

int optimiser::NueRecoOptimiser::EfficiencyNeutrinoMaker(optimiser::MetricHolder& Metric, std::string& sample){


  //Get the metric branch 
  std::string branchname = Metric.GetBranchName();
  std::vector<float> * branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  //Get the partner metric branch
  std::string partnername = Metric.GetPartnerMetricName();
  std::vector<float>* partner_vals;
  err = GetBranchPointer(partnername,sample,partner_vals);
  if(err){return 1;}


  std::string NeutrinoRecoEBeforeName = "NeutrinoRecoEBefore";
  std::string NeutrinoRecoEBeforeExtraName = "NeutrinoRecoEBeforeExtra"; 
  std::string NeutrinoRecoEAfterName = "NeutrinoRecoEAfter";
  std::string NeutrinoRecoEAfterExtraName = "NeutrinoRecoEAfterExtra"; 

  if(!Metric.CheckHistogram(NeutrinoRecoEBeforeName,sample)){Metric.SetHistogram(NeutrinoRecoEBeforeName,50,0,1000,sample);}
  if(!Metric.CheckHistogram(NeutrinoRecoEBeforeExtraName,sample)){Metric.SetHistogram(NeutrinoRecoEBeforeExtraName,50,0,1000,sample);}
  if(!Metric.CheckHistogram(NeutrinoRecoEAfterName,sample)){Metric.SetHistogram(NeutrinoRecoEAfterName,50,0,1000,sample);}
  if(!Metric.CheckHistogram(NeutrinoRecoEAfterExtraName,sample)){Metric.SetHistogram(NeutrinoRecoEAfterExtraName,50,0,1000,sample);}
 
  // for(auto const&partner_val: *partner_vals){
  //   Metric.GetHistogram(NeutrinoRecoEBeforeName,sample)->Fill(partner_val);
  // }

  std::vector<float> alreadyadded;
  int i=-1;
  for(auto const&branch_val: *branch_vals){
    Metric.GetHistogram(NeutrinoRecoEBeforeName,sample)->Fill(branch_val);
    ++i;
     // if(std::find(partner_vals->begin(),partner_vals->end(),branch_val) == partner_vals->end()){
    //   std::cerr << "why does your efficiency after not exist you silly billy" << std::endl;
    //   return 1;
    // }
    if(partner_vals->at(i) != 1){
      continue;
    }

    if(std::find(alreadyadded.begin(),alreadyadded.end(),branch_val) != alreadyadded.end()){
      Metric.GetHistogram(NeutrinoRecoEAfterExtraName,sample)->Fill(branch_val);
      continue;
    }

    Metric.GetHistogram(NeutrinoRecoEAfterName,sample)->Fill(branch_val);
    alreadyadded.push_back(branch_val);
  }
    


  return 0;
}

int optimiser::NueRecoOptimiser::MVAAnalysisAllsCut(optimiser::MetricHolder& Metric, std::string& sample, TString& MVAMethod){

  //Get the metric branch
  std::string branchname = Metric.GetBranchName();
  std::vector<float> * branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}


  //Get the nu reco energy  
  std::string nu_reco_energy = "nu_reco_energy";
  std::vector<float>* nu_reco_energy_vals;
  err = GetBranchPointer(nu_reco_energy,sample,nu_reco_energy_vals);
  if(err){return 1;}

  //Get the nu truth energy  
  std::string nu_truth_energy = "nu_truth_energy";
  std::vector<float> * nu_truth_energy_vals;
  err = GetBranchPointer(nu_truth_energy,sample,nu_truth_energy_vals);
  if(err){return 1;}

  //Get the nu truth energy  
  std::string nu_E = "nu_E";
  std::vector<float> * nu_E_vals;
  err = GetBranchPointer(nu_E,sample,nu_E_vals);
  if(err){return 1;}

  std::map<float,std::vector<float> > NeutrinoE;
  std::map<float,float > MaxNeutrinoE;
  std::map<float, int >  MaxNeutrinoE_iter;
  
  
  //Re initialise the maps
  for(int n_iter=0; n_iter<nu_E_vals->size(); ++n_iter){
    MaxNeutrinoE[nu_E_vals->at(n_iter)]      = -99999;
    MaxNeutrinoE_iter[nu_E_vals->at(n_iter)] = -99999;   
  }

  std::string NeutrinoTrueEBeforeName = "NeutrinoTrueEBefore" + (std::string) MVAMethod; 
  std::string NeutrinoRecoEBeforeName = "NeutrinoRecoEBeforeName" +(std::string) MVAMethod;
  std::string NeutrinoRecoEBeforeExtraName = "NeutrinoRecoEBeforeExtraName" +(std::string) MVAMethod;
  std::string NeutrinoRecoEAfterName = "NeutrinoRecoEAfterName" +(std::string) MVAMethod;
  std::string NeutrinoTrueEAfterName = "NeutrinoTrueEAfterName" +(std::string) MVAMethod;
  std::string NeutrinoRecoEAfterExtraName = "NeutrinoRecoEAfterExtraName" +(std::string) MVAMethod;
  std::string NeutrinoTrueEAfterExtraName = "NeutrinoTrueEAfterExtraName" +(std::string) MVAMethod;

  if(!Metric.CheckHistogram(NeutrinoTrueEBeforeName,sample)){Metric.SetHistogram(NeutrinoTrueEBeforeName,30,0,3000,sample);}
  if(!Metric.CheckHistogram(NeutrinoRecoEBeforeName,sample)){Metric.SetHistogram(NeutrinoRecoEBeforeName,30,0,3000,sample);}
  if(!Metric.CheckHistogram(NeutrinoRecoEBeforeExtraName,sample)){Metric.SetHistogram(NeutrinoRecoEBeforeExtraName,30,0,3000,sample);}
  if(!Metric.CheckHistogram(NeutrinoRecoEAfterName,sample)){Metric.SetHistogram(NeutrinoRecoEAfterName,30,0,3000,sample);}
  if(!Metric.CheckHistogram(NeutrinoTrueEAfterName,sample)){Metric.SetHistogram(NeutrinoTrueEAfterName,30,0,3000,sample);}
  if(!Metric.CheckHistogram(NeutrinoRecoEAfterExtraName,sample)){Metric.SetHistogram(NeutrinoRecoEAfterExtraName,30,0,3000,sample);}
  if(!Metric.CheckHistogram(NeutrinoTrueEAfterExtraName,sample)){Metric.SetHistogram(NeutrinoTrueEAfterExtraName,30,0,3000,sample);}
  
  
  //Loop over the neutrinos and get their true energy
  for(int n_iter=0; n_iter<nu_E_vals->size(); ++n_iter){
    Metric.GetHistogram(NeutrinoTrueEBeforeName,sample)->Fill(nu_truth_energy_vals->at(n_iter)*1000);
  }

  for(int n_iter=0; n_iter<nu_reco_energy_vals->size(); ++n_iter){
    //check if the true neutrino is already be identified 
    if(MaxNeutrinoE[nu_truth_energy_vals->at(n_iter)] < branch_vals->at(n_iter)){
      MaxNeutrinoE[nu_truth_energy_vals->at(n_iter)] = branch_vals->at(n_iter);
      MaxNeutrinoE_iter[nu_truth_energy_vals->at(n_iter)] = n_iter;
    }
  }
  //Add the max reco values to the efficiency and the others to the additional graph
  for(int n_iter=0; n_iter<nu_reco_energy_vals->size(); ++n_iter){
    if(MaxNeutrinoE[nu_truth_energy_vals->at(n_iter)] == branch_vals->at(n_iter)){
      Metric.GetHistogram(NeutrinoRecoEBeforeName,sample)->Fill(branch_vals->at(n_iter));
    }
    else{
      Metric.GetHistogram(NeutrinoRecoEBeforeExtraName,sample)->Fill(branch_vals->at(n_iter)); 
    }
  }

  //Perform the BDT Cut
  
  //Loop over the neutrinos
  for(int n_iter=0; n_iter<nu_reco_energy_vals->size(); ++n_iter){

    //Perform the BDT Cut
    if(reader->EvaluateMVA(MVAMethod) < 0){continue;}

    //Save the Neutrino E
    NeutrinoE[nu_truth_energy_vals->at(n_iter)].push_back(branch_vals->at(n_iter));
  }

  //Identify if this neutrino is the largest
  for(auto const& trueNeutrino: NeutrinoE){
    float MaxNeutrinoEval      = -99999;
    int   MaxNeutrinoEval_iter = -99999; 
    //Identify the remianing reco neutrino with the max energy
    for(int n_iter=0; n_iter<(trueNeutrino.second.size()); ++n_iter){
      if(MaxNeutrinoEval < trueNeutrino.second.at(n_iter)){
	MaxNeutrinoEval = trueNeutrino.second.at(n_iter);
	MaxNeutrinoEval_iter = n_iter;
      }
    }
    for(int n_iter=0; n_iter<(trueNeutrino.second.size()); ++n_iter){ 
      //Add the biggest to account for the efficiency.
      if(n_iter == MaxNeutrinoEval_iter){
	Metric.GetHistogram(NeutrinoRecoEAfterName,sample)->Fill(trueNeutrino.second.at(n_iter));
	Metric.GetHistogram(NeutrinoTrueEAfterName,sample)->Fill(trueNeutrino.first*1000);
      }
      //Add the extra reco
      else{
	Metric.GetHistogram(NeutrinoRecoEAfterExtraName,sample)->Fill(trueNeutrino.second.at(n_iter));
	Metric.GetHistogram(NeutrinoTrueEAfterExtraName,sample)->Fill(trueNeutrino.first*1000);
      }
    }
  }
return 0;
}

  
