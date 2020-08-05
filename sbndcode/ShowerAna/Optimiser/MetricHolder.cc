//Framework includes
#include "MetricHolder.hh"

//Root includes
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TGraph2D.h"
#include "TProfile.h" 
#include "TEfficiency.h"

//C++ includes
#include <iostream>

void optimiser::MetricHolder::Perform2DCutFinder(){

  std::cout << "Perform2DCutFinder" << std::endl;

  //  SignalHistPartner->Scale(1/TotalSigPOT);
  //  BackgroundHistPartner->Scale(1/TotalBKPOT);
  
  int max_ind = -999;
  int maxsigbk_ind = -999;
  int n_x = SignalHistPartner->GetNbinsX();
  int n_y = SignalHistPartner->GetNbinsY();

  double bktot_err = 0;
  double sigtot_err = 0;
  double bktot = BackgroundHistPartner->IntegralAndError(1, n_x+1,1,n_y+1, bktot_err);
  double sigtot = SignalHistPartner->IntegralAndError(1, n_x+1,1,n_y+1, sigtot_err);
  double maxeffpur = -999;
  double maxsigbk = -999;
  double efficiency[n_x*n_y];
  double purity[n_x*n_y];
  double bkRej[n_x*n_y];
  double effpur[n_x*n_y];
  double effpurErr[n_x*n_y];
  double sigbk[n_x*n_y];
  double sigbkErr[n_x*n_y];
  double xval[n_x*n_y];
  double yval[n_x*n_y];

  int maxeffpur_i;
  int maxeffpur_j;
  int maxsigbk_i;
  int maxsigbk_j;

  double maxeffpur_err = 0;
  double maxsigbkErr = 0;
  double maxefficiency_effpur_err = 0;
  double maxpurity_err = 0;
  double maxefficiency_effbk_err = 0;
  double maxbackgnd_err = 0;

  //Get the values.
  for(int j=1; j<n_y+1; ++j){
    for(int i=1; i<n_x+1; ++i){
      
      int ipj = i-1 + n_x*(j-1);

      double presig_err = 0;
      double prebk_err = 0;
      double presig = 0;
      double prebk =  0;
      
      if(fLessThan && fPartnerLessThan && fAND){
	presig = SignalHistPartner->IntegralAndError(1, i, 1, j, presig_err);
	prebk  = BackgroundHistPartner->IntegralAndError(1, i, 1, j, prebk_err);
      }
      else if(!fLessThan && fPartnerLessThan && fAND){
	presig = SignalHistPartner->IntegralAndError(i,n_x+1, 1, j, presig_err);
	prebk  = BackgroundHistPartner->IntegralAndError(i,n_x+1, 1, j, prebk_err);
      }
      else if(fLessThan && !fPartnerLessThan && fAND){
        presig = SignalHistPartner->IntegralAndError(1,i, j, n_y+1, presig_err);
        prebk  = BackgroundHistPartner->IntegralAndError(1,i,j, n_y+1, prebk_err);
      }
      else if(!fLessThan && !fPartnerLessThan && fAND){
        presig = SignalHistPartner->IntegralAndError(i,n_x+1, j,n_y+1, presig_err);
        prebk  = BackgroundHistPartner->IntegralAndError(i,n_x+1,j,n_y+1, prebk_err);
      }
      else if(fLessThan && fPartnerLessThan && !fAND){
	double presig_err_two=0;
	presig = SignalHistPartner->IntegralAndError(1, i, j+1, n_y+1, presig_err) + SignalHistPartner->IntegralAndError(1, n_x+1, 1, j, presig_err_two);
	presig_err = TMath::Sqrt(presig_err*presig_err + presig_err_two*presig_err_two);
	double prebk_err_two=0;
	prebk  = BackgroundHistPartner->IntegralAndError(1, i, j+1, n_y+1, prebk_err) + BackgroundHistPartner->IntegralAndError(1, n_x+1, 1, j, prebk_err_two);
	prebk_err = TMath::Sqrt(prebk_err*prebk_err + prebk_err_two*prebk_err_two);
      }
      else if(!fLessThan && fPartnerLessThan && !fAND){
	double presig_err_two=0;
	presig = SignalHistPartner->IntegralAndError(i, n_x+1, j+1, n_y+1, presig_err) + SignalHistPartner->IntegralAndError(1, n_x+1, 1, j, presig_err_two);
	presig_err = TMath::Sqrt(presig_err*presig_err + presig_err_two*presig_err_two);
	double prebk_err_two=0;
	prebk  = BackgroundHistPartner->IntegralAndError(i, n_x+1, j+1, n_y+1, prebk_err) + BackgroundHistPartner->IntegralAndError(1, n_x+1, 1, j, prebk_err_two);
	prebk_err = TMath::Sqrt(prebk_err*prebk_err + prebk_err_two*prebk_err_two);
      }
      else if(fLessThan && !fPartnerLessThan && !fAND){
	double presig_err_two=0;
	presig = SignalHistPartner->IntegralAndError(1,i, 1, j-1, presig_err) + SignalHistPartner->IntegralAndError(1, n_x+1,j,n_y+1, presig_err_two);
	presig_err = TMath::Sqrt(presig_err*presig_err + presig_err_two*presig_err_two);
	double prebk_err_two=0;
	prebk  = BackgroundHistPartner->IntegralAndError(1,i,1, j-1, prebk_err) + BackgroundHistPartner->IntegralAndError(1, n_x+1,j,n_y+1, prebk_err_two);
	prebk_err = TMath::Sqrt(prebk_err*prebk_err + prebk_err_two*prebk_err_two);
      }
      else if(!fLessThan && !fPartnerLessThan && !fAND){
	double presig_err_two=0;
	presig = SignalHistPartner->IntegralAndError(i,n_x+1,1, j-1, presig_err) + SignalHistPartner->IntegralAndError(1, n_x+1,j,n_y+1, presig_err_two);
	presig_err = TMath::Sqrt(presig_err*presig_err + presig_err_two*presig_err_two);
	double prebk_err_two=0;
	prebk  = BackgroundHistPartner->IntegralAndError(i,n_x+1,1, j-1, prebk_err) + BackgroundHistPartner->IntegralAndError(1, n_x+1,j,n_y+1, prebk_err_two);
	prebk_err = TMath::Sqrt(prebk_err*prebk_err + prebk_err_two*prebk_err_two);
      }



      double efficencyErr = 0;
      double bkRejErr     = 0;
      double purityErr    = 0; 
      double effpurErr    = 0;
      double sigbkErr     = 0;

      efficiency[ipj] = presig / sigtot;
      efficencyErr = efficiency[ipj] * TMath::Sqrt(((presig_err / presig) * (presig_err / presig)) + (sigtot_err / sigtot) * (sigtot_err / sigtot));
      
      bkRej[ipj] = 1 - prebk / bktot;
      bkRejErr = bkRej[ipj] * TMath::Sqrt((prebk_err / prebk) * (prebk_err / prebk) + (bktot_err / bktot) * (bktot_err / bktot));
      
      if (presig + prebk != 0) {
	purity[ipj] = presig / (presig + prebk);
	double purityErr_dem = TMath::Sqrt((presig_err * presig_err) + (prebk_err * prebk_err));
	purityErr = purity[ipj] * TMath::Sqrt(((presig_err / presig) * (presig_err / presig)) + (purityErr_dem / (presig + prebk)) * (purityErr_dem / (presig + prebk)));
      } 
      else {
      purity[ipj] = 0;
      purityErr = 0;
      }
      
      effpur[ipj] = efficiency[ipj] * purity[ipj];
      if (efficiency[ipj] == 0 || purity[ipj] == 0) {
	effpurErr = 0;
      } 
      else {
	effpurErr = effpur[ipj] * TMath::Sqrt(((efficencyErr / efficiency[ipj]) * (efficencyErr / efficiency[ipj])) + (purityErr / purity[ipj]) * (purityErr / purity[ipj]));
      }
      
      sigbk[ipj] = efficiency[ipj] * bkRej[ipj];
      if (efficiency[ipj] == 0 || bkRej[ipj] == 0) {
	sigbkErr = 0;
      } 
      else {
	sigbkErr = sigbk[ipj] * TMath::Sqrt(((efficencyErr / efficiency[ipj]) * (efficencyErr / efficiency[ipj])) + (bkRejErr / bkRej[ipj]) * (bkRejErr / bkRej[ipj]));
      }


      xval[ipj] = ((TAxis*)SignalHistPartner->GetXaxis())->GetBinCenter(i);
      yval[ipj] = ((TAxis*)SignalHistPartner->GetYaxis())->GetBinCenter(j);

      if(efficiency[ipj]>1){continue;}


      //Keep the largest value.
      if (effpur[ipj] > maxeffpur) {
	max_ind = ipj;
	maxeffpur = effpur[ipj];
	maxeffpur_err = sigbkErr;
	maxefficiency_effpur_err = efficencyErr;
	maxpurity_err = purityErr;
	maxeffpur_i = i;
	maxeffpur_j = j;
      }
      if (sigbk[ipj] > maxsigbk) {
	maxsigbk = sigbk[ipj];
	maxsigbk_ind = ipj;
	maxsigbkErr = sigbkErr;
	maxefficiency_effbk_err = efficencyErr;
	maxbackgnd_err = bkRejErr;
	maxsigbk_i = i;
	maxsigbk_j = j;
      }
    }
  }

  //Draw the curves
  TGraph2D *efficency_graph  = new TGraph2D("eff2D","eff2D",n_x*n_y,xval,yval,efficiency);
  TGraph2D *background_graph = new TGraph2D("bkrej2D","bjreg2D",n_x*n_y,xval,yval,bkRej);
  TGraph2D *purity_graph     = new TGraph2D("purity2D","purity2D",n_x*n_y,xval,yval,purity);  
  TGraph2D *effbk_graph      = new TGraph2D("effbk2D","effbk2D",n_x*n_y,xval,yval,sigbk);
  TGraph2D *effpur_graph     = new TGraph2D("effpur2D","effpur2D",n_x*n_y,xval,yval,effpur);

  //Draw the histograms 
  SignalHistPartner->Write();
  BackgroundHistPartner->Write();

  auto signal_canvas = new TCanvas("signal2D_canvas", "signal2D_canvas", 600, 400);
  SignalHistPartner->Draw("colz");
  signal_canvas->Write();

  auto bkgrd_canvas = new TCanvas("background2D_canvas", "background2D_canvas", 600, 400);
  BackgroundHistPartner->Draw("colz");
  bkgrd_canvas->Write();

  //Draw the profiles
  TProfile* SignalProfileX = SignalHistPartner->ProfileX();
  TProfile* SignalProfileY = SignalHistPartner->ProfileY();
  TProfile* BackgroundProfileX = BackgroundHistPartner->ProfileX();
  TProfile* BackgroundProfileY = BackgroundHistPartner->ProfileY();

  TGraphErrors* RMS = new TGraphErrors();
  for(int i=0;i<n_x;++i){
    RMS->SetPoint(RMS->GetN(),SignalProfileX->GetXaxis()->GetBinCenter(i),TMath::Sqrt(SignalProfileX->GetBinEntries(i))*SignalProfileX->GetBinError(i));
    RMS->SetPointError(RMS->GetN()-1,SignalProfileX->GetXaxis()->GetBinWidth(i)/2,TMath::Sqrt(SignalProfileX->GetBinEntries(i))*SignalProfileX->GetBinError(i)/TMath::Sqrt(2*(SignalProfileX->GetBinEntries(i)-1)));
  }
  RMS->SetName("RMS");
  RMS->Write();

  SignalProfileX->GetYaxis()->SetTitle("(True Energy - Reco Energy)/True Energy");
  SignalProfileX->SetTitle("");
  SignalProfileX->GetYaxis()->SetRangeUser(0,1);

  BackgroundProfileX->SetTitle("");


  SignalProfileX->SetName("Signal");
  BackgroundProfileX->SetName("Background");


  BackgroundProfileX->SetLineColor(kRed);
  BackgroundProfileY->SetLineColor(kRed);

  auto xprofile_canvas = new TCanvas("xprofile_canvas", "xprofile_canvas", 600, 400);
  SignalProfileX->Draw();
  BackgroundProfileX->Draw("SAME");
  xprofile_canvas->BuildLegend();
  xprofile_canvas->Write();

  auto yprofile_canvas = new TCanvas("yprofile_canvas", "yprofile_canvas", 600, 400);
  SignalProfileY->Draw();
  BackgroundProfileY->Draw("SAME");
  yprofile_canvas->BuildLegend();
  yprofile_canvas->Write();

  SignalProfileX->Write();
  SignalProfileY->Write();
  BackgroundProfileX->Write();
  BackgroundProfileY->Write();

  efficency_graph->Write();
  background_graph->Write();
  purity_graph->Write();
  effbk_graph->Write();
  effpur_graph->Write();
  
  auto efficency_canvas = new TCanvas("eff2D_canvas", "eff2D_canvas", 600, 400);
  efficency_graph->Draw("colz");
  efficency_canvas->Write();

  auto background_canvas = new TCanvas("bk2D_canvas", "bk2D_canvas", 600, 400);
  background_graph->Draw("colz");
  background_canvas->Write();

  auto purity_canvas = new TCanvas("pur2D_canvas", "pur2D_canvas", 600, 400);
  purity_graph->Draw("colz");
  purity_canvas->Write();

  auto effbk_canvas = new TCanvas("effbk2D_canvas", "effbk2D_canvas", 600, 400);
  effbk_graph->Draw("colz");
  effbk_canvas->Write();

  auto effpur_canvas = new TCanvas("effpur2D_canvas", "effpur2D_canvas", 600, 400);
  effpur_graph->Draw("colz");
  effpur_canvas->Write();

  delete effpur_canvas;
  delete effbk_canvas;
  delete purity_canvas;
  delete background_canvas;
  delete efficency_canvas;
  delete efficency_graph;
  delete background_graph;
  delete purity_graph;
  delete effbk_graph;
  delete effpur_graph;
  delete xprofile_canvas;
  delete yprofile_canvas;

  //Make the 1D projections
  if(fPartnerLessThan){
    Perform1DCutFinder(SignalHistPartner->ProjectionX("_pxsig",1,maxsigbk_j,"e"), BackgroundHistPartner->ProjectionX("_pxbk",1,maxsigbk_j ,"e")); 
  }
  else{
    Perform1DCutFinder(SignalHistPartner->ProjectionX("_pxsig",maxsigbk_j,n_y+1,"e"), BackgroundHistPartner->ProjectionX("_pxbk",maxsigbk_j,n_y+1 ,"e"));
  }
  if(fLessThan){
    Perform1DCutFinder(SignalHistPartner->ProjectionY("_pysig",1,maxsigbk_i,"e"), BackgroundHistPartner->ProjectionY("_pjbk",1,maxsigbk_i,"e")); 
  }
  else{
    Perform1DCutFinder(SignalHistPartner->ProjectionY("_pysig",maxsigbk_i,n_x+1,"e"), BackgroundHistPartner->ProjectionY("_pjbk",maxsigbk_i,n_x+1,"e")); 
  }

  std::cout << "Best Cut is at x: " << xval[max_ind] << " y: " << yval[max_ind] << std::endl;
  std::cout << "Efficiency: " << efficiency[max_ind] << " +- " << maxefficiency_effpur_err << std::endl;
  std::cout << "Purity    : " << purity[max_ind] << " +- " << maxpurity_err << std::endl;
  std::cout << "Eff*Purity: " << effpur[max_ind] << " +- " << maxeffpur_err << std::endl;

  std::cout << "Best sig*bkrej Cut is at x: " << xval[maxsigbk_ind] << " y: " << yval[maxsigbk_ind] << std::endl;
  std::cout << "Efficiency:     " << efficiency[maxsigbk_ind] << " +- " << maxefficiency_effbk_err << std::endl;
  std::cout << "Background Rej: " << bkRej[maxsigbk_ind] << " +- " << maxbackgnd_err << std::endl;
  std::cout << "Eff*Purity:     " << sigbk[maxsigbk_ind] << " +- " << maxsigbkErr << std::endl;


  return;
}

void optimiser::MetricHolder::Perform1DCutFinder(){
  //  Perform1DCutFinder(SignalHist,BackgroundHist,SignalHistNotNorm,BackgroundHistNotNorm);
  Perform1DCutFinder(SignalHist,BackgroundHist);
}

void optimiser::MetricHolder::Perform1DCutFinder(TH1D* signalhist, TH1D* backgroundhist){

  // TH1::AddDirectory(kFALSE);
  //  signalhist->Scale(1/TotalSigPOT);
  //  backgroundhist->Scale(1/TotalBKPOT);
  
  //Find the min and max start positions
  float minsig = signalhist->GetBinCenter(0);
  float maxsig = signalhist->GetBinCenter(signalhist->GetNbinsX());
  float minbk = backgroundhist->GetBinCenter(0);
  float maxbk = backgroundhist->GetBinCenter(backgroundhist->GetNbinsX());

  if (minsig != minbk) {
    std::cout << "min point does not match" << std::endl;
    return;
  }
  if (maxbk != maxbk) {
    std::cout << "max point does not match" << std::endl;
    return;
  }
  if (signalhist->GetNbinsX() != backgroundhist->GetNbinsX()) {
    std::cout << "bin numbers sig " << signalhist->GetNbinsX() << " and background " << backgroundhist->GetNbinsX() << "do not match" << std::endl;
    return;
  }

  int max_ind = -999;
  int maxsigbk_ind = -999;
  int n = signalhist->GetNbinsX();
  double bktot_err = 0;
  double bktot = backgroundhist->IntegralAndError(1, n+1, bktot_err);
  double sigtot_err = 0;
  double sigtot = signalhist->IntegralAndError(1, n+1, sigtot_err);
  double maxeffpur = -999;
  double maxsigbk = -999;
  double efficiency[n];
  double efficencyErr[n];
  double purity[n];
  double purityErr[n];
  double bkRej[n];
  double bkRejErr[n];
  double effpur[n];
  double effpurErr[n];
  double sigbk[n];
  double sigbkErr[n];
  double xval[n];
  double xvalErr[n];
  double maxpurity = -99;
  int maxpurity_iter = -99;
  //Get the values.
  for (int i = 1; i < n+1; ++i) {

    double presig_err = 0;
    double prebk_err = 0;
    double presig = 0;
    double prebk  = 0;

    if(fLessThan){
      presig = signalhist->IntegralAndError(1, i, presig_err);
      prebk  = backgroundhist->IntegralAndError(1, i, prebk_err);
    }
    else{
      presig = signalhist->IntegralAndError(i, n+1, presig_err);
      prebk  = backgroundhist->IntegralAndError(i, n+1, prebk_err);
    }

    //    std::cout << "presig: " << presig << " prebk: " << prebk << std::endl;

    efficiency[i-1] = (presig / sigtot);
    efficencyErr[i-1] = efficiency[i-1] * TMath::Sqrt(((presig_err / presig) * (presig_err / presig)) + (sigtot_err / sigtot) * (sigtot_err / sigtot));

    bkRej[i-1] = 1 - ((prebk / bktot));
    bkRejErr[i-1] = bkRej[i-1] * TMath::Sqrt((prebk_err / prebk) * (prebk_err / prebk) + (bktot_err / bktot) * (bktot_err / bktot));

    if (presig + prebk != 0) {
      purity[i-1] = presig / (presig + prebk);
      double purityErr_dem = TMath::Sqrt((presig_err * presig_err) + (prebk_err * prebk_err));
      purityErr[i-1] = purity[i-1] * TMath::Sqrt(((presig_err / presig) * (presig_err / presig)) + (purityErr_dem / (presig + prebk)) * (purityErr_dem / (presig + prebk)));
    } else {
      purity[i-1] = 0;
      purityErr[i-1] = 0;
    }


    effpur[i-1] = efficiency[i-1] * purity[i-1];
    if (efficiency[i-1] == 0 || purity[i-1] == 0) {
      effpurErr[i-1] = 0;
    } else {
      effpurErr[i-1] = effpur[i-1] * TMath::Sqrt(((efficencyErr[i-1] / efficiency[i-1]) * (efficencyErr[i-1] / efficiency[i-1])) + (purityErr[i-1] / purity[i-1]) * (purityErr[i-1] / purity[i-1]));
    }

    if(purity[i-1]+ effpur[i-1] >  maxpurity){
      maxpurity = purity[i-1]+ effpur[i-1];
      maxpurity_iter = i-1;
    }

    sigbk[i-1] = efficiency[i-1] * bkRej[i-1];
    if (efficiency[i-1] == 0 || bkRej[i-1] == 0) {
      sigbkErr[i-1] = 0;
    } else {
      sigbkErr[i-1] = sigbk[i-1] * TMath::Sqrt(((efficencyErr[i-1] / efficiency[i-1]) * (efficencyErr[i-1] / efficiency[i-1])) + (bkRejErr[i-1] / bkRej[i-1]) * (bkRejErr[i-1] / bkRej[i-1]));
    }

    xval[i-1] = signalhist->GetBinCenter(i);
    xvalErr[i-1] = signalhist->GetBinWidth(i)/2;

    //    std::cout << "xval " << xval[i-1] << " efficiency: " << efficiency[i-1] << " purity: " << purity[i-1] << " eff * pur "
    //          << effpur[i-1] << "background Rejection: " << bkRej[i-1] << " and sig*bkrej: " << sigbk[i-1] << std::endl;

    if(efficiency[i-1]>1){continue;}

    //Keep the largest value.
    if (effpur[i-1] > maxeffpur) {
      max_ind = i-1;
      maxeffpur = effpur[i-1];
    }
    if (sigbk[i-1] > maxsigbk) {
      maxsigbk = sigbk[i-1];
      maxsigbk_ind = i-1;
    }
  }

  //Draw the curves
  auto effpur_canvas = new TCanvas("effpur_canvas", "effpur_canvas", 600, 400);
  auto sigbk_canvas = new TCanvas("sigbk_canvas", "sigbk_canvas", 600, 400);

  TGraphErrors* efficiencyPlot = new TGraphErrors(n, xval, efficiency, xvalErr, efficencyErr);
  efficiencyPlot->GetYaxis()->SetTitle("Efficiency");
  efficiencyPlot->SetTitle("Efficiency");
  efficiencyPlot->SetLineColor(kBlue);
  efficiencyPlot->Draw("a4");
  efficiencyPlot->Write();

  TGraphErrors* bkRejPlot = new TGraphErrors(n, xval, bkRej, xvalErr, bkRejErr);
  bkRejPlot->GetYaxis()->SetTitle("Background Rejection");
  bkRejPlot->SetTitle("Background Rejection");
  bkRejPlot->SetLineColor(kRed);
  bkRejPlot->Draw("a4");
  bkRejPlot->Write();

  TGraphErrors* purityPlot = new TGraphErrors(n, xval, purity, xvalErr, purityErr);
  purityPlot->GetYaxis()->SetTitle("Purity");
  purityPlot->SetTitle("Purity");
  purityPlot->SetLineColor(kRed);
  purityPlot->Draw("a4");
  purityPlot->Write();

  TGraphErrors* effxpurPlot = new TGraphErrors(n, xval, effpur, xvalErr, effpurErr);
  effxpurPlot->GetYaxis()->SetTitle("Effxpur");
  effxpurPlot->SetTitle("Efficiency x Purity");
  effxpurPlot->SetLineColor(kBlack);
  effxpurPlot->Draw("a4");
  effxpurPlot->Write();

  TGraphErrors* sigbkPlot = new TGraphErrors(n, xval, sigbk, xvalErr, sigbkErr);
  sigbkPlot->GetYaxis()->SetTitle("Efficiency x Background Rejection");
  sigbkPlot->SetTitle("Efficiency x Background Rejection");
  sigbkPlot->SetLineColor(kBlack);
  sigbkPlot->Draw("a4");
  sigbkPlot->Write();

  TGraph* roccurve = new TGraph(n, efficiency,bkRej);
  std::cout << "writing roc curve" << std::endl;
  roccurve->GetYaxis()->SetTitle("Efficiency");
  roccurve->GetXaxis()->SetTitle("Background Rejection");
  roccurve->Write();

  //Make the graph to present
  effpur_canvas->cd();
  auto mg1 = new TMultiGraph("mg1", "mg1");
  mg1->Add(efficiencyPlot);
  mg1->Add(purityPlot);
  mg1->Add(effxpurPlot);
  mg1->GetYaxis()->SetRangeUser(0, 1);
  mg1->Write();
  mg1->Draw("AP");
  effpur_canvas->BuildLegend();
  effpur_canvas->Write();

  purityPlot->SetLineColor(kRed);
  effxpurPlot->SetLineColor(kBlack);

  //Make the graph to present
  sigbk_canvas->cd();
  auto mg2 = new TMultiGraph("mg2", "mg2");
  mg2->Add(efficiencyPlot);
  mg2->Add(bkRejPlot);
  mg2->Add(sigbkPlot);
  //  mg2->Add(purityPlot);
  // mg2->Add(effxpurPlot);

  auto mg3 = new TMultiGraph("mg3", "mg3");
  mg3->Add(purityPlot);  
  mg3->Add(effxpurPlot); 


  mg2->GetYaxis()->SetRangeUser(0, 1);
  mg2->Write();

  mg2->Draw("AP");
  sigbk_canvas->BuildLegend(0.6, 0.4, 0.9, 0.6);
  sigbk_canvas->Write();

  signalhist->Scale(1/signalhist->Integral());
  backgroundhist->Scale(1/backgroundhist->Integral());

  double maxScale = TMath::Max(backgroundhist->GetBinContent(backgroundhist->GetMaximumBin()),
			       signalhist->GetBinContent(signalhist->GetMaximumBin()));

  signalhist->SetLineColor(kBlue);
  signalhist->SetFillStyle(3003);
  signalhist->SetFillColor(6);
  signalhist->Scale(1. / maxScale);
  signalhist->SetTitle("Signal");
  signalhist->Write();
  backgroundhist->SetLineColor(kRed);
  backgroundhist->SetFillStyle(3003);
  backgroundhist->SetFillColor(42);
  backgroundhist->Scale(1. / maxScale);
  backgroundhist->SetTitle("Background");
  backgroundhist->Write();

  bkRejPlot->GetXaxis()->SetLabelOffset(999);
  bkRejPlot->GetXaxis()->SetLabelSize(0);

  auto full_canvas = new TCanvas("full_canvas", "full_canvas", 600, 400);
  full_canvas->Draw();
  Float_t sf = 0.3;
  //  TPad *p3 = new TPad("p2","p2",0,1-sf,1,1);
  TPad *p2 = new TPad("p2","p2",0,0,1,sf);
  TPad *p1 = new TPad("p1","p1",0,sf,1,1);
  p2->SetTopMargin(0.02);
  //  p1->SetTopMargin(0.01);

  p1->SetBottomMargin(0.03);
  p2->SetBottomMargin(0.3);

  p1->cd();
  //  p1->Draw();
  mg2->SetTitle("");
  mg2->GetXaxis()->SetLabelOffset(999);
  mg2->GetXaxis()->SetLabelSize(0);
  Double_t upTM = p1->GetTopMargin();
  Double_t upBM = p1->GetBottomMargin();
  Double_t lowTM = p2->GetTopMargin();
  Double_t lowBM = p2->GetBottomMargin();
  Double_t ratio = ( (upBM-(1-upTM))*(1-sf) ) / ( (lowBM-(1-lowTM))*sf ) ;
  float textsize = mg2->GetYaxis()->GetLabelSize()*ratio; 
  mg2->Draw("AP");
  signalhist->Draw("HIST SAME e");
  backgroundhist->Draw("HIST SAME e");
  TLine* lfull_1 = new TLine(xval[maxsigbk_ind], 0, xval[maxsigbk_ind], 1);
  lfull_1->SetLineColor(kBlack);
  p1->BuildLegend();
  lfull_1->Draw("same");
  gPad->SetTickx(2);
  full_canvas->cd();
  p1->Draw();


  // p3->cd();
  // bkRejPlot->GetXaxis()->SetLabelOffset(999);
  // bkRejPlot->GetXaxis()->SetLabelSize(0);
  // bkRejPlot->GetYaxis()->SetLabelSize(textsize);
  // bkRejPlot->Draw("AP");
  // full_canvas->cd();
  // p3->Draw("same");
  
  //  full_canvas->cd();
  p2->cd();
  //  p2->cd();
  mg3->SetTitle("");
  mg3->GetYaxis()->SetLabelSize(0.08);
  mg3->GetXaxis()->SetLabelSize(textsize);
  mg3->GetXaxis()->SetTitleSize(textsize);
  mg3->GetYaxis()->SetRangeUser(0,maxpurity+0.001);
  mg3->Draw("AP");
  p2->BuildLegend();

  TLine* lfull_2 = new TLine(xval[max_ind], 0, xval[max_ind], maxpurity+0.001);
  lfull_2->SetLineColor(kBlack);
  lfull_2->Draw("same");
  full_canvas->cd();
  p2->Draw("same");

  
 
  full_canvas->Write();

  delete p2;
  delete p1;
  delete full_canvas;
  delete lfull_2;
  delete lfull_1;

  auto sigbk_cut_canvas = new TCanvas("sigbk_cut_canvas", "sigbk_cut_canvas", 600, 400);
  sigbk_cut_canvas->cd();
  mg2->GetHistogram()->SetTitle(Form("Cut at %0.1f " + units + " with efficiency %0.2f and background rejection %0.2f", xval[maxsigbk_ind], efficiency[maxsigbk_ind], bkRej[maxsigbk_ind]));
  mg2->GetXaxis()->SetTitle(axisTitle);
  mg2->Draw("AP");
  signalhist->Draw("HIST SAME e");
  backgroundhist->Draw("HIST SAME e");
  sigbk_cut_canvas->BuildLegend(0.6, 0.4, 0.9, 0.6);
  TLine* l1 = new TLine(xval[maxsigbk_ind], 0, xval[maxsigbk_ind], 1);
  l1->SetLineColor(kBlack);
  l1->Draw("same");

  sigbk_cut_canvas->Write();

  //Draw the overlayed samples with the cut.
  auto cut_canvas = new TCanvas("cut_canvas", "cut_canvas", 600, 400);
  cut_canvas->cd();
  mg1->GetHistogram()->SetTitle(Form("Cut at %0.1f " + units + " with efficiency %0.2f and purity %0.2f", xval[max_ind], efficiency[max_ind], purity[max_ind]));
  mg1->GetXaxis()->SetTitle(axisTitle);
  mg1->Draw("AP");
  signalhist->Draw("HIST SAME");
  backgroundhist->Draw("HIST SAME");

  cut_canvas->Update();
  cut_canvas->BuildLegend();

  TLine* l = new TLine(xval[max_ind], 0, xval[max_ind], 1);
  l->SetLineColor(kBlack);
  l->Draw("same");

  cut_canvas->Write();

  std::cout << "Best Cut is at: " << xval[max_ind] << std::endl;
  std::cout << "Efficiency: " << efficiency[max_ind] << " +- " << efficencyErr[max_ind] << std::endl;
  std::cout << "Purity    : " << purity[max_ind] << " +- " << purityErr[max_ind] << std::endl;
  std::cout << "Eff*Purity: " << effpur[max_ind] << " +- " << effpurErr[max_ind] << std::endl;

  std::cout << "Best sig*bkrej Cut is at: " << xval[maxsigbk_ind] << std::endl;
  std::cout << "Efficiency:     " << efficiency[maxsigbk_ind] << " +- " << efficencyErr[maxsigbk_ind] << std::endl;
  std::cout << "Background Rej: " << bkRej[maxsigbk_ind] << " +- " << bkRejErr[maxsigbk_ind] << std::endl;
  std::cout << "Eff*Purity:     " << sigbk[maxsigbk_ind] << " +- " << sigbkErr[maxsigbk_ind] << std::endl;


  delete efficiencyPlot;
  delete bkRejPlot;
  delete purityPlot;
  delete effxpurPlot;
  delete sigbkPlot;
  delete roccurve;

  delete mg1;
  delete mg2;

  delete effpur_canvas;
  delete sigbk_canvas;
  delete sigbk_cut_canvas;
  delete cut_canvas;
}


void optimiser::MetricHolder::Perform1DCutFinder(TH1D* signalhist, TH1D* backgroundhist,
						 TH1D* signalhistnotnorm, TH1D* backgroundhistnotnorm){

  // TH1::AddDirectory(kFALSE);
  //  signalhist->Scale(1/TotalSigPOT);
  //  backgroundhist->Scale(1/TotalBKPOT);
  
  //Find the min and max start positions
  float minsig = signalhist->GetBinCenter(0);
  float maxsig = signalhist->GetBinCenter(signalhist->GetNbinsX());
  float minbk = backgroundhist->GetBinCenter(0);
  float maxbk = backgroundhist->GetBinCenter(backgroundhist->GetNbinsX());

  if (minsig != minbk) {
    std::cout << "min point does not match" << std::endl;
    return;
  }
  if (maxbk != maxbk) {
    std::cout << "max point does not match" << std::endl;
    return;
  }
  if (signalhist->GetNbinsX() != backgroundhist->GetNbinsX()) {
    std::cout << "bin numbers sig " << signalhist->GetNbinsX() << " and background " << backgroundhist->GetNbinsX() << "do not match" << std::endl;
    return;
  }

  int max_ind = -999;
  int maxsigbk_ind = -999;
  int n = signalhist->GetNbinsX();
  double bktot_err = 0;
  double bktot = backgroundhist->IntegralAndError(1, n+1, bktot_err);
  double sigtot_err = 0;
  double sigtot = signalhist->IntegralAndError(1, n+1, sigtot_err);
  double maxeffpur = -999;
  double maxsigbk = -999;
  double efficiency[n];
  double efficencyErr[n];
  double purity[n];
  double purityErr[n];
  double bkRej[n];
  double bkRejErr[n];
  double effpur[n];
  double effpurErr[n];
  double sigbk[n];
  double sigbkErr[n];
  double xval[n];
  double xvalErr[n];

  //Get the values.
  for (int i = 1; i < n+1; ++i) {

    double presig_err = 0;
    double prebk_err = 0;
    double presig = 0;
    double prebk  = 0;

    if(fLessThan){
      presig = signalhist->IntegralAndError(1, i, presig_err);
      prebk  = backgroundhist->IntegralAndError(1, i, prebk_err);
    }
    else{
      presig = signalhist->IntegralAndError(i, n+1, presig_err);
      prebk  = backgroundhist->IntegralAndError(i, n+1, prebk_err);
    }

    //    std::cout << "presig: " << presig << " prebk: " << prebk << std::endl;

    efficiency[i-1] = presig / sigtot;
    efficencyErr[i-1] = efficiency[i-1] * TMath::Sqrt(((presig_err / presig) * (presig_err / presig)) + (sigtot_err / sigtot) * (sigtot_err / sigtot));

    bkRej[i-1] = 1 - prebk / bktot;
    bkRejErr[i-1] = bkRej[i-1] * TMath::Sqrt((prebk_err / prebk) * (prebk_err / prebk) + (bktot_err / bktot) * (bktot_err / bktot));

    if (presig + prebk != 0) {
      purity[i-1] = presig / (presig + prebk);
      double purityErr_dem = TMath::Sqrt((presig_err * presig_err) + (prebk_err * prebk_err));
      purityErr[i-1] = purity[i-1] * TMath::Sqrt(((presig_err / presig) * (presig_err / presig)) + (purityErr_dem / (presig + prebk)) * (purityErr_dem / (presig + prebk)));
    } else {
      purity[i-1] = 0;
      purityErr[i-1] = 0;
    }

    effpur[i-1] = efficiency[i-1] * purity[i-1];
    if (efficiency[i-1] == 0 || purity[i-1] == 0) {
      effpurErr[i-1] = 0;
    } else {
      effpurErr[i-1] = effpur[i-1] * TMath::Sqrt(((efficencyErr[i-1] / efficiency[i-1]) * (efficencyErr[i-1] / efficiency[i-1])) + (purityErr[i-1] / purity[i-1]) * (purityErr[i-1] / purity[i-1]));
    }

    sigbk[i-1] = efficiency[i-1] * bkRej[i-1];
    if (efficiency[i-1] == 0 || bkRej[i-1] == 0) {
      sigbkErr[i-1] = 0;
    } else {
      sigbkErr[i-1] = sigbk[i-1] * TMath::Sqrt(((efficencyErr[i-1] / efficiency[i-1]) * (efficencyErr[i-1] / efficiency[i-1])) + (bkRejErr[i-1] / bkRej[i-1]) * (bkRejErr[i-1] / bkRej[i-1]));
    }

    xval[i-1] = signalhist->GetBinCenter(i);
    xvalErr[i-1] = signalhist->GetBinWidth(i)/2;

    //    std::cout << "xval " << xval[i-1] << " efficiency: " << efficiency[i-1] << " purity: " << purity[i-1] << " eff * pur "
    //          << effpur[i-1] << "background Rejection: " << bkRej[i-1] << " and sig*bkrej: " << sigbk[i-1] << std::endl;

    //Keep the largest value.
    if (effpur[i-1] > maxeffpur) {
      max_ind = i-1;
      maxeffpur = effpur[i-1];
    }
    if (sigbk[i-1] > maxsigbk) {
      maxsigbk = sigbk[i-1];
      maxsigbk_ind = i-1;
    }
  }


  int max_ind_notnorm = -999;
  int maxsigbk_notnorm_ind_notnorm = -999;
  double bktot_notnorm_err_notnorm = 0;
  double bktot_notnorm = backgroundhist->IntegralAndError(1, n+1, bktot_notnorm_err_notnorm);
  double sigtot_notnorm_err_notnorm = 0;
  double sigtot_notnorm = signalhist->IntegralAndError(1, n+1, sigtot_notnorm_err_notnorm);
  double maxeffpur_notnorm_notnorm = -999;
  double maxsigbk_notnorm_notnorm = -999;
  double efficiency_notnorm[n];
  double efficencyErr_notnorm[n];
  double purity_notnomr[n];
  double purity_notnomrErr[n];
  double bkRej_notnorm[n];
  double bkRej_notnormErr[n];
  double effpur_notnorm[n];
  double effpur_notnormErr[n];
  double sigbk_notnorm[n];
  double sigbk_notnormErr[n];
  double xval_notnorm[n];
  double xval_notnormErr[n];

  //Get the values.
  for (int i = 1; i < n+1; ++i) {

    double presig_err = 0;
    double prebk_err = 0;
    double presig = 0;
    double prebk  = 0;

    if(fLessThan){
      presig = signalhist->IntegralAndError(1, i, presig_err);
      prebk  = backgroundhist->IntegralAndError(1, i, prebk_err);
    }
    else{
      presig = signalhist->IntegralAndError(i, n+1, presig_err);
      prebk  = backgroundhist->IntegralAndError(i, n+1, prebk_err);
    }

    //    std::cout << "presig: " << presig << " prebk: " << prebk << std::endl;

    efficiency_notnorm[i-1] = (presig / sigtot_notnorm);
    efficencyErr_notnorm[i-1] = efficiency_notnorm[i-1] * TMath::Sqrt(((presig_err / presig) * (presig_err / presig)) + (sigtot_notnorm_err_notnorm / sigtot_notnorm) * (sigtot_notnorm_err_notnorm / sigtot_notnorm));

    bkRej_notnorm[i-1] = 1 - ((prebk / bktot_notnorm));
    bkRej_notnormErr[i-1] = bkRej_notnorm[i-1] * TMath::Sqrt((prebk_err / prebk) * (prebk_err / prebk) + (bktot_notnorm_err_notnorm / bktot_notnorm) * (bktot_notnorm_err_notnorm / bktot_notnorm));

    if (presig + prebk != 0) {
      purity_notnomr[i-1] = presig / (presig + prebk);
      double purity_notnomrErr_dem = TMath::Sqrt((presig_err * presig_err) + (prebk_err * prebk_err));
      purity_notnomrErr[i-1] = purity_notnomr[i-1] * TMath::Sqrt(((presig_err / presig) * (presig_err / presig)) + (purity_notnomrErr_dem / (presig + prebk)) * (purity_notnomrErr_dem / (presig + prebk)));
    } else {
      purity_notnomr[i-1] = 0;
      purity_notnomrErr[i-1] = 0;
    }

    effpur_notnorm[i-1] = efficiency_notnorm[i-1] * purity_notnomr[i-1];
    if (efficiency_notnorm[i-1] == 0 || purity_notnomr[i-1] == 0) {
      effpur_notnormErr[i-1] = 0;
    } else {
      effpur_notnormErr[i-1] = effpur_notnorm[i-1] * TMath::Sqrt(((efficencyErr_notnorm[i-1] / efficiency_notnorm[i-1]) * (efficencyErr_notnorm[i-1] / efficiency_notnorm[i-1])) + (purity_notnomrErr[i-1] / purity_notnomr[i-1]) * (purity_notnomrErr[i-1] / purity_notnomr[i-1]));
    }

    sigbk_notnorm[i-1] = efficiency_notnorm[i-1] * bkRej_notnorm[i-1];
    if (efficiency_notnorm[i-1] == 0 || bkRej_notnorm[i-1] == 0) {
      sigbk_notnormErr[i-1] = 0;
    } else {
      sigbk_notnormErr[i-1] = sigbk_notnorm[i-1] * TMath::Sqrt(((efficencyErr_notnorm[i-1] / efficiency_notnorm[i-1]) * (efficencyErr_notnorm[i-1] / efficiency_notnorm[i-1])) + (bkRej_notnormErr[i-1] / bkRej_notnorm[i-1]) * (bkRej_notnormErr[i-1] / bkRej_notnorm[i-1]));
    }

    xval_notnorm[i-1] = signalhist->GetBinCenter(i);
    xval_notnormErr[i-1] = signalhist->GetBinWidth(i)/2;

    //    std::cout << "xval_notnorm " << xval_notnorm[i-1] << " efficiency_notnorm: " << efficiency_notnorm[i-1] << " purity_notnomr: " << purity_notnomr[i-1] << " eff * pur "
    //          << effpur_notnorm[i-1] << "background Rejection: " << bkRej_notnorm[i-1] << " and sig*bkrej: " << sigbk_notnorm[i-1] << std::endl;

    //Keep the largest value.
    if (effpur_notnorm[i-1] > maxeffpur_notnorm_notnorm) {
      max_ind_notnorm = i-1;
      maxeffpur_notnorm_notnorm = effpur_notnorm[i-1];
    }
    if (sigbk_notnorm[i-1] > maxsigbk_notnorm_notnorm) {
      maxsigbk_notnorm_notnorm = sigbk_notnorm[i-1];
      maxsigbk_notnorm_ind_notnorm = i-1;
    }
  }


  TGraphErrors* efficiencyPlotNotNorm = new TGraphErrors(n, xval, efficiency_notnorm, xvalErr, efficencyErr_notnorm);
  efficiencyPlotNotNorm->GetYaxis()->SetTitle("Efficency");
  efficiencyPlotNotNorm->SetTitle("Efficency");
  efficiencyPlotNotNorm->SetLineColor(kBlue);
  efficiencyPlotNotNorm->Draw("a4");
  efficiencyPlotNotNorm->Write();

  TGraphErrors* bkRejPlotNotNorm = new TGraphErrors(n, xval, bkRej_notnorm, xvalErr, xval_notnormErr);
  bkRejPlotNotNorm->GetYaxis()->SetTitle("Background Rejection");
  bkRejPlotNotNorm->SetTitle("Background Rejection");
  bkRejPlotNotNorm->SetLineColor(kRed);
  bkRejPlotNotNorm->Draw("a4");
  bkRejPlotNotNorm->Write();

  std::cout << "Writing roccuve" << std::endl;
  TGraph* roccurve = new TGraph(n, efficiency_notnorm,bkRej_notnorm);
  roccurve->SetName("ROCCURVE");
  roccurve->GetYaxis()->SetTitle("Efficiency");
  roccurve->GetXaxis()->SetTitle("Background Rejection");
  roccurve->Write();



  //Draw the curves
  auto effpur_canvas = new TCanvas("effpur_canvas", "effpur_canvas", 600, 400);
  auto sigbk_canvas = new TCanvas("sigbk_canvas", "sigbk_canvas", 600, 400);

  TGraphErrors* efficiencyPlot = new TGraphErrors(n, xval, efficiency, xvalErr, efficencyErr);
  efficiencyPlot->GetYaxis()->SetTitle("Efficency");
  efficiencyPlot->SetTitle("Efficency");
  efficiencyPlot->SetLineColor(kOrange-1);
  efficiencyPlot->SetLineStyle(2);

  efficiencyPlot->Draw("a4");
  efficiencyPlot->Write();

  TGraphErrors* bkRejPlot = new TGraphErrors(n, xval, bkRej, xvalErr, bkRejErr);
  bkRejPlot->GetYaxis()->SetTitle("Background Rejection");
  bkRejPlot->SetTitle("Background Rejection");
  bkRejPlot->SetLineColor(kCyan-4);
  bkRejPlot->SetLineStyle(2);
  bkRejPlot->Draw("a4");
  bkRejPlot->Write();

  TGraphErrors* purityPlot = new TGraphErrors(n, xval, purity, xvalErr, purityErr);
  purityPlot->GetYaxis()->SetTitle("Purity");
  purityPlot->SetTitle("Purity");
  purityPlot->SetLineColor(kRed);
  purityPlot->Draw("a4");
  purityPlot->Write();

  TGraphErrors* effxpurPlot = new TGraphErrors(n, xval, effpur, xvalErr, effpurErr);
  effxpurPlot->GetYaxis()->SetTitle("Effxpur");
  effxpurPlot->SetTitle("Eff*pur");
  effxpurPlot->SetLineColor(kBlack);
  effxpurPlot->Draw("a4");
  effxpurPlot->Write();

  TGraphErrors* sigbkPlot = new TGraphErrors(n, xval, sigbk, xvalErr, sigbkErr);
  sigbkPlot->GetYaxis()->SetTitle("Signal*Background Rejection");
  sigbkPlot->SetTitle("Signal*Background Rejection");
  sigbkPlot->SetLineColor(kBlack);
  sigbkPlot->Draw("a4");
  sigbkPlot->Write();


  //Make the graph to present
  effpur_canvas->cd();
  auto mg1 = new TMultiGraph("mg1", "mg1");
  mg1->Add(efficiencyPlot);
  mg1->Add(purityPlot);
  mg1->Add(effxpurPlot);
  mg1->GetYaxis()->SetRangeUser(0, 1);
  mg1->Write();
  mg1->Draw("AP");
  effpur_canvas->BuildLegend();
  effpur_canvas->Write();

  purityPlot->SetLineColor(kMagenta-1);
  effxpurPlot->SetLineStyle(2);

  //Make the graph to present
  sigbk_canvas->cd();
  auto mg2 = new TMultiGraph("mg2", "mg2");
  mg2->Add(efficiencyPlot);
  mg2->Add(bkRejPlot);
  mg2->Add(sigbkPlot);
  mg2->Add(purityPlot);
  mg2->Add(effxpurPlot);
  mg2->Add(efficiencyPlotNotNorm);
  mg2->Add(bkRejPlotNotNorm);


  mg2->GetYaxis()->SetRangeUser(0, 1);
  mg2->Write();

  mg2->Draw("AP");
  sigbk_canvas->BuildLegend(0.6, 0.4, 0.9, 0.6);
  sigbk_canvas->Write();

  double maxScale = TMath::Max(backgroundhist->GetBinContent(backgroundhist->GetMaximumBin()),
			       signalhist->GetBinContent(signalhist->GetMaximumBin()));

  signalhist->SetLineColor(kBlue);
  signalhist->SetFillStyle(3003);
  signalhist->SetFillColor(6);
  signalhist->Scale(1. / maxScale);
  signalhist->SetTitle("Signal");
  signalhist->Write();
  backgroundhist->SetLineColor(kRed);
  backgroundhist->SetFillStyle(3003);
  backgroundhist->SetFillColor(42);
  backgroundhist->Scale(1. / maxScale);
  backgroundhist->SetTitle("Background");
  backgroundhist->Write();


  auto sigbk_cut_canvas = new TCanvas("sigbk_cut_canvas", "sigbk_cut_canvas", 600, 400);
  sigbk_cut_canvas->cd();
  mg2->GetHistogram()->SetTitle(Form("Cut at %0.1f " + units + " with efficiency %0.2f and background rejection %0.2f", xval[maxsigbk_ind], efficiency[maxsigbk_ind], bkRej[maxsigbk_ind]));
  mg2->GetXaxis()->SetTitle(axisTitle);
  mg2->Draw("AP");
  signalhist->Draw("HIST SAME e");
  backgroundhist->Draw("HIST SAME e");
  sigbk_cut_canvas->BuildLegend(0.6, 0.4, 0.9, 0.6);
  TLine* l1 = new TLine(xval[maxsigbk_ind], 0, xval[maxsigbk_ind], 1);
  l1->SetLineColor(kBlack);
  l1->Draw("same");

  sigbk_cut_canvas->Write();

  //Draw the overlayed samples with the cut.
  auto cut_canvas = new TCanvas("cut_canvas", "cut_canvas", 600, 400);
  cut_canvas->cd();
  mg1->GetHistogram()->SetTitle(Form("Cut at %0.1f " + units + " with efficiency %0.2f and purity %0.2f", xval[max_ind], efficiency[max_ind], purity[max_ind]));
  mg1->GetXaxis()->SetTitle(axisTitle);
  mg1->Draw("AP");
  signalhist->Draw("HIST SAME");
  backgroundhist->Draw("HIST SAME");

  cut_canvas->Update();
  cut_canvas->BuildLegend();

  TLine* l = new TLine(xval[max_ind], 0, xval[max_ind], 1);
  l->SetLineColor(kBlack);
  l->Draw("same");

  cut_canvas->Write();

  std::cout << "Best Cut is at: " << xval[max_ind] << std::endl;
  std::cout << "Efficiency: " << efficiency[max_ind] << " +- " << efficencyErr[max_ind] << std::endl;
  std::cout << "Purity    : " << purity[max_ind] << " +- " << purityErr[max_ind] << std::endl;
  std::cout << "Eff*Purity: " << effpur[max_ind] << " +- " << effpurErr[max_ind] << std::endl;

  std::cout << "Best sig*bkrej Cut is at: " << xval[maxsigbk_ind] << std::endl;
  std::cout << "Efficiency:     " << efficiency[maxsigbk_ind] << " +- " << efficencyErr[maxsigbk_ind] << std::endl;
  std::cout << "Background Rej: " << bkRej[maxsigbk_ind] << " +- " << bkRejErr[maxsigbk_ind] << std::endl;
  std::cout << "Eff*Purity:     " << sigbk[maxsigbk_ind] << " +- " << sigbkErr[maxsigbk_ind] << std::endl;

  delete efficiencyPlotNotNorm;
  delete bkRejPlotNotNorm;
  delete efficiencyPlot;
  delete bkRejPlot;
  delete purityPlot;
  delete effxpurPlot;
  delete sigbkPlot;
  delete roccurve;

  delete mg1;
  delete mg2;

  delete effpur_canvas;
  delete sigbk_canvas;
  delete sigbk_cut_canvas;
  delete cut_canvas;
}



//This analysis the histograms are already integrated so no need for that. They also start off 2D.
void optimiser::MetricHolder::MakeStandardNumShowerAnalysisGraphs(float totalsig, float totalbk){

  //  SignalHistPartner->Scale(1/TotalSigPOT);
  //  BackgroundHistPartner->Scale(1/TotalBKPOT);

  TH2D* EfficiencyHist     = (TH2D*) SignalHistPartner->Clone();
  TH2D* PurityHist         = (TH2D*) SignalHistPartner->Clone(); 
  TH2D* BackgroundRejHist  = (TH2D*) BackgroundHistPartner->Clone();
  TH2D* EffBkRejHist       = (TH2D*) SignalHistPartner->Clone();
  TH2D* EffPurHist         = (TH2D*) SignalHistPartner->Clone();

  EfficiencyHist->SetName("Efficiency");
  PurityHist->SetName("Purity");
  BackgroundRejHist->SetName("Background Rejection");
  EffBkRejHist->SetName("Efficiency*Background Rejection");
  EffPurHist->SetName("Efficiency*Purity");
  EfficiencyHist->SetTitle("Efficiency");
  PurityHist->SetTitle("Purity");
  BackgroundRejHist->SetTitle("Background Rejection");
  EffBkRejHist->SetTitle("Efficiency*Background Rejection");
  EffPurHist->SetTitle("Efficiency*Purity");


  int n_x = SignalHistPartner->GetNbinsX();
  int n_y = SignalHistPartner->GetNbinsY();

  //  double totalsig = SignalHistPartner->GetEntries()/(SignalHistPartner->GetNbinsX()*SignalHistPartner->GetNbinsY());
  // double totalbk  = BackgroundHistPartner->GetEntries()/(SignalHistPartner->GetNbinsX()*SignalHistPartner->GetNbinsY());

  double max_effbk = -999;
  int max_effbk_i  = -999;
  int max_effbk_j  = -999;

  double max_effpur = -999;
  int max_effpur_i  = -999;
  int max_effpur_j  = -999;


  for(int i=1; i<n_x+1; ++i) {
    for(int j=1; j<n_y+1; ++j){ 

      int bin = EfficiencyHist->GetBin(i,j);

      EfficiencyHist->SetBinContent(i,j,SignalHistPartner->GetBinContent(i,j)/totalsig);
      EfficiencyHist->SetBinError(bin,EfficiencyHist->GetBinError(bin)/totalsig);

      if(SignalHistPartner->GetBinContent(i,j)+ BackgroundHistPartner->GetBinContent(i,j) > 0){
	PurityHist->SetBinContent(i,j,SignalHistPartner->GetBinContent(i,j)/(SignalHistPartner->GetBinContent(i,j)+ BackgroundHistPartner->GetBinContent(i,j)));
	double PurityDemErr = TMath::Sqrt(SignalHistPartner->GetBinError(bin)*SignalHistPartner->GetBinError(bin) + BackgroundHistPartner->GetBinError(bin)*BackgroundHistPartner->GetBinError(bin));
	double PurityError = (SignalHistPartner->GetBinContent(i,j)/(SignalHistPartner->GetBinContent(i,j)+ BackgroundHistPartner->GetBinContent(i,j)))*TMath::Sqrt((SignalHistPartner->GetBinError(bin)/SignalHistPartner->GetBinContent(i,j))*(SignalHistPartner->GetBinError(bin)/SignalHistPartner->GetBinContent(i,j)) + (PurityDemErr/(SignalHistPartner->GetBinContent(i,j)+ BackgroundHistPartner->GetBinContent(i,j)))*(PurityDemErr/(SignalHistPartner->GetBinContent(i,j)+ BackgroundHistPartner->GetBinContent(i,j))));
	PurityHist->SetBinError(bin,PurityError); 
      }
      else{
	PurityHist->SetBinContent(i,j,0);
	PurityHist->SetBinError(bin,0);
      }

      BackgroundRejHist->SetBinContent(i,j,1-(BackgroundHistPartner->GetBinContent(i,j)/totalbk));
      BackgroundRejHist->SetBinError(bin,BackgroundRejHist->GetBinError(bin)/totalbk);

      EffBkRejHist->SetBinContent(i,j,EfficiencyHist->GetBinContent(i,j)*BackgroundRejHist->GetBinContent(i,j));
      double EffBkRejErr = EfficiencyHist->GetBinContent(i,j)*BackgroundRejHist->GetBinContent(i,j)*TMath::Sqrt((EfficiencyHist->GetBinError(bin)/EfficiencyHist->GetBinContent(i,j))*(EfficiencyHist->GetBinError(bin)/EfficiencyHist->GetBinContent(i,j)) + (BackgroundRejHist->GetBinError(bin)/BackgroundRejHist->GetBinContent(i,j))*(BackgroundRejHist->GetBinError(bin)/BackgroundRejHist->GetBinContent(i,j)));
      EffBkRejHist->SetBinError(bin,EffBkRejErr);

      EffPurHist->SetBinContent(i,j,EfficiencyHist->GetBinContent(i,j)*PurityHist->GetBinContent(i,j));
      double EffPurErr = EfficiencyHist->GetBinContent(i,j)*PurityHist->GetBinContent(i,j)*TMath::Sqrt((EfficiencyHist->GetBinError(bin)/EfficiencyHist->GetBinContent(i,j))*(EfficiencyHist->GetBinError(bin)/EfficiencyHist->GetBinContent(i,j)) + (PurityHist->GetBinError(bin)/PurityHist->GetBinContent(i,j))*(PurityHist->GetBinError(bin)/PurityHist->GetBinContent(i,j)));
      EffPurHist->SetBinError(bin,EffPurErr);

      //Store the max values 
      if(max_effbk < EffBkRejHist->GetBinContent(i,j)){
	max_effbk = EffBkRejHist->GetBinContent(i,j);
	max_effbk_i = i;
	max_effbk_j = j;
      } 
      if(max_effpur < EffPurHist->GetBinContent(i,j)){
	max_effpur = EffPurHist->GetBinContent(i,j);
	max_effpur_i = i;
	max_effpur_j = j;
      } 
    }
  }
  
  auto efficency_canvas = new TCanvas("eff2D_canvas", "eff2D_canvas", 600, 400);
  EfficiencyHist->Draw("colz");
  efficency_canvas->Write();

  auto background_canvas = new TCanvas("bk2D_canvas", "bk2D_canvas", 600, 400);
  BackgroundRejHist->Draw("colz");
  background_canvas->Write();

  auto purity_canvas = new TCanvas("pur2D_canvas", "pur2D_canvas", 600, 400);
  PurityHist->Draw("colz");
  purity_canvas->Write();

  auto effbk_canvas = new TCanvas("effbk2D_canvas", "effbk2D_canvas", 600, 400);
  EffBkRejHist->Draw("colz");
  effbk_canvas->Write();

  auto effpur_canvas = new TCanvas("effpur2D_canvas", "effpur2D_canvas", 600, 400);
  EffPurHist->Draw("colz");
  effpur_canvas->Write();

  //Make the 1D Graphs
  TH1D * Efficiency1D_sigbk    =  EfficiencyHist->ProjectionX("_pxsig",max_effbk_j,max_effbk_j,"e");
  Efficiency1D_sigbk->SetMarkerColor(kBlue);
  Efficiency1D_sigbk->SetLineColor(kBlue);
  Efficiency1D_sigbk->GetYaxis()->SetRangeUser(0, 1);
  
  TH1D * BackgroundRej1D_sigbk =  BackgroundRejHist->ProjectionX("_pxbk",max_effbk_j,max_effbk_j,"e");
  BackgroundRej1D_sigbk->SetMarkerColor(kRed);
  BackgroundRej1D_sigbk->SetLineColor(kRed);

  TH1D * EffBkRejHist1D_sigbk  =  EffBkRejHist->ProjectionX("_pxeffbk",max_effbk_j,max_effbk_j,"e");
  EffBkRejHist1D_sigbk->SetMarkerColor(kBlack);
  EffBkRejHist1D_sigbk->SetLineColor(kBlack);

  TH1D * PurityHist1D_sigbk  =  PurityHist->ProjectionX("_pxeffpure",max_effbk_j,max_effbk_j,"e");
  PurityHist1D_sigbk->SetMarkerColor(kRed);
  PurityHist1D_sigbk->SetLineColor(kRed);

  TH1D * EffPurHist1D_sigbk  =  EffPurHist->ProjectionX("_pxeffpur",max_effbk_j,max_effbk_j,"e");
  EffPurHist1D_sigbk->SetMarkerColor(kBlack);
  EffPurHist1D_sigbk->SetLineColor(kBlack);

  Efficiency1D_sigbk->SetTitle("Efficiency");
  BackgroundRej1D_sigbk->SetTitle("Background Rejection");
  EffBkRejHist1D_sigbk->SetTitle("Efficiency x Background Rejection");
  PurityHist1D_sigbk->SetTitle("Purity");
  EffPurHist1D_sigbk->SetTitle("Efficiency x Purity");
  
 
  auto effbk1D_canvas = new TCanvas("effbk1D_canvas", "effbk1D_canvas", 600, 400);
  Efficiency1D_sigbk->Draw("P");
  BackgroundRej1D_sigbk->Draw("P SAME");
  EffBkRejHist1D_sigbk->Draw("P SAME");
  effbk1D_canvas->BuildLegend();
  //  TLine* l1 = new TLine(((TAxis*)SignalHistPartner->GetYaxis())->GetBinCenter(max_effbk_j), 0, ((TAxis*)SignalHistPartner->GetYaxis())->GetBinCenter(max_effbk_j), 1);
  double max_effbk_1d = EffBkRejHist1D_sigbk->GetMaximumBin();
  //  TLine* l1 = new TLine(((TAxis*)SignalHistPartner->GetYaxis())->GetBinCenter(max_effbk_1d), 0, ((TAxis*)SignalHistPartner->GetYaxis())->GetBinCenter(max_effbk_1d), 1); 
  TLine* l1 = new TLine(EffBkRejHist1D_sigbk->GetBinCenter(max_effbk_1d),0,EffBkRejHist1D_sigbk->GetBinCenter(max_effbk_1d),1);

  l1->SetLineColor(kBlack);
  l1->Draw();
  effbk1D_canvas->Write();


  double max_effpur_1d = EffPurHist1D_sigbk->GetMaximumBin();


  auto full_canvas = new TCanvas("full_canvas", "full_canvas", 600, 400);
  full_canvas->Draw();
  Float_t sf = 0.3;
  TPad *p2 = new TPad("p2","p2",0,0,1,sf);
  TPad *p1 = new TPad("p1","p1",0,sf,1,1);
  p2->SetTopMargin(0.02);
  p1->SetBottomMargin(0.03);
  p2->SetBottomMargin(0.3);

  p1->cd();
  Efficiency1D_sigbk->GetXaxis()->SetLabelOffset(999);
  Efficiency1D_sigbk->GetXaxis()->SetLabelSize(0);
  Efficiency1D_sigbk->GetXaxis()->SetRangeUser(0,800);

  Double_t upTM = p1->GetTopMargin();
  Double_t upBM = p1->GetBottomMargin();
  Double_t lowTM = p2->GetTopMargin();
  Double_t lowBM = p2->GetBottomMargin();
  Double_t ratio = ( (upBM-(1-upTM))*(1-sf) ) / ( (lowBM-(1-lowTM))*sf ) ;
  float textsize = Efficiency1D_sigbk->GetYaxis()->GetLabelSize()*ratio; 
  Efficiency1D_sigbk->Draw("P");
  BackgroundRej1D_sigbk->Draw("P SAME");
  EffBkRejHist1D_sigbk->Draw("P SAME");
  p1->BuildLegend();
  l1->Draw("same");
  gPad->SetTickx(2);
  full_canvas->cd();
  p1->Draw();

  double maxpurity = PurityHist1D_sigbk->GetBinContent(PurityHist1D_sigbk->GetMaximumBin());
  p2->cd();
  PurityHist1D_sigbk->GetYaxis()->SetLabelSize(0.08);
  PurityHist1D_sigbk->GetXaxis()->SetLabelSize(textsize);
  PurityHist1D_sigbk->GetXaxis()->SetTitleSize(textsize);
  PurityHist1D_sigbk->GetYaxis()->SetRangeUser(0,maxpurity+0.001);
  PurityHist1D_sigbk->GetXaxis()->SetRangeUser(0,800);

  PurityHist1D_sigbk->Draw("P");
  EffPurHist1D_sigbk->Draw("P SAME");
  p2->BuildLegend();

  TLine* lfull_2 = new TLine(EffPurHist1D_sigbk->GetBinCenter(max_effpur_1d), 0, EffPurHist1D_sigbk->GetBinCenter(max_effpur_1d), maxpurity+0.001);
  lfull_2->SetLineColor(kBlack);
  lfull_2->Draw("same");
  full_canvas->cd();
  p2->Draw("same");

  
 
  full_canvas->Write();

  delete p2;
  delete p1;
  delete full_canvas;
  delete lfull_2;

  

 

  std::cout << "Best Cut is at Energy: " << ((TAxis*)SignalHistPartner->GetXaxis())->GetBinCenter(max_effpur_i) << " Num Showers: " << ((TAxis*)SignalHistPartner->GetYaxis())->GetBinCenter(max_effpur_j) << std::endl;
  std::cout << "Efficiency: " << EfficiencyHist->GetBinContent(max_effpur_i,max_effpur_j) << " - " << EfficiencyHist->GetBinErrorLow(max_effpur_i,max_effpur_j) << " + " << EfficiencyHist->GetBinErrorUp(max_effpur_i,max_effpur_j)  << std::endl;
  std::cout << "Purity: " << PurityHist->GetBinContent(max_effpur_i,max_effpur_j) << " - " << PurityHist->GetBinErrorLow(max_effpur_i,max_effpur_j) << " + " << PurityHist->GetBinErrorUp(max_effpur_i,max_effpur_j)  << std::endl;
  std::cout << "Eff*Purity: " << EffPurHist->GetBinContent(max_effpur_i,max_effpur_j) << " - " << EffPurHist->GetBinErrorLow(max_effpur_i,max_effpur_j) << " + " << EffPurHist->GetBinErrorUp(max_effpur_i,max_effpur_j)  << std::endl;

  std::cout << "Best Cut Background Rej is at Energy: " << ((TAxis*)SignalHistPartner->GetXaxis())->GetBinCenter(max_effbk_i) << " Num Showers: " << ((TAxis*)SignalHistPartner->GetYaxis())->GetBinCenter(max_effbk_j) << std::endl;
  std::cout << "Efficiency: " << EfficiencyHist->GetBinContent(max_effbk_i,max_effbk_j) << " - " << EfficiencyHist->GetBinErrorLow(max_effbk_i,max_effbk_j) << " + " << EfficiencyHist->GetBinErrorUp(max_effbk_i,max_effbk_j)  << std::endl;
  std::cout << "Background Rejection: " << BackgroundRejHist->GetBinContent(max_effbk_i,max_effbk_j) << " - " << BackgroundRejHist->GetBinErrorLow(max_effbk_i,max_effbk_j) << " + " << BackgroundRejHist->GetBinErrorUp(max_effbk_i,max_effbk_j)  << std::endl;
  std::cout << "Eff*Purity: " << EffBkRejHist->GetBinContent(max_effbk_i,max_effbk_j) << " - " << EffBkRejHist->GetBinErrorLow(max_effbk_i,max_effbk_j) << " + " << EffBkRejHist->GetBinErrorUp(max_effbk_i,max_effbk_j)  << std::endl;

  delete efficency_canvas;
  delete background_canvas;
  delete purity_canvas;
  delete effbk_canvas;
  delete effpur_canvas;

  return;
}

void optimiser::MetricHolder::Plot1D(){
  Plot1D(SignalHist,BackgroundHist);
}


void optimiser::MetricHolder::Plot1D(TH1D* signalhist, TH1D* backgroundhist){

  double maxScale = TMath::Max(backgroundhist->GetBinContent(backgroundhist->GetMaximumBin()),
			       signalhist->GetBinContent(signalhist->GetMaximumBin()));

  signalhist->SetLineColor(kBlue);
  signalhist->SetFillStyle(3003);
  signalhist->SetFillColor(6);
  signalhist->Scale(1. / maxScale);
  signalhist->SetTitle("Signal");
  signalhist->Write();
  backgroundhist->SetLineColor(kRed);
  backgroundhist->SetFillStyle(3003);
  backgroundhist->SetFillColor(42);
  backgroundhist->Scale(1. / maxScale);
  backgroundhist->SetTitle("Background");
  backgroundhist->Write();
  
  auto single_canvas = new TCanvas("single_canvas", "single_canvas", 600, 400);
  single_canvas->cd();
  signalhist->Draw();
  backgroundhist->Draw("SAME");
  single_canvas->BuildLegend(0.6, 0.4, 0.9, 0.6);
  single_canvas->Write();
}

void optimiser::MetricHolder::Plot2D(){
  Plot2D(SignalHistPartner,BackgroundHistPartner);
}

void optimiser::MetricHolder::Plot2D(TH2D* signalhist, TH2D* backgroundhist){

  auto signal_canvas = new TCanvas("signal_canvas", "signal_canvas", 600, 400);
  signalhist->Draw("colz");
  signal_canvas->Write();

  auto background_canvas = new TCanvas("background_canvas", "background_canvas", 600, 400);
  backgroundhist->Draw("colz");
  background_canvas->Write();
  

}

void optimiser::MetricHolder::MakeEfficiencyPlots(TH1D* postcut, TH1D* precut, TH1D* extra){

  TString eff_name = "eff_name";
  TString extraname = "extra_name";

  std::cout << "PostCut Size: " << postcut->GetEntries() << " PreCut Size; " << precut->GetEntries() << " Extra: " << extra->GetEntries() << std::endl;

  float maxSignal = precut->GetBinContent(precut->GetMaximumBin());
  
  if(TEfficiency::CheckConsistency(*postcut,*precut)){
    //Make the eff plot.
    TEfficiency* reco_eff = new TEfficiency(*postcut,*precut);
    //    reco_eff->SetTitle(title);
    //reco_eff->SetName(eff_name);
    reco_eff->SetTitle("title"); 
    reco_eff->SetName("eff_name"); 
    reco_eff->Write();

    //    postcut->Scale(1./maxSignal);
    //    precut->Scale(1./maxSignal);
      
    postcut->SetLineColor(kBlue);
    postcut->SetFillStyle(3003);
    postcut->SetFillColor(6);
    //    postcut->SetTitle(postcutname);
    postcut->SetTitle("postcutname"); 
    postcut->Write();
    precut->SetLineColor(kRed);
    precut->SetFillStyle(3003);
    precut->SetFillColor(42);
    precut->SetTitle("precutname");
    //    precut->SetTitle(precutname);
    precut->Write();

    TString canavs_name = eff_name + "_canvas";
    auto cut_canvas = new TCanvas(canavs_name, canavs_name, 600, 400);
    precut->Draw("HIST SAME");
    cut_canvas->cd();
    reco_eff->Draw("same");
    postcut->Draw("HIST SAME");
    cut_canvas->BuildLegend(0.6, 0.4, 0.9, 0.6);
    cut_canvas->Write();
        
    //    postcut->Scale(maxSignal);
    //    precut->Scale(maxSignal);

    std::cout << "Total Efficiency: " << (float) postcut->GetEntries() / (float) precut->GetEntries() << std::endl; 
    
  }

  if(TEfficiency::CheckConsistency(*extra,*precut)){
    //Make the eff plot.
    TEfficiency* reco_eff = new TEfficiency(*extra,*precut);
    //    reco_eff->SetTitle(title);
    //reco_eff->SetName(eff_name);
    reco_eff->SetTitle("reco eff");
    reco_eff->SetName("reco eff");
    
    reco_eff->Write();

    //    extra->Scale(1./maxSignal);
    //    precut->Scale(1./maxSignal);
      
    extra->SetLineColor(kBlue);
    extra->SetFillStyle(3003);
    extra->SetFillColor(6);
    extra->SetTitle(extraname);
    extra->Write();
    precut->SetLineColor(kRed);
    precut->SetFillStyle(3003);
    precut->SetFillColor(42);
    //    precut->SetTitle(precutname);
    precut->SetTitle("precutname");
    precut->Write();

    TString canavs_name = eff_name + "_canvas_extra";
    auto cut_canvas = new TCanvas(canavs_name, canavs_name, 600, 400);
    precut->Draw("HIST SAME");
    cut_canvas->cd();
    reco_eff->Draw("same");
    extra->Draw("HIST SAME");
    cut_canvas->BuildLegend(0.6, 0.4, 0.9, 0.6);
    cut_canvas->Write();
    
      
    //    extra->Scale(maxSignal);
    //  precut->Scale(maxSignal);
  }



}


void optimiser::MetricHolder::MakeEfficiencyPlots(TH1D* postcut_eff, TH1D* precut_eff, TH1D* extra_eff, TH1D* postcut_bk_en, TH1D* precut_bk, TH1D* extra_bk,
						  TH1D* postcut_eff_norm, TH1D* precut_eff_norm, TH1D* postcut_bk_en_norm, TH1D* precut_bk_norm){

  postcut_eff_norm->SaveAs("eff.C");
  postcut_bk_en_norm->SaveAs("bk.C");

  TString eff_name = "eff_name";
  TString extra_effname = "extra_eff_name";
  TString bk_name = "bk_name";
  TString extra_bkname = "extra_bk_name";
  
  TH1D* postcut_bk = (TH1D*) precut_bk->Clone();
  TH1D* After_bk  =  (TH1D*) postcut_bk_en->Clone();
  After_bk->Scale(-1);
  std::cout << "testone" << std::endl;
  postcut_bk->Add(After_bk);

  std::cout << "testtwo" << std::endl;

  TH1D* postcut_bk_norm = (TH1D*) precut_bk_norm->Clone();
  TH1D* After_bk_norm  =  (TH1D*) postcut_bk_en_norm->Clone();
  After_bk_norm->Scale(-1);
  postcut_bk_norm->Add(After_bk_norm);
  std::cout << "testthree" << std::endl;

  
  TH1D* purity = (TH1D*) postcut_bk_en->Clone();
  purity->Add(postcut_eff);
  std::cout << "testfour" << std::endl;


  TEfficiency* pruity_teff;
  if(TEfficiency::CheckConsistency(*postcut_eff,*purity)){
    //Make the bk plot.
    pruity_teff = new TEfficiency(*postcut_eff,*purity);
    //    pruity_teff->SetTitle(title);
    //pruity_teff->SetName(bk_name);
    pruity_teff->SetTitle("title"); 
    pruity_teff->SetName("pure_name"); 
    pruity_teff->Write();
    
  }

  std::cout << "entries: " << postcut_eff_norm->GetEntries() << " " << precut_eff_norm->GetEntries() << " " << postcut_bk_norm->GetEntries() << " " << precut_bk_norm->GetEntries() << std::endl;

  TEfficiency* reco_eff_norm;
  if(TEfficiency::CheckConsistency(*postcut_eff_norm,*precut_eff_norm)){
    reco_eff_norm = new TEfficiency(*postcut_eff_norm,*precut_eff_norm);
  }

  TEfficiency* reco_bk_norm;
  if(TEfficiency::CheckConsistency(*postcut_bk_norm,*precut_bk_norm)){
    reco_bk_norm = new TEfficiency(*postcut_bk_norm,*precut_bk_norm);
  }



              
  std::cout << "Postcut_Eff Size: " << postcut_eff->GetEntries() << " Precut_Eff Size; " << precut_eff->GetEntries() << " Extra_Eff: " << extra_eff->GetEntries() << std::endl;
  std::cout << "Postcut_Eff Size int: " << postcut_eff->Integral() << " Precut_Eff Size; " << precut_eff->Integral() << " Extra_Eff: " << extra_eff->Integral() << std::endl;


  float maxSignal = precut_eff->GetBinContent(precut_eff->GetMaximumBin());
  
  TEfficiency* reco_eff;
  if(TEfficiency::CheckConsistency(*postcut_eff,*precut_eff)){
    //Make the eff plot.
    reco_eff = new TEfficiency(*postcut_eff,*precut_eff);
    //    reco_eff->SetTitle(title);
    //reco_eff->SetName(eff_name);
    reco_eff->SetTitle("title"); 
    reco_eff->SetName("eff_name"); 
    reco_eff->Write();

    //    postcut_eff->Scale(1./maxSignal);
    //    precut_eff->Scale(1./maxSignal);
      
    postcut_eff->SetLineColor(kBlue);
    postcut_eff->SetFillStyle(3003);
    postcut_eff->SetFillColor(6);
    //    postcut_eff->SetTitle(postcut_effname);
    postcut_eff->SetTitle("postcut_effname"); 
    postcut_eff->Write();
    precut_eff->SetLineColor(kRed);
    precut_eff->SetFillStyle(3003);
    precut_eff->SetFillColor(42);
    precut_eff->SetTitle("precut_effname");
    //    precut_eff->SetTitle(precut_effname);
    precut_eff->Write();

    TString canavs_name = eff_name + "_canvas";
    auto cut_canvas = new TCanvas(canavs_name, canavs_name, 600, 400);
    precut_eff->Draw("HIST SAME");
    cut_canvas->cd();
    reco_eff->Draw("same");
    postcut_eff->Draw("HIST SAME");
    cut_canvas->BuildLegend(0.6, 0.4, 0.9, 0.6);
    cut_canvas->Write();
        
    //    postcut_eff->Scale(maxSignal);
    //    precut_eff->Scale(maxSignal);

    std::cout << "Total Efficiency: " << (float) postcut_eff->Integral() / (float) precut_eff->Integral() << std::endl; 
    
  }

  TEfficiency* reco_eff_extra;
  if(TEfficiency::CheckConsistency(*extra_eff,*precut_eff)){
    //Make the eff plot.
    reco_eff_extra = new TEfficiency(*extra_eff,*precut_eff);
    //    reco_eff_extra->SetTitle(title);
    //reco_eff_extra->SetName(eff_name);
    reco_eff_extra->SetTitle("reco eff");
    reco_eff_extra->SetName("reco eff");
    
    reco_eff_extra->Write();

    //    extra_eff->Scale(1./maxSignal);
    //    precut_eff->Scale(1./maxSignal);
      
    extra_eff->SetLineColor(kBlue);
    extra_eff->SetFillStyle(3003);
    extra_eff->SetFillColor(6);
    extra_eff->SetTitle(extra_effname);
    extra_eff->Write();
    precut_eff->SetLineColor(kRed);
    precut_eff->SetFillStyle(3003);
    precut_eff->SetFillColor(42);
    //    precut_eff->SetTitle(precut_effname);
    precut_eff->SetTitle("precut_effname");
    precut_eff->Write();

    TString canavs_name = eff_name + "_canvas_extra_eff";
    auto cut_canvas = new TCanvas(canavs_name, canavs_name, 600, 400);
    precut_eff->Draw("HIST SAME");
    cut_canvas->cd();
    reco_eff_extra->Draw("same");
    extra_eff->Draw("HIST SAME");
    cut_canvas->BuildLegend(0.6, 0.4, 0.9, 0.6);
    cut_canvas->Write();
    
      
    //    extra_eff->Scale(maxSignal);
    //    precut_eff->Scale(maxSignal);
  }


  std::cout << "Postcut_Eff Size: " << postcut_eff->GetEntries() << " Precut_Eff Size; " << precut_eff->GetEntries() << " Extra_Eff: " << extra_eff->GetEntries() << std::endl;
  
  maxSignal = precut_bk->GetBinContent(precut_bk->GetMaximumBin());

  TEfficiency* reco_bk; 
  if(TEfficiency::CheckConsistency(*postcut_bk,*precut_bk)){
    //Make the bk plot.
    reco_bk = new TEfficiency(*postcut_bk,*precut_bk);

    for(int i=0; i<postcut_bk->GetNbinsX(); ++i){
    }

    //    reco_bk->SetTitle(title);
    //reco_bk->SetName(bk_name);
    reco_bk->SetTitle("title"); 
    reco_bk->SetName("bk_name"); 
    reco_bk->Write();

    //    postcut_bk->Scale(1./maxSignal);
    //    precut_bk->Scale(1./maxSignal);
      
    postcut_bk->SetLineColor(kBlue);
    postcut_bk->SetFillStyle(3003);
    postcut_bk->SetFillColor(6);
    //    postcut_bk->SetTitle(postcut_bkname);
    postcut_bk->SetTitle("postcut_bkname"); 
    postcut_bk->Write();
    precut_bk->SetLineColor(kRed);
    precut_bk->SetFillStyle(3003);
    precut_bk->SetFillColor(42);
    precut_bk->SetTitle("precut_bkname");
    //    precut_bk->SetTitle(precut_bkname);
    precut_bk->Write();

    TString canavs_name = bk_name + "_canvas";
    auto cut_canvas = new TCanvas(canavs_name, canavs_name, 600, 400);
    precut_bk->Draw("HIST SAME");
    cut_canvas->cd();
    reco_bk->Draw("same");
    postcut_bk->Draw("HIST SAME");
    cut_canvas->BuildLegend(0.6, 0.4, 0.9, 0.6);
    cut_canvas->Write();
        
    //    postcut_bk->Scale(maxSignal);
    //    precut_bk->Scale(maxSignal);

    std::cout << "Total Bkiciency: " << (float) postcut_bk_en->Integral() / (float) precut_bk->Integral() << std::endl; 

    std::cout << "Total Purity: " << (float) postcut_eff->Integral() / ((float) postcut_bk_en->Integral() + (float) postcut_eff->Integral()) << std::endl; 

    
  }

  TEfficiency *reco_bk_extra;
  if(TEfficiency::CheckConsistency(*extra_bk,*precut_bk)){
    //Make the bk plot.
    reco_bk_extra = new TEfficiency(*extra_bk,*precut_bk);
    //    reco_bk_extra->SetTitle(title);
    //reco_bk_extra->SetName(bk_name);
    reco_bk_extra->SetTitle("reco bk");
    reco_bk_extra->SetName("reco bk");
    
    reco_bk_extra->Write();

    //    extra_bk->Scale(1./maxSignal);
    //    precut_bk->Scale(1./maxSignal);
      
    extra_bk->SetLineColor(kBlue);
    extra_bk->SetFillStyle(3003);
    extra_bk->SetFillColor(6);
    extra_bk->SetTitle(extra_bkname);
    extra_bk->Write();
    precut_bk->SetLineColor(kRed);
    precut_bk->SetFillStyle(3003);
    precut_bk->SetFillColor(42);
    //    precut_bk->SetTitle(precut_bkname);
    precut_bk->SetTitle("precut_bkname");
    precut_bk->Write();

    TString canavs_name = bk_name + "_canvas_extra_bk";
    auto cut_canvas = new TCanvas(canavs_name, canavs_name, 600, 400);
    precut_bk->Draw("HIST SAME");
    cut_canvas->cd();
    reco_bk_extra->Draw("same");
    extra_bk->Draw("HIST SAME");
    cut_canvas->BuildLegend(0.6, 0.4, 0.9, 0.6);
    cut_canvas->Write();
    
      
    //    extra_bk->Scale(maxSignal);
    //    precut_bk->Scale(maxSignal);
  }


  pruity_teff->SetFillColor(0);
  
  precut_bk->SetFillColor(kRed);
  precut_bk->SetLineColor(kRed);

  precut_eff->SetFillColor(kBlue);
  precut_eff->SetLineColor(kBlue);

  reco_eff->SetLineColor(kOrange-1);
  reco_bk->SetLineColor(kCyan-3);

  reco_bk->SetFillColor(0);
  reco_eff->SetFillColor(0);

  precut_bk->SetName("Preselection Background Distribution");
  precut_eff->SetName("Preselection Signal Distribution");

  reco_bk->SetName("Normalised Background Rejection");
  reco_eff->SetName("Normalised Efficiency");

  reco_bk->SetTitle("Normalised Background Rejection");
  reco_eff->SetTitle("Normalised Efficiency");


  pruity_teff->SetName("Purity");

  precut_bk->SetTitle("Preselection Background Distribution");
  precut_eff->SetTitle("Preselection Signal Distribution");

  //  reco_bk->SetTitle("Background Rejection");
  //  reco_eff->SetTitle("Efficiency");

  pruity_teff->SetTitle("Purity");

  reco_bk_norm->Write();
  reco_eff_norm->Write();


  reco_eff_norm->SetLineColor(kBlue);
  reco_bk_norm->SetLineColor(kRed);


  reco_eff->SetLineStyle(2);
  reco_bk->SetLineStyle(2);


  reco_eff_norm->SetName("Efficiency");
  reco_bk_norm->SetName("Background Rejection");

  reco_eff_norm->SetFillColor(0);
  reco_bk_norm->SetFillColor(0);


  precut_bk->Scale(1/precut_bk->Integral());
  precut_eff->Scale(1/precut_eff->Integral());
  
  maxSignal = TMath::Max(precut_bk->GetBinContent(precut_bk->GetMaximumBin()),precut_eff->GetBinContent(precut_eff->GetMaximumBin()));

  precut_bk->Scale(1/maxSignal);
  precut_eff->Scale(1/maxSignal);

  precut_bk->GetYaxis()->SetRangeUser(0,1.1);

  TString canavs_name = "effbk_canvas";
  auto cut_canvas = new TCanvas(canavs_name, canavs_name, 600, 400);
  cut_canvas->cd();
  precut_bk->Draw("HIST");
  precut_eff->Draw("HIST SAME");
  reco_bk->Draw("same");
  reco_eff->Draw("same");
  pruity_teff->Draw("same");
  reco_eff_norm->Draw("same");
  reco_bk_norm->Draw("same");

  cut_canvas->BuildLegend(0.6, 0.4, 0.9, 0.6);
  cut_canvas->Write();


  precut_bk->Scale(maxSignal);
  precut_eff->Scale(maxSignal);
 

  // TString canavs_name_extra = "effbk_canvas_extra";
  // auto cut_canvas_extra = new TCanvas(canavs_name_extra, canavs_name_extra, 600, 400);
  // cut_canvas_extra->cd();
  // precut_bk_extra->Draw("HIST");
  // precut_eff_extra->Draw("HIST SAME");
  // reco_bk_extra->Draw("same");
  // reco_eff_extra->Draw("same");
  // cut_canvas_extra->BuildLegend(0.6, 0.4, 0.9, 0.6);
  //  cut_canvas_extra->Write();

    


}



