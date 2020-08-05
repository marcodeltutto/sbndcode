
KSTest(){

  std::ofstream myfile;
  myfile.open ("histcomp.txt");
  myfile << "0";
  myfile.close();
  
  //Get the histogram for the signal training
  auto trainingfile = new TFile("mcp2electrons/showervalidationGraphs_test.root");
  TDirectory* dir_train = signalfile->GetDirectory("ana");
  TH1D* signalhist_train = (TH1D*) dir_train->Get("MetricTree");     
  TH1D* backgroundhist_train   = (TH1D*) dir_train->Get("MetricTree");     

  //Get the histogram for the background
  auto validationfile = new TFile("mcp2photons/showervalidationGraphs_test.root");
  TDirectory* dir_val = backgroundfile->GetDirectory("ana");
  TH1D* signalhist_val = (TH1D*) dir_val->Get("MetricTree");     
  TH1D* backgroundhist_val = (TH1D*) dir_val->Get("MetricTree");     
    
  if (signalhist_train == NULL || backgroundhist_train == NULL || signalhist_val == NULL || backgroundhist_val == NULL) {
    std::cout << "Histogram is Null, Returning" << std::endl;
    return;
  }

  float KS_sig = signalhist_val.KolmogorovTest(signalhist_train,"X");
  float KS_bk = backgroundhist_val.KolmogorovTest(backgroundhist_train,"X");

  std::ofstream myfile;
  myfile.open ("histcomp.txt");
  myfile << KS_sig*KS_bk;
  myfile.close();

    
}
