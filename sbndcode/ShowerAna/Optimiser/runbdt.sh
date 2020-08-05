echo "hello"

source /grid/fermiapp/products/sbnd/setup_sbnd.sh; 
setup mrb;
source /sbnd/app/users/dbarker/larsoft_sbnd_shower/localProducts_larsoft_v08_40_00_e19_prof/setup;
mrbslp


make

#Tune the MVA 
#./MakeSelectionEfficiencyPlots -s ../Signal/showervalidationGraphs_test_train.root -b ../Background/showervalidationGraphs_test_train.root --MVA

#Run over the training 
./MakeSelectionEfficiencyPlots -s ../Signal/showervalidationGraphs_test_train.root -b ../Background/showervalidationGraphs_test_train.root
mv CutFile.root CutFileTrain.root

#Run over the validation 
./MakeSelectionEfficiencyPlots -s ../Signal/showervalidationGraphs_test_val.root -b ../Background/showervalidationGraphs_test_val.root
mv CutFile.root CutFileVal.root

#Run the KSScore
root -q ./KSTest.C
