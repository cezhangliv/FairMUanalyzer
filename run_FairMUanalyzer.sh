#./run_FairMUanalyzer [input.root] [output_prefix] [RunN] [Tgt] 
#./run_FairMUanalyzer root/event123.root custom/output/path/myresult [RunN] [Tgt] [mf_flag]
# output to custom/output/path/myresult.root


#./run_FairMUanalyzer root/single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1.root result/test 100 1

#./batch_run root/single_muon_interaction_1_CDbugfix11July25_run8_21July25
#./batch_run root/single_muon_interaction_0_CDbugfix11July25_run8_21July25

: '
./run_FairMUanalyzer root/single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1.root \
	result/FairMUanalyzer_18July25_single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1 -1 1

./run_FairMUanalyzer root/single_muon_interaction_0_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1.root \
	result/FairMUanalyzer_18July25_single_muon_interaction_0_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1 -1 0

./run_FairMUanalyzer root/single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1.root \
	result/FairMUanalyzer_18July25_single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF0 -1 1 0

./run_FairMUanalyzer root/single_muon_interaction_0_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1.root \
	result/FairMUanalyzer_18July25_single_muon_interaction_0_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF0 -1 0 0

./run_FairMUanalyzer root/single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1.root \
	result/FairMUanalyzer_18July25_single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1_NoTightTrackCutTgt2 -1 1
'

##### Milosz

: '
./run_FairMUanalyzer /afs/cern.ch/work/m/mzdybal/public/run8_reconstruction.root \
	result/FairMUanalyzer_06Aug25_Milosz_run8_reconstruction_1_CDbugfix11July25_MF1 -1 1

./run_FairMUanalyzer /afs/cern.ch/work/m/mzdybal/public/run8_reconstruction.root \
	result/FairMUanalyzer_06Aug25_Milosz_run8_reconstruction_0_CDbugfix11July25_MF1 -1 0

./run_FairMUanalyzer /afs/cern.ch/work/m/mzdybal/public/run17_interaction1_reconstruction.root \
	result/FairMUanalyzer_06Aug25_Milosz_run17_reconstruction_1_CDbugfix11July25_MF1 -1 1

./run_FairMUanalyzer /afs/cern.ch/work/m/mzdybal/public/run17_interaction0_reconstruction.root \
	result/FairMUanalyzer_06Aug25_Milosz_run17_reconstruction_0_CDbugfix11July25_MF1 -1 0
'

### 21 Aug 25 test new run24

#./run_FairMUanalyzer root/single_muon_interaction_0_CDbugfix11July25_run24_muedaq04-1753473533_MF1.root \
#	result/FairMUanalyzer_21Aug25_run24_reconstruction_0_CDbugfix11July25_MF1 -1 0
#
#./run_FairMUanalyzer root/single_muon_interaction_1_CDbugfix11July25_run24_muedaq04-1753473533_MF1.root \
#	result/FairMUanalyzer_21Aug25_run24_reconstruction_1_CDbugfix11July25_MF1 -1 1

#./batch_run root/single_muon_interaction_0_CDbugfix11July25_run24_muedaq04
#./batch_run root/single_muon_interaction_1_CDbugfix11July25_run24_muedaq04

### 29 Aug 25 GA discussion test, approved plots

./run_FairMUanalyzer100 root/single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1.root \
	result/FairMUanalyzer_29Aug25_GAtest100_single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1 -1 1


#./run_FairMUanalyzer101 root/single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1.root result/FairMUanalyzer_29Aug25_GAtest101_single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1 -1 1
#./run_FairMUanalyzer201 root/single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1.root result/FairMUanalyzer_29Aug25_GAtest201_single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1 -1 1