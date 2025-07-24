#./run_FairMUanalyzer [input.root] [output_prefix] [RunN] [Tgt] 
#./run_FairMUanalyzer root/event123.root custom/output/path/myresult [RunN] [Tgt] [mf_flag]
# output to custom/output/path/myresult.root


#./run_FairMUanalyzer root/single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1.root result/test 100 1

#./batch_run root/single_muon_interaction_1_CDbugfix11July25_run8_21July25
#./batch_run root/single_muon_interaction_0_CDbugfix11July25_run8_21July25

./run_FairMUanalyzer root/single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1.root \
	result/FairMUanalyzer_18July25_single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1 -1 1

./run_FairMUanalyzer root/single_muon_interaction_0_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1.root \
	result/FairMUanalyzer_18July25_single_muon_interaction_0_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1 -1 0

./run_FairMUanalyzer root/single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1.root \
	result/FairMUanalyzer_18July25_single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF0 -1 1 0

./run_FairMUanalyzer root/single_muon_interaction_0_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1.root \
	result/FairMUanalyzer_18July25_single_muon_interaction_0_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF0 -1 0 0



