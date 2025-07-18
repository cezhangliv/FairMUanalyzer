#./run_FairMUanalyzer [input.root] [output_prefix] [RunN] [Tgt] 
#./run_FairMUanalyzer root/event123.root custom/output/path/myresult [RunN] [Tgt] 
# output to custom/output/path/myresult.root

./run_FairMUanalyzer root/single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1.root \
	result/FairMUanalyzer_ -1 1
