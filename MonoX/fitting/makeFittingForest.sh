echo monojet_zmm_0tag monojet_wmn_0tag monojet_wen_0tag monojet_tme_0tag monojet_tem_0tag monojet_zee_0tag | xargs -n 1 -P 8 python makeFittingForest.py --region
#echo monojet_zmm_1tag monojet_wmn_1tag monojet_wen_1tag monojet_tme_1tag monojet_tem_1tag monojet_zee_1tag monojet_signal_1tag | xargs -n 1 -P 8 python makeFittingForest.py --region
#echo monojet_zmm_2tag monojet_wmn_2tag monojet_wen_2tag monojet_tme_2tag monojet_tem_2tag monojet_zee_2tag monojet_signal_2tag | xargs -n 1 -P 8 python makeFittingForest.py --region
#echo zmm wmn wen tm te zee signal | xargs -n 1 -P 8 python makeFittingForest.py --ddt True --region 
#echo zmm_fail wmn_fail wen_fail tm_fail te_fail zee_fail signal_fail | xargs -n 1 -P 8 python makeFittingForest.py --ddt True --region
#echo signal signal_fail | xargs -n 1 -P 8 python makeFittingForest.py --ddt True --region

#hadd /uscms_data/d3/jmartine/panda/v_8026_0_4/limits/fitting_mass3/fittingForest_mass3.root /uscms_data/d3/jmartine/panda/v_8026_0_4/limits/fitting_mass3/*.root
#echo wmn_fail wen_fail zmm_fail zee_fail pho_fail te_fail tm_fail | xargs -n 1 -P 8 python makeFittingForest.py --input /data/t3home000/mcremone/lpc/jorgem/panda/v_8026_0_4/flat/control/ --output /data/t3home000/mcremone/lpc/jorgem/panda/v_8026_0_4/flat/control/ --region

#echo signal | xargs -n 1 -P 8 python makeFittingForest.py --input /data/t3home000/mcremone/lpc/jorgem/panda/v_8026_0_4/flat/signal/ --output /data/t3home000/mcremone/lpc/jorgem/panda/v_8026_0_4/flat/signal/ --var True --region
