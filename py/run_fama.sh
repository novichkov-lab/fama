#!/bin/sh
#python fama.py --config config.ini --project project.ini --sample test_sample --end pe1

# ENIGMA FW306 (prepilot) sediment, all samples together
#~ python3 fama.py -c config.ini -p project_FW306_nitrogen9.ini
#~ python3 fama.py -c config.ini -p project_FW306_universal1.ini

# ENIGMA isolates from water, Enterobacteria
#~ python3 fama.py -c config.ini -p project_water_isolates_nitrogen9.ini --prot
#~ python3 fama.py -c config.ini -p project_water_isolates_universal1.ini --prot

# ENIGMA pilot EB106
# python3 fama.py -c config.ini -p project_EB106-01_nitrogen9.ini
#python3 fama.py -c config.ini -p project_EB106-01_universal1.ini

#python3 fama.py -c config.ini -p project_EB106-02_nitrogen9.ini
#python3 fama.py -c config.ini -p project_EB106-02_universal1.ini

#python3 fama.py -c config.ini -p project_EB106-03_nitrogen9.ini
#python3 fama.py -c config.ini -p project_EB106-03_universal1.ini

#python3 fama.py -c config.ini -p project_EB106-04_nitrogen9.ini
#python3 fama.py -c config.ini -p project_EB106-04_universal1.ini

#python3 fama.py -c config.ini -p project_EB106-05_nitrogen9.ini
#python3 fama.py -c config.ini -p project_EB106-05_universal1.ini

#python3 fama.py -c config.ini -p project_EB106-06_20_nitrogen9.ini
#python3 fama.py -c config.ini -p project_EB106-06_20_universal1.ini

# ENIGMA pilot EB271
#~ python3 fama.py -c config.ini -p project_EB271-ZV-D103_universal1.ini
#~ python3 fama.py -c config.ini -p project_EB271-ZV-D103_nitrogen9.ini

#~ python3 fama.py -c config.ini -p project_EB271-ZV-D126_universal1.ini
#~ python3 fama.py -c config.ini -p project_EB271-ZV-D126_nitrogen9.ini

#~ python3 fama.py -c config.ini -p project_EB271-ZV-D142_universal1.ini
#~ python3 fama.py -c config.ini -p project_EB271-ZV-D142_nitrogen9.ini

#python3 fama.py -c config.ini -p project_EB271-ZV-D194_universal1.ini
#python3 fama.py -c config.ini -p project_EB271-ZV-D194_nitrogen9.ini

#~ python3 fama.py -c config.ini -p project_EB271-ZV-D217_universal1.ini
#~ python3 fama.py -c config.ini -p project_EB271-ZV-D217_nitrogen9.ini

#~ python3 fama.py -c config.ini -p project_EB271-ZV-D243_universal1.ini
#~ python3 fama.py -c config.ini -p project_EB271-ZV-D243_nitrogen9.ini

#~ python3 fama.py -c config.ini -p project_EB271-ZC-D286_universal1.ini
#~ python3 fama.py -c config.ini -p project_EB271-ZC-D286_nitrogen9.ini

#~ python3 fama.py -c config.ini -p project_EB271-ZC-D309_universal1.ini
#~ python3 fama.py -c config.ini -p project_EB271-ZC-D309_nitrogen9.ini

#~ python3 fama.py -c config.ini -p project_EB271-ZC-D331_universal1.ini
#~ python3 fama.py -c config.ini -p project_EB271-ZC-D331_nitrogen9.ini

#~ python3 fama.py -c config.ini -p project_EB271-ZC-D353_universal1.ini
#~ python3 fama.py -c config.ini -p project_EB271-ZC-D353_nitrogen9.ini

#~ python3 fama.py -c config.ini -p project_EB271-ZS-D377_universal1.ini
#~ python3 fama.py -c config.ini -p project_EB271-ZS-D377_nitrogen9.ini

#~ python3 fama.py -c config.ini -p project_EB271-ZS-D399_universal1.ini
#~ python3 fama.py -c config.ini -p project_EB271-ZS-D399_nitrogen9.ini

#~ python3 fama.py -c config.ini -p project_EB271-ZS-D421_universal1.ini
#~ python3 fama.py -c config.ini -p project_EB271-ZS-D421_nitrogen9.ini

# python3 fama.py -c config.ini -p project_EB271-ZS-D444_universal1.ini
# python3 fama.py -c config.ini -p project_EB271-ZS-D444_nitrogen9.ini

# Test groundwater metagenomes
# This project is for wrong FW301 fastq file. Use fw106 single project instead
#python3 fama_se.py --config config.ini --project project_fw106_fw301_nitrogen9.ini --sample fw106
#File deleted python3 fama_se.py --config config.ini --project project_fw106_fw301_nitrogen9.ini --sample fw301

#~ python3 fama.py -c config.ini -p project_fw301_nitrogen9.ini
#~ python3 fama.py -c config.ini -p project_fw301_universal1.ini

# Benchmark
#python3 fama.py -c config.ini -p project_benchmark_nitrogen9.ini
#python3 fama.py -c config.ini -p project_benchmark_universal1.ini

# 12wells groundwater samples
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample DP16D_1 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample DP16D_1 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample DP16D_2 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample DP16D_2 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample DP16D_3 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample DP16D_3 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample DP16D_4 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample DP16D_4 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample DP16D_5 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample DP16D_5 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-021_1 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-021_1 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-021_2 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-021_2 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-021_3 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-021_3 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-021_4 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-021_4 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-021_5 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-021_5 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-104_1 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-104_1 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-104_2 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-104_2 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-104_3 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-104_3 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-104_4 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-104_4 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-104_5 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-104_5 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-106_1 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-106_1 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-106_2 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-106_2 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-106_3 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-106_3 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-106_4 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-106_4 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-106_5 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-106_5 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-215_1 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-215_1 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-215_2 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-215_2 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-215_3 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-215_3 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-215_4 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-215_4 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-215_5 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-215_5 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-300_1 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-300_1 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-300_2 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-300_2 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-300_3 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-300_3 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-300_4 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-300_4 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-300_5 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-300_5 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-301_1 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-301_1 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-301_2 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-301_2 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-301_3 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-301_3 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-301_4 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-301_4 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-301_5 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-301_5 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-305_1 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-305_1 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-305_2 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-305_2 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-305_3 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-305_3 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-305_4 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-305_4 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-305_5 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-305_5 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-602-26_1 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-602-26_1 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-602-26_2 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-602-26_2 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-602-26_3 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-602-26_3 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-602-26_4 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-602-26_4 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-602-26_5 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample FW-602-26_5 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample GW-199_1 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample GW-199_1 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample GW-199_2 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample GW-199_2 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample GW-199_3 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample GW-199_3 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample GW-199_4 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample GW-199_4 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample GW-199_5 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample GW-199_5 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample GW-715_1 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample GW-715_1 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample GW-715_2 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample GW-715_2 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample GW-715_3 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample GW-715_3 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample GW-715_4 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample GW-715_4 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample GW-715_5 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample GW-715_5 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample GW928_1 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample GW928_1 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample GW928_2 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample GW928_2 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample GW928_3 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample GW928_3 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample GW928_4 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample GW928_4 --end pe2
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample GW928_5 --end pe1
#python fama.py --config config.ini --project project_12wells_nitrogen.ini --sample GW928_5 --end pe2

# Cellulose adapted consortia samples, passages B1F01-B1F14
#python fama.py --config config.ini --project project_B1F_cazy_t.ini --sample B1F01 --end pe1
#python fama.py --config config.ini --project project_B1F_cazy_t.ini --sample B1F01 --end pe2
#python fama.py --config config.ini --project project_B1F_cazy_t.ini --sample B1F02 --end pe1
#python fama.py --config config.ini --project project_B1F_cazy_t.ini --sample B1F02 --end pe2
#python fama.py --config config.ini --project project_B1F_cazy_t.ini --sample B1F03 --end pe1
#python fama.py --config config.ini --project project_B1F_cazy_t.ini --sample B1F03 --end pe2
#python fama.py --config config.ini --project project_B1F_cazy_t.ini --sample B1F04 --end pe1
#python fama.py --config config.ini --project project_B1F_cazy_t.ini --sample B1F04 --end pe2
#python fama.py --config config.ini --project project_B1F_cazy_t.ini --sample B1F05 --end pe1
#python fama.py --config config.ini --project project_B1F_cazy_t.ini --sample B1F05 --end pe2
#python fama.py --config config.ini --project project_B1F_cazy_t.ini --sample B1F06 --end pe1
#python fama.py --config config.ini --project project_B1F_cazy_t.ini --sample B1F06 --end pe2
#python fama.py --config config.ini --project project_B1F_cazy_t.ini --sample B1F07 --end pe1
#python fama.py --config config.ini --project project_B1F_cazy_t.ini --sample B1F07 --end pe2
#python fama.py --config config.ini --project project_B1F_cazy_t.ini --sample B1F08 --end pe1
#python fama.py --config config.ini --project project_B1F_cazy_t.ini --sample B1F08 --end pe2
#python fama.py --config config.ini --project project_B1F_cazy_t.ini --sample B1F09 --end pe1
#python fama.py --config config.ini --project project_B1F_cazy_t.ini --sample B1F09 --end pe2
#python fama.py --config config.ini --project project_B1F_cazy_t.ini --sample B1F10 --end pe1
#python fama.py --config config.ini --project project_B1F_cazy_t.ini --sample B1F10 --end pe2
#python fama.py --config config.ini --project project_B1F_cazy_t.ini --sample B1F11 --end pe1
#python fama.py --config config.ini --project project_B1F_cazy_t.ini --sample B1F11 --end pe2
#python fama.py --config config.ini --project project_B1F_cazy_t.ini --sample B1F13 --end pe1
#python fama.py --config config.ini --project project_B1F_cazy_t.ini --sample B1F13 --end pe2
#python fama.py --config config.ini --project project_B1F_cazy_t.ini --sample B1F14 --end pe1
#python fama.py --config config.ini --project project_B1F_cazy_t.ini --sample B1F14 --end pe2


# EB271 samples Torben's trimmed reads
#python fama.py --config config.ini --project project_EB271_nitrogen8_t.ini --sample HL1G --end pe1
#python fama.py --config config.ini --project project_EB271_nitrogen8_t.ini --sample HL1G --end pe2
#python fama.py --config config.ini --project project_EB271_nitrogen8_t.ini --sample HL1H --end pe1
#python fama.py --config config.ini --project project_EB271_nitrogen8_t.ini --sample HL1H --end pe2
#python fama.py --config config.ini --project project_EB271_nitrogen8_t.ini --sample HL1I --end pe1
#python fama.py --config config.ini --project project_EB271_nitrogen8_t.ini --sample HL1I --end pe2
#python fama.py --config config.ini --project project_EB271_nitrogen8_t.ini --sample HL1J --end pe1
#python fama.py --config config.ini --project project_EB271_nitrogen8_t.ini --sample HL1J --end pe2
#python fama.py --config config.ini --project project_EB271_nitrogen8_t.ini --sample HL1K --end pe1
#python fama.py --config config.ini --project project_EB271_nitrogen8_t.ini --sample HL1K --end pe2
#python fama.py --config config.ini --project project_EB271_nitrogen8_t.ini --sample HL1L --end pe1
#python fama.py --config config.ini --project project_EB271_nitrogen8_t.ini --sample HL1L --end pe2
#python fama.py --config config.ini --project project_EB271_nitrogen8_t.ini --sample HL1M --end pe1
#python fama.py --config config.ini --project project_EB271_nitrogen8_t.ini --sample HL1M --end pe2
#python fama.py --config config.ini --project project_EB271_nitrogen8_t.ini --sample HL1N --end pe1
#python fama.py --config config.ini --project project_EB271_nitrogen8_t.ini --sample HL1N --end pe2
#python fama.py --config config.ini --project project_EB271_sulfate_t.ini --sample HL1G --end pe1
#python fama.py --config config.ini --project project_EB271_sulfate_t.ini --sample HL1G --end pe2
#python fama.py --config config.ini --project project_EB271_sulfate_t.ini --sample HL1H --end pe1
#python fama.py --config config.ini --project project_EB271_sulfate_t.ini --sample HL1H --end pe2
#python fama.py --config config.ini --project project_EB271_sulfate_t.ini --sample HL1I --end pe1
#python fama.py --config config.ini --project project_EB271_sulfate_t.ini --sample HL1I --end pe2
#python fama.py --config config.ini --project project_EB271_sulfate_t.ini --sample HL1J --end pe1
#python fama.py --config config.ini --project project_EB271_sulfate_t.ini --sample HL1J --end pe2
#python fama.py --config config.ini --project project_EB271_sulfate_t.ini --sample HL1K --end pe1
#python fama.py --config config.ini --project project_EB271_sulfate_t.ini --sample HL1K --end pe2
#python fama.py --config config.ini --project project_EB271_sulfate_t.ini --sample HL1L --end pe1
#python fama.py --config config.ini --project project_EB271_sulfate_t.ini --sample HL1L --end pe2
#python fama.py --config config.ini --project project_EB271_sulfate_t.ini --sample HL1M --end pe1
#python fama.py --config config.ini --project project_EB271_sulfate_t.ini --sample HL1M --end pe2
#python fama.py --config config.ini --project project_EB271_sulfate_t.ini --sample HL1N --end pe1
#python fama.py --config config.ini --project project_EB271_sulfate_t.ini --sample HL1N --end pe2

# EB271 samples
#python fama.py --config config.ini --project project_EB271_nitrogen.ini --sample HL1G --end pe1
#python fama.py --config config.ini --project project_EB271_nitrogen.ini --sample HL1G --end pe2
#python fama.py --config config.ini --project project_EB271_nitrogen.ini --sample HL1H --end pe1
#python fama.py --config config.ini --project project_EB271_nitrogen.ini --sample HL1H --end pe2
#python fama.py --config config.ini --project project_EB271_nitrogen.ini --sample HL1I --end pe1
#python fama.py --config config.ini --project project_EB271_nitrogen.ini --sample HL1I --end pe2
#python fama.py --config config.ini --project project_EB271_nitrogen.ini --sample HL1J --end pe1
#python fama.py --config config.ini --project project_EB271_nitrogen.ini --sample HL1J --end pe2
#python fama.py --config config.ini --project project_EB271_nitrogen.ini --sample HL1K --end pe1
#python fama.py --config config.ini --project project_EB271_nitrogen.ini --sample HL1K --end pe2
#python fama.py --config config.ini --project project_EB271_nitrogen.ini --sample HL1L --end pe1
#python fama.py --config config.ini --project project_EB271_nitrogen.ini --sample HL1L --end pe2
#python fama.py --config config.ini --project project_EB271_nitrogen.ini --sample HL1M --end pe1
#python fama.py --config config.ini --project project_EB271_nitrogen.ini --sample HL1M --end pe2
#python fama.py --config config.ini --project project_EB271_nitrogen.ini --sample HL1N --end pe1
#python fama.py --config config.ini --project project_EB271_nitrogen.ini --sample HL1N --end pe2

#python fama.py --config config.ini --project project_EB271_sulfur.ini --sample HL1G --end pe1
#python fama.py --config config.ini --project project_EB271_sulfur.ini --sample HL1G --end pe2
#python fama.py --config config.ini --project project_EB271_sulfur.ini --sample HL1H --end pe1
#python fama.py --config config.ini --project project_EB271_sulfur.ini --sample HL1H --end pe2
#python fama.py --config config.ini --project project_EB271_sulfur.ini --sample HL1I --end pe1
#python fama.py --config config.ini --project project_EB271_sulfur.ini --sample HL1I --end pe2
#python fama.py --config config.ini --project project_EB271_sulfur.ini --sample HL1J --end pe1
#python fama.py --config config.ini --project project_EB271_sulfur.ini --sample HL1J --end pe2
#python fama.py --config config.ini --project project_EB271_sulfur.ini --sample HL1K --end pe1
#python fama.py --config config.ini --project project_EB271_sulfur.ini --sample HL1K --end pe2
#python fama.py --config config.ini --project project_EB271_sulfur.ini --sample HL1L --end pe1
#python fama.py --config config.ini --project project_EB271_sulfur.ini --sample HL1L --end pe2
#python fama.py --config config.ini --project project_EB271_sulfur.ini --sample HL1M --end pe1
#python fama.py --config config.ini --project project_EB271_sulfur.ini --sample HL1M --end pe2
#python fama.py --config config.ini --project project_EB271_sulfur.ini --sample HL1N --end pe1
#python fama.py --config config.ini --project project_EB271_sulfur.ini --sample HL1N --end pe2

# Hans's samples
#python fama.py --config config.ini --project project_Hans_nitrogen_t.ini --sample HL1A --end pe1
#python fama.py --config config.ini --project project_Hans_nitrogen_t.ini --sample HL1A --end pe2
#python fama.py --config config.ini --project project_Hans_nitrogen_t.ini --sample HL1B --end pe1
#python fama.py --config config.ini --project project_Hans_nitrogen_t.ini --sample HL1B --end pe2
#python fama.py --config config.ini --project project_Hans_nitrogen_t.ini --sample HL1C --end pe1
#python fama.py --config config.ini --project project_Hans_nitrogen_t.ini --sample HL1C --end pe2
#python fama.py --config config.ini --project project_Hans_nitrogen_t.ini --sample HL1D --end pe1
#python fama.py --config config.ini --project project_Hans_nitrogen_t.ini --sample HL1D --end pe2
#python fama.py --config config.ini --project project_Hans_nitrogen_t.ini --sample HL1E --end pe1
#python fama.py --config config.ini --project project_Hans_nitrogen_t.ini --sample HL1E --end pe2
#python fama.py --config config.ini --project project_Hans_nitrogen_t.ini --sample HL1F --end pe1
#python fama.py --config config.ini --project project_Hans_nitrogen_t.ini --sample HL1F --end pe2

#python fama.py --config config.ini --project project_Hans_sulfate_t.ini --sample HL1A --end pe1
#python fama.py --config config.ini --project project_Hans_sulfate_t.ini --sample HL1A --end pe2
#python fama.py --config config.ini --project project_Hans_sulfate_t.ini --sample HL1B --end pe1
#python fama.py --config config.ini --project project_Hans_sulfate_t.ini --sample HL1B --end pe2
#python fama.py --config config.ini --project project_Hans_sulfate_t.ini --sample HL1C --end pe1
#python fama.py --config config.ini --project project_Hans_sulfate_t.ini --sample HL1C --end pe2
#python fama.py --config config.ini --project project_Hans_sulfate_t.ini --sample HL1D --end pe1
#python fama.py --config config.ini --project project_Hans_sulfate_t.ini --sample HL1D --end pe2
#python fama.py --config config.ini --project project_Hans_sulfate_t.ini --sample HL1E --end pe1
#python fama.py --config config.ini --project project_Hans_sulfate_t.ini --sample HL1E --end pe2
#python fama.py --config config.ini --project project_Hans_sulfate_t.ini --sample HL1F --end pe1
#python fama.py --config config.ini --project project_Hans_sulfate_t.ini --sample HL1F --end pe2


#python fama.py --config config.ini --project project_Hans_nitrogen.ini --sample HL1A --end pe1
#python fama.py --config config.ini --project project_Hans_nitrogen.ini --sample HL1A --end pe2
#python fama.py --config config.ini --project project_Hans_nitrogen.ini --sample HL1B --end pe1
#python fama.py --config config.ini --project project_Hans_nitrogen.ini --sample HL1B --end pe2
#python fama.py --config config.ini --project project_Hans_nitrogen.ini --sample HL1C --end pe1
#python fama.py --config config.ini --project project_Hans_nitrogen.ini --sample HL1C --end pe2
#python fama.py --config config.ini --project project_Hans_nitrogen.ini --sample HL1D --end pe1
#python fama.py --config config.ini --project project_Hans_nitrogen.ini --sample HL1D --end pe2
#python fama.py --config config.ini --project project_Hans_nitrogen.ini --sample HL1E --end pe1
#python fama.py --config config.ini --project project_Hans_nitrogen.ini --sample HL1E --end pe2
#python fama.py --config config.ini --project project_Hans_nitrogen.ini --sample HL1F --end pe1
#python fama.py --config config.ini --project project_Hans_nitrogen.ini --sample HL1F --end pe2

#python fama.py --config config.ini --project project_Hans_sulfate.ini --sample HL1A --end pe1
#python fama.py --config config.ini --project project_Hans_sulfate.ini --sample HL1A --end pe2
#python fama.py --config config.ini --project project_Hans_sulfate.ini --sample HL1B --end pe1
#python fama.py --config config.ini --project project_Hans_sulfate.ini --sample HL1B --end pe2
#python fama.py --config config.ini --project project_Hans_sulfate.ini --sample HL1C --end pe1
#python fama.py --config config.ini --project project_Hans_sulfate.ini --sample HL1C --end pe2
#python fama.py --config config.ini --project project_Hans_sulfate.ini --sample HL1D --end pe1
#python fama.py --config config.ini --project project_Hans_sulfate.ini --sample HL1D --end pe2
#python fama.py --config config.ini --project project_Hans_sulfate.ini --sample HL1E --end pe1
#python fama.py --config config.ini --project project_Hans_sulfate.ini --sample HL1E --end pe2
#python fama.py --config config.ini --project project_Hans_sulfate.ini --sample HL1F --end pe1
#python fama.py --config config.ini --project project_Hans_sulfate.ini --sample HL1F --end pe2

# FW 306 samples, Trimmomatic

#python fama.py --config config.ini --project project_FW306_cazy1_t.ini --sample sample1 --end pe1
#python fama.py --config config.ini --project project_FW306_cazy1_t.ini --sample sample2 --end pe1
#python fama.py --config config.ini --project project_FW306_cazy1_t.ini --sample sample3 --end pe1
#python fama.py --config config.ini --project project_FW306_cazy1_t.ini --sample sample4 --end pe1
#python fama.py --config config.ini --project project_FW306_cazy1_t.ini --sample sample5 --end pe1
#python fama.py --config config.ini --project project_FW306_cazy1_t.ini --sample sample6 --end pe1
#python fama.py --config config.ini --project project_FW306_cazy1_t.ini --sample sample1 --end pe2
#python fama.py --config config.ini --project project_FW306_cazy1_t.ini --sample sample2 --end pe2
#python fama.py --config config.ini --project project_FW306_cazy1_t.ini --sample sample3 --end pe2
#python fama.py --config config.ini --project project_FW306_cazy1_t.ini --sample sample4 --end pe2
#python fama.py --config config.ini --project project_FW306_cazy1_t.ini --sample sample5 --end pe2
#python fama.py --config config.ini --project project_FW306_cazy1_t.ini --sample sample6 --end pe2

#python fama.py --config config.ini --project project_FW306_nitrogen8_t.ini --sample sample1 --end pe1
#python fama.py --config config.ini --project project_FW306_nitrogen8_t.ini --sample sample2 --end pe1
#python fama.py --config config.ini --project project_FW306_nitrogen8_t.ini --sample sample3 --end pe1
#python fama.py --config config.ini --project project_FW306_nitrogen8_t.ini --sample sample4 --end pe1
#python fama.py --config config.ini --project project_FW306_nitrogen8_t.ini --sample sample5 --end pe1
#python fama.py --config config.ini --project project_FW306_nitrogen8_t.ini --sample sample6 --end pe1
#python fama.py --config config.ini --project project_FW306_nitrogen8_t.ini --sample sample1 --end pe2
#python fama.py --config config.ini --project project_FW306_nitrogen8_t.ini --sample sample2 --end pe2
#python fama.py --config config.ini --project project_FW306_nitrogen8_t.ini --sample sample3 --end pe2
#python fama.py --config config.ini --project project_FW306_nitrogen8_t.ini --sample sample4 --end pe2
#python fama.py --config config.ini --project project_FW306_nitrogen8_t.ini --sample sample5 --end pe2
#python fama.py --config config.ini --project project_FW306_nitrogen8_t.ini --sample sample6 --end pe2

#python fama.py --config config.ini --project project_FW306_nitrogen_t.ini --sample sample1 --end pe1
#python fama.py --config config.ini --project project_FW306_nitrogen_t.ini --sample sample2 --end pe1
#python fama.py --config config.ini --project project_FW306_nitrogen_t.ini --sample sample3 --end pe1
#python fama.py --config config.ini --project project_FW306_nitrogen_t.ini --sample sample4 --end pe1
#python fama.py --config config.ini --project project_FW306_nitrogen_t.ini --sample sample5 --end pe1
#python fama.py --config config.ini --project project_FW306_nitrogen_t.ini --sample sample6 --end pe1

#python fama.py --config config.ini --project project_FW306_nitrogen_t.ini --sample sample1 --end pe2
#python fama.py --config config.ini --project project_FW306_nitrogen_t.ini --sample sample2 --end pe2
#python fama.py --config config.ini --project project_FW306_nitrogen_t.ini --sample sample3 --end pe2
#python fama.py --config config.ini --project project_FW306_nitrogen_t.ini --sample sample4 --end pe2
#python fama.py --config config.ini --project project_FW306_nitrogen_t.ini --sample sample5 --end pe2
#python fama.py --config config.ini --project project_FW306_nitrogen_t.ini --sample sample6 --end pe2

#python fama.py --config config.ini --project project_FW306_sulfate_t.ini --sample sample1 --end pe1
#python fama.py --config config.ini --project project_FW306_sulfate_t.ini --sample sample2 --end pe1
#python fama.py --config config.ini --project project_FW306_sulfate_t.ini --sample sample3 --end pe1
#python fama.py --config config.ini --project project_FW306_sulfate_t.ini --sample sample4 --end pe1
#python fama.py --config config.ini --project project_FW306_sulfate_t.ini --sample sample5 --end pe1
#python fama.py --config config.ini --project project_FW306_sulfate_t.ini --sample sample6 --end pe1

#python fama.py --config config.ini --project project_FW306_sulfate_t.ini --sample sample1 --end pe2
#python fama.py --config config.ini --project project_FW306_sulfate_t.ini --sample sample2 --end pe2
#python fama.py --config config.ini --project project_FW306_sulfate_t.ini --sample sample3 --end pe2
#python fama.py --config config.ini --project project_FW306_sulfate_t.ini --sample sample4 --end pe2
#python fama.py --config config.ini --project project_FW306_sulfate_t.ini --sample sample5 --end pe2
#python fama.py --config config.ini --project project_FW306_sulfate_t.ini --sample sample6 --end pe2



# FW 306 samples, deduplicated

#python fama.py --config config.ini --project project_fw306_sediment_nitrogen.ini --sample sample1 --end pe1
#python fama.py --config config.ini --project project_fw306_sediment_nitrogen.ini --sample sample2 --end pe1
#python fama.py --config config.ini --project project_fw306_sediment_nitrogen.ini --sample sample3 --end pe1
#python fama.py --config config.ini --project project_fw306_sediment_nitrogen.ini --sample sample4 --end pe1
#python fama.py --config config.ini --project project_fw306_sediment_nitrogen.ini --sample sample5 --end pe1
#python fama.py --config config.ini --project project_fw306_sediment_nitrogen.ini --sample sample6 --end pe1

#python fama.py --config config.ini --project project_fw306_sediment_nitrogen.ini --sample sample1 --end pe2
#python fama.py --config config.ini --project project_fw306_sediment_nitrogen.ini --sample sample2 --end pe2
#python fama.py --config config.ini --project project_fw306_sediment_nitrogen.ini --sample sample3 --end pe2
#python fama.py --config config.ini --project project_fw306_sediment_nitrogen.ini --sample sample4 --end pe2
#python fama.py --config config.ini --project project_fw306_sediment_nitrogen.ini --sample sample5 --end pe2
#python fama.py --config config.ini --project project_fw306_sediment_nitrogen.ini --sample sample6 --end pe2

#python fama.py --config config.ini --project project_fw306_sediment_sulfate.ini --sample sample1 --end pe1
#python fama.py --config config.ini --project project_fw306_sediment_sulfate.ini --sample sample2 --end pe1
#python fama.py --config config.ini --project project_fw306_sediment_sulfate.ini --sample sample3 --end pe1
#python fama.py --config config.ini --project project_fw306_sediment_sulfate.ini --sample sample4 --end pe1
#python fama.py --config config.ini --project project_fw306_sediment_sulfate.ini --sample sample5 --end pe1
#python fama.py --config config.ini --project project_fw306_sediment_sulfate.ini --sample sample6 --end pe1

#python fama.py --config config.ini --project project_fw306_sediment_sulfate.ini --sample sample1 --end pe2
#python fama.py --config config.ini --project project_fw306_sediment_sulfate.ini --sample sample2 --end pe2
#python fama.py --config config.ini --project project_fw306_sediment_sulfate.ini --sample sample3 --end pe2
#python fama.py --config config.ini --project project_fw306_sediment_sulfate.ini --sample sample4 --end pe2
#python fama.py --config config.ini --project project_fw306_sediment_sulfate.ini --sample sample5 --end pe2
#python fama.py --config config.ini --project project_fw306_sediment_sulfate.ini --sample sample6 --end pe2

