#!/bin/sh

python3 fama.py -c config_test.ini -p /mnt/data2/ENIGMA/genomes/internal/Fama/internal_protein_univ.ini --prot
python3 fama.py -c config_test.ini -p /mnt/data2/ENIGMA/genomes/internal/Fama/internal_protein_nitr.ini --prot
python3 fama.py -c config_test.ini -p /mnt/data2/ENIGMA/genomes/internal/Fama/internal_protein_cazy.ini --prot

#python3 fama.py -c config_test.ini -p /mnt/data2/Bacter/Acidovorax/comparative_genomics/proteins/acidovorax_nitrogen_project2.ini --prot

#CAMI testing
#python3 fama.py -c config_test.ini -p ../p/project_cami_universal1.ini
#python3 fama.py -c config_test.ini -p ../p/project_cami_nitrogen10.ini
#python3 fama_prepare.py -c config_test.ini -i /mnt/data2/Metagenomics/benchmark/cami/CAMI_medium/fastq.list -r universal_v1.4 -p cami_toy_med_universal_v.1.4
#python3 fama.py -c config_test.ini -p /mnt/data2/Metagenomics/benchmark/cami/CAMI_high/project_universal_v1.4_cami.ini
#python3 fama_assembly.py -c config_test.ini -p /mnt/data2/Metagenomics/benchmark/cami/CAMI_high/project_universal_v1.4_cami.ini --coassembly
#python3 fama.py -c config_test.ini -p /mnt/data2/Metagenomics/benchmark/cami/CAMI_high/ref_data/fama_cami_high_ref_universal.ini --prot
#python3 fama.py -c config_test.ini -p /mnt/data2/Metagenomics/benchmark/cami/CAMI_medium/project_universal_v1.4.ini

# Test project
# python3 fama.py -c config_test.ini -p /mnt/data2/SEED/fama/p/project_FW3062M_universal1.4.ini
#python3 fama.py -c config_test.ini -p /mnt/data2/SEED/fama/p/project_FW3062M_universal1.ini

# ENIGMA isolates, non-pseudomonad, nitrogen v.10
#python3 fama.py -c config_test.ini -p ../p/project_protein_enigma_nonpseudomonas_nitrogen10.ini --prot
#python3 fama.py -c config_test.ini -p ../p/project_protein_ht_ok_isolates_nitrogen10.ini --prot
#python3 fama.py -c config_test.ini -p ../p/project_protein_enigma_isolates_nitrogen10.ini --prot
#python3 fama.py -c config_test.ini -p ../p/project_protein_ht1_isolates_nitrogen10.ini -s GW821-FHT02F09 --prot

# Protein projects
# python3 fama.py -c config_test.ini -p /mnt/data2/SEED/fama/p/project_protein_ht1_isolates_universal1.4.ini --prot
# python3 fama.py -c config_test.ini -p /mnt/data2/SEED/fama/p/project_protein_enigma_nonpseudomonas_universal1.4.ini --prot
# python3 fama.py -c config_test.ini -p /mnt/data2/SEED/fama/p/project_protein_pseudo_isolates_universal_v1.4.ini --prot


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

