#!/bin/sh
#python parse_diamond_result.py --config config.ini --col nitrogen_v7 --infile /mnt/data2/SEED/fama/t/data/test_diamond_output.txt --fastq /mnt/data2/FEBA/1/deduped/4701_nodup_PE1.fastq.trimmed.paired.gz
python parse_bgr_search_result.py --config config.ini --project project.ini --sample test_sample --end pe1
