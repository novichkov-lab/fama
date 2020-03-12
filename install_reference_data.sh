#!/bin/sh
script_dir="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"
REFDIR=$script_dir/refdata
if [ -f "$script_dir/refdata/__READY__" ]; then
    echo "Reference data initialization is skipped because it was already prepared"
else
    echo "Reference data initialization"
    if [ -d "$script_dir/refdata" ]; then
        rm -r $script_dir/refdata/*
    else
        mkdir $script_dir/refdata
    fi
    cd $REFDIR
    echo "downloading nitrogen cycle database: http://iseq.lbl.gov/mydocs/fama_downloads/fama_nitrogen_v10.tar.gz"
    curl -LJO -q http://iseq.lbl.gov/mydocs/fama_downloads/fama_nitrogen_v10.tar.gz
    tar xvf fama_nitrogen_v10.tar.gz
    rm fama_nitrogen_v10.tar.gz
    diamond makedb --in fama_nitrogen-cycle_classification_db_v.10.0.faa --db fama_nitrogen-cycle_classification_db_v.10.0
    diamond makedb --in fama_nitrogen-cycle_preselection_db_v.10.0.faa --db fama_nitrogen-cycle_preselection_db_v.10.0
    rm fama_nitrogen-cycle_classification_db_v.10.0.faa
    rm fama_nitrogen-cycle_preselection_db_v.10.0.faa

    echo "downloading Microbe Censys data: http://iseq.lbl.gov/mydocs/fama_downloads/microbecensus_data.tar.gz"
    curl -LJO -q http://iseq.lbl.gov/mydocs/fama_downloads/microbecensus_data.tar.gz
    tar xvf microbecensus_data.tar.gz
    rm microbecensus_data.tar.gz
    diamond makedb --in seqs.fa --db seqs
    rm seqs.fa
    
    echo "downloading taxonomy database: http://iseq.lbl.gov/mydocs/fama_downloads/fama1_taxonomy.tar.gz"
    curl -LJO -q http://iseq.lbl.gov/mydocs/fama_downloads/fama1_taxonomy.tar.gz
    tar xvf fama1_taxonomy.tar.gz
    rm fama1_taxonomy.tar.gz
fi

cd $script_dir
cp config.ini.template config.ini
echo "taxonomy_file = ${REFDIR}/fama_taxonomy.tsv" >> config.ini
echo "microbecensus_data = ${REFDIR}" >> config.ini
echo "" >> config.ini
echo "[nitrogen_v10]" >> config.ini
echo "functions_file = ${REFDIR}/fama_nitrogen-cycle_v.10.0_functions_thresholds.tsv" >> config.ini
echo "taxonomy_file = ${REFDIR}/fama_nitrogen-cycle_v.10.0_taxonomy.tsv" >> config.ini
echo "proteins_list_file = ${REFDIR}/fama_nitrogen-cycle_v.10.0_proteins.txt" >> config.ini
echo "reference_diamond_db = ${REFDIR}/preselection_db_nr100_nitrogen10.dmnd" >> config.ini
echo "reference_db_size = 18388755" >> config.ini
echo "background_diamond_db = ${REFDIR}/classification_db_nr100_nitrogen10.dmnd" >> config.ini
echo "background_db_size = 64090197" >> config.ini

if [ -s "$script_dir/refdata/fama_taxonomy.tsv" ] ; then
  echo "DATA DOWNLOADED SUCCESSFULLY"
  touch $REFDIR/__READY__
else
  echo "Init failed"
  exit 1
fi
