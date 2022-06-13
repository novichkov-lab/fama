#!/bin/sh
script_dir="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"

REFDIR=$1
if ! [ -d "$REFDIR" ]; then
	echo "Reference data directory $REFDIR doesn't exist"
	exit
fi

if [ -f "$REFDIR/fama/__READY__" ]; then
    echo "Reference data generation is skipped because it was already prepared"
else
    echo "Reference data initialization"
    if [ -d "$REFDIR/fama" ]; then
        rm -r $REFDIR/fama/*
    else
        mkdir $REFDIR/fama
    fi
    cd $REFDIR
    echo "downloading Fama reference database: http://iseq.lbl.gov/mydocs/fama_downloads/cgcms_fama_ref.tar.gz"
    curl -LJO -q http://iseq.lbl.gov/mydocs/fama_downloads/cgcms_fama_ref.tar.gz
    tar xvf cgcms_fama_ref.tar.gz
    rm cgcms_fama_ref.tar.gz
	cd "$REFDIR/fama"

    echo "downloading Microbe Census data: http://iseq.lbl.gov/mydocs/fama_downloads/microbecensus_data.tar.gz"
    curl -LJO -q http://iseq.lbl.gov/mydocs/fama_downloads/microbecensus_data.tar.gz
    tar xvf microbecensus_data.tar.gz
    rm microbecensus_data.tar.gz
    diamond makedb --in seqs.fa --db seqs
    rm seqs.fa
    
    echo "downloading taxonomy database: http://iseq.lbl.gov/mydocs/fama_downloads/fama1_taxonomy.tar.gz"
    curl -LJO -q http://iseq.lbl.gov/mydocs/fama_downloads/fama_taxonomy.tar.gz
    tar xvf fama_taxonomy.tar.gz
    rm fama_taxonomy.tar.gz

    echo "Making DIAMOND databases: http://iseq.lbl.gov/mydocs/fama_downloads/fama1_taxonomy.tar.gz"
	
	cd "$REFDIR/fama/nitrogen11"
    diamond makedb --in preselection_db_nr100.faa --db preselection_db_nr100
    diamond makedb --in classification_db_nr100.faa --db classification_db_nr100
    rm preselection_db_nr100.faa
    rm classification_db_nr100.faa

	cd "$REFDIR/fama/universal1.4"
    diamond makedb --in preselection_db_nr100.faa --db preselection_db_nr100
    diamond makedb --in classification_db_nr100.faa --db classification_db_nr100
    rm preselection_db_nr100.faa
    rm classification_db_nr100.faa
	
	cd "$REFDIR/fama/cazy2"
    diamond makedb --in preselection_database.faa --db preselection_database
    diamond makedb --in classification_database.faa --db classification_database.faa
    rm preselection_database.faa
    rm classification_database.faa

fi

echo "Generating config.ini"
cd $script_dir
cp config.ini.template config_cgcms.ini
echo "taxonomy_file = ${REFDIR}/fama/fama_taxonomy.tsv" >> config.ini
echo "microbecensus_data = ${REFDIR}/fama" >> config.ini
echo "" >> config.ini
echo "[nitrogen_v11]" >> config.ini
echo "functions_file = ${REFDIR}/fama/nitrogen11/fama_nitrogen-cycle_v.11.0_functions_thresholds.tsv" >> config.ini
echo "taxonomy_file = ${REFDIR}/fama/nitrogen11/fama_nitrogen-cycle_v.11.0_taxonomy.tsv" >> config.ini
echo "proteins_list_file = ${REFDIR}/fama/nitrogen11/fama_nitrogen-cycle_v.11.0_proteins.txt" >> config.ini
echo "reference_diamond_db = ${REFDIR}/fama/nitrogen11/preselection_db_nr100.dmnd" >> config.ini
echo "reference_db_size = 18388755" >> config.ini
echo "background_diamond_db = ${REFDIR}/fama/nitrogen11/classification_db_nr100.dmnd" >> config.ini
echo "background_db_size = 64090197" >> config.ini

echo "" >> config.ini
echo "[universal_v1.4]" >> config.ini
echo "functions_file = ${REFDIR}/fama/universal1.4/fama_function_thresholds_v.1.4.txt" >> config.ini
echo "taxonomy_file = ${REFDIR}/fama/universal1.4/fama_universal_taxonomy_v.1.4.txt" >> config.ini
echo "proteins_list_file = ${REFDIR}/fama/universal1.4/fama_universal_v.1.4.txt" >> config.ini
echo "reference_diamond_db = ${REFDIR}/fama/universal1.4/preselection_db_nr100.dmnd" >> config.ini
echo "reference_db_size = 31895938" >> config.ini
echo "background_diamond_db = ${REFDIR}/fama/universal1.4/classification_db_nr100.dmnd" >> config.ini
echo "background_db_size = 49177580" >> config.ini

echo "" >> config.ini
echo "[cazy_v2]" >> config.ini
echo "functions_file = ${REFDIR}/fama/cazy2/collection_functions.txt" >> config.ini
echo "taxonomy_file = ${REFDIR}/fama/cazy2/collection_taxonomy.txt" >> config.ini
echo "proteins_list_file = ${REFDIR}/fama/cazy2/final_gene_list.txt" >> config.ini
echo "reference_diamond_db = ${REFDIR}/fama/cazy2/preselection_database.dmnd" >> config.ini
echo "reference_db_size = 238821857" >> config.ini
echo "background_diamond_db = ${REFDIR}/fama/cazy2/classification_database.dmnd" >> config.ini
echo "background_db_size = 1308452213" >> config.ini

if [ -s "${REFDIR}/fama/fama_taxonomy.tsv" ] ; then
  echo "DATA DOWNLOADED SUCCESSFULLY"
  touch $REFDIR/fama/__READY__
else
  echo "Init failed"
  exit 1
fi
