###############################################################################
# Prepare analysis folder
###############################################################################

## Create needed folders
mkdir scripts
mkdir data
mkdir results
mkdir figures

## Move downloaded data to data
sudo mv /home/remote/Descargas/*.tar data/

## Decompress files
for file in data/*.tar
do
  pref="${file//_RAW.tar/}"
  echo $pref
  mkdir $pref
  tar -xvf ${pref}_RAW.tar  -C $pref/
  gunzip -d $pref/*.idat.gz
done

## Copy 450K bad probes from epimutations folder
cp ../Epimutations_analysis/data/HM450.hg19.manifest.rds data/

