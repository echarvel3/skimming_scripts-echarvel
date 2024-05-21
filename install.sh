# Installation guide for Genome Skim processing pipelines

#**Important Note on Conda:** Most of the tools would be installed through the conda distribution within an activated conda environment in the working shell. 

#* Whenever you are using any of the pipelines, make sure that you change the env name in the `conda_source.sh` [script](https://github.com/smirarab/skimming_scripts/blob/master/conda_source.sh). Refer to `CONDAENV=GSkim4` in the script where `GSkim4` should be changed to your environment's name corresponding to all the tool installations.
#* In the example below, we would change the `conda_source.sh` [script](https://github.com/smirarab/skimming_scripts/blob/master/conda_source.sh) to say `CONDAENV=tutorial`

#iThe `conda_source.sh` [script](https://github.com/smirarab/skimming_scripts/blob/master/conda_source.sh) is meant to allow users to edit the name of their conda environment to easily switch between various configurations. A sample conda environment configuration can be seen [here](https://github.com/smirarab/skimming_scripts/blob/master/Obsolete/environment.yml). 

echo "Please enter:"
read -p "Name of New Conda Environment: -> " conda_environment

conda create --name ${conda_environment} python=3.8 --yes --quiet
conda activate ${conda_environment}

echo "Adding Channels to Conda Enviroment:(${conda_environment})..."
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels https://conda.anaconda.org/gurobi

echo "Installing SKMER..."
conda install skmer==3.2.1 --yes --quiet
skmer -h || echo "ERROR: SKMER INSTALLATION FAILED." && exit 1 

jellyfish -h || echo "Installing Jellyfish..." && conda install jellyfish --yes --quiet
seqtk -h || echo "Installing SEQTK..." && conda install seqtk --yes --quiet
mash -h || echo "Installing MASH..." && conda install mash --yes --quiet

echo "Installing GUROBI..."
conda install gurobi --yes --quiet

echo "Finished Setting up Conda Environment: (${conda_environment})"
echo "==================="

echo "Installing RESPECT..."
git clone https://github.com/shahab-sarmashghi/RESPECT.git
cd RESPECT
python setup.py install
cd ..

#echo "Installing CONSULT-II..."
#git clone https://github.com/bo1929/CONSULT-II
#cd CONSULT-II/
#make all
#wget https://ter-trees.ucsd.edu/data/consult/CONSULT-II/library-v030-WoL18G.tar.gz
#cd ..

echo "Installing KRANK..."
git clone https://github.com/bo1929/KRANK.git
make -C ./KRANK/
cd ./KRANK/
echo "Downloading Bacterial Decontamination Library..."
wget https://ter-trees.ucsd.edu/data/krank/lib_reps_adpt-k29_w35_h13_b16_s8.tar.gz
tar -zxf ./lib_reps_adpt-k29_w35_h13_b16_s8.tar.gz
rm ./lib_reps_adpt-k29_w35_h13_b16_s8.tar.gz
cd ../

### Installing kraken2
#git clone https://github.com/DerrickWood/kraken2
#cd kraken2/
#./install_kraken2.sh ./
### Installing relevant library for human contaminant removal using kraken2
#./kraken2-build --download-taxonomy --db krakenlib
#./kraken2-build --download-library human --db krakenlib/
#./kraken2-build --build --db krakenlib/
#cd ..

