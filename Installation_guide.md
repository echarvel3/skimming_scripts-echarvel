# Installation guide for Genome Skim processing pipelines

**Important Note on Conda:** Most of the tools would be installed through the conda distribution within an activated conda environment in the working shell. 

**Important note on Gurobi** : To be able to run Gurobi and RESPECT, you will need to create an academic license through this [link](https://www.gurobi.com/documentation/9.1/quickstart_mac/obtaining_a_grb_license.html). Download the license and transfer it to whatever computer or cluster you use. When Gurobi is run, it looks for the license file in the defualt locations. On the expiration of the license, you will have to repeat the procedure after removing the `gurobi.lic` file from its default location.

## Installation Instructions

Cloning the repository 
```
git clone https://github.com/echarvel3/skimming_scripts-echarvel.git
cd skimming_scripts-echarvel
```
### Quick Installation
To install the pipeline, run the script called **install.sh**, which should prompt you for a conda enviroment name, and should download and install all the tools correctly. Otherwise, you can install with the instructions below.
```
bash ./install.sh
```
### Manual Installation 

You can **either**...
_Quickly create conda environmet from spec file._ 
!RECOMMENDED!
```
conda create --name "${environment_name}" --file ./Obsolete/respect-spec-file.txt
```
**or** create from scratch.
```
conda create --name "${environment_name}" python=3.7
conda activate "${environment_name}"

### Order matters here
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels https://conda.anaconda.org/gurobi

### Install Skmer

conda install skmer==3.2.1
skmer -h

### The following tools should ideally be installed along with skmer. 
### If not, you can always run this command to install them separately.

conda install jellyfish seqtk mash 

### Run this command to install gurobi solver for respect

conda install gurobi 
```
Then install **RESPECT** (requires gurobi license!).
```
### Install RESPECT
git clone https://github.com/shahab-sarmashghi/RESPECT.git
cd RESPECT
python setup.py install
cd ..
```
```
### BBMap has been made available as a part of the repository and you can use it directly when we clone the repository later
### You can always run the code below to download BBMap separately, if you want to
### wget -O bbmap.tar.gz https://sourceforge.net/projects/bbmap/files/BBMap_39.01.tar.gz/download
### tar xvfz bbmap.tar.gz
### rm bbmap.tar.gz

### FastME has been made available as a part of the repository and you can use it directly when we clone the repository later
### You can also install FastME (to get backbone trees) separately using the following commands
### wget http://www.atgc-montpellier.fr/download/sources/fastme/fastme-2.1.5.tar.gz
### tar xvfz fastme-2.1.5.tar.gz
### chmod +x fastme-2.1.5/binaries/fastme-2.1.5-linux64
### Change "linux64" at the end if using other platforms (linux32 or windows).
### ./fastme-2.1.5/binaries/fastme-2.1.5-linux64 -h
```
Installing KRANK (for contamination removal)
```
git clone https://github.com/bo1929/KRANK.git
make -C KRANK
cd ./KRANK/
wget https://ter-trees.ucsd.edu/data/krank/lib_reps_adpt-k29_w35_h13_b16_s8.tar.gz
tar -zxf ./https://ter-trees.ucsd.edu/data/krank/lib_reps_adpt-k29_w35_h13_b16_s8.tar.gz
```
Running the pipeline on sample data
```
gunzip ./test/skims/*
bash ../fast_skims_pipeline2.sh -i ./test/skims/ -t "${threads}"

```

