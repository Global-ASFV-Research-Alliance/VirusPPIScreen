# VirusPPIScreen

Scripts to run genome-wide virus-virus PPI screens with AlphaFold v2 

## Installation and Dependencies
You will need a machine that can run bash scripts 
1. install [Miniconda](https://docs.anaconda.com/free/miniconda/).
1. clone this repo and 'cd' into it
    ```bash
    git clone https://github.com/Global-ASFV-Research-Alliance/VirusPPIScreen.git 
    cd ./VirusPPIScreen 
    ```
1. create conda environment and activate it
    ```bash
    conda env create --file=VirusPPIScreen_env.yml
    conda activate VirusPPIScreen
    ```
1. download software dependencies with scripts/download_dependencies.sh
Please visit [hhsuite](https://mmseqs.com/hhsuite/) and choose the appropriate compiled 
binary filename for your system
    ```bash
    chmod +x scripts/download_dependencies.sh
    ./scripts/download_dependencies.sh hhsuite-linux-sse2.tar.gz
    ```

## Test the installation
examples/test_pmsa.sh generates a genome wide pMSAs for the ASFV_test_genome.fa 
using the test database example_viral_db.fasta and the test AlphaFold outputs
AlphaFold_MSAs/

    ```bash
    chmod +x examples/test_pmsa.sh
    ./examples/test_pmsa.sh
    ```
As the ASFV_test_genome.fa has two proteins, this script will produce one heterodimer
and two homodimer pMSAs in examples/example_pmsa_out/pmsas/
This scripts should take about 10-20 minutes to complete and will use two CPU cores 

## Author 
Jacob A. Fenster
