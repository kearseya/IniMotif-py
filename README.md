# IniMotif-py

Beta POC code for https://github.com/chengl7/IniMotif

## IniMotif in python

*Alex Kearsey 
*Dr Lu Cheng
*Joel Southgate
*Anna Price
*WebLogo source code

Python scripts for extracting TF motif from SELEXseq data.

## Setup
The recommended way to setup your computer for running Inimotif is to install the official distribution of Python 3.7. You can download the official Python 3.7 distribution [here](https://www.python.org/downloads/release/python-375). Once you've installed Python, open the terminal/prompt and enter the following command:
```
pip install -r requirements.txt
```
This will install the required python packages that Inimotif needs to run. Note that Inimotif may work for versions of Python 3 prior to 3.7, but this is untested. For Mac you may need to update Tcl/Tk to use the GUI for versions of Python prior to 3.7.2. See [here](https://www.python.org/download/mac/tcltk/) for more information.
<br/>
<br/>
Alternatively, you can also run Inimotif using [Anaconda](https://www.anaconda.com/distribution/). To setup a conda environment and install the requirements:
```
conda create -y --name inimotif python==3.7
conda install -fyq --name inimotif -c conda-forge --file requirements.txt
conda activate inimotif
```


## Running IniMotif
You can download Inimotif [here](https://github.com/kearseya/IniMotif-py/archive/master.zip). Open the HTMLmaker.py in python3 and enter data into the pop up window.

* On Mac/Linux <br />
```bash
$ python3 HTMLmaker.py
```
* On Windows <br />
```bash
> python3 HTMLmaker.py
```
<br />
Press the "Enter" button and exit the GUI window. If the GUI does not show, there will be command line input prompts.
<br />
Figures are saved in a folder with the name of the analysis identidier and split into subfolders for the logos, hamming distance, position bias, and kmer frequency.

## Variable inputs:

**Analysis identifier** - The output name for the results page and the folder that the figures are organised in <br />
**Path to directory** - The path to the directory the data files are <br />
**Number of runs** - How many files are being analysed <br />
**Reverse compliment** - Should the kmers have their reverse compliment also counted <br />
**Logo format** - What should the Y-axis be (bits or frequency) for the logos generated <br />
**Min k** - Minimum length of kmers that should be counted <br />
**Max k** - Maximum length of kmers that should be counted <br />
**Experiment type** - Specify what experiment type has the data files come from <br />
**Multiple rounds in the same file** - Are multiple rounds from the experiment contained in the same file (*NB* the 5' and 3' barcodes must be specified to distinguish between the runs) <br />
**Known barcode** - If the 5' and 3' barcodes are known, they can be entered and removed from kmer counting (*this is due to giving a false reading if counted*) <br />
**Start round** - What round/run number does the analysis start from (*with large numbers of rounds in eg. a SELEX experiment, excluding the first few rounds can be faster for analysis*) <br />
**Number of motifs** - How many motifs should be calculated (*generation of more than 1 motif will remove kmers used for previous logos for their generation*) <br />
**Allowed hamming distance** - At what hamming distance value (compared to the consensus kmer) should kmers be used for the making of the position weight matrix <br />

## Tutorial

[ChIP](https://github.com/kearseya/IniMotif-py/blob/master/tutorial/ChIPseq%20tutorial.md) [SELEX](https://github.com/kearseya/IniMotif-py/blob/master/tutorial/SELEXseq%20tutorial.md) [Masker](https://github.com/kearseya/IniMotif-py/blob/master/tutorial/Masker%20tutorial.md) 
