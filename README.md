# Sternberg Lab crRNA Tool
This tool is to help find and evaluate crRNA spacers for the VchINTEGRATE system

## Installation
This library requires a working installation of bowtie2 in the user's path. 
See [here](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2) for information on installing bowtie2. 

Python dependencies can be installed with `pip install -r requirements.txt` run in this directory. 

## Usage
There are two primary functions: spacer generation, and spacer evaluation

#### Spacer generation
This function will generate a number of spacers per given region of the genome. It can be set to target a set of genes, intergenic regions, or custom regions. Modify the variables in `spacer_gen.py` and run the file. A csv will output in the specified directory containing valid spacers targeting the region with minimal off-target potential, as well as information about potential offtargets. 

Be aware that each region can take 2-5 minutes on an ordinary personal computer, so it is not currently practical to run this code against thousands of genes in one run. 

#### Spacer evalution
This function can be used to check for potential off-target activity of a given spacer or set of spacers. Modify the variables in `spacer_eval.py` and run the file. A csv will output with information on off-target potential for each provided spacer, and, for spacers in which potential off-targets were found, a detailed text file with information about each is also generated. 