<h1 align="center">LIT-SpaceDock</h1>  

[![Generic badge](https://img.shields.io/badge/version-1.0.0-blue.svg)](https://shields.io/)  

<p align="center">
<img src="https://github.com/litfsindt/spacedock/blob/dc7f6994b584a5a7657549023677f6bcebe78eaf/docs/images/SpaceDock_illus.png" width="600" />
</p>

## Introduction

## Requirements
1. /!\ Python 3.9+
2. A Linux operating system
3. GOLD (https://www.ccdc.cam.ac.uk/solutions/software/gold/) docking poses of Enamine building blocks

## Installation
SpaceDock package :
- spacedock_launcher.py 
- SpaceDock scripts (reagent.py and reaction.py) in SpaceDock/
- Mol2 parser (mol2parser.py) in SpaceDock/
- tag_table_REAL.tsv in datas/ to annote Enamine REAL Space building blocks

There is no external libraries (e.g. NumPy, SciPy, ...).
Downloading the folder and owning Python3.9 is sufficient to run the scripts

## Usage

``` bash
$ python3.9 SpaceDock_launcher.py -pbb docked_BBs/ -o DRD3 -r "1 10"
```
Outputs :
- DRD3.mol2
- DRD3.tsv
For more option, use the SpaceDock_launcher.py --help command

##### Contacts
François Sindt: f.sindt'[at]'unistra.fr  
Didier Rognan, PhD: rognan'[at]'unistra.fr

## References
