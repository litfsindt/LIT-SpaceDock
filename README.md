<h1 align="center">LIT-SpaceDock</h1>  

[![Generic badge](https://img.shields.io/badge/version-1.0.0-blue.svg)](https://shields.io/)  

<p align="center">
<img src="https://github.com/litfsindt/spacedock/blob/dc7f6994b584a5a7657549023677f6bcebe78eaf/docs/images/SpaceDock_illus.png" width="600" />
</p>

## Introduction
SpaceDock is a straightforward and rapid computational approach designed to screen virtually 'ultra-large' synthetically accessible chemical libraries. It first requires docking commercially available chemical reagents into the target of interest, assembling them in the target's cavity based on standard organic chemistry reactions, thereby proposing multibillion compound libraries achievable in one or two synthetic steps.

## Requirements
1. /!\ Python 3.9+
2. A Linux operating system
3. Enamine building blocks poses docked with GOLD*

*https://www.ccdc.cam.ac.uk/solutions/software/gold/

## Installation
SpaceDock package :
- spacedock_launcher.py 
- SpaceDock scripts (reagent.py and reaction.py) in SpaceDock/
- Mol2 parser (mol2parser.py) in SpaceDock/
- tag_table_REAL.tsv in datas/ to annote Enamine REAL Space building blocks

Downloading the folder and having Python 3.9 installed is sufficient to run the scripts. 
No external libraries (e.g., NumPy, SciPy, ...) are required.

## Usage

``` bash
$ python3.9 SpaceDock_launcher.py -pbb docked_BBs/ -o DRD3 -r "1 10"
```
This command will initiate SpaceDock using the path containing building block poses (-pbb: docked_BBs/) and using reactions 1 and 10 (-r "1 10").
Outputs :
- DRD3_1.mol2
- DRD3_1.tsv
- DRD3_10.mol2
- DRD3_10.tsv

It produces two different outputs for each selected reaction.

Avalaible reactions (-r) : 
- all : All possible reactions
- 1 : Amide
- 10 : Sulfonamide
- 27 : Benzoxazole

For more option, use the SpaceDock_launcher.py --help command

## Contacts
François Sindt: f.sindt'[at]'unistra.fr  
Didier Rognan, PhD: rognan'[at]'unistra.fr

## References
- Hartenfeller M, et al. A collection of robust organic synthesis reactions for in silico molecule design. J Chem Inf Model 51, 3093-3098 (2011).

