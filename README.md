<h1 align="center">LIT-SpaceDock</h1>  

[![Generic badge](https://img.shields.io/badge/version-1.0.0-blue.svg)](https://shields.io/)  

<p align="center">
<img src="https://github.com/litfsindt/spacedock/blob/dc7f6994b584a5a7657549023677f6bcebe78eaf/docs/images/SpaceDock_illus.png" width="600" />
</p>

## Introduction
SpaceDock is a simple and fast computational approach to screen virtualy "ultra-large" synthetic accessible chemical library. It first requires docking commercially available chemical reagents in the target of interest, in order to assemble them in the target's cavity 3D coordinates according to standard organic chemistry reactions to propose multibillion compound libraries in one or two synthetic steps.

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

Downloading the folder and owning Python3.9 is sufficient to run the scripts. </br>
There is no needed external libraries (e.g. NumPy, SciPy, ...).

## Usage

``` bash
$ python3.9 SpaceDock_launcher.py -pbb docked_BBs/ -o DRD3 -r "1 10"
```
Path containing building blocks docking poses (-pbb : docked_BBs/) and reactions 1 and 10 selected.
Outputs :
- DRD3_1.mol2
- DRD3_1.tsv
- DRD3_10.mol2
- DRD3_10.tsv

It produces two differents outputs per selected reactions.
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
- 
