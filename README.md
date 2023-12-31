<h1 align="center">LIT-SpaceDock</h1>  

[![Generic badge](https://img.shields.io/badge/version-1.0.0-blue.svg)](https://shields.io/)  

<p align="center"><img src="https://github.com/litfsindt/LIT-SpaceDock/blob/master/img/SpaceDock_illus.png" width="600" /></p>


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
$ python3.9 SpaceDock_launcher.py -pbb docked_BBs_DRD3/ -o DRD3 -r 1 10
```
This command will initiate SpaceDock using the path containing building block poses (-pbb docked_BBs/) and using reactions 1 and 10 (-r 1 10).

#### Outputs :
- DRD3_1.mol2
- DRD3_1.tsv
- DRD3_10.mol2
- DRD3_10.tsv

It produces two different outputs for each selected reaction :
- mol2 file : Multiple combinations of two building blocks,, referred to "SpaceDock poses"
- tsv file : A tab-separated file containing the file names of the two building block poses for each combination.

#### Avalaible reactions (-r) : Version 1.0.0
- all : All possible reactions
- 1 : Amide
- 10 : Sulfonamide
- 27 : Benzoxazole

For more option, use the SpaceDock_launcher.py --help command

#### Filtering based on optional IChem IFP similarity
``` bash
$ python3.9 SpaceDock_launcher.py -pbb docked_BBs_DRD3/ -o DRD3 -r 1 10 --ichem IChem/ichem_DRD3.conf
```
This command allows the *on the fly* filtering of combinations based on IChem* interactions fingerprint similarity. 
A reference ligand is required, and the ichem.conf file (examples provided in the IChem/ folder) must be filled out.
*Executable and a license key is available in the IChem/ folder.

## Energy-minimization with SZYBKI
Once the "SpaceDock poses" are generated, they can be subject to a brief energy minimization using SZYBKI* under protein constraints to relax the newly formed bond.
``` bash
$ szybki -protein DRD3_3PBL_protein.mol2 -in DRD3_1.mol2 -out DRD3_1_minimized.mol2 -prefix DRD3_1_minimized -ff mmff94 -optGeometry cart -exact_vdw -heavy_rms
```
*https://www.eyesopen.com/szybki

## Examples
Certain building blocks have been already docked with GOLD in the Dopamine D3 Receptor (located in the docked_BBs_DRD3/ folder) and in the Estrogen Beta Receptor (docked_BBs_ERB/). These are ready for use in SpaceDock example combinations (with -pbb option).

#### Structures
Dopamine D3 Receptor (PDB : 3PBL) and Estrogen Beta Receptor (PDB : 1QKM) structures used for this work are available in structures/ folder.

## Contacts
François Sindt: f.sindt'[at]'unistra.fr

Didier Rognan, PhD: rognan'[at]'unistra.fr

## References
- Hartenfeller M, et al. A collection of robust organic synthesis reactions for in silico molecule design. J Chem Inf Model 51, 3093-3098 (2011).
- Da Silva F, Desaphy J, Rognan D. IChem: A Versatile Toolkit for Detecting, Comparing, and Predicting Protein-Ligand Interactions. ChemMedChem 13, 507-510 (2018).
- SZBYKI v2.4.0.0: OpenEye, Cadence Molecular Sciences, Santa Fe, NM, U.S.A., www.eyesopen.com

