Big_data_Project
================

# Automated Code to Build and search a Chemistry database using Mongo DB.

Script ideas from Matt Swain https://github.com/mcs07/mongodb-chemistry
The scripts generates database and adds fingerprints and perform search.

Required tools
 * RDKit
 * Pymongo
 * Mongo DB

To build a database
 * db_build.py  
 * addfps.py
 * query.py


A small dataset is given in data directory to test 

## Usage

python db_build.py -h

python db_build.py --db moltest--i benzodiazepine.smi --tag chembl_id --fpname mfp1 morgan

python addfps.py -h

python addfps.py --db moltest --fpSize 1024 --fpname mfp2 morgan

python query.py -h

python query.py --db chembldb --smi 'CC1=NN=C2N1C3=C(C=C(C=C3)Cl)C(=NC2)C4=CC=CC=C4' --fpSize 512 --fpname mfp1 --t 0.5 --tag chembl_id morgan --radius 2


Check the Readme.pdf for usage or email abhik1368@gmail.com 
