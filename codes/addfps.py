# coding=utf-8
###############################################################################################
### Adding fingerprints to mongodb database with RDKit chemical fingerprints
### This code is similar to db_build.py it is only used if we already have a database of 
### molecule and we only need to generate different fingerprints for the molecules.
###
### Abhik Seal
### 1 Nov 2014
##############################################################################################

"""
usage: addfps.py [-h] --db DB [--fpSize FPSIZE] --fpname FPNAME
                 {morgan,rdkfp,rdmaccs} ...

Generate fingerprints for a Given database in MongoDB

optional arguments:
  -h, --help            show this help message and exit
  --db DB               Input Database Name
  --fpSize FPSIZE       Length of the fingerprints
  --fpname FPNAME       Name of the fp ex:mfp1,mfp2 .. etc

subcommands:
  valid subcommands

  {morgan,rdkfp,rdmaccs}
                        additional help
    morgan              Generate Morgan type fingerprints
    rdkfp               Generate RDKFingerprint
    rdmaccs             Generate MACCS Keys

    example : python addfps.py --db chemtest --fpSize 1024 --fpname mfp2 morgan

"""

import uuid
import pymongo
from bson.binary import Binary
from rdkit import Chem
from rdkit.Chem import AllChem
import os.path
import argparse
from rdkit.Chem import MACCSkeys

def add_fps(fpname,**param):
    """
    Generate fingerprints for every molecule in the database.
    Each fingerprint is stored as an array of "on" bits, along with the total count of this array.

    :param param: add multiple parameters when morgan or rdkfp is used
    :param fpname: name of the fingerprint by which the fingerprints will be stored in mongo db

    """
    min =param.get('minPath')
    max =param.get('maxPath')
    nbitsh=param.get('nBits')
    rad=param.get('rad')
    fpsize=param.get('fpSize')

    if rad is not None :
        for molecule in db.molecules.find({str(fpname): {'$exists': False}}, timeout=False):
            rdmol = Chem.Mol(molecule['rdmol'])
            mfp = list(AllChem.GetMorganFingerprintAsBitVect(rdmol, int(rad) , nBits=int(fpsize)).GetOnBits())
            molecule[str(fpname)] = {'bits': mfp, 'count': len(mfp)}

            db.molecules.save(molecule)
        print ("fingerprints %s done ..." % fpname)
    elif nbitsh is not None:
        for molecule in db.molecules.find({str(fpname): {'$exists': False}}, timeout=False):
            rdmol = Chem.Mol(molecule['rdmol'])
            mfp = list(Chem.RDKFingerprint(rdmol, minPath=int(min),maxPath=int(max),fpSize=fpsize, nBitsPerHash=int(nbitsh)).GetOnBits())
            molecule[str(fpname)] = {'bits': mfp, 'count': len(mfp)}

            db.molecules.save(molecule)
        print ("fingerprints %s done ..." % fpname)
    else:
        for molecule in db.molecules.find({str(fpname): {'$exists': False}}, timeout=False):
            rdmol = Chem.Mol(molecule['rdmol'])
            mfp = list(MACCSkeys.GenMACCSKeys(rdmol).GetOnBits())
            molecule[str(fpname)] = {'bits': mfp, 'count': len(mfp)}

            db.molecules.save(molecule)
        print ("fingerprints %s done ..." % fpname)

def count_fps(fpname):
    """
    Build collection containing total counts of all occurrences of each fingerprint bit.
    This creates a collection that stores counts of the number of times each individual bit occurs in the molecules
    collection. This is used for an efficiency shortcut. The resulting documents have an _id that corresponds to the
    fingerprint bit and a single count field. (e.g. { "_id" : 511, "count" : 148 })

    :param fpname: name of the fingerprint by which the fingerprints will be stored in mongo db

    """
    for fp in [str(fpname)]:
        db['{}_counts'.format(fp)].drop()
        counts = {}
        for molecule in db.molecules.find({fp: {'$exists': True}}, timeout=False):
            for bit in molecule[fp]['bits']:
                counts[bit] = counts.get(bit, 0) + 1
        for k, v in counts.items():
            db['{}_counts'.format(fp)].insert({'_id': k, 'count': v})


def ensure_indices(fpname):
    """
    Build index on relevant fields.

    :param fpname: name of the fingerprint by which the fingerprints will be stored in mongo db

    """
    fpb=str(fpname)+'.bits'
    fpc=str(fpname)+'.count'
    print ("Building Indices...")
    db.molecules.ensure_index(fpb)
    db.molecules.ensure_index(fpc)


if __name__ == '__main__':

    """
       Using command line arguments

    """
    parser = argparse.ArgumentParser(description="Generate fingerprints for a Given database in MongoDB")
    parser.add_argument("--db", required=True, help="Input Database Name")
    parser.add_argument("--fpSize",type=int,help="Length of the fingerprints" )
    parser.add_argument("--fpname" , required=True,help=" Name of the fp ex:mfp1,mfp2 .. etc")

    subparsers = parser.add_subparsers(title='subcommands', description='valid subcommands',
                                         help='additional help')

    morgan_p = subparsers.add_parser('morgan', help="Generate Morgan type fingerprints")
    morgan_p.set_defaults(which='morgan')
    morgan_p.add_argument("--radius",type=int,help="Radius for morgan fingerprints")

    rdkfinger_p=subparsers.add_parser('rdkfp',help="Generate RDKFingerprint")
    rdkfinger_p.set_defaults(which='rdkfp')
    rdkfinger_p.add_argument('--minPath' , help="minimum number of bonds to include in the subgraphs Default 1.")
    rdkfinger_p.add_argument('--maxPath',help="maximum number of bonds to include in the subgraphs Default 7.")
    rdkfinger_p.add_argument('--nBitsPerHash',help="number of bits to set per path Defaults 2.")

    rdmaccs_p=subparsers.add_parser('rdmaccs',help="Generate MACCS Keys")
    rdmaccs_p.set_defaults(which='rdmaccs')

    args = vars(parser.parse_args())
    #name=args['db']

    dbconn="pymongo.MongoClient('localhost', 27017)."+ args['db']
    db = eval(dbconn)
    db.admin.command('shardCollection', args['db']+'.'+args['fpname'], key={'_id': hashed})
    if args['which'] == 'morgan':
        if args['fpSize'] is None and args['radius'] is None:
            add_fps(args['fpname'],fpSize=512,rad=2)
        elif args['fpSize'] is None:
            add_fps(args['fpname'],fpSize=512,rad=args['radius'])
        elif args['radius'] is None :
            add_fps(args['fpname'],fpSize=args['fpSize'],rad=2)
        else:
            add_fps(args['fpname'],fpSize=args['fpSize'],rad=args['radius'])

    elif args['which'] == 'rdkfp':
        if args['fpSize'] is None and args['minPath'] is None and args['maxPath'] is None and args['nBitsPerHash'] is None :
            add_fps(args['fpname'],fpSize=512,minPath=1,maxPath=7,nBitsh=2)
        elif args['minPath'] is None and args['maxPath'] is None and args['nBitsPerHash'] is None :
            add_fps(args['fpname'],fpSize=args['fpSize'],minPath=1,maxPath=7,nBitsh=2)
        else:
            add_fps(args['fpname'],fpSize=args['fpSize'],minPath=args['minPath'],maxPath=args['maxPath'],nBitsh=args['nBitsPerHash'])
    elif args['which'] == 'rdmaccs':
        add_fps(args['fpname'])

    count_fps(args['fpname'])
    ensure_indices(args['fpname'])
