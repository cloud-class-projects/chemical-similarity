##############################################################################
### Searching the mongodb database with RDKit chemical fingerprints
### This performs searching with a single query in smiles format
###
### Abhik Seal
### 1 November 2014
### Code taken from Matt Swain's github https://github.com/mcs07/mongodb-chemistry
###
##############################################################################

"""
usage: query.py [-h] --db DB --smi SMI [--fpSize FPSIZE] --fpname FPNAME
                [--t T] [--tag TAG]
                {morgan,rdkfp,rdmaccs} ...

Search MongoDB database

optional arguments:
  -h, --help            show this help message and exit
  --db DB               Input Database Name
  --smi SMI             Enter the smiles string
  --fpSize FPSIZE       Length of the fingerprints
  --fpname FPNAME       Name of the fp ex:mfp1,mfp2 .. etc
  --t T                 Similarity threshold
  --tag TAG             tag name in the original database Ex:chembl_id

subcommands:
  valid subcommands

  {morgan,rdkfp,rdmaccs}
                        additional help
    morgan              Generate Morgan type fingerprints
    rdkfp               Generate RDKFingerprint
    rdmaccs             Generate MACCS Keys

  example : python monQuery.py --db chemtest --smi 'CC1=NN=C2N1C3=C(C=C(C=C3)Cl)C(=NC2)C4=CC=CC=C4' --fpSize 512
            --fpname mfp1 --t 0.8 --tag chembl_id morgan --radius 2

"""

import pymongo
from rdkit import Chem
from rdkit.Chem import AllChem
import time
import numpy as np
import datetime
from rdkit.Chem import MACCSkeys
import argparse


def similarity(qfp,con,fpname,tag ,t=0.7):
    """Perform a similarity search using aggregation with new operators from MongoDB 2.6.

    :rtype : object
    :param qfp: The query fingerprint
    :param threshold: The tanimoto threshold
    """
    db = eval(con)
    qn = qfp[1]                           # Number of bits in query fingerprint
    qmin = int(qn * t)               # Minimum number of bits in results fingerprints
    qmax = int(qn / t)               # Maximum number of bits in results fingerprints
    ncommon = qn - qmin + 1                 # Number of fingerprint bits in which at least 1 must be in common
    mc='db.'+fpname+'_counts.find'
    fpc=fpname+'.count'
    fpb=fpname+'.bits'

    reqbits = [count['_id'] for count in eval(mc)({'_id': {'$in': qfp[0]}}).sort('count', 1).limit(ncommon)]
    #bit1 = qfp[:ncommon]
    print "start Aggregate .. "
    aggregate = [
        {'$match': {str(fpc): {'$gte': qmin, '$lte': qmax}, str(fpb): {'$in': reqbits}}},
        {'$project': {
            'tanimoto': {'$let':
            {
                'vars': {'common': {'$size': {'$setIntersection': [str('$'+fpb), qfp[0]]}}},
                'in': {'$divide': ['$$common', {'$subtract': [{'$add': [qn, str('$'+fpc)]}, '$$common']}]}
            }},
            tag: 1
        }},
        {'$match':  {'tanimoto': {'$gte': t}}}
    ]

    response = db.molecules.aggregate(aggregate)
    print "response done .. "
    # for result in response['result']:
    #     print '%s : %s : %s' % (result['tanimoto'] * 100, result['chembl_id'], result['smiles'])
    return response['result']


def genfp(mol,t,con,**param):
    """Generate fingerprints for every molecule in the database.

       Each fingerprint is stored as an array of "on" bits, along with the total count of this array.
    """
    q=[]
    report = []
    times = []
    counts = []
    min =param.get('minPath')
    max =param.get('maxPath')
    nbitsh=param.get('nBitsPerHash')
    rad=param.get('radius')
    fpsize=param.get('fpSize')
    fpname=param.get('fpname')
    tag=param.get('tag')
    if rad is not None :
        mfp = list(AllChem.GetMorganFingerprintAsBitVect(mol, int(rad), nBits=int(fpsize)).GetOnBits())
        q.append(mfp)
        q.append(len(mfp))
        print ("fingerprints %s done ..." % fpname)
    elif nbitsh is not None:
        mfp = list(Chem.RDKFingerprint(mol, minPath=int(min),maxPath=int(max),fpSize=fpsize, nBitsPerHash=int(nbitsh)).GetOnBits())
        q.append(mfp)
        q.append(len(mfp))
    else:
        mfp = list(MACCSkeys.GenMACCSKeys(mol).GetOnBits())
        q.append(mfp)
        q.append(len(mfp))
    start = time.time()
    results = similarity(q,con,fpname,tag,t)
    end=time.time()
    counts.append(len(results))
    report.append('Results ({})'.format(len(results)))
    times.append(end - start)
    # Produce a report of the results
    report.append('Counts: {}'.format(' '.join(str(c) for c in counts)))
    report.append('Times: {}'.format(' '.join(str(t) for t in times)))
    print "Hits: ", len(results),"Time : ",end-start
    for r in results:
        report.append('{}: {}'.format(r['tanimoto'], r['chembl_id']))
        print report[-1]


def main():

    """
          Searching the database with a smiles string Using command line arguments

    """
    parser = argparse.ArgumentParser(description="Search MongoDB database")
    parser.add_argument("--db", required=True, help="Input Database Name")
    parser.add_argument("--smi",required=True,help="Enter the smiles string")
    parser.add_argument("--fpSize",type=int,help="Length of the fingerprints" )
    parser.add_argument("--fpname" , required=True,help=" Name of the fp ex:mfp1,mfp2 .. etc")
    parser.add_argument("--t",type=float,help="Similarity threshold")
    parser.add_argument("--tag",help="tag name in the original database Ex:chembl_id")
    subparsers = parser.add_subparsers(title='subcommands', description='valid subcommands',
        help='additional help')

    morgan_p = subparsers.add_parser('morgan', help="Generate Morgan type fingerprints")
    morgan_p.set_defaults(which='morgan')
    morgan_p.add_argument("--radius",required=True,type=int,help="Radius for morgan fingerprints")

    rdkfinger_p=subparsers.add_parser('rdkfp',help="Generate RDKFingerprint")
    rdkfinger_p.set_defaults(which='rdkfp')
    rdkfinger_p.add_argument('--minPath' , help="minimum number of bonds to include in the subgraphs Default 1.")
    rdkfinger_p.add_argument('--maxPath',help="maximum number of bonds to include in the subgraphs Default 7.")
    rdkfinger_p.add_argument('--nBitsPerHash',help="number of bits to set per path Defaults 2.")

    rdmaccs_p=subparsers.add_parser('rdmaccs',help="Generate MACCS Keys")
    rdmaccs_p.set_defaults(which='rdmaccs')

    args = vars(parser.parse_args())
    con="pymongo.MongoClient('localhost', 27017)."+ args['db']

    mol = Chem.MolFromSmiles(str(args['smi']))
    t=args['t']

    if args['which'] == 'morgan':

        genfp(mol,t,con,fpname=args['fpname'],fpSize=args['fpSize'],radius=args['radius'],tag=args['tag'])

    elif args['which'] == 'rdkfp':
        genfp(mol,t,con,fpname=args['fpname'],fpSize=args['fpSize'],minPath=args['minPath'],maxPath=args['maxPath'],nBitsPerHash=args['nBitsPerHash'],tag=args['tag'])
    elif args['which'] == 'rdmaccs':
        genfp(mol,t,con,fpname=args['fpname'],tag=args['tag'])

if __name__ == '__main__':
    main()