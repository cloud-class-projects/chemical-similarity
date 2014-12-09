# coding=utf-8

###############################################################################################
### Generate text result of chemical search using 100,250,500 and 1000 molecules using
### at different similarity thresholds 0.65,0.75,0.85,0.95
### __author__ = 'abhikseal'
### 1 Nov 2014
### usage : python profile.py
##############################################################################################

import random
import time
import numpy as np
import pymongo


def choose_random(N):
    """Choose 1000 random molecules from the database.

    :param N : Integer  Number of random molecules

    """
    db = pymongo.MongoClient('localhost',27020).chembldb
    # Get all CHEMBL IDs
    db.molecules.ensure_index('chembl_id')
    chembl_ids = [m['chembl_id'] for m in db.molecules.find().sort('chembl_id')]
    print len(chembl_ids)
    random.seed(201405291515)
    rands = random.sample(chembl_ids, N)
    return(rands)

def similarity(qfp,db, threshold=0.8):
    """
    Perform similarity search using the aggregate function

    :param qfp:  query fingerprint
    :param db: A mongo database connection
    :param threshold: integer similarity threshold
    :return: return aggregated results similarity scores , chembl_ids

    """
    qn = len(qfp)                           # Number of bits in query fingerprint
    qmin = int(qn * threshold)               # Minimum number of bits in results fingerprints
    qmax = int(qn / threshold)               # Maximum number of bits in results fingerprints
    ncommon = qn - qmin + 1                 # Number of fingerprint bits in which at least 1 must be in common
    reqbits = [count['_id'] for count in db.mfp2_counts.find({'_id': {'$in': qfp}}).sort('count', 1).limit(ncommon)]
    aggregate = [
        {'$match': {'mfp2.count': {'$gte': qmin, '$lte': qmax}, 'mfp2.bits': {'$in': reqbits}}},
        {'$project': {
            'tanimoto': {'$let': {
                'vars': {'common': {'$size': {'$setIntersection': ['$mfp2.bits', qfp]}}},
                'in': {'$divide': ['$$common', {'$subtract': [{'$add': [qn, '$mfp2.count']}, '$$common']}]}
            }},
            'smiles': 1,
            'chembl_id': 1
        }},
        {'$match':  {'tanimoto': {'$gte': threshold}}}
    ]
    response = db.molecules.aggregate(aggregate)
    return response['result']

def profile_mongodb(chembl_ids):

    """
    Do multiple molecules queries and generates report on mean , median and 95th percentile timings
    :param chembl_ids: set of chembl_ids against which similarity search will be performed

    """
    db = pymongo.MongoClient('localhost',27020).chembldb
    size=len(chembl_ids)
    for threshold in [0.65,0.75,0.85,0.95]:
        report = []
        times = []
        counts = []

        # Perform the queries
        for chembl_id in chembl_ids:
            print chembl_id
            qmol = db.molecules.find_one({'chembl_id': chembl_id})
            report.append('Query: {}'.format(qmol['chembl_id']))
            start = time.time()
            results = similarity(qmol['mfp2']['bits'],db,threshold)
            end = time.time()
            counts.append(len(results))
            report.append('Results ({})'.format(len(results)))
            print report[-1]
            for r in results:
                report.append('{}: {}'.format(r['tanimoto'], r['chembl_id']))
                print report[-1]
            times.append(end - start)

        # Produce a report of the results
        report.append('Counts: {}'.format(' '.join(str(c) for c in counts)))
        report.append('Mean: {}'.format(np.mean(counts)))
        report.append('Median: {}'.format(np.median(counts)))
        report.append('Times: {}'.format(' '.join(str(t) for t in times)))
        report.append('Mean: {}'.format(np.mean(times)))
        report.append('Median: {}'.format(np.median(times)))
        report.append('95th Percentile: {}'.format(np.percentile(times, 95)))
        report = '\n'.join(report)
        #print report
        with open('results-{}-{}--rdmaccs.txt'.format(threshold,size), 'w') as f:
            f.write(report)

if __name__ == '__main__':

    # number of molecules to run multiple molecule queries

    n=[100,250,500,1000]
    for i in n:
        rands=choose_random(N=i)
        profile_mongodb(rands)
        print "Computation done for .. %d " %i

