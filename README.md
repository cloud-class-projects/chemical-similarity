Big_data_Project
================

# Automated Code to Build and search a Chemistry database using Mongo DB by Sharding

Script ideas from Matt Swain https://github.com/mcs07/mongodb-chemistry
The scripts generates database and adds fingerprints and perform search.For configuring Mongo DB across multiple machine Cloudmesh Python API is used . The vmstart.py code starts the ubuntu machines as shard machines in India grid and write, executes the filescript.sh script given below in the box.

### Required tools
 * RDKit
 * Pymongo
 * Mongo DB
 * cloudmesh

### Setting up the machines

Config servers are lie the brains of the cluster: they hold all of the metadata about which servers hold what data.Thus, they must be set up first and the data they hold is extremely important. Each config server should be on a separate physical machine .For the config server and mongos router separatelya VMs are started and Mongo DB is installed on it. The config servers must be started before any of the mongos processes, as mongos pulls its configuration from them. Config servers are standalone mongod processes, and it is started as the same way as normal mongod processes. For production purposes it is recomended to use 3 config servers and for testing puporses one config server is fine to go with. Also we start the shard servers on each of the VMs using the command given below. After shard servers,config servers and and router server is ready we copy the ssh keys to each of the machines for secure data transfer between machines.

```
mongod --configsvr --dbpath=/data/configdb
mongod --shardsvr --dbpath=/data/shardb
```
Also when we start up config servers, we do not use the --replSet option: config servers are not members of a replica set. The --configsvr option indicates to the mongod that VM is used as a config server. It is not strictly required, as all it does is change the default port mongod  which listens on to 27019 and the default data directory to /data/configdb. Once we have the shard servers and confif servers running we start the mongos router server which needs to know where the config servers are so mongos is started with the --configdb option:

```
mongos --configdb config-1:27019 > -f /var/lib/mongos.conf
```
mongos runs on port 27017. Note that it does not need a data directory. Adding shards to the mongos is easy once in the shell the command below adds three shards,

```
> sh.addShard("server-1:27017,server-2:27017,server-3:27017")
```
Mongo DB won’t distribute your data automatically until we explicitly tell which database and collections to shard. We explicitly shard the chembldb database,molecules collection and mfp_counts collection. On the mongos shell the database and collection are sharded by,

```
>db.runCommand({enablesharding: "chembldb"});
>sh.shardCollection("chembldb.molecules", { "_id": "hashed" } )
>sh.shardCollection("chembldb.mfp1_counts", { "_id": "hashed" } )  
```

### To build a database of fingerprints 
To build a chemical database of fingerprints use the db\_build.py program in the codes folder. Before using db\_build program one needs to install the Pymongo and RDKit chemical toolkit . In ubuntu machines it is installed using 
```
pip install pymongo
sudo apt-get install python-rdkit librdkit1 rdkit-data
```
db\_build.py is a command line argument program where one can input .sdf,.sdf.gz and .smi format, the pattern of fingerprint, fingerprint length and fingerprint tag name . Currently morgan type, RDKFingerprint and rdkit maccs keys are supported.

#### Usage of db_build.py
```
$ python db_build.py -h
usage: db_build.py [-h] --i I --db DB [--tag TAG] [--fpSize FPSIZE]
                   [--fpname FPNAME]
                   {morgan,rdkfp,rdmaccs} ...
python db_build.py --db moltest--i benzodiazepine.smi --tag chembl_id --fpname mfp1 morgan

## To generate morgan type fingerprints with default parameters and fpSize 1024  
$ python dbbuild.py --db moltest --i benzodiazepine.smi --tag chembl_id 
                    --fpSize 1024 --fpname mfp2 morgan    
                    
## To generate morgan type fingerprints with fpSize 1024 and radius of 4
$ python db_build.py --db moltest --i benzodiazepine.smi --tag chembl_id 
                    --fpSize 1024 --fpname mfp3 morgan --radius 4

## For this projet we used chembl.sdf which has ~1.3 million molecules
## To generate morgan type fingerprints with radius 2 into chembldb database.
python db_build.py --db chembldb --i chembl.sdf --tag chembl_id fpname mfp1 morgan

```
Once the database is built with one fingerprint addfps.py can be called to generate more fingerprints over the molecular data. This way multiple fingerprints can be generated and stored. The code for addfps.py is stored in codes folder. The code is almost similar to db\_buildy but here you need to specify the database for which you will generate the fingerprint. Currently addfps.py supports morgan type, RDKFingerprint and maccs keys fingerprints. To execute it, 
```
python addfps.py -h
python addfps.py --db moltest --fpSize 1024 --fpname mfp2 morgan
```
Now to see everything is working fine we should go to Mongos machine and login using the admin and check wether the chembldb is created along with molecules and mfp1_counts collections.
```
## Go to the mongo terminal and check chembldb is created
> show dbs
admin       (empty)
chembldb    11.15GB
local       0.078GB
> use chembldb
# Shows the mfp_1 counts are generated.
> show collections 
mfp1_counts
mfp2_counts
molecules
system.indexes
```

The query to the mongo database can be done using the Query.py script in the codes folder. In Query script is a command line program where user should give smiles string(smi) , database name to search(db) , the tag name of ids which is given for generation of original database(tag) , size of the fingerprint(fpsize) and its parameters for generation of similar type of fingerprint and fingperprint name (fpname) as given in the search database (ex: mfp1). Below shows the script how it is executed. 
MongoDB version 2.6 introduced some new aggregation features that may have better performance. Of particular interest is the $setIntersection operator, which is exactly what we need to calculate the number of bits in common between two fingerprints.

```
$ python query.py -h
usage: query.py [-h] --db DB --smi SMI [--fpSize FPSIZE] --fpname FPNAME
                   [--t T] [--tag TAG]
                   {morgan,rdkfp,rdmaccs} ...

## Example showing how the query is executed
$ python query.py --db chembldb --smi 'CC1=NN=C2N1C3=C(C=C(C=C3)Cl)C(=NC2)C4=CC=CC=C4' --fpSize 512 --fpname mfp1 --t 0.8 --tag chembl_id morgan --radius 2

fingerprints mfp1 done ...
start Aggregate .. 
response done .. 
Hits:  6 Time :  0.0532460212708
1.0: 450819
0.952380952381: 19002642
1.0: 2118
0.952380952381: 178274
0.813953488372: 12562523
0.818181818182: 21489341

```
