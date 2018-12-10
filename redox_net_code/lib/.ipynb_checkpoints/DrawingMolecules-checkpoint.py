#utility library
import re
import sys
from urllib.request import urlopen
from urllib.error import HTTPError
import os
# scientific python libraries
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
#custom thermochemistry modules
from lib.thermodynamic_constants import *
from lib.thermodynamic_new import *
from lib.parseOrcaGOHA import *
from lib.seedCalcs import read_input_file_dict_with_rxn
#other: warning options
pd.options.mode.chained_assignment = None  # default='warn'
from IPython.display import Image, SVG
from rdkit import Chem
from rdkit.Chem import AllChem,Draw,Descriptors
from rdkit.Chem.Draw import IPythonConsole, ReactionToImage
# import imolecule
import scipy
import json


def keggID_to_SMILES(keggID, ROOTDIR):
    '''Takes a KEGG ID. Returns a SMILE string.
    '''
    #print keggID
    #if a specific smiles is specified in the input_file_specs.txt file, use it!

    outFileMOL = os.path.join(ROOTDIR, 'urlTemp.mol')
    outFile = open(outFileMOL, 'wb')

    urlStr = 'http://rest.kegg.jp/get/cpd:C%05d/mol' % int(keggID[1:].strip())
    strTemp = urlopen(urlStr).read()
    outFile.write(strTemp)
    outFile.close()
    outFileSMILES = os.path.join(ROOTDIR, 'temp.smi')
    os.system('babel -imol ' + outFileMOL + ' -osmi ' + outFileSMILES)
    smiles = open(outFileSMILES,'rb').readlines()[0].strip()
    return smiles

def keggID_to_SMILES_mol(keggID):
    #Receives a KEGG ID and returns a SMILES string.
    #read in json file
    with open('Kegg_Canonical_Final.json', 'r') as fp:
        kegg_dict = json.load(fp)
    #check if we have Kegg ID
    if keggID in kegg_dict:
        smiles = kegg_dict[keggID]
    else:
        try:
            outFileMOL = 'urlTemp.mol'
            outFile = open(outFileMOL, 'wb')
            urlStr = 'http://rest.kegg.jp/get/cpd:C%05d/mol' % int(keggID[1:].strip())
            strTemp = urlopen(urlStr).read()
            outFile.write(strTemp)
            outFile.close()
            outFileSMILES = 'temp.smi'
            os.system('babel -imol ' + outFileMOL + ' -osmi ' + outFileSMILES)
            smiles = open(outFileSMILES,'rb').readlines()[0].strip()
            kegg_dict[keggID] = smiles
            #write value updated dictionary into file
            with open('Kegg_Dict.json', 'w') as fp:
                json.dump(kegg_dict, fp)
        except HTTPError:
            return 'NA'
    return smiles

def getSmiles(ID):
    dirName = '/Users/eudoraolsen/Dropbox/Free_Energy_ Work_Eudora/smilesFiles'
    smilesFile = ID+'.smiles'
    smilesStr = open(os.path.join(dirName, smilesFile), 'rb').readlines()[0].strip()
    return smilesStr

def drawSmilesRxn(rxnString):
    '''give a kegg ID and a directory
    go find the smiles strings
    '''
    subs1 = [c.strip() for c in rxnString.split('=')[0].split('+')]
    subs = [ re.search('C\d\d\d\d\d', s).group() for s in subs1]

    prods1 = [c.strip() for c in rxnString.split('=')[1].split('+')]
    prods = [ re.search('C\d\d\d\d\d', p).group() for p in prods1]

    rxnStr2 = '.'.join([keggID_to_SMILES_mol(ID) for ID in subs]) + '>>' + '.'.join([keggID_to_SMILES_mol(ID) for ID in prods])
    #print rxnStr2

    #print rxnStr2
    #return AllChem.ReactionFromSmarts(rxnStr2)
    return rxnStr2

