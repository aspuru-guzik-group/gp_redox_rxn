from rdkit import Chem
from rdkit.Chem import AllChem,Draw,Descriptors
import json
import pandas as pd
from lib.DrawingMolecules import drawSmilesRxn
import numpy as np
from rdkit.Chem.Draw import IPythonConsole, ReactionToImage
import os
import re
from lib.Reaction_Library_Mother import *


def match_with_smiles(rxn_dict, kegg2smiles):

    match_dict = {}
    
    for key in rxn_dict.keys():
        #if not empty list
        if len(rxn_dict[key]):
            match_list = []
            for prod in rxn_dict[key]:
                smiles = Chem.MolToSmiles(prod[0], isomericSmiles = True)
                #SUBSTITUTE with comparison against KEGG SMILES
                match_prod = [k for k in kegg2smiles.keys() if kegg2smiles[k].lower() == smiles.lower()]
                match_list += match_prod
            if len(match_list):
                match_dict[key] = match_list

    return match_dict

def match_with_substructure(rxn_dict, match_dict, kegg2smiles):
    #find all substrates that generated a product. 
    key_with_prod = sort([key for key in rxn_dict.keys() if len(rxn_dict[key])])
    #find substrates that did not match using the smiles matching function
    #These are potentially false negatives. 
    no_match_keys = set(key_with_prod) - set(match_dict.keys())

    match_SS_dict = {}
    for key in no_match_keys: #for all the substrates who didn't match using smiles
        match_list = []
        for prod in rxn_dict[key]: #for all of their products. 
            m = prod[0]
            #do substructure matching. 
            match_prod = [k for k in kegg2smiles.keys() if AreTwoMoleculesSame(m, Chem.MolFromSmiles(kegg2smiles[k]))]
            match_list += match_prod
        #if list not empty
        if len(match_list):
            match_SS_dict[key] = match_list
    return match_SS_dict

def consider_delta_stoich(match_dict, match_SS_dict, kegg2smiles):
    #Write down how this works:
    new_match_SS_dict = {}
    for key in match_SS_dict.keys(): #the rxn matches using rdkit substructres (no chirality)
        match_list = []
        for product in match_SS_dict[key]:
            #we compare the number of stereocenters of the substrate and its product.  
            if delta_stereochem(kegg2smiles, match_SS_dict, key, product) != 0:
                match_list += [product]
        if len(match_list):       
            new_match_SS_dict[key] = match_list
    
    match_dict_final = dict(match_dict.items() + new_match_SS_dict.items())
    return match_dict_final

#This helps you get rid of false positives. How? 
def delta_stereochem(ks, match_ald_SS_dict, key1, key2):
    return len(re.findall('\@+',ks[key1])) - len(re.findall('\@+',ks[key2]))


def AreTwoMoleculesSame(mol_obj1, mol_obj2, use_chirality = False):
    return( mol_obj1.HasSubstructMatch(mol_obj2, useChirality=use_chirality) 
           and mol_obj2.HasSubstructMatch(mol_obj1, useChirality=use_chirality))

def get_subs(x):
    return x.Reaction.split(' = ')[0]

def get_prod(x):
    return x.Reaction.split(' = ')[1]

def rxn_list_from_dict(match_dict):
    rxn_list = []
    for key in match_dict:
        for x in match_dict[key]:
            rxn = key + ' = ' + x
            rxn_list.append(rxn)
    rxn_list = list(set(rxn_list))
    return rxn_list

def draw_subset(rxn_list):
    for rxn in rxn_list:
        print(rxn)
        smarts = drawSmilesRxn(rxn)
        print(smarts)
        display(AllChem.ReactionFromSmarts(str(smarts)))

def categorize_rxn(row):
    
    from lib.Reaction_Library_Mother import G2_function, G3_function

    rxn = row.rxn
    sub = rxn.split('>>')[0]
    prod = rxn.split('>>')[1]
    # Apply reaction_type to substrate
    #G1_prods = G1_function(sub) #list of G1 products 
    G2_prods = G2_function(sub) #list of G2 products
    G3_prods = G3_function(sub) #list of G3 products
    #G4_prods = G4_function(sub) #list of G4 products 
    
    # if len(G1_prods) > 0:
    #     for x in G1_prods:
    #         if x == prod:
    #             return 1
    if len(G2_prods) > 0:
        for x in G2_prods:
            if x == prod:
                return 2
    if len(G3_prods) > 0:
        for x in G3_prods:
            if x == prod:
                return 3
    # if len(G4_prods) > 0:
    #     for x in G4_prods:
    #         if x == prod:
    #             return 4

def categorize_rxn_ox(row):
    rxn = row.rxn
    sub = rxn.split('>>')[0]
    prod = rxn.split('>>')[1]
    # Apply reaction_type to substrate
    # G1_prods = G1_function_ox(sub) #list of G1 products 
    G2_prods = G2_function_ox(sub) #list of G2 products
    G3_prods = G3_function_ox(sub) #list of G3 products
    # G4_prods = G4_function_ox(sub) #list of G4 products 
    
    # if len(G1_prods) > 0:
    #     for x in G1_prods:
    #         if x == prod:
    #             return 1
    if len(G2_prods) > 0:
        for x in G2_prods:
            if x == prod:
                return 2
    if len(G3_prods) > 0:
        for x in G3_prods:
            if x == prod:
                return 3
    # if len(G4_prods) > 0:
    #     for x in G4_prods:
    #         if x == prod:
    #             return 4


def count_carbons(x):
    smiles = x.smiles
    new_smiles = smiles.replace('Cl', '')
    numC = new_smiles.count('C') + smiles.count('c')
    return numC

def reverse_reaction(rxn):
    sub = rxn.split('>>')[0]
    prod = rxn.split('>>')[1]
    new_rxn = prod + '>>' + sub
    return new_rxn

def update_rxn_file(rxn_network_file, rxn_list):
    with open(rxn_network_file, "w") as outfile:
        for rxn in rxn_list:
            outfile.write(rxn +'\n')

#Reaction Strings

#Reduction

#G1 Rxn:

G1_carb_acid = AllChem.ReactionFromSmarts('[CX3:1](=O)[OX2H1]>>[CX3H1:1](=O)')

#G2 Rxn:

#aldehydes
G2_ald = AllChem.ReactionFromSmarts('[CX3H1:2](=O)[#6:1] >>[#6:1][CX4H2:2][OX2H1]')
#Ketones
G2_ket = AllChem.ReactionFromSmarts('[#6:1][CX3:2](=O)[#6:3] >>[#6:1][CX4H1:2]([#6:3])[OX2H1]')
#6 carbon, 6 membered rings
c6r6_smarts = '[C:6][C:5]1[O:10][C:1][C:2][C:3][C:4]1'
c6o_smarts = '[C:6][C:5]([O:10])[C:4][C:3][C:2][C:1]'
sugar_6c6r_rxn = c6r6_smarts + '>>' + c6o_smarts
G2_6c6r = AllChem.ReactionFromSmarts(sugar_6c6r_rxn)
#6 carbon, 5 membered rings
c6r5_smarts = '[C:6][C:5]1[O:10][C:2](O)([C:1])[C:3][C:4]1'
c6o_smarts = '[C:6][C:5]([O:10])[C:4][C:3][C:2](O)[C:1]'
sugar_6c5r_rxn = c6r5_smarts + '>>' + c6o_smarts
G2_6c5r = AllChem.ReactionFromSmarts(sugar_6c5r_rxn)
#5 carbon, 6 membered rings
c5r6_smarts = '[C:5]1[O:10][C:1](O)[C:2][C:3][C:4]1'
c5o_smarts = '[C:5]([O:10])[C:4][C:3][C:2][C:1](O)'
sugar_5c6r_rxn = c5r6_smarts + '>>' + c5o_smarts
G2_5c6r = AllChem.ReactionFromSmarts(sugar_5c6r_rxn)
#5 carbon, 5 membered rings
c5r5_smarts = '[C:5]1[O:10][C:2]([C:1])[C:3][C:4]1'
c5o_smarts = '[C:5]([O:10])[C:4][C:3][C:2][C:1]'
sugar_5c5r_rxn = c5r5_smarts + '>>' + c5o_smarts
G2_5c5r = AllChem.ReactionFromSmarts(sugar_5c5r_rxn)

#G3 Rxn:

#aldehydes
#G3_ald = AllChem.ReactionFromSmarts('[CX3H1:2](=O)[#6:1]>>[#6:1][CX4H2:2][NX3;H2,H1;!$(NC=O)]')
G3_ald = AllChem.ReactionFromSmarts('[#6:1][CX3H1:2](=O)>>[#6:1][CX4H2:2][NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])]')

#Ketones
G3_ket = AllChem.ReactionFromSmarts('[#6:1][CX3:2](=O)[#6:3]>>[C:1][CX4H1:2]([C:3])[NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])]')

#G4 Rxn:

# end hydroxyl
G4_end_hydr = AllChem.ReactionFromSmarts('[CX4H2:2][OX2H1]>>[CX4H3:2]')
# mid hydroxyl
G4_mid_hydr = AllChem.ReactionFromSmarts('[#6:1][#6H1:2]([#6:3])[OX2H1]>>[#6:1][#6H1:2][#6:3]')
# arom hydroxyl
G4_arom_hydr = AllChem.ReactionFromSmarts('[c:1][OX2H1]>>[cH1:1]')

#Oxidation

#G1 Rxn:

G1_carb_acid_ox = AllChem.ReactionFromSmarts('[CX3H1:1](=O)>>[CX3:1](=O)[OX2H1]')

#G2 Rxn:

#aldehydes
G2_ald_ox = AllChem.ReactionFromSmarts('[#6:1][CX4H2:2][OX2H1] >> [CX3H1:2](=O)[#6:1]')
#Ketones
G2_ket_ox = AllChem.ReactionFromSmarts('[#6:1][CX4H1:2]([#6:3])[OX2H1] >> [#6:1][CX3H0:2](=O)[#6:3]')

#G3 Rxn:

#aldehydes
#G3_ald_ox = AllChem.ReactionFromSmarts('[#6:1][CX4H2:2][NX3;H2,H1;!$(NC=O)]>>[CX3H1:2](=O)[#6:1]')
G3_ald_ox = AllChem.ReactionFromSmarts('[#6:1][CX4H2:2][NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])]>>[#6:1][CX3H1:2](=O)')

#Ketones
#G3_ket_ox = AllChem.ReactionFromSmarts('[C:1][CX4H1:2]([C:3])[NX3;H2,H1;!$(NC=O)]>>[#6:1][CX3H0:2](=O)[#6:3]')
G3_ket_ox = AllChem.ReactionFromSmarts('[C:1][CX4H1:2]([C:3])[NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])]>>[#6:1][CX3H0:2](=O)[#6:3]')

#G4 Rxn:

# end hydroxyl
G4_end_hydr_ox = AllChem.ReactionFromSmarts('[CX4H3:2]>>[CX4H2:2][OX2H1]')
# mid hydroxyl
G4_mid_hydr_ox = AllChem.ReactionFromSmarts('[#6:1][#6H2:2][#6:3]>>[#6:1][#6H1:2]([#6:3])[OX2H1]')
# arom hydroxyl
G4_arom_hydr_ox = AllChem.ReactionFromSmarts('[cH1:1]>>[c:1][OX2H1]')
