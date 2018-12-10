from rdkit import Chem
import json
from rdkit.Chem import AllChem,Draw,Descriptors
import copy
import pandas as pd
import re
import matplotlib.pyplot as plt
import lib.Reaction_Library as RL 


### Reduction Functions: 

def G1_function(metab_smiles):
    reactant = Chem.MolFromSmiles(metab_smiles)
    prod_list = []
    prod = RL.G1_carb_acid.RunReactants([reactant])
    if len(prod) > 0:
        for x in prod:
            smiles = Chem.MolToSmiles(x[0])
            prod_list.append(smiles)
    return prod_list


def G2_function(metab_smiles):
    reactant = Chem.MolFromSmiles(metab_smiles)
    prod_list = []
    
    #Aldehyde reduction
    prod_1 = RL.G2_ald.RunReactants([reactant])
    if len(prod_1) > 0:
        for x in prod_1:
            smiles = Chem.MolToSmiles(x[0])
            prod_list.append(smiles)
    #Ketone reduction        
    prod_2 = RL.G2_ket.RunReactants([reactant])
    if len(prod_2) > 0:
        for x in prod_2:
            smiles = Chem.MolToSmiles(x[0])
            prod_list.append(smiles)

    #6 carbon ring reactions... 
    #prod_3 = RL.G2_6c6r.RunReactants([reactant])
    #if len(prod_3) > 0:
    #    for x in prod_3:
    #        smiles = Chem.MolToSmiles(x[0])
    #        prod_list.append(smiles)
    #prod_4 = RL.G2_6c5r.RunReactants([reactant])
    #if len(prod_4) > 0:
    #    for x in prod_4:
    #        smiles = Chem.MolToSmiles(x[0])
    #        prod_list.append(smiles)
    #prod_5 = RL.G2_5c6r.RunReactants([reactant])
    #if len(prod_5) > 0:
    #    for x in prod_5:
    #        smiles = Chem.MolToSmiles(x[0])
    #        prod_list.append(smiles)
    #prod_6 = RL.G2_5c5r.RunReactants([reactant])
    #if len(prod_6) > 0:
    #    for x in prod_6:
    #        smiles = Chem.MolToSmiles(x[0])
    #        prod_list.append(smiles)
    
    return prod_list


def G3_function(metab_smiles):
    reactant = Chem.MolFromSmiles(metab_smiles)
    prod_list = []
    
    #Aldehyde reduction
    prod_1 = RL.G3_ald.RunReactants([reactant])
    if len(prod_1) > 0:
        for x in prod_1:
            smiles = Chem.MolToSmiles(x[0])
            prod_list.append(smiles)
    #Ketone reduction        
    prod_2 = RL.G3_ket.RunReactants([reactant])
    if len(prod_2) > 0:
        for x in prod_2:
            smiles = Chem.MolToSmiles(x[0])
            prod_list.append(smiles)
    return prod_list


def G4_function(metab_smiles):
    reactant = Chem.MolFromSmiles(metab_smiles)
    prod_list = []
    # Hydroxyl in the end of a molecule
    prod_1 = RL.G4_end_hydr.RunReactants([reactant])
    if len(prod_1) > 0:
        for x in prod_1:
            smiles = Chem.MolToSmiles(x[0])
            prod_list.append(smiles)
    # Hydroxyl in the middle of a molecule. 
    prod_2 = RL.G4_mid_hydr.RunReactants([reactant])
    if len(prod_2) > 0:
        for x in prod_2:
            smiles = Chem.MolToSmiles(x[0])
            prod_list.append(smiles)
    prod_3 = RL.G4_arom_hydr.RunReactants([reactant])
    if len(prod_3) > 0:
        for x in prod_3:
            smiles = Chem.MolToSmiles(x[0])
            prod_list.append(smiles)
    return prod_list


### Oxidation Functions: 


def G1_function_ox(metab_smiles):
    reactant = Chem.MolFromSmiles(metab_smiles)
    prod_list = []
    prod = RL.G1_carb_acid_ox.RunReactants([reactant])
    if len(prod) > 0:
        for x in prod:
            smiles = Chem.MolToSmiles(x[0])
            prod_list.append(smiles)
    return prod_list


def G2_function_ox(metab_smiles):
    reactant = Chem.MolFromSmiles(metab_smiles)
    prod_list = []
    
    #Aldehyde reduction
    prod_1 = RL.G2_ald_ox.RunReactants([reactant])
    if len(prod_1) > 0:
        for x in prod_1:
            smiles = Chem.MolToSmiles(x[0])
            prod_list.append(smiles)
    #Ketone reduction        
    prod_2 = RL.G2_ket_ox.RunReactants([reactant])
    if len(prod_2) > 0:
        for x in prod_2:
            smiles = Chem.MolToSmiles(x[0])
            prod_list.append(smiles)
    return prod_list


def G3_function_ox(metab_smiles):
    reactant = Chem.MolFromSmiles(metab_smiles)
    prod_list = []
    
    #Aldehyde reduction
    prod_1 = RL.G3_ald_ox.RunReactants([reactant])
    if len(prod_1) > 0:
        for x in prod_1:
            smiles = Chem.MolToSmiles(x[0])
            prod_list.append(smiles)
    #Ketone reduction        
    prod_2 = RL.G3_ket_ox.RunReactants([reactant])
    if len(prod_2) > 0:
        for x in prod_2:
            smiles = Chem.MolToSmiles(x[0])
            prod_list.append(smiles)
    return prod_list


def G4_function_ox(metab_smiles):
    reactant = Chem.MolFromSmiles(metab_smiles)
    prod_list = []
    # Hydroxyl in the end of a molecule
    prod_1 = RL.G4_end_hydr_ox.RunReactants([reactant])
    if len(prod_1) > 0:
        for x in prod_1:
            smiles = Chem.MolToSmiles(x[0])
            prod_list.append(smiles)
    # Hydroxyl in the middle of a molecule. 
    prod_2 = RL.G4_mid_hydr_ox.RunReactants([reactant])
    if len(prod_2) > 0:
        for x in prod_2:
            smiles = Chem.MolToSmiles(x[0])
            prod_list.append(smiles)
    prod_3 = RL.G4_arom_hydr_ox.RunReactants([reactant])
    if len(prod_3) > 0:
        for x in prod_3:
            smiles = Chem.MolToSmiles(x[0])
            prod_list.append(smiles)
    return prod_list


# Old Reduce metabolite: 

def reduce_metabolite_old(metabolite_smiles):
    #empty rxn_list
    rxn_list = [ ]
    #single starting substrate
    substrate_list = [ metabolite_smiles ]
    # initialize this flag to 1
    num_new_reactions = 1

    while (num_new_reactions): #keep going as long as you have new reactions:
        new_reactions = [ ] # this is the list of new reactions in the latest iteration
        substrate_list_temp = []
        #For every substrate in substrate list:
        for sub in substrate_list:
            # For every reaction_type: (looping through all G1 - G4 functions):

            #G1. Apply reaction_type to substrate
            G1_prods = G1_function(sub) #list of G1 products 
            #G2
            G2_prods = G2_function(sub) #list of G2 products
            #G3
            #G3_prods = G3_function(sub) #list of G3 products
            #G4
            G4_prods = G4_function(sub) #list of G4 products 

            #All prods: 
            # GAll_prods = list(set(G1_prods + G2_prods + G3_prods + G4_prods)) #unique new products.
            GAll_prods = list(set(G1_prods + G2_prods + G4_prods)) #unique new products. 

            
            substrate_list_temp += GAll_prods
            #2. Store resulting (uniqe) reactions in new_reactions list
            for prod in GAll_prods: 
                new_reactions.append( sub+'>>'+prod )

            rxn_list += list(set(new_reactions))

        #After reducing all current substrates, update substrate list. 
        substrate_list = list(set(substrate_list_temp))

        if len(GAll_prods)==0:
            num_new_reactions = 0
    
    return list(set(rxn_list))


def reduce_metabolite(metabolite_smiles, Master_Sub_List = [], do_G3 = False):
	rxn_list = [ ] #empty rxn_list
	substrate_list = [ metabolite_smiles ] #single starting substrate
	num_new_reactions = 1     # initialize this flag to 1

	while (num_new_reactions): #keep going as long as you have new reactions:
		new_reactions = [ ] # this is the list of new reactions in the latest iteration
		substrate_list_temp = []
		# Remove substrates that are already in the Master Substrate List: 
		substrate_list = list( set(substrate_list) - set(Master_Sub_List) )
		Master_Sub_List += substrate_list # --> shouldn't you turn this into a set? 
        
        #For every substrate in substrate list:
		for sub in substrate_list:
            # For every reaction_type: (looping through all G1 - G4 functions):
			#G1_prods = G1_function(sub) # list of G1 products 
			G2_prods = G2_function(sub) # list of G2 products
			G3_prods = G3_function(sub) # list of G4 products 
            
			# if do_G3:
			# 	G3_prods = G3_function(sub) # list of G3 products
			# 	GAll_prods = list(set(G1_prods + G2_prods + G3_prods + G4_prods))
			# else:
			GAll_prods = list(set(G2_prods + G3_prods)) #unique new products.
             
			substrate_list_temp += GAll_prods # this puts together the products of all substrates? 
            
            # Store resulting (uniqe) reactions in new_reactions list
			for prod in GAll_prods: 
				new_reactions.append( sub+'>>'+prod )
			rxn_list += list(set(new_reactions))

        #After reducing all current substrates, update substrate list. 
		substrate_list = list(set(substrate_list_temp))
		if len(substrate_list)==0:
			num_new_reactions = 0
	
	return list(set(rxn_list)), Master_Sub_List


def oxidize_metabolite_old(metabolite_smiles):
    #empty rxn_list
    rxn_list = [ ]
    #single starting substrate
    substrate_list = [ metabolite_smiles ]
    # initialize this flag to 1
    num_new_reactions = 1

    while (num_new_reactions): #keep going as long as you have new reactions:
        new_reactions = [ ] # this is the list of new reactions in the latest iteration
        substrate_list_temp = []
        #For every substrate in substrate list:
        for sub in substrate_list:
            # For every reaction_type: (looping through all G1 - G4 functions):

            #G1. Apply reaction_type to substrate
            G1_prods = G1_function_ox(sub) #list of G1 products 
            #G2
            G2_prods = G2_function_ox(sub) #list of G2 products
            #G3
            # G3_prods = G3_function_ox(sub) #list of G3 products
            #G4
            G4_prods = G4_function_ox(sub) #list of G4 products 

            #All prods: 
            GAll_prods = list(set(G1_prods + G2_prods + G4_prods)) #unique new products. 
            
            substrate_list_temp += GAll_prods
            #2. Store resulting (uniqe) reactions in new_reactions list
            for prod in GAll_prods: 
                new_reactions.append( sub+'>>'+prod )
            

            rxn_list += list(set(new_reactions))
        #print 'GAll:' + str(GAll_prods)
        #print 'new:' + str(new_reactions)
        #for rxn in new_reactions:
        #    display(AllChem.ReactionFromSmarts(rxn))
        #After reducing all current substrates, update substrate list. 
        substrate_list = list(set(substrate_list_temp))
        #print substrate_list
        if len(GAll_prods)==0:
            num_new_reactions = 0
    return list(set(rxn_list))


def oxidize_metabolite(metabolite_smiles, Master_Sub_List = None, rxn_network_file = None, rxn_limit = None):
    # empty rxn_list
    rxn_list = [ ]
    # single starting substrate
    substrate_list = [ metabolite_smiles ]
    # initialize this flag to 1
    num_new_reactions = 1
    # num reaction status
    len_status = 0
    
    while (num_new_reactions): #keep going as long as you have new reactions:
        new_reactions = [ ] # this is the list of new reactions in the latest iteration
        substrate_list_temp = []
        
        # Remove substrates that are already in the Master Substrate List: 
        # print 'Substrate List before filtering:', substrate_list
        substrate_list = list(set(substrate_list) - set(Master_Sub_List))
        # print 'Substrate List after filtering:', substrate_list
        Master_Sub_List += substrate_list
        # print 'Updated MasterSubList', Master_Sub_List, '\n'
        
        #For every substrate in substrate list, perform all 4 types of oxidations: 
        for sub in substrate_list:
            # For every reaction_type: (looping through all G1 - G4 functions):

            #G1. Apply reaction_type to substrate
            # G1_prods = G1_function_ox(sub) #list of G1 products 
            #G2
            G2_prods = G2_function_ox(sub) #list of G2 products
            #G3
            G3_prods = G3_function_ox(sub) #list of G3 products
            #G4
            # G4_prods = G4_function_ox(sub) #list of G4 products 

            # These are all the products from all reactions acting on current substrate
            GAll_prods = list(set(G2_prods + G3_prods)) #unique new products. 
            
            substrate_list_temp += GAll_prods
            #2. Store resulting (uniqe) reactions in new_reactions list
            for prod in GAll_prods: 
                new_reactions.append( sub+'>>'+prod )

            # Update the list of reactions
            rxn_list += new_reactions
            # Get rid of duplicates --> 
            rxn_list = list(set(rxn_list))
              
        #After reducing all current substrates, update substrate list. 
        substrate_list = list(set(substrate_list_temp))
        #print substrate_list
        if len(substrate_list)==0:
            num_new_reactions = 0
            
        # If you only want to run up to a certain number of reactions: 
        if rxn_limit:
                if len(rxn_list) > rxn_limit:
                    return list(set(rxn_list)), Master_Sub_List 
                
        if len(rxn_list) > len_status:
            # print 'number of reactions', len(rxn_list)
            len_status = len(rxn_list) + 1000
            
    return list(set(rxn_list)), Master_Sub_List


def visualize_reaction_old(metabolite_smiles):
    rxn_list = reduce_metabolite_old(metabolite_smiles)
    rxn_list.sort(key = len, reverse = True)
    print("Starting metabolite: " + metabolite_smiles)
    print("Number of reactions: " + str(len(rxn_list)))
    for rxn in rxn_list:
        print(rxn)
        display(AllChem.ReactionFromSmarts(rxn))

        
def visualize_reaction_new(metabolite_smiles):
    rxn_list = reduce_metabolite(metabolite_smiles)
    rxn_list.sort(key = len, reverse = True)
    print("Starting metabolite: " + metabolite_smiles)
    print("Number of reactions: " + str(len(rxn_list)))
    for rxn in rxn_list:
        print(rxn)
        display(AllChem.ReactionFromSmarts(rxn))

        
def visualize_reaction_ox_old(metabolite_smiles):
    rxn_list = oxidize_metabolite_old(metabolite_smiles)
    rxn_list.sort(key = len)
    print("Starting metabolite: " + metabolite_smiles)
    print("Number of reactions: " + str(len(rxn_list)))
    for rxn in rxn_list:
        print(rxn)
        display(AllChem.ReactionFromSmarts(rxn))
