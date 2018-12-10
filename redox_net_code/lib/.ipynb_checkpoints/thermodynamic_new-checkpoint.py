import pandas as pd
import sys
import numpy as np
from lib.thermodynamic_constants import *
import os
sys.path.append(os.path.join(os.getcwd(), 'Modules'))
import lib.util
from scipy.stats import scoreatpercentile
import matplotlib.pyplot as plt
import random

#global variable
#ATP_ADP_GFE = 208.44  #kcal/mol
#NAD_NADH_GFE = 14.96  #kcal/mol
#NADP_NADPH_GFE = 15.63 #kcal/mol

ATP_ADP_GFE = 0  #kcal/mol
NAD_NADH_GFE = 0  #kcal/mol
NADP_NADPH_GFE = 0 #kcal/mol

def H_rotational(x, T):
    #T = 298.15
    Hrot = 1.5 * R_kCal * T
    return Hrot


def S_rotational(x, T):
    T = 298.15
    qrotRoomTemp = x.QRot
    qrot = qrotRoomTemp * (T / 298.15) ** (1.5)
    Srot = R_kCal * (np.log(qrot) + 1.5)
    return Srot


def S_translational(x, T):
    # Translational
    mass = x.Mass
    Strans = R_kCal * (np.log((2.0 * np.pi * amutokg * mass * kb * T)
                              ** (1.5) * kb * T / 101325.0 / (hJ ** 3.0)) + 2.5)
    return Strans


def H_translational(x, T):
    # Translational
    Htrans = 1.5 * R_kCal * T
    return Htrans


def H_vibrational(x, T):
    # Vibrational
    freqs = x.Freq
    hVibArray = []
    ZPVE = np.sum(freqs)
    ZPVE = 0.5 * hJ * float(100) * speedc * na * ZPVE / (kcaljou)
    hvib = 0
    for i in xrange(len(freqs)):
        freqsi = freqs[i]
        ui = freqsi * float(100) * speedc * hJ / (kb * T)
        expui = np.exp(ui)
        hvib = hvib + freqsi * \
            float(100) * speedc * np.exp(-ui) / (1.0 - np.exp(-ui))
    hvib = hvib * hJ * na / kcaljou + ZPVE
    hVibArray.append(hvib)
    Hvib = sum(hVibArray)

    return Hvib


def S_vibrational(x, T):
    # Vibrational
    freqs = x.Freq
    sVibArray = []
    ZPVE = sum(freqs)
    ZPVE = 0.5 * hJ * float(100) * speedc * na * ZPVE / (kcaljou)
    svib = 0
    for i in xrange(len(freqs)):
        freqsi = freqs[i]
        ui = freqsi * float(100) * speedc * hJ / (kb * T)
        expui = np.exp(ui)
        svib = svib + ui / (expui - 1.0) - log(1.0 - 1.0 / expui)
    svib = float(1000) * svib * R_kCal
    sVibArray.append(svib)
    Svib = sum(sVibArray)

    return Svib


def H_total(x, T, SPE = 1, vib = True):
    '''
    Inputs:
        SPE --> Flag, determines whether to use alternative 
        single point energy calculation. 
            1 --> alternative (e.g. DLPNO), SPE calculation
            0 --> SPE taken from e.g. B3LYP energy minimization
    '''
    # Total
    if SPE=='DLPNO': #SPE = 1,  Take SPE from alternative electronic structure method
        Htot = x.HTrans + x.HVib + x.HRot + x.DLPNO_SPE + R_kCal * T
    elif SPE=='OMEGA':
        Htot = x.HTrans + x.HVib + x.HRot + x.OMEGA_SPE + R_kCal * T
    elif SPE=='SPE':
        Htot = x.HTrans + x.HVib + x.HRot + x.SPE + R_kCal * T
    else: #SPE = 0, Take SPE from same energy minimization, Harmonic Analysis method. 
        Htot = x.HTrans + x.HVib + x.HRot + x.EMin + R_kCal * T 
        if not vib:
            Htot = x.HTrans + x.HRot + x.EMin + R_kCal * T

    return Htot


def S_total(x, T, vib = True):
    # Total
    Stot = x.STrans + x.SVib + x.SRot
    if not vib:
        Stot = x.STrans + x.SRot    
    return Stot


def G_total(x, T):
    # Total
    Htot = x.HTot
    Stot = x.STot
    G = Htot - T * Stot
    return G


def AlbertyTransform(x, T, pH, mu):
    G_transf = array_transform_ab_initio(np.array([x.GTot]), np.array([x.numberH]),
                                         np.array([x.charge]), default_nMg, pH,
                                         default_pMg, default_I, T, mu)
    return G_transf


def combine_conformers(x, T, Alberty):
    #T = 298.15
    if Alberty: 
        G_array = np.array(x['G_app'].values)
        G_combo = -R_kCal * T * util.log_sum_exp(G_array / (-R_kCal * T))
    
    else: 
        G_array = x['GTot'].values
        G_combo = -R_kCal * T * util.log_sum_exp(G_array / (-R_kCal * T))

    return pd.Series([G_combo], index=["G_avg"])


def combine_protonation(x, T):
    #T = 298.15
    G_array = x['G_avg'].values
    G_combo = -R_kCal * T * util.log_sum_exp(G_array / (-R_kCal * T))
    return pd.Series([G_combo], index=["G_f"])


def G_single_molecule(df, molStr, T, pH, mu, SPE = 'EMin', vib = True):
    '''
    WARNING: you are using merged as a global variable
    inputs: 
        --> mu --> proton solvation potential
        --> SPE --> Flag, determines whether to use alternative 
        single point energy calculation. 
    '''
    #Filter data by molecule string. 
    df_mol = df[df.Molecule == molStr]
    # Translational
    df_mol['HTrans'] = df_mol.apply(H_translational, axis=1, args=[T])
    df_mol['STrans'] = df_mol.apply(S_translational, axis=1, args=[T])
    # Rotational
    df_mol['HRot'] = df_mol.apply(H_rotational, axis=1, args=[T])
    df_mol['SRot'] = df_mol.apply(S_rotational, axis=1, args=[T])
    # Vibrational
    df_mol['HVib'] = df_mol.apply(H_vibrational, axis=1, args=[T])
    df_mol['SVib'] = df_mol.apply(S_rotational, axis=1, args=[T])
    # Total
    df_mol['HTot'] = df_mol.apply(H_total, axis=1, args=[T, SPE, vib])
    df_mol['STot'] = df_mol.apply(S_total, axis=1, args=[T, vib])
    df_mol['GTot'] = df_mol.apply(G_total, axis=1, args=[T])
    # Transformed
    df_mol['G_app'] = df_mol.apply(AlbertyTransform, axis=1, args=[T, pH, mu])
    return df_mol

def average_conformers(df_mol,T, Alberty = True):
    # Average conformers
    grouped1 = df_mol.groupby(["Molecule", "charge"])
    grouped1 = grouped1.apply(combine_conformers, T, Alberty)
    grouped1 = grouped1.reset_index()
    return grouped1

def average_protonation(grouped1,T):
    # Average protonation states
    grouped2 = grouped1.groupby('Molecule')
    grouped2 = grouped2.apply(combine_protonation, (T))
    grouped2 = grouped2.reset_index()
    return grouped2

def filter_data_IQR(df, mol, z, Alberty = True, numStd= 1, iqr_factor = 1.349):   
    data = df[(df['charge'] == z) & 
                (df['Molecule'] == mol)]
    
    if Alberty: #Use transformed Gibbs Energies. 
        #"Filter using the Alberty transformed Gibbs energies"
        tempMedian = data.median().G_app
        tempIQR = scoreatpercentile(data.G_app.values,75) - scoreatpercentile(data.G_app.values,25)
        tempStd = iqr_factor*tempIQR
        #and filter
        dataFiltered = data[(np.abs(data.G_app-tempMedian)
                             <=numStd*tempStd) ]
        #append to empty data frame

    else: #Use untransformed Gibbs Energies
        #"Filter using the Alberty transformed Gibbs energies"
        tempMedian = data.median().GTot
        tempIQR = scoreatpercentile(data.GTot.values,75) - scoreatpercentile(data.GTot.values,25)
        tempStd = iqr_factor*tempIQR
        #and filter
        dataFiltered = data[(np.abs(data.GTot-tempMedian)
                             <=numStd*tempStd) ]
        #append to empty data frame
    
    return dataFiltered

def sample_conformers(data, sampleSize):
    '''
    input: df --> dataframe with thermodynamic data 
    for single molecule/charge pair
    '''
    rows = random.sample(data.index, sampleSize)
    sampledData = data.ix[rows]
    return sampledData

def get_substrates(rxnList):
    subs = [s.strip() for s in (rxnList.split('=')[0]).split('+')]
    return subs

def get_products(rxnList):
    prod = [p.strip() for p in (rxnList.split('=')[1]).split('+')]
    return prod

def Gf_reactants(reacts, df, T, pH, mu, sampleSize, SPE = 'EMin', vib = True, Alberty = True):
    G_f_r = 0
    for r in reacts:
        #Get raw thermodynamic data for this reactant
        df_mol = G_single_molecule(df, r, T, pH, mu, SPE, vib)
        #Get all protononation states
        zList = set(df_mol[df_mol['Molecule']==r].charge.values)
        #Data frame where we'll store dG_f of each protonation state
        df_dGProt_all = pd.DataFrame()
        for z in zList: #For each protonation state
            #Filter data
            dataFiltered = filter_data_IQR(df_mol,r,z,Alberty) #option: with or without Alberty transform
            #subsample data
            dataSampled = sample_conformers(dataFiltered, sampleSize)
            #Average over conformers and store in df_dGProt_all
            df_dGProt_temp = average_conformers(dataSampled, T, Alberty)
            df_dGProt_all = df_dGProt_all.append(df_dGProt_temp)
        #average over all protonation states for this reactant
        df_dGMol = average_protonation(df_dGProt_all, T)
        #update substrate formation energy
        G_f_r += df_dGMol.G_f.values[0]
    return G_f_r

def getGr(rxn_str, df, T, pH, mu, sampleSize, SPE = 'EMin', vib = True, Alberty = True): 
    '''
    rxn_str --> The reaction string, with numW
    df --> The data frame with raw electronic structure data. 
    T, ph, mu, sampleSize, SPE --> self explanatoryd
    '''

    subs = get_substrates(rxn_str)
    prods = get_products(rxn_str)
    #adding waters
    numWats = subs[0].split('_')[1]
    G_f_other = 0

    #ATP/ADP reaction
    if 'C00002' in rxn_str and 'C00008' in rxn_str:
        if 'C00002_' + numWats in subs:
            G_f_other = ATP_ADP_GFE
            subs.remove('C00002_' + numWats)
            prods.remove('C00008_' + numWats)
        else:
            G_f_other = -ATP_ADP_GFE
            subs.remove('C00008_' + numWats)
            prods.remove('C00002_' + numWats)
        

    #NADPH oxidoreductase reaction
    if 'C00005' in rxn_str and 'C00006' in rxn_str:
        if 'C00005_' + numWats in subs:
            G_f_other = -NADP_NADPH_GFE
            subs.remove('C00005_' + numWats)
            prods.remove('C00006_' + numWats)
        else:
            G_f_other = NAD_NADH_GFE
            subs.remove('C00006_' + numWats)
            prods.remove('C00005_' + numWats)
        
    
    #NADH oxidoreductase reaction
    if 'C00003' in rxn_str and 'C00004' in rxn_str:
        if 'C00003_' + numWats in subs:
            G_f_other = NAD_NADH_GFE
            subs.remove('C00003_' + numWats)
            prods.remove('C00004_' + numWats)
        else:
            G_f_other = -NAD_NADH_GFE
            subs.remove('C00004_' + numWats)
            prods.remove('C00003_' + numWats)
        


    G_f_s = Gf_reactants(subs, df, T, pH, mu, sampleSize, SPE, vib, Alberty)
    G_f_p = Gf_reactants(prods, df, T, pH, mu, sampleSize, SPE, vib, Alberty)
    G_r = G_f_p - G_f_s + G_f_other
    return (G_f_s, G_f_p, G_r)

def getGr_pandas(x, rxn_str, dfAll, mu, sampleSize, numIter, SPE = 'EMin', vib = True):
    '''
    x --> data frame with experimental values for a single reaction
    rxn_str --> The reaction string, with numW
    df --> The data frame with raw electronic structure data. 
    T, ph, mu, sampleSize, SPE --> self explanatoryd
    '''
    #Set some variables
    Gr_array = np.zeros(numIter)
    #Get values from data frame with experimental data
    T = x.Temp #temperature
    pH = x.PH #PH
    GR_exp = x.GR_exp #Experimental Gibbs reaction energy. 
    
    #Perform G_R estimate many times, with random subsampling
    print(rxn_str, T, pH)
    for n in range(numIter):
        Gfs, Gfp, Gr = getGr(rxn_str, dfAll, T, pH, mu,
                             sampleSize, SPE, vib) 
        Gr_array[n] = Gr
    
    #Get summary statistics of the estimates (deviations from experiment, etc.)
    Gr_median = np.median(Gr_array)
    Gr_std = np.std(Gr_array)
    Gr_err = Gr_median - GR_exp
    return Gr_err

def get_data_s(df, s):
    df.head()
    data = df[df['s'] == s].Gr.values
    return data

def box_plot_with_points_rxn(df):
    fig, ax = plt.subplots(figsize = (8,5))
    temp = df[['Arxn','s','Gr']].boxplot(by=['s'],ax=ax, return_type = 'axes')
    
    xticks = temp['Gr'].xaxis.get_ticklabels()
    counter = 1
    for tick in xticks:
        s = int(float(tick.get_text()))
        data = get_data_s(df, s)
        plt.plot(np.repeat(counter, len(data)), data, 'ok')
        counter += 1
    
    plt.suptitle('')
    fig.autofmt_xdate()
    plt.ylabel('$G_r$ Estimate')
    plt.xlabel('Sample Size')

def G_mols_in_rxn(rxn_str, dfAll, T, pH, mu, SPE = 'EMin', vib = True):
    df_mol = pd.DataFrame()
    subs = get_substrates(rxn_str)
    prods = get_products(rxn_str)
    for s in subs:
        df_mol_temp = G_single_molecule(dfAll, s, T, pH, mu, SPE, vib)
        df_mol = df_mol.append(df_mol_temp)
    for p in prods:
        df_mol_temp = G_single_molecule(dfAll, p, T, pH, mu, SPE, vib)
        df_mol = df_mol.append(df_mol_temp)
    return df_mol

#Convert equilibrium constant to delta G.
def k2DeltaG(T,K):
    return -R_kCal*T*np.log(K)

def pandask2DeltaG(row):
    return -R_kCal*row.Temp*np.log(row.Kprime)

def string_to_stringwH20(rxnString, numW):
    substrates  = rxnString.split('=')[0].split('+')    
    products = rxnString.split('=')[1].split('+')
    if len(substrates) == len(products):
            rxnWat = [(sub.strip(), str(numW)) for sub in substrates]
            prodWat = [(prod.strip(), str(numW)) for prod in products]
            molListWat = rxnWat + prodWat
    elif len(substrates) == 1:
            rxnWat = [(sub.strip(), str(2*numW)) for sub in substrates]
            prodWat = [(prod.strip(), str(numW)) for prod in products]
            molListWat = rxnWat + prodWat
    else:
            rxnWat = [(sub.strip(), str(numW)) for sub in substrates]
            prodWat = [(prod.strip(), str(2*numW)) for prod in products]
            molListWat = rxnWat + prodWat
    
    prod_str = ' + '.join([ mol+'_'+wat+'w' for (mol, wat) in prodWat])
    rxn_str = ' + '.join([ mol+'_'+wat+'w' for (mol, wat) in rxnWat])
    return ' = '.join([rxn_str, prod_str])
