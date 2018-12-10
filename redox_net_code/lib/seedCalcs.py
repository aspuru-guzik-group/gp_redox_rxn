#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ==== Libraries ====
import os


# ==== Functions ====
def replace_allnew(astr, bstr, dict):
    with open(astr, 'r') as afile:
        with open(bstr, 'w') as bfile:
            data = afile.read()
            for key, value in dict.items():
                data = data.replace(key, str(value))
            bfile.write(data)
    return data


def getWaters(afile):
    nWaters = ""
    with open(afile, 'r') as bfile:
        nWaters = int(bfile.read())
    return nWaters


def getMultCharge(afile):

    with open(afile, 'r') as afile:
        text = afile.read()
        startloc = text.rfind("$molecule\n") + len("$molecule\n")
        good = [int(d) for d in text[startloc:].split("\n")[0].split()]
    return good[0], good[1]

def read_input_file_dict(start_dir, paramFile = 'input_file_specs.txt'):
    calcDict = {}
    #What happens if we want a single point energy and need a different input file? 
    parameter_file = os.path.join(start_dir, paramFile)

    with open(parameter_file, 'rb') as inpFile:
        for line in inpFile:
            key_string = '{#'+line.split(' = ')[0].strip()+'}'
            value_string = line.split(' = ')[-1].strip()
            calcDict[key_string] = value_string
    if calcDict['{#IMP_SOLV}'] == 'SMD':
        calcDict['{#IMP_SOLV}'] = "%cosmo smd true \nend"
        calcDict['{#IMP_SOLV_NAME}'] = "SMD"
    else: 
        calcDict['{#IMP_SOLV}'] = "%cosmo epsilon 80.4\nrefrac 1.33\nend"
        calcDict['{#IMP_SOLV_NAME}'] = "COSMO"

    return calcDict

def read_input_file_dict_with_rxn(paramFile = 'input_file_specs.txt'):
    calcDict = {}
    #What happens if we want a single point energy and need a different input file? 
    rxn_flag = 0
    param_flag = 0
    calcDict['WORKDIR'] = '/'.join(paramFile.split('/')[:-1])
    with open(paramFile, 'rb') as inpFile:
        for line in inpFile:
            if '$RXN_LIST' in line:
                rxn_flag = 1
                calcDict['RXN_LIST'] = []
                continue
            if '$PARAMETERS' in line: 
                rxn_flag = 0
                param_flag = 1 
                continue
            if rxn_flag:
                if 'END' in line:
                    rxn_flag = 0
                    continue
                else:
                    calcDict['RXN_LIST'].append(line.strip())
            if param_flag:
                if "END" in line:
                    param_flag = 0
                else:
                    key_string = '{#'+line.split(' = ')[0].strip()+'}'
                    value_string = line.split(' = ')[-1].strip()
                    calcDict[key_string] = value_string

    #Implicit solvation details
    if calcDict['{#IMP_SOLV}'] == 'SMD':
        calcDict['{#IMP_SOLV}'] = "%cosmo smd true \nend"
        calcDict['{#IMP_SOLV_NAME}'] = "SMD"
    else: 
        calcDict['{#IMP_SOLV}'] = "%cosmo epsilon 80.4\nrefrac 1.33\nend"
        calcDict['{#IMP_SOLV_NAME}'] = "COSMO"

    #Dispersion correction details: 
    if calcDict['{#DISPERSION}'] == 'NONE':
        calcDict['{#DISPERSION}'] = ''

    #Basis set and input template details
    if calcDict['{#BASIS}'] == 'P6_31G':
        calcDict['{#BASIS}'] = '_6_31G'
        calcDict['{#TEMPLATE_INP}'] = 'template_input_pople.inp'
    elif calcDict['{#BASIS}'] == 'def2_SVP':
        calcDict['{#BASIS}'] = 'def2-SVP'
    elif calcDict['{#BASIS}'] == 'def2_TZVP':
        calcDict['{#BASIS}'] = 'def2-TZVP'       

    return calcDict

def read_input_file_dict_SPE(paramFile = 'input_file_specs.txt'):
    calcDict = {}
    #What happens if we want a single point energy and need a different input file? 
    rxn_flag = 0
    param_flag = 0
    calcDict['WORKDIR'] = '/'.join(paramFile.split('/')[:-1])
    with open(paramFile, 'rb') as inpFile:
        for line in inpFile:
            if '$RXN_LIST' in line:
                rxn_flag = 1
                calcDict['RXN_LIST'] = []
                continue
            if '$PARAMETERS' in line: 
                rxn_flag = 0
                param_flag = 1 
                continue
            if rxn_flag:
                if 'END' in line:
                    rxn_flag = 0
                    continue
                else:
                    calcDict['RXN_LIST'].append(line.strip())
            if param_flag:
                if "END" in line:
                    param_flag = 0
                else:
                    key_string = '{#'+line.split(' = ')[0].strip()+'}'
                    value_string = line.split(' = ')[-1].strip()
                    calcDict[key_string] = value_string

    #Implicit solvation details
    if calcDict['{#IMP_SOLV}'] == 'SMD':
        calcDict['{#IMP_SOLV}'] = "%cosmo smd true \nend\n"
        calcDict['{#IMP_SOLV_NAME}'] = "SMD"
    elif calcDict['{#IMP_SOLV}'] == 'NONE':
        calcDict['{#IMP_SOLV}'] = ''
        calcDict['{#IMP_SOLV_NAME}'] = "NONE"
    else: 
        calcDict['{#IMP_SOLV}'] = "%cosmo epsilon 80.4\nrefrac 1.33\nend\n"
        calcDict['{#IMP_SOLV_NAME}'] = "COSMO"

    #Dispersion correction details: 
    if calcDict['{#DISPERSION}'] == 'NONE':
        calcDict['{#DISPERSION}'] = ''

    #Basis set and input template details
    if calcDict['{#BASIS}'] == 'P6_31G':
        calcDict['{#BASIS}'] = '_6_31G'
    elif calcDict['{#BASIS}'] == 'def2_SVP':
        calcDict['{#BASIS}'] = 'def2-SVP'
    elif calcDict['{#BASIS}'] == 'def2_TZVP':
        calcDict['{#BASIS}'] = 'def2-TZVP'
    elif calcDict['{#BASIS}'] == 'cc_pVTZ_cc_pVTZ_C':
        calcDict['{#BASIS}'] = 'cc-pVTZ cc-pVTZ/C'
    elif calcDict['{#BASIS}'] == 'cc_pVDZ_cc_pVDZ_C':
        calcDict['{#BASIS}'] = 'cc-pVDZ cc-pVDZ/C'     

    if calcDict['{#METHOD}'] == 'DLPNO_CCSDT':
        calcDict['{#METHOD}'] = 'DLPNO-CCSD(T)'
        calcDict['{#MISCELLANEOUS}'] = "%mdci\nTCutPNO 3.33e-7\nTCutPairs 1e-4\nTCutMKN 1e-3\nend\n"
        calcDict['{#MEMORY}'] = '13000' 
    elif calcDict['{#METHOD}'] == 'DLPNO_CCSD':
        calcDict['{#METHOD}'] = 'DLPNO-CCSD'
        calcDict['{#MISCELLANEOUS}'] = "%mdci\nTCutPNO 3.33e-7\nTCutPairs 1e-4\nTCutMKN 1e-3\nend\n"
        calcDict['{#MEMORY}'] = '13000'
    else:
        calcDict['{#MISCELLANEOUS}'] = ''
        calcDict['{#MEMORY}'] = '8000'


    return calcDict

def seedMagic(start_dir, calcDict, template_inp='template_input.inp', template_job='template_job.sl', water_file='numWat.txt', charge_file='solute.qcinp'):
    
    '''
    READ IN INPUTS FROM A FILE!
    '''

    template_inp = calcDict['{#TEMPLATE_INP}']
    template_job = calcDict['{#TEMPLATE_JOB}']

    #required files for magical seeding. 
    charge_file = "solute.qcinp"
    geom_file = "solute.xyz"
    water_file = 'numWat.txt'
    req_files = [geom_file, water_file, charge_file]
    
    for root, dir, files in os.walk(start_dir):
        if set(req_files).issubset(set(files)):
            print("Seeding %s" % root)
            # get required info
            calcDict['{#NWATERS}'] = getWaters(os.path.join(root, water_file))
            calcDict['{#CHARGE}'], calcDict['{#MULT}'] = getMultCharge(
                os.path.join(root, charge_file))
            calcDict['{#MOLNAME}'] = root.split('/')[-2]
            calcDict['{#CONFNUM}'] = int(root.split('/')[-1].split('_')[1])
            calcDict['{#PROTNUM}'] = int(root.split('/')[-1].split('_')[2])
            calcDict['{#TITLE}'] = "molecule %s, conformer # %d and protonation = %d " % (
                calcDict['{#MOLNAME}'],
                calcDict['{#CONFNUM}'], calcDict['{#PROTNUM}'])
            calcDict['{#STARTDIR}'] = root
            
            if not calcDict['{#JOB_TYPE}'] == 'SPE':
                calcDict["{#NAME}"] = "%s_%s_%s" % (
                    calcDict['{#METHOD}'].split('/')[0].split()[0], calcDict['{#MOLNAME}'], calcDict['{#IMP_SOLV_NAME}'].split('(')[0])
                calcDict['{#XYZFILE}'] = geom_file
            
            else:
                calcDict["{#NAME}"] = "%s_%s_%s_SPE" % (
                    calcDict['{#METHOD}'].split('-')[0], calcDict['{#MOLNAME}'], calcDict['{#IMP_SOLV_NAME}'].split('(')[0])
            
                optName = "%s_%s_%s" % (
                    calcDict['{#GO_METHOD}'].split('-')[0], calcDict['{#MOLNAME}'], calcDict['{#GO_IMP_SOLV}'].split('(')[0])
                opt_geom_file = optName + '_opt.xyz'
                calcDict['{#XYZFILE}'] = opt_geom_file
            
            calcDict["{#FULLPATH}"] = os.path.join(root, calcDict["{#NAME}"])

            new_job = os.path.join(root, "%s.sl" % calcDict["{#NAME}"])
            new_inp = os.path.join(root, "%s.inp" % calcDict["{#NAME}"])
           
            # create new pbs and inp files
            
            replace_allnew(template_job, new_job, calcDict)
            replace_allnew(template_inp, new_inp, calcDict)

