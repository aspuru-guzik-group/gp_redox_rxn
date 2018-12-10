'''
Parses an Orca file that performed a geometry optimization and harmonic analysis. 
'''
import sys
import os
import csv
#sys.path.append('./Modules')
sys.path.append('/n/home09/ajinich/deltaG/mother/code/ian_ben_adrian/Modules')
from lib.thermodynamic_constants import *
import fnmatch


def getFormationEnergy(filename):
    '''parse terms relevant to calculating formation free energy
    '''
    freqs = []
    Ee = []
    readingFreqs = False
    with open(filename, 'r') as f:
        lines = f.readlines()
        for i in range(len(lines)):
            if 'VIBRATIONAL FREQUENCIES' in lines[i - 3]:
                readingFreqs = True
            if readingFreqs is True:
                freqs.append(float(lines[i].lstrip().split()[1]))
                if '------------' in lines[i + 3]:
                    readingFreqs = False
            if "FINAL SINGLE POINT ENERGY" in lines[i]:

                Ee.append(float(lines[i].split()[-1].strip()))
            if "Total Mass" in lines[i]:
                mass = float(lines[i].strip().split()[3])
            if "qrot =" in lines[i]:
                qrotRoomTemp = float(lines[i].strip().split()[2])
            if "The first frequency considered to be a vibration is" in lines[i]:
                firstVibration = int(lines[i].strip().split()[9]) - 1
        EeMin = Eh_to_kcalpermol(Ee[-1])
        freqs = [f for f in freqs[firstVibration:] if f>0]

    return freqs, qrotRoomTemp, mass, EeMin

# explore a directory tree to get relevant files:

def get_spe(file_name):
    spe = 'NONE'
    with open(file_name, 'rb') as f_in:
        for line in f_in:
            if "FINAL SINGLE POINT ENERGY" in line:
                spe = float(line.split()[-1])
    print (spe)
    spe=Eh_to_kcalpermol(spe)
    return spe

def get_file_list(root_dir, match_string):
    file_list = []
    for root, dirs, files in os.walk(root_dir):
        for file in os.listdir(root):
            if fnmatch.fnmatch(file, match_string):
                file_list.append(os.path.join(root, file))
    return file_list


def get_mol_and_confo(path_str):
    mol_ind = -3
    confo_ind = -2
    mol = path_str.split('/')[mol_ind]
    confo = path_str.split('/')[confo_ind]
    return mol, confo


def filter_crashed_files(root_dir, file_list):
    successful_file = []
    crashed_file = []
    match_str = 'ORCA TERMINATED NORMALLY'
    for afile in file_list:
        file_str = os.path.join(root_dir, afile)
        SUCCESS = 0
        with open(file_str, 'rb') as fin:
            for line in fin:
                # print line
                if match_str in line:
                    print (afile, 'SUCCESS')
                    successful_file.append(afile)
                    break
            else:
                print (afile, 'CRASHED')
                crashed_file.append(afile)
    return successful_file, crashed_file


def getHydrogens(soluteFile):
    """Get number of hydrogen atoms
    The # hydrogens is used in when applying Alberty's Legendre Transform.
    """
    with open(soluteFile, 'r') as f:
        lines = f.readlines()[2:]
        nH = len([line.split()[0] for line in lines if line.split()[0] == 'H'])
    return nH

def getCharge(filename):
	"""Get charge of species
	"""
	flag=False
	with open(filename,'r') as f:
		for line in f.readlines():	
			if flag:	
				charge=int(line.split()[0])
				break
			if 'molecule' in line:
				flag=True
	return charge


# Get the file list
#root_dir = '/media/TOSHIBA EXT/backup/batches_trimmed_09_13_13_GORDON/batches_trimmed_01/'
'''
root_dir = '/n/aspuru_lab2/ajinich/semi_emp/molecules/'

#match_string = 'orca_go_ha_cosmo_*.out'
match_string ="zindo*.out"

file_list = get_file_list(root_dir, match_string)

file_out = 'zindo.csv'
file_crashed = 'zindoCRASHED.csv'
success, crashed = filter_crashed_files(root_dir, file_list)

soluteFileStr = 'solute.xyz'
chargeFileStr='solute.qcinp'

# Write the set of crashed files.
with open(file_crashed, 'wb') as fout:
    for afile in crashed:
        fout.write(afile + '\n')


print len(crashed), 'crashed'
print len(success), 'success'
print len(file_list), 'total'

with open(file_out, 'wb') as fout:
    writer = csv.writer(fout, delimiter='\t')
    writer.writerow(
        ['Molecule', 'conformer', 'proto_num', 'numberH', 'charge','EMin', 'Mass', 'QRot', 'Freq'])
    for afile in success:
        path_str = os.path.join(root_dir, afile)
        mol, confo = get_mol_and_confo(path_str)
        print mol, confo
        proto_number = confo.split('_')[-2]
        #number of hydrogens. 
        soluteFile = os.path.join('/'.join(path_str.split('/')[:-1]), soluteFileStr)
        nH = getHydrogens(soluteFile)
        #charge
        chargeFile = os.path.join('/'.join(path_str.split('/')[:-1]), chargeFileStr)
        charge = getCharge(chargeFile)
        
        freqs, qrot, mass, EeMin = getFormationEnergy(path_str)
        freq_str = '_'.join([str(f) for f in freqs])
        writer.writerow(
            [mol, confo, proto_number, nH, charge, EeMin, mass, qrot, freq_str])
'''
