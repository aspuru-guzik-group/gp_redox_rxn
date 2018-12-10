import csv
import numpy as np

# R_kCal = 1.986e-3 # kCal/(K*mol)
R_kCal = 1.987191683e-3  # kCal/(mol*K)
R = 8.31e-3  # kJ/(K*mol)
F = 96.48  # kC/mol
J_per_cal = 4.184
default_T = 298.15  # K
default_I = 0.25  # mM
default_pH = 7.0
default_c0 = 1  # M
default_pMg = 10
default_nMg = 0  # added AJG
default_RT = R * default_T
default_c_mid = 1e-3  # M
default_c_range = (1e-6, 1e-2)  # M
dG0_f_Mg = -455.3  # kJ/mol, formation energy of Mg2+
# mu_proton = -241.3560 # kcal/mol, chemical potential in Water,  estimated via DFT. (03/07/2013 --> refer to protonPotential.txt)
# mu_proton = -360 #kcal/mol # Fitted value 10 explicit waters only-
# Adrian Jinich
mu_proton_cosmo = -268.61  # kcal/mol #4 waters + hydronium COSMO - Ian Dunn
mu_proton_experimental = -265.9
mu_proton = mu_proton_cosmo
# mu_proton = -450 #kcal/mol, fitted solvation free energy, for PEP and 3PG, 2PG
#mu_proton = 0

# Units from Thomas Markovich.
JtoeV = 6.241509 * 1e18
hbar = 1.054571726 * 1e-34 * JtoeV  # planck constant/(2 pi) units are eV*s
h = 6.62606957 * 1e-34 * JtoeV     # planck constant - units are eV*s
speedc = float(299792458)             # speed of light - units are m/s
kb = 1.3806488 * 1e-23           # boltzmann constant -units are J/K
# Define conversion factors
mtocm = 100         # 100 cm in one meter
# convert time to length by using speedc (speed of light)
stocm = speedc * mtocm
fstos = 1e-15                   # by definition 10^15 fs per second
# E(eV) = h(eV*s)*v(1/s) so eVtoHz = 1/h(eV*s) = 2.41804*10^14
eVtoHz = 1 / h
HztoeV = 1 / eVtoHz               # transform 1/s to energy eV
stoeVm1 = 1 / HztoeV              # transfom seconds to 1/energy 1/eV
# E(J) = h(J*cm)*v(1/cm) so Jtocmm1 = 1/h(J*cm) = 5.03445*10^22
Jtocmm1 = 1 / (h / JtoeV * stocm)
# E(eV) = h(eV*cm)*v(1/cm) 1/(h c * 100 cm/m) =  8065.3947633434
eVtocmm1 = 1 / (h * stocm)
angularMomentuminAtomicUnits = 1 / 1.054571726e-34
auPerWaveNumber = 4.5563e-6  # converts wavenumbers to au.
auPereV = 0.03674932        # Converts eV to Hartree
kcaljou = 4.184e3       # Converts joules to kcals
# hJ = 6.62617363e-34       # Plank's constant in joules
hJ = 6.62606957e-34
# na = 6.02204531e23        # Avogadro's number from Thomas.
na = 6.02214129e23
bh = hJ / (300 * kb)
hartreeToJoule = 4.359744 * 10e-18
angstromInMeters = 10e-10
psToSeconds = 1e-12
gramToKilogram = 1. / 1000.  # Converts grams to kilograms
# This is amber velocity in meters/second
amberVel = angstromInMeters / (20.455 * psToSeconds)
amutokg = 1.660538921e-27

symbol_d_G = "&Delta;G"
symbol_d_G0 = "&Delta;G&deg;"
symbol_d_G_prime = "&Delta;G'"
symbol_d_G0_prime = "&Delta;G'&deg;"

symbol_dr_G = "&Delta;<sub>r</sub>G"
symbol_dr_G0 = "&Delta;<sub>r</sub>G&deg;"
symbol_dr_G_prime = "&Delta;<sub>r</sub>G'"
symbol_dr_G0_prime = "&Delta;<sub>r</sub>G'&deg;"
symbol_dr_Gc_prime = "&Delta;<sub>r</sub>G'<sup>c</sup>"

symbol_df_G = "&Delta;<sub>f</sub>G"
symbol_df_G0 = "&Delta;<sub>f</sub>G&deg;"
symbol_df_G_prime = "&Delta;<sub>f</sub>G'"
symbol_df_G0_prime = "&Delta;<sub>f</sub>G'&deg;"


def Eh_to_kcalpermol(Eh_energy):
    return Eh_energy * 627.5095


def debye_huckel(I):
    return (2.91482 * np.sqrt(I)) / (1 + 1.6 * np.sqrt(I))


def correction_function(nH, z, nMg, pH, pMg, I, T):
    """
        nH and z - are the species parameters (can be vectors)
        pH and I - are the conditions, must be scalars
        returns the correction element used in the transform function

    Returns:
        The correction, in units of kJ/mol.
    """
    DH = debye_huckel(I)
    print(nH, DH, z, pH)
# return nMg * (R*T*np.log(10)*pMg - dG0_f_Mg) + nH * (R*T*np.log(10)*pH +
# DH) - (z**2) * DH
    return nMg * (R_kCal * T * np.log(10) * pMg - dG0_f_Mg) + nH * (R_kCal * T * np.log(10) * pH + DH) - (z ** 2) * DH


def correction_function_ab_initio(nH, z, nMg, pH, pMg, I, T, mu = mu_proton):
    """
        nH and z - are the species parameters (can be vectors)
        pH and I - are the conditions, must be scalars
        returns the correction element used in the transform function

    Returns:
        The correction, in units of kCal/mol.
        mu_proton
        Addition of a term for the chemical potential of proton, with respect to electrons and nucleii in vacuum
        based on DFT estimates of hydronium and proton in 20 Water cluster.
    """
    DH = debye_huckel(I)
    #print nH, DH, z, pH
# return nMg * (R*T*np.log(10)*pMg - dG0_f_Mg) + nH * (R*T*np.log(10)*pH +
# DH) - (z**2) * DH
    return nMg * (R_kCal * T * np.log(10) * pMg - dG0_f_Mg) + nH * (R_kCal * T * np.log(10) * pH + DH - mu) - (z ** 2) * DH


def transform(dG0, nH, z, nMg, pH, pMg, I, T):
    # return dG0 + correction_function(nH, z, nMg, pH, pMg, I, T)
    return dG0 + correction_function(nH, z, nMg, pH, pMg, I, T)


def array_transform(dG0, nH, z, nMg, pH, pMg, I, T):
    """
        dG0, nH and z - are the species parameters (can be vectors)
        pH and I - are the conditions, must be scalars
        returns the transformed gibbs energy: dG0'
    """
    from util import log_sum_exp
    ddG0 = correction_function(nH, z, nMg, pH, pMg, I, T)
    print("correction term original")
    print(ddG0)
    dG0_tag = dG0 + ddG0
#    return -R * T * log_sum_exp(dG0_tag / (-R*T))
    return -R_kCal * T * log_sum_exp(dG0_tag / (-R_kCal * T))


def array_transform_ab_initio(dG0, nH, z, nMg, pH, pMg, I, T, mu = mu_proton):
    """
        dG0, nH and z - are the species parameters (can be vectors)
        pH and I - are the conditions, must be scalars
        returns the transformed gibbs energy: dG0'
    """
    from util import log_sum_exp
    ddG0 = correction_function_ab_initio(nH, z, nMg, pH, pMg, I, T, mu)
    #print "correction term ab initio"
    #print ddG0
    dG0_tag = dG0 + ddG0
#    return -R * T * log_sum_exp(dG0_tag / (-R*T))
    return -R_kCal * T * log_sum_exp(dG0_tag / (-R_kCal * T))


class RedoxCarrier(object):

    # dG0 =  -E'*F * deltaE - R*T*ln(10)*pH * deltaH
    # Where:
    # F = 96.48 # kC/mol
    #    R*T*ln(10) = 5.7 kJ/mol
    #    deltaE - change in e-
    #    deltaH - change in H+
    #    pH - the conditions in which the E' was measured

    def __init__(self, cid_ox, cid_red, nH_ox, nH_red, z_ox, z_red, E0_prime, pH, ref):
        self.cid_ox = cid_ox
        self.cid_red = cid_red
        self.nH_ox = nH_ox
        self.nH_red = nH_red
        self.z_ox = z_ox
        self.z_red = z_red
        self.E0_prime = E0_prime
        self.pH = pH
        self.ref = ref
        self.delta_H = nH_red - nH_ox
        # difference in no. of electrons
        self.delta_e = (nH_red - nH_ox) - (z_red - z_ox)
        self.ddG0_prime = -E0_prime * F * self.delta_e

        # this calculation is not correct, one must use the reverse Lagendre transform
        # in order to convert G' to G.
        self.ddG0 = self.ddG0_prime - R * \
            default_T * np.log(10) * pH * self.delta_H


class RedoxCarriers(dict):

    def __init__(self):
        for row in csv.DictReader(open('../data/thermodynamics/redox.csv', 'r')):
            name = row['name']
            cid_ox = int(row['CID_ox'])
            cid_red = int(row['CID_red'])
            nH_ox = int(row['nH_ox'])
            z_ox = int(row['charge_ox'])
            nH_red = int(row['nH_red'])
            z_red = int(row['charge_red'])
            E0_prime = float(row["E'0"])
            pH = float(row['pH'])
            ref = row['ref']
            self[name] = RedoxCarrier(cid_ox, cid_red, nH_ox, nH_red,
                                      z_ox, z_red, E0_prime, pH, ref)

    def GetAllCids(self):
        cids_ox = [val.cid_ox for val in self.values()]
        cids_red = [val.cid_red for val in self.values()]
        return set(cids_ox + cids_red)
