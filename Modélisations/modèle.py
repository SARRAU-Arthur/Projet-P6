# Bibliothèques locales

from constantes import *
from base_données import *
from code_david_louapre import *

# Bibliothèques externes

from scipy.integrate import quad
from scipy.interpolate import interp1d
from os import system, name

import numpy as np
import matplotlib.pyplot as plt
import warnings
import time

# Signatures et implémentation de fonctions
    
def fonction_mathématique_corps_noir(T):
    """ Renvoie la fonction de luminance d'un corps noir de température T 
        et en fonction de la longueur d'onde lambda_ """
    
    def valeur_luminance_corps_noir (lambda_):
        """ Renvoie l'image d'une valeur lamba_ donnée à travers la fonction de luminance 
        d'un corps noir de température T """
        term_1 = C_1 / (lambda_ ** 5)
        term_2 = np.exp(C_2 / (T * lambda_)) - 1
        return term_1 / term_2
    
    return valeur_luminance_corps_noir

def fonction_mathématique_facteur_terre_soleil():
    flux_émis_corps_noir_Soleil = C_S * T_S
    # r = flux_émis_corps_noir_Soleil[0] * (1 - a_T) * (R_T * R_S / (2 * d_TS ** 2)) ** 2 / 16
    p = (R_T / (2 * d_TS)) ** 2
    surface_soleil = 4 * np.pi * R_S ** 2
    surface_terre = 4 * np.pi * R_T ** 2
    densité_flux_soleil = flux_émis_corps_noir_Soleil * surface_soleil
    # print(f"{densité_flux_soleil = :.2e}")
    phi_abs_terre = (1 - a_T) * densité_flux_soleil * p
    # print(f"{phi_abs_terre = :.2e}")
    flux_surfacique_terre_sol = phi_abs_terre / surface_terre
    # print(f"{flux_surfacique_terre_sol = :.2e}")
    facteur = flux_surfacique_terre_sol / flux_émis_corps_noir_Soleil
    return facteur

def spectre_luminance_corps_noir (intervalle, T, corps):
    """ Représentation graphique luminance corps noir en fonction de la longuer d'onde """
    y = fonction_mathématique_corps_noir(T)
    plt.xscale('log')
    plt.plot(intervalle, y(intervalle), label = f'T = {T} K')
    plt.xlabel("Longueur d'onde (en m)")
    plt.ylabel("Luminance spectrale (en kg.m^-1.s^-3)")
    plt.legend()
    plt.title(f"Spectre de luminance corps noir: {corps}")
    plt.grid(True)
    plt.show()
    return None

def produit_de_fonctions(fonction1, fonction2):
    """ Réalise un produit de fonctions et retourne un objet de type <function> """
    
    def fonction_produit(x):
        """ Réalise un produit de deux images d'une valeur x à travers deux fonctions """
        return fonction1(x) * fonction2(x)
    
    return fonction_produit

def intégrale(fonction, borne_inf, borne_sup): 
    return quad(fonction, borne_inf, borne_sup, full_output = 1)

def spectre_transmission_CO2 (données_abscisses, données_ordonnées):
    """ Représentation graphique spectre transmission CO2 en fonction de la longueur d'onde """
    plt.plot(données_abscisses, données_ordonnées * 100)
    plt.xlabel("Longueur d'onde (en m)")
    plt.ylabel("Transmittance du CO2 (en %)")
    plt.grid(True)
    plt.show()
    return None

def éléments_graphe(grandeur):
    plt.xlabel("Altitude (en m)")
    plt.legend()
    plt.axvline(x = z_trop, color = 'grey', linestyle = '--', linewidth = 0.85, label = str(z_trop))
    plt.text(z_trop + 100, np.mean(grandeur), 'Troposphère', rotation = 'vertical')
    plt.axvline(x = z_strat1, color = 'grey', linestyle = '--', linewidth = 0.85, label = z_strat1)
    plt.text(z_strat1 + 100, np.mean(grandeur), 'Stratosphère 1', rotation = 'vertical')
    plt.axvline(x = z_strat2, color = 'grey', linestyle = '--', linewidth = 0.85, label = z_strat2)
    plt.text(z_strat2 + 100, np.mean(grandeur), 'Stratosphère 2', rotation = 'vertical')
    plt.axvline(x = z_meso, color = 'grey', linestyle = '--', linewidth = 0.85, label = z_meso)
    plt.text(z_meso + 100, np.mean(grandeur), 'Mésosphère 1', rotation = 'vertical')
    plt.grid(True)
    return None

def profile_P_fonction_altitude():
    z = np.linspace(0, z_meso - 1)
    fonction_P = fonction_mathématique_pression_altitude()
    P =  np.array([fonction_P(z_i) for z_i in z])
    plt.ylabel("Pression (en Pa)")
    plt.plot(z, pressure(z), marker = '*', label = 'Atmosphère isotherme')
    plt.plot(z, P, marker = '+', label = 'Atmosphère Standard 1976')
    éléments_graphe(P)
    plt.show()
    return None

def profile_T_fonction_altitude():
    z = np.linspace(0, z_meso - 1)
    fonction_T = fonction_mathématique_température_altitude()
    T = np.array([fonction_T(z_i) for z_i in z])
    plt.xlabel("Température (en K)")
    plt.plot(T, z, marker = '+')
    plt.ylabel("Altitude (en m)")
    plt.axhline(y = z_trop, color = 'grey', linestyle = '--', linewidth = 0.85, label = z_trop)
    plt.text(np.mean(T), z_trop + 100, 'Troposphère')
    plt.axhline(y = z_strat1, color = 'grey', linestyle = '--', linewidth = 0.85, label = z_strat1)
    plt.text(np.mean(T), z_strat1 + 100, 'Stratosphère 1')
    plt.axhline(y = z_strat2, color = 'grey', linestyle = '--', linewidth = 0.85, label = z_strat2)
    plt.text(np.mean(T), z_strat2 + 100, 'Stratosphère 2')
    plt.axhline(y = z_meso, color = 'grey', linestyle = '--', linewidth = 0.85, label = z_meso)
    plt.text(np.mean(T), z_meso + 100, 'Mésosphère 1')
    plt.grid(True)
    plt.show()
    return None

def profile_n_fonction_altitude():
    z = np.linspace(0, z_meso - 1)
    fonction_n = fonction_mathématique_quantité_matière_altitude()
    n =  np.array([fonction_n(z_i) for z_i in z])
    plt.ylabel("Densité particulaire volumique (en molec.m^-3)")
    plt.plot(z, air_number_density(z), marker = '*', label = 'Atmosphère isotherme')
    plt.plot(z, n, marker = '+', label = 'Atmosphère Standard 1976')
    éléments_graphe(n)
    plt.show()
    return None

def fonction_transmittance_vers_absorbance(transmittance):
    """ Réalise une fonction absorbance qui retourne un objet de type <function> """
    
    def absorbance(x):
        """ Réalise l'opération 1 - transmittance pour une valeur x donnée"""
        return 1 - transmittance(x)
    
    return absorbance

def évaluer_fonction_interpolation (fonction, valeur):
    """ Obtenir l'image d'une valeur à travers une fonction de classe <interp1D> """
    return (fonction.__call__(valeur)).tolist()
 
def fonction_mathématique_interpolation (longueur_onde, taux_CO2):
    """ Retourne l'absorbance et la transmittance du CO2 en fonction de 
    la longueur d'onde (en m) sous forme d'un objet de type <function> 
    qui est mathématiquement continu """
    
    def valeur_interpolation (longueur_onde, taux_CO2):
        """ Retourne la transmittance et l'absorbance du CO2 à une longueur d'onde donnée, 
        sous la forme d'un objet de classe <interp1D> """
        return interp1d(longueur_onde, taux_CO2, kind = 'linear', 
                        bounds_error = False, fill_value = (1, 1))
    
    transmittance = lambda x: évaluer_fonction_interpolation(valeur_interpolation(longueur_onde, taux_CO2), x) \
                    if (min(longueur_onde) <= x) and (x <= max(longueur_onde)) \
                    else 1
    absorbance = fonction_transmittance_vers_absorbance(transmittance)
    return transmittance, absorbance

def fonction_mathématique_interpolation_transmittance (longueur_onde, taux_CO2):
    """ Retourne la transmittance du CO2 en fonction de la longueur d'onde (en m)
    sous forme d'un objet de classe <function> qui est mathématiquement continue """
    return fonction_mathématique_interpolation (longueur_onde, taux_CO2)[0]

def fonction_mathématique_interpolation_absorbance (longueur_onde, taux_CO2):
    """ Retourne l'absorbance du CO2 en fonction de la longueur d'onde (en m)
    sous forme d'un objet de classe <function> qui est mathématiquement continue """
    return fonction_mathématique_interpolation (longueur_onde, taux_CO2)[1]

def affichage_physique (paramètre, M_0):
    """ Affichage de résultat avec valeur, unité et incertitude associées """
    print(f"Exitance totale {paramètre}: M_0 = {M_0[0] :.2e} ± {M_0[1] :.2e} W.m^-2 = kg.s^-3.K^-4")

# Programme principal

""" Paramétrage de l'affichage """
system('cls' if name == 'nt' else 'clear')
warnings.filterwarnings('ignore')

""" Chargement données selon le modèle choisi """
longueur_onde, taux_CO2 = chargement_données()
# intervalle = np.linspace(min(longueur_onde), max(longueur_onde))
intervalle = np.linspace(1E-6, 20E-6, 200)
# intervalle = np.linspace(1E-10, 1E10)

""" Transformation vers fonctions mathématiques continues """
CO2_transmittance = fonction_mathématique_interpolation_transmittance(longueur_onde, taux_CO2)
CO2_absorbance = fonction_mathématique_interpolation_absorbance(longueur_onde, taux_CO2)

""" Affichage spectre tranmission CO2 """
# plt.plot(longueur_onde, taux_CO2 * 100, c = 'b', label = 'Avant interpolation')
# plt.plot(intervalle, [CO2_transmittance(i) * 100 for i in intervalle], 
#             marker = '+', c = 'g', label = 'Après interpolation')
# plt.xlabel("Longueur d'onde (en m)")
# plt.ylabel("Transmittance du CO2 (en %)")
# plt.grid(True)
# plt.legend()
# plt.show()

""" Affichage spectres luminance spectrale corps noir """
# intervalle = np.linspace(1E-7, 100E-6, 10000)
# spectre_luminance_corps_noir(intervalle, T_S, 'Soleil')
# intervalle = np.linspace(1E-6, 1000E-6, 10000)
# spectre_luminance_corps_noir(intervalle, T_T, 'Terre')

""" Profiles T, P, n en fonction de l'altitude """
# profile_T_fonction_altitude()
# profile_P_fonction_altitude()
# profile_n_fonction_altitude()

""" Calculs différents flux par système sous hypothèse de corps noirs """
T = T_S
print(f'> Système Soleil: T = {T} K \n')
luminance_corps_noir_Soleil = fonction_mathématique_corps_noir(T)

flux_émis_corps_noir_Soleil = (C_S * T ** 4, 0.00)
affichage_physique('théorique', flux_émis_corps_noir_Soleil)

facteur = fonction_mathématique_facteur_terre_soleil()

intégrande_absorbance_Soleil = produit_de_fonctions(CO2_absorbance, luminance_corps_noir_Soleil)
intégrande_transmittance_Soleil = produit_de_fonctions(CO2_transmittance, luminance_corps_noir_Soleil)

start_time = time.time()
M_0_absorbance_Soleil = intégrale(intégrande_absorbance_Soleil, min(longueur_onde), max(longueur_onde))
end_time = time.time()
execution_time = end_time - start_time
affichage_physique('absorbance',  (facteur * M_0_absorbance_Soleil[0], M_0_absorbance_Soleil[1])) # Flèche 1
print(f"Temps d'exécution: {execution_time} secondes")

start_time = time.time()
M_0_transmittance_Soleil = intégrale(intégrande_transmittance_Soleil, min(longueur_onde), max(longueur_onde))
end_time = time.time()
M_0_transmittance_Soleil = (facteur * M_0_transmittance_Soleil[0] * (1 - a_T), M_0_transmittance_Soleil[1]) 
execution_time = end_time - start_time
affichage_physique('transmittance', M_0_transmittance_Soleil) # Flèche 2
print(f"Temps d'exécution: {execution_time} secondes")

T = T_T
print(f'> Système Terre: T = {T} K \n')
luminance_corps_noir_Terre = fonction_mathématique_corps_noir(T)

flux_émis_corps_noir_Terre = (C_S * T ** 4, 0.00)
affichage_physique('théorique', flux_émis_corps_noir_Terre)

intégrande_absorbance_Terre = produit_de_fonctions(CO2_absorbance, luminance_corps_noir_Terre)
intégrande_transmittance_Terre = produit_de_fonctions(CO2_transmittance, luminance_corps_noir_Terre)

start_time = time.time()
M_0_absorbance_Terre = intégrale(intégrande_absorbance_Terre, min(longueur_onde), max(longueur_onde))
end_time = time.time()
execution_time = end_time - start_time
affichage_physique('absorbance', M_0_absorbance_Terre) # Flèche 3
print(f"Temps d'exécution: {execution_time} secondes")

start_time = time.time()
M_0_transmittance_Terre = intégrale(intégrande_transmittance_Terre, min(longueur_onde), max(longueur_onde))
end_time = time.time()
execution_time = end_time - start_time
affichage_physique('transmittance', M_0_transmittance_Terre) # Flèche 4
print(f"Temps d'exécution: {execution_time} secondes")

print(f'\n> Système Atmosphère: T = {T} K \n')

M_0_atmosphère = M_0_absorbance_Terre[0] / 2
affichage_physique('transmittance & absorbance', M_0_atmosphère) # Flèches 5, 6

""" Bilan thermique """
émis = flux_émis_corps_noir_Soleil
affichage_physique('Émis', émis)

absorbé = M_0_transmittance_Terre[0] + M_0_transmittance_Soleil[0]
affichage_physique('Absorbé', absorbé)