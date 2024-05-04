# Bibliothèques locales

from constantes import *
from base_données import *

# Bibliothèques externes

from scipy.integrate import quad
from scipy.interpolate import interp1d
from os import system, name

import numpy as np
import matplotlib.pyplot as plt
import warnings 

# Signatures et implémentation de fonctions
    
def fonction_mathématique_corps_noir():
    """ Renvoie la fonction de luminance d'un corps noir de température T 
        et en fonction de la longueur d'onde lambda_ """
    
    def valeur_luminance_corps_noir (lambda_, T = T):
        """ Renvoie l'image d'une valeur lamba_ donnée à travers la fonction de luminance 
        d'un corps noir de température T """
        term_1 = np.pi * (C_1 * (10 ** 6) ** 4) / lambda_ ** 5
        term_2 = np.exp(C_2 * (10 ** 6) / (T * lambda_)) - 1
        return term_1 / term_2
    
    return valeur_luminance_corps_noir()

def produit_de_fonctions(fonction1, fonction2):
    """ Réalise un produit de fonctions et retourne un objet de type <function> """
    
    def fonction_produit(x):
        """ Réalise un produit de deux images d'une valeur x à travers deux fonctions """
        return fonction1(x) * fonction2(x)
    
    return fonction_produit

def fonction_transmittance_vers_absorbance(transmittance):
    """ Réalise une fonction abssorbance retourne un objet de type <function> """
    
    def absorbance(x):
        """ Réalise l'opération 1 - transmittance(x) """
        return 1 - transmittance(x)
    
    return absorbance

def intégrale(fonction, borne_inf, borne_sup): 
    return quad(fonction, borne_inf, borne_sup, limit = 10 ** 7, full_output = 1)

def spectre_transmission_CO2 (données_abscisses, données_ordonnées):
    """ Représentation graphique spectre transmission CO2 en fonction de la longueur d'onde """
    plt.plot(données_abscisses, données_ordonnées)
    plt.title("Taux de transmission du CO2 en fonction de la longueur d'onde")
    plt.xlabel("Longueur d'onde (en m)")
    plt.ylabel("Taux de transmission du CO2 (en %)")
    plt.grid(True)
    plt.show()
    return None

def spectre_luminance_corps_noir (abscisses, fonction):
    """ Représentation graphique luminance corps noir en fonction de la longuer d'onde """
    affichage_fonction_continue(abscisses, fonction)
    plt.title("")
    plt.xlabel("Longueur d'onde (en m)")
    plt.ylabel('Luminance spectrale (en kg.m^-1.s^-3)')
    plt.grid(True)
    plt.show()
    return None

def affichage_fonction_continue (tableau_absicces, fonction):
    """ Affiche le graphique d'une fonction parfaitement continue """
    x = np.logspace(min(tableau_absicces), max(tableau_absicces))
    plt.plot(x, fonction(x))
    return None
        
def fonction_mathématique_interpolation (longueur_onde, taux_CO2):
    """ Retourne l'absorbance et la transmittance du CO2 en fonction de 
    la longueur d'onde (en m) sous forme d'un objet de classe <function> 
    qui est mathématiquement continue """
    
    def valeur_interpolation (longueur_onde, taux_CO2):
        """ Retourne la transmittance du CO2 à une longueur d'onde donnée, 
        sous la forme d'un objet de classe <interp1D> """
        return interp1d(longueur_onde, taux_CO2, kind = 'linear', 
                        bounds_error = False, fill_value = (100, 100))
    
    transmittance = lambda x: évaluer_fonction_interpolation(valeur_interpolation(longueur_onde, taux_CO2), x) \
                    if min(longueur_onde) <= x <= max(longueur_onde) \
                    else 1.0
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
                    
def évaluer_fonction_interpolation (fonction, valeur):
    """ Obtenir l'image d'une valeur à travers une fonction de classe <interp1D> """
    return (fonction.__call__(valeur)).tolist()

def affichage_physique (paramètre, M_0, T):
    """ Affichage de résultat avec valeur, unité et incertitude associées """
    print(f"Exitance totale {paramètre}: M_0({T} K) = {M_0[0]} ± {M_0[1]} W.m^-2 = kg.s^-3.K^-4")
        
# Programme principal

system('cls' if name == 'nt' else 'clear')
warnings.filterwarnings('ignore')

""" Chargement données selon le modèle choisi """
longueur_onde, taux_CO2 = chargement_données()

""" Transformation vers fonctions mathématiques continues """
CO2_transmittance = fonction_mathématique_interpolation_transmittance(longueur_onde, taux_CO2)
CO2_absorbance = fonction_mathématique_interpolation_absorbance(longueur_onde, taux_CO2)

""" Système Terre"""
T = T_0

flux_émis_corps_noir_Terre = (C_S * T ** 4, 0.00)
affichage_physique('théorique', flux_émis_corps_noir_Terre, T)

luminance_corps_noir_Soleil = fonction_mathématique_corps_noir()
intégrande_absorbance_Soleil = produit_de_fonctions(CO2_absorbance, luminance_corps_noir_Soleil)
M_0_sans_CO2 = intégrale(intégrande_absorbance_Soleil, 0, np.inf)
affichage_physique('Soleil', flux_émis_corps_noir_Terre, T)

""" Système Soleil"""
T = T_S
flux_émis_corps_noir_Soleil = (C_S * T ** 4, 0.00)
affichage_physique('théorique', flux_émis_corps_noir_Soleil, T)

""" Flux sans CO2 """
intégrande_sans_taux_CO2 = fonction_mathématique_corps_noir()
M_0_sans_CO2 = intégrale(intégrande_sans_taux_CO2, 0, np.inf)
affichage_physique('sans CO2', M_0_sans_CO2, T)

""" Flux avec CO2 """
intégrande_avec_taux_CO2 = produit_de_fonctions(intégrande_sans_taux_CO2, CO2_transmittance)
M_0_avec_CO2 = intégrale(intégrande_avec_taux_CO2, 0, np.inf)
affichage_physique('avec CO2', M_0_avec_CO2, T)

""" Affichage spectre tranmission CO2 """
# spectre_transmission_CO2(longueur_onde, taux_CO2)

""" Affichage spectre luminance spectrale corps noir """
# spectre_luminance_corps_noir(longueur_onde, fonction_mathématique_corps_noir())