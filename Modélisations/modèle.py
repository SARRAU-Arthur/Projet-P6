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
    
def fonction_mathématique_corps_noir(T):
    """ Renvoie la fonction de luminance d'un corps noir de température T 
        et en fonction de la longueur d'onde lambda_ """
    
    def valeur_luminance_corps_noir (lambda_):
        """ Renvoie l'image d'une valeur lamba_ donnée à travers la fonction de luminance 
        d'un corps noir de température T """
        term_1 = C_1 / lambda_ ** 5
        term_2 = np.exp(C_2 / (T * lambda_)) - 1
        print(term_1, term_2)
        return term_1 / term_2
    
    return valeur_luminance_corps_noir

def spectre_luminance_corps_noir ():
    """ Représentation graphique luminance corps noir en fonction de la longuer d'onde """
    x = np.linspace(0, 10)
    y = fonction_mathématique_corps_noir(T_S)
    plt.plot(x, y(x))
    y = fonction_mathématique_corps_noir(T_T)
    plt.plot(x, y(x))
    plt.xlabel("Longueur d'onde (en m)")
    plt.ylabel('Luminance spectrale (en kg.m^-1.s^-3)')
    plt.grid(True)
    plt.show()
    return None

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

def profile_fonction_altitude():
    z = np.linspace(0, z_meso - 1)
    f = fonction_mathématique_température_altitude()
    plt.xlabel("Altitude atmosphérique (en m)")
    plt.scatter(z, f(z), marker = '+')
    plt.grid(True)
    plt.show()
    return None

def fonction_mathématique_pression_altitude():
  
    def pression_altitude(z):
        a, b, P_base = fonction_mathématiques_coefficients()
        liste_z = np.array([0, z_trop, z_strat1, z_strat2])
        indice =  np.argmin([np.abs(z - élément_liste) for élément_liste in liste_z])
        z_base = liste_z[indice - 1] if z < liste_z[indice] and indice > 0 else liste_z[indice]
        return np.piecewise(z, 
                            [(z_trop <= z) & (z < z_strat1),
                            (z_trop > z) & (z >= z_strat1)],
                            [P_base(z) * np.exp((-M * g * (z_trop - z)) / (R * T_T)),
                            P_base(z) * (a(z) * z / (z_base * a(z) + b(z))) ** (-M * g / (R * a(z)))])
        
    return pression_altitude
        
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
                    
def évaluer_fonction_interpolation (fonction, valeur):
    """ Obtenir l'image d'une valeur à travers une fonction de classe <interp1D> """
    return (fonction.__call__(valeur)).tolist()

def affichage_physique (paramètre, M_0):
    """ Affichage de résultat avec valeur, unité et incertitude associées """
    print(f"Exitance totale {paramètre}: M_0 = {M_0[0]} ± {M_0[1]} W.m^-2 = kg.s^-3.K^-4")
        
# Programme principal

""" Affichage """
system('cls' if name == 'nt' else 'clear')
warnings.filterwarnings('ignore')

""" Chargement données selon le modèle choisi """
longueur_onde, taux_CO2 = chargement_données()

""" Transformation vers fonctions mathématiques continues """
CO2_transmittance = fonction_mathématique_interpolation_transmittance(longueur_onde, taux_CO2)
CO2_absorbance = fonction_mathématique_interpolation_absorbance(longueur_onde, taux_CO2)

""" Affichage spectre tranmission CO2 """
# spectre_transmission_CO2(longueur_onde, taux_CO2)

""" Affichage spectre luminance spectrale corps noir """
spectre_luminance_corps_noir()

""" Calculs différents flux par système sous hypothèse de corps noirs """
# Système Terre
T = T_T
print(f'> Système Terre: T = {T} K \n')
luminance_corps_noir_Terre = fonction_mathématique_corps_noir(T)

flux_émis_corps_noir_Terre = (C_S * T ** 4, 0.00)
affichage_physique('théorique', flux_émis_corps_noir_Terre)

plt.plot(longueur_onde, CO2_absorbance(longueur_onde))
plt.plot(longueur_onde, CO2_transmittance(longueur_onde))
plt.show()
intégrande_absorbance_Terre = produit_de_fonctions(CO2_absorbance, luminance_corps_noir_Terre)
M_0_absorbance_Terre = intégrale(intégrande_absorbance_Terre, 0, np.inf)
affichage_physique('absorbance', M_0_absorbance_Terre) # Flèche 3

intégrande_transmittance_Terre = produit_de_fonctions(CO2_transmittance, luminance_corps_noir_Terre)
M_0_transmittance_Terre = intégrale(intégrande_transmittance_Terre, 0, np.inf)
affichage_physique('transmittance', M_0_transmittance_Terre) # Flèche 4

# Système Soleil
T = T_S
print(f'\n > Système Soleil: T = {T} K \n')
luminance_corps_noir_Soleil = fonction_mathématique_corps_noir(T)

flux_émis_corps_noir_Soleil = (C_S * T ** 4, 0.00)
affichage_physique('théorique', flux_émis_corps_noir_Soleil)

intégrande_absorbance_Soleil = produit_de_fonctions(CO2_absorbance, luminance_corps_noir_Soleil)
M_0_absorbance_Soleil = intégrale(intégrande_absorbance_Soleil, 0, np.inf)
affichage_physique('absorbance', M_0_absorbance_Soleil) # Flèche 1

intégrande_transmittance_Soleil = produit_de_fonctions(CO2_transmittance, luminance_corps_noir_Soleil)
M_0_transmittance_Soleil = intégrale(intégrande_transmittance_Soleil, 0, np.inf) * (1 - a_T)
affichage_physique('transmittance', M_0_transmittance_Soleil) # Flèche 2

# Système Atmosphère
print(f'\n > Système Atmosphère: T = {T} K \n')

M_0_atmosphère = M_0_absorbance_Terre[0] / 2
affichage_physique('transmittance & absorbance', M_0_atmosphère) # Flèches 5, 6

""" Bilan thermique """
tolérance = 10

émis = flux_émis_corps_noir_Soleil
affichage_physique('Émis', émis)

absorbé = M_0_transmittance_Terre[0] + M_0_transmittance_Soleil[0]
affichage_physique('Absorbé', absorbé)

print(f'Conclusion modèle: {np.abs(émis - absorbé) < tolérance}')