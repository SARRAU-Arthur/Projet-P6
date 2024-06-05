from hapi import *
from constantes import *
from scipy.integrate import quad

import numpy as np

def chemin_acces(langue, lettre, extension):
    """ Raccourci pour aller chercher le fichier au format choisi 
    dans un dossier potentiellement différent de l'actuel """
    return f'../{langue}/{lettre}.{extension}'

def chargement_données_NIST():
    """ Chargement données de la base donnée NIST en deux tableaux: 
    taux transmission CO2 (en %) et longueur d'onde (en m) """
    data = np.loadtxt(chemin_acces('Bases de données','CO2 Absorption NIST','csv'), delimiter = ';')
    taux_CO2 = []
    nombre_onde = 1E-2 * data[:-1,0] # NaN dernière ligne, on exclu pour éviter erreurs à l'exécution
    for i in range(0,np.size(data[:,0]) - 1):
        taux_CO2.append(np.mean(data[i,1:5]))
    return 1 / nombre_onde, taux_CO2

def chargement_données_HITRAN_complet_z_constant():
    """ Chargement données de la base donnée depuis le site HITRAN.org
    et redistribution dans un fichier csv en deux colonnes:
    taux transmission CO2 (pas en %) et nombre d'onde (en m) """
    db_begin('data') # Chargement données depuis site
    fetch('CO2', 2, 1, 500, 2100) # Accès aux données: (numéro molécule CO2 = 2) entre 500 et 2100 cm^-1
    nombre_onde, transmittance = transmittanceSpectrum(nombre_onde, k_abs)
    nom_fichier = chemin_acces('Bases de données','CO2 Absorption z constant HITRAN','csv')
    with open(nom_fichier, mode = 'w', newline = '') as fichier_csv:
        for lignes in range(0,len(nombre_onde)):
            fichier_csv.write(str(nombre_onde[lignes]) + ' ; ' + str(transmittance[lignes]) + '\n')
    return None

def fonction_mathématiques_coefficients():
    
    def coefficient_a(z):
        a = [-6.5E-3, 0, 1E-3, 2.8E-3]
        return np.piecewise(z,
                [z < z_trop,
                    (z_trop <= z) & (z < z_strat1),
                    (z_strat1 <= z) & (z < z_strat2),
                    (z_strat2 <= z) & (z < z_meso)],
                [a[0],
                    a[1],
                    a[2],
                    a[3]])
    
    def coefficient_b(z):
        # b = [288, 219, 199, 138.5]
        b = [288.15, 216.65, 216.65, 228.65]
        return np.piecewise(z,
                [z < z_trop,
                    (z_trop <= z) & (z < z_strat1),
                    (z_strat1 <= z) & (z < z_strat2),
                    (z_strat2 <= z) & (z < z_meso)],
                [b[0],
                    b[1],
                    b[2],
                    b[3]])
        
    def coefficient_P_0(z):
        P = [P_0, P_trop, P_strat1, P_strat2]
        return np.piecewise(z,
                [z < z_trop,
                    (z_trop <= z) & (z < z_strat1),
                    (z_strat1 <= z) & (z < z_strat2),
                    (z_strat2 <= z) & (z < z_meso)],
                [P[0],
                    P[1],
                    P[2],
                    P[3]])
            
    return coefficient_a, coefficient_b, coefficient_P_0

def fonction_mathématique_température_altitude():
    
    def température_altitude(z):
        a, b, _ = fonction_mathématiques_coefficients()
        z_base = fonction_mathématique_altitude_base()
        return a(z) * (z - z_base(z)) + b(z)
    
    return température_altitude

def fonction_mathématique_altitude_base():
    
    def altitude_base(z):
        liste_z = np.array([0, z_trop, z_strat1, z_strat2])
        indice = np.argmin([np.abs(z - élément_liste) for élément_liste in liste_z])
        z_base = liste_z[indice - 1] if z < liste_z[indice] and indice > 0 else liste_z[indice]
        return z_base
    
    return altitude_base

def fonction_mathématique_pression_altitude():
  
    def pression_altitude(z):
        a, b, P_base = fonction_mathématiques_coefficients()
        z_0 = fonction_mathématique_altitude_base()
        T = fonction_mathématique_température_altitude()
        return P_base(z) * np.exp(- (M * g * (z - z_trop)) / (R * T(z))) if (z_trop <= z) & (z < z_strat1) \
                            else P_base(z) * ((a(z) * z + b(z)) / (z_0(z) * a(z) + b(z))) ** (-M * g / (R * a(z)))
    
    return pression_altitude

def fonction_mathématique_quantité_matière_altitude():

    def quantité_matière_altitude(z):
        T = fonction_mathématique_température_altitude()  
        P = fonction_mathématique_pression_altitude()
        return P(z) / (k_B * T(z))
    
    return quantité_matière_altitude

def chargement_données_HITRAN_complet_fonction_z():
    """ Chargement données de la base donnée depuis le site HITRAN.org
    et redistribution dans un fichier csv en deux colonnes:
    taux transmission CO2 (pas en %) et nombre d'onde (en m) """
    db_begin('data') # Chargement données depuis site
    fetch('CO2', 2, 1, 500, 2100) # Accès aux données: (numéro molécule CO2 = 2) entre 500 et 2100 cm^-1
    nombre_onde, k_abs = absorptionCoefficient_Lorentz(SourceTables = 'CO2', 
                                         Diluent = {'air': 1.0})
    intégrale_densité_moléculaire = quad(fonction_mathématique_quantité_matière_altitude(), 0, z_meso, \
                                    limit = 10 ** 7, full_output = 1)[0]
    transmittance = 1 - (k_abs * 10E-4 * CO2_fraction) * intégrale_densité_moléculaire 
    nom_fichier = chemin_acces('Bases de données', 'CO2 Absorption fonction z HITRAN', 'csv')
    with open(nom_fichier, mode = 'w', newline = '') as fichier_csv:
        for lignes in range(0, len(nombre_onde)):
            fichier_csv.write(str(nombre_onde[lignes]) + ' ; ' + str(transmittance[lignes]) + '\n')
    return nombre_onde, k_abs

# def chargement_données_HITRAN_complet_fonction_z_k_abs():
#     """ Chargement données de la base donnée depuis le site HITRAN.org
#     et redistribution dans un fichier csv en deux colonnes:
#     taux transmission CO2 (pas en %) et nombre d'onde (en m) """
#     db_begin('data') # Chargement données depuis site
#     fetch('CO2', 2, 1, 500, 2100) # Accès aux données: (numéro molécule CO2 = 2) entre 500 et 2100 cm^-1
#     nombre_onde, k_abs = absorptionCoefficient_Lorentz(SourceTables = 'CO2', 
#                                          Diluent = {'air': 1.0})
#     nom_fichier = chemin_acces('Bases de données', 'CO2 Absorption fonction z HITRAN k_abs', 'csv')
#     with open(nom_fichier, mode = 'w', newline = '') as fichier_csv:
#         for lignes in range(0, len(nombre_onde)):
#             fichier_csv.write(str(nombre_onde[lignes]) + ' ; ' + str(k_abs[lignes]) + '\n')
#     return None

def chargement_données_HITRAN(nom_fichier):
    """ Chargement données de la base donnée HITRAN en un tableau à 2 colonnes:
    taux transmission CO2 (en %) et longueur d'onde (en m). Préciser dans <nom_fichier>
    si l'on souhaite utiliser les données avec z constant ou en fonction de z """
    data = np.loadtxt(chemin_acces('Bases de données', nom_fichier, 'csv'), delimiter = ';')
    taux_CO2 = []
    nombre_onde = 1E-2 * data[:,0]
    for i in range(0,np.size(data[:,0])):
        taux_CO2.append(data[i,1] * 100)
    return 1 / nombre_onde, taux_CO2

def chargement_données():
    return chargement_données_HITRAN('CO2 Absorption fonction z HITRAN')