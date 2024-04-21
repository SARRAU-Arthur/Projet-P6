from hapi import *

import matplotlib.pyplot as plt
import numpy as np
import csv

def chemin_acces(langue, lettre, extension):
    """ Raccourci pour aller chercher le fichier au format choisi 
    dans un dossier potentiellement différent de l'actuel """
    return f'../{langue}/{lettre}.{extension}'

def chargement_données_NIST():
    """ Chargement données de la base donnée en deux tableaux: 
    taux transmission CO2 (en %) et longueur d'onde (en m) """
    data = np.loadtxt(chemin_acces('Bases de données','CO2 Absorption NIST','csv'), ';')
    taux_CO2 = []
    nombre_onde = 1E-2 * data[:-1,0] # Nan dernière ligne, on exclu pour éviter erreurs à l'exécution
    for i in range(0,np.size(data[:,0]) - 1):
        taux_CO2.append(np.mean(data[i,1:5]) * 100)
    return 1 / nombre_onde, taux_CO2

db_begin('data')

fetch('CO2', 2, 1, 500, 2100)

k, coef = absorptionCoefficient_Lorentz(SourceTables = 'CO2', Diluent = {'air':1.0})

nom_fichier = "CO2 Absorption HAPI.csv"
with open(nom_fichier, mode='w', newline='') as fichier_csv:
    writer = csv.writer(fichier_csv)
    writer.writerow(coef)

# largeur = 10
# hauteur = len(coef) // largeur
# print(np.reshape(coef, (hauteur, largeur)))

nu=1/k*10000
# plt.plot(nu,coef)
# plt.show()