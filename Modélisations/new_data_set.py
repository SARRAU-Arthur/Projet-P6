from hapi import *

import matplotlib.pyplot as plt
import numpy as np
import csv

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