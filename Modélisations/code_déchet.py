def fonction_taux_transmission_CO2(longueur_onde, taux_CO2):
    fonction_par_morceaux = []
    for i in range(len(longueur_onde)-1):
        fonction_par_morceaux.append(lambda x, 
                                     borne_inf = longueur_onde[i + 1], 
                                     borne_sup = longueur_onde[i], 
                                     y = taux_CO2[i]:                                         
                                     0 if borne_inf <= x <= borne_sup else y)
        
    def fonction(x):
        for f in fonction_par_morceaux:
            if f(x) != 0:
                return f(x)
            else:
                pass
    
    return fonction