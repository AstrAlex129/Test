#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 15:51:06 2026

@author: skraburskia
"""

# Bibliothèque
import numpy as np 
import matplotlib.pyplot as plt

#=========================================
# Equation Spatio-temporelle
#=========================================

#=========================================
# Fonction initiale
#=========================================

def F_spatiotemp_euler(epsilon : list, t : list, rho : list, R : list, v : list ) -> list : 
    R0 = R[0]
    t0 = t[0]
    S = np.zeros(len(t))
    F = np.zeros(len(t))
    L = np.zeros(len(t))
    
    depsilon = epsilon[1]-epsilon[0]
    dt = t[1] -t[0]
    
    S[1:-1] = epsilon[1:-1]*rho[1:-1]*((R[2:]-R[:-2])/(2*dt)) - (2*R0*rho[1:-1]*v[1:-1])/(epsilon[1:-1]*R[1:-1]*t0) - (rho[1:-1]*R0)/(t0*R[1:-1])*((v[2:]-v[:-2])/(2*depsilon))
    
    F[1:-1] = -(v[1:-1]*R0)/(t0*R[1:-1]) * (rho[2:]-rho[:-2])/(2*depsilon)
    
    L[1:-1] = F[1:-1] + S[1:-1]
    
    return L

def F_spatiotemp_continuite(epsilon : list, t : list, rho : list, R : list, v : list, K : float, gamma : float, rho0 : float, U : list) -> list : 
    R0 = R[0]
    t0 = t[0]
    S = np.zeros(len(t))
    F = np.zeros(len(t))
    L = np.zeros(len(t))
    
    depsilon = epsilon[1]-epsilon[0]
    dt = t[1] -t[0]
    
    S[1:-1] = -(gamma*rho[1:-1]*K*rho0**(gamma-1))/(rho[1:-1]*R[1:-1]) * ((rho[2:] - rho[:-2])/(2*depsilon)) - (rho[1:-1]*rho0)/R[1:-1] * ((U[2:] - U[:-2])/(2*depsilon)) + (rho[1:-1]*rho0*R0*epsilon[1:-1]*v[1:-1])/t0 * ((R[2:] - R[:-2])/(2*dt))
    
    F[1:-1] = -(rho[1:-1]*rho0*R0**2*v[1:-1])/(t0**2*R[1:-1])*((v[2:] - v[:-2])/(2*depsilon))
    
    L[1:-1] = F[1:-1] + S[1:-1]
    
    return L
    
    
def F_spatiotemp_poisson(epsilon : list, rho : list, U : list, rho0 : float):
    pass

#=========================================
# Initialisation Runge-Kutta 3
#=========================================

def RK4_spatiotemp_euler(h : float,f, Y : list, n : float, epsilon : list, t : list, rho : list, R : list, v : list ) -> list[list]:
    rho_approx = []
    # a voir quelle vont être les conditions initiales. dans Y il doit y avoir les variable rho et cie pour les utiliser en suite
    # on veut calculer les 3 premières valeurs de rho ainsi que les autres variable dans S non ?
    t = 0
    while n>= 0 :
        
        rho_approx.append(Y)
        
        k1 = f(t,Y)
        k2 = f(t+h/2, Y+(h*k1)/2,epsilon, rho, R, v)
        k3 = f(t+h/2, Y +(h*k2)/2,epsilon, rho, R, v)
        k4 = f(t+h, Y + h*k3,epsilon, rho, R, v)
        Y = Y + h/6*(k1 +2*k2+ 2*k3 + k4)
        
        t = t+h
        n -= 1
    
    return rho_approx # doit contenir list[list,list,list]

def RK4_spatiotemp_continuite(h : float ,n : float,f, Y : list, epsilon : list, t : list, rho : list, R : list, v : list, K : float, gamma : float, rho0 : float, U : list) -> list[list]: # même commentaire que RK3 euleur
    rho_approx = []
    # a voir quelle vont être les conditions initiales. dans Y il doit y avoir les variable rho et cie pour les utiliser en suite
    # on veut calculer les 3 premières valeurs de rho ainsi que les autres variable dans S non ?
    t = 0
    while n>= 0 :
        
        rho_approx.append(Y)
        
        k1 = f(t,Y)
        k2 = f(t+h/2, Y+(h*k1)/2)
        k3 = f(t+h/2, Y +(h*k2)/2)
        k4 = f(t+h, Y + h*k3)
        Y = Y + h/6*(k1 +2*k2+ 2*k3 + k4)
        
        t = t+h
        n -= 1
    
    return rho_approx # doit contenir list[list,list,list]

def RK4_spatiotemp_poisson(h,T,f,n):
    pass
    

 #=========================================
 # Adam-B 3
 #=========================================   
    
def AB3_spatiotemp_euler(f,n : float, epsilon : list, t : list, rho : list, R : list, v : list ):
    rho_mat = np.zeros((n,len(t)))
    R_mat = np.zeros((n,len(t)))
    v_mat = np.zeros((n,len(t)))
    
    F_mat = np.zeros((n,len(t)))
    
    depsilon = epsilon[1]-epsilon[0]
    dt = t[1] -t[0]
    
    
    
    # avec RK3, on obtient les 3 termes précédent n. Seulement pour n >= 2
    
    
    
    
    
    
    
    pass

def AB3_spatiotemp_continuite(f):
    pass

def AB3_spatiotemp_poisson(f):
    pass
    
    
    
    
    
    
    
    
    
    
    
    