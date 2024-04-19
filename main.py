import matplotlib
import numpy as np
import sympy as sym
from Helpers import identifier, isCharacter
import math
from numpy import matrix, array
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, show, xlabel, ylabel, legend, title
import scipy.optimize as opt
from GPII import *
from math import sqrt


plt.rcParams["text.usetex"] = True
tex_fonts = {
    # Use LaTeX to write all text
    "text.usetex": True,
    "font.family": "serif",
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 10,
    "font.size": 10,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8
}

plt.rcParams.update(tex_fonts)
plt.rc('text', usetex=True)







def gauss(term):
    ids = identifier(term)
    symbols = []
    for str1 in ids:
        symbols.append(sym.sympify(str1))
    termSymbol = sym.sympify(term)
    values = []
    for ident in ids:
        exec("values.append(" + ident + ")")

    derivatives = []
    i = 0
    while i < len(symbols):
        r = sym.diff(termSymbol, symbols[i])
        j = 0
        while j < len(symbols):
            # exec('r.evalf(subs={symbols[j]: ' + values[j] + '})')
            r = r.evalf(subs={symbols[j]: values[j]})
            j += 1
        derivatives.append(r.evalf())
        i += 1
    i = 0
    while i < len(derivatives):
        exec("derivatives[i] *= sigma_" + ids[i])
        i = i + 1
    res = 0
    for z in derivatives:
        res += z ** 2
    return math.sqrt(res)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



m1 = matrix("""
1 2 3;
4 5 6


""")


#Gitter

# 1. Messung
L1 = 0.675
g = 1005 #linien pro cm
x1 = matrix("""
8.5;
17;
25.9;
34.9;
44.5;
54.8;
66
""")# von minus bis plus, in cm
sigma_x = 4/1000 #mm
# 2. Messung
L2 = 0.37
x2 = matrix ("""
4.5;
9.4;
14.2;
19.3;
24.5;
30.1;
36.4
""")

#3.Messung
L3 = 0.443
x3 = matrix ("""
5.7;
14.2;
17.1;
23.2;
29.4;
36.2;
43.6
""")

# Einzelspalt
lamda = 633/(10**9)
b_literatur = 80/(10**6)
#1. Messung mit Lineal
L1 = 0.867
x1 = array([1.3, 2.7, 4.1, 5.5, 6.9])# wieder doppelte abstände in cm


L2 = 0.659
x2 = array([1.05, 2.05, 3.1, 4.15, 5.25])

L3 = 0.466
x3 = array([0.8, 1.5, 2.4, 3, 3.7, 4.5])
# 2. Photodiode
# 0.5 cm abstände
L4 = 0.376
U = array([0.02, 0.1, 0.19, 0.21, 0.1, 0.01, 0.2, 0.73, 0.15, 0.92, 0.32, 0.16, 2.0, 6.7, 13.8, 13.8, 13.8, 13.8, 13.8, 11.5, 4.87, 1.19, 0.07, 0.53, 1.09, 1.06, 0.57, 0.1, 0.03, 0.2, 0.34, 0.29, 0.12, 0]) # in Volt
#3 Linse
#Bildgröße
B = 0.0015 #m
#Gegenstandsgröße
G = 0.000080
#Bildweite
bw = 1.31
#Gegenstandsweite
gw = 0.12

#Spektrum
L_gitter = 1.25 #m

L_linse = 1.335

L_spalt = 1.45

#Fresnelbeugung Lochblende

L_linsfres = 0.995# zur wand

W_Fresnel = matrix("""
2 48.5;
3 78.2;
4 83.2;
5 86.2;
6 88.6

""")# anzahl ordnungen, entfernung schirm zu blende

