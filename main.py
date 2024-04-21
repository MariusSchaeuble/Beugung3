import matplotlib
import numpy
import numpy as np
import sympy as sym
from Helpers import identifier, isCharacter
import math
from numpy import matrix, array, mean, std, max, linspace
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, show, xlabel, ylabel, legend, title, savefig
import scipy.optimize as opt
from GPII import *
from math import sqrt
pi = math.pi


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


matplotlib.rc('xtick', labelsize=20)
matplotlib.rc('ytick', labelsize=20)





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
sigma_lamda = 0
b_literatur = 80/(10**6)
#1. Messung mit Lineal
L1 = 0.867
x1 = array([1.3, 2.7, 4.1, 5.5, 6.9])# wieder doppelte abstände in cm


L2 = 0.659
x2 = array([1.05, 2.05, 3.1, 4.15, 5.25])

L3 = 0.466
x3 = array([0.8, 1.5, 2.4, 3, 3.7, 4.5])

ordnung5 = array([1, 2, 3, 4, 5])
ordnung6 = array([1, 2, 3, 4, 5, 6])
alpha1 = x1/(200*L1)
alpha2 = x2/(200*L2)
alpha3 = x3/(200*L3)

def linear(x, m):
    return m*x



plot(ordnung5, alpha1, 'x', label='L = 86cm')
optimizedParameters, s = opt.curve_fit(linear, ordnung5, alpha1)
m_std = np.sqrt(np.diag(s))
plt.plot(ordnung5, linear(ordnung5, *optimizedParameters), label="fit")
xlabel('Ordnung des Minimums', fontsize=20)
ylabel('Winkel in rad', fontsize=20)
legend(fontsize=20)
plt.tight_layout()
savefig('geo1.png')
show()


m1 = optimizedParameters[0]
sigma_m1 = s[0]

b1_geo = lamda/m1
sigma_b1_geo = gauss("lamda/m1")


plot(ordnung5, alpha2, 'x', label='L = 66cm')
optimizedParameters, s = opt.curve_fit(linear, ordnung5, alpha2)
m_std = np.sqrt(np.diag(s))
plt.plot(ordnung5, linear(ordnung5, *optimizedParameters), label="fit")
xlabel('Ordnung des Minimums', fontsize=20)
ylabel('Winkel in rad', fontsize=20)
legend(fontsize=20)
plt.tight_layout()
savefig('geo2.png')
show()

m2 = optimizedParameters[0]
sigma_m2 = s[0]

b2_geo = lamda/m2
sigma_b2_geo = gauss("lamda/m2")



plot(ordnung6, alpha3, 'x', label='L = 47cm')
optimizedParameters, s = opt.curve_fit(linear, ordnung6, alpha3)
m_std = np.sqrt(np.diag(s))
plt.plot(ordnung6, linear(ordnung6, *optimizedParameters), label="fit")
xlabel('Ordnung des Minimums', fontsize=20)
ylabel('Winkel in rad', fontsize=20)
legend(fontsize=20)
plt.tight_layout()
savefig('geo3.png')
show()

m3 = optimizedParameters[0]
sigma_m3 = s[0]

b3_geo = lamda/m3
sigma_b3_geo = gauss("lamda/m2")

b_geo = mean([b1_geo, b2_geo, b3_geo])
sigma_b_geo = std([b1_geo, b2_geo, b3_geo])

print(latexTable(roundCol(array([b_geo]), array([sigma_b_geo]), "\mu m", 6)))



# 2. Photodiode
# 0.5 cm abstände
L4 = 0.376
U = array([0.02, 0.1, 0.19, 0.21, 0.1, 0.01, 0.2, 0.73, 0.15, 0.92, 0.32, 0.16, 2.0, 6.7, 13.8, 13.8, 13.8, 13.8, 13.8, 11.5, 4.87, 1.19, 0.07, 0.53, 1.09, 1.06, 0.57, 0.1, 0.03, 0.2, 0.34, 0.29, 0.12, 0]) # in Volt
U_fit = array([0.02, 0.1, 0.19, 0.21, 0.1, 0.01, 0.2, 0.73, 0.15, 0.92, 0.32, 0.16, 2.0, 6.7, 11.5, 4.87, 1.19, 0.07, 0.53, 1.09, 1.06, 0.57, 0.1, 0.03, 0.2, 0.34, 0.29, 0.12, 0]) # in Volt

l = len(U)
space = numpy.ones(l)
for i in range(l):
    space[i] = 0.005*(i)

space_fit = np.ones(len(U_fit))
space_fit[0:14] = space[0:14]
space_fit[14:] = space[19:]

def sincM(x, I0, b, x0):
    #return I0*(np.sin(math.pi*b*lamda*(x-x0)/L4)/(math.pi*lamda*b*(x-x0)/L4))**2
    return I0*np.sinc(b*(x-x0))**2


plot(space, U, 'x', label='Messwerte')
optimizedParameters, s = opt.curve_fit(sincM, space_fit, U_fit, maxfev=10000)
plt.plot(space, sincM(space, *optimizedParameters), label="fit")
xlabel('Verschiebung in m', fontsize=20)
ylabel('Spannung in V', fontsize=20)
legend(fontsize=20)
plt.tight_layout()
savefig('photo.png')
show()

freq = optimizedParameters[1]
b_photo = lamda*L4*freq/pi





#3 Linse
#Bildgröße
B = 0.0015 #m
sigma_B = 0.001
#Gegenstandsgröße
G = 0.000080
#Bildweite
bw = 1.31
sigma_bw = 0.03
#Gegenstandsweite
gw = 0.12
sigma_gw = 0.03


b_linse = B*gw/bw
sigma_b_linse = gauss("B*gw/bw")
print(latexTable(roundCol(array([b_linse]), array([sigma_b_linse]), "\mu m", 6)))

#Spektrum
L_gitter = 1.25 #m
sigma_L_gitter = 0.03


L_linse = 1.335

L_spalt = 1.45

spaltabstand = 1/(g*100)
sigma_spaltabstand = 0
strich_rot = 0.165
sigma_strich_rot = 0.01
strich_blau = 0.106
sigma_strich_blau = sigma_strich_rot

lamda_rot = spaltabstand*0.5*strich_rot/L_gitter
sigma_lamda_rot = gauss("spaltabstand*0.5*strich_rot/L_gitter")

lamda_blau = spaltabstand*0.5*strich_blau/L_gitter
sigma_lamda_blau = gauss("spaltabstand*0.5*strich_blau/L_gitter")

print(latexTable(roundCol(array([lamda_blau]), array([sigma_lamda_blau]), "n m", 9)))
print(latexTable(roundCol(array([lamda_rot]), array([sigma_lamda_rot]), "n m", 9)))


#Fresnelbeugung Lochblende

L_linsfres = 0.995# zur wand
L_linsfres -= 0.025#brennnweite

W_Fresnel = matrix("""
2 48.5;
3 78.2;
4 83.2;
5 86.2;
6 88.6

""")# anzahl ordnungen, entfernung schirm zu blende

d1 = numpy.ones(5)
for i in range(5):
    d1[i] = L_linsfres - 0.01*W_Fresnel[i, 1]

distance = np.ones(5)
for i in range(5):
    distance[i] = 1/d1[i] + 1/(0.01*W_Fresnel[i, 1])


ordnung5 = array([2, 3, 4, 5, 6])
plot(ordnung5, distance, 'x', label='Messwerte')
optimizedParameters, s = opt.curve_fit(linear, ordnung5, distance, maxfev=10000)
plt.plot(ordnung5, linear(ordnung5, *optimizedParameters), label="fit")
xlabel('Ordnung', fontsize=20)
ylabel('1/a + 1/b in 1/m', fontsize=20)
legend(fontsize=20)
plt.tight_layout()
savefig('fresnel.png')
show()

steig = optimizedParameters[0]
sigma_steig = s[0]
r_fres = sqrt(2*lamda/steig)
sigma_r_fres = gauss("sqrt(2*lamda/steig)")


print(latexTable(roundCol(array([r_fres]), array([sigma_r_fres]), "\mu m", 6)))