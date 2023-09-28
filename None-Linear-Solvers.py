from typing import List

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from pandas import DataFrame

f = lambda x: np.power(x,3)+10*np.cos(2*x)+np.log(11+x) #funksjon f(x) fra oppgave 1

fprime = lambda x: 3*np.power(x,2)-20*np.sin(2*x)+(1/(11+x))

xpts = np.linspace(-1, 2, 1000) # Setting x range for plot

#plot
plt.plot(xpts, f(xpts), label ='f(x)')
plt.xlabel("x")
plt.ylabel("f(x)")
plt.axhline(0, linewidth=1, color='r', linestyle='-')
plt.grid(color='grey', linestyle='--', linewidth=0.1)
plt.title("plot 1.1")
plt.show()




def bisection_method(f, a, b, tol=0.5e-4):
    '''Tar inn funksjon f, start intervall [a,b] og ønsket nøyaktighet (stop kriterie).
    tester om det finnes røtter innenfor ingervall med skjærings sentingen f(a)xf(b) < 0.
    gitt det finnes røtter halverer den inervallet til (a-b)/2 > tol. Returnerer pandas DataFrame
    med iterasjoner og næyaktighet per iterasjon'''

    if f(a)*f(b) > 0:
        print('No roots on this intervall. Try a new guess')
    else:

        approx_root = []
        accuracy = []
        a = a
        b = b
        while (b - a) / 2.0 > tol:
            c = (a + b) / 2.0
            accuracy.append(abs(f(c)))
            approx_root.append(c)
            if f(a) * f(c) < 0:
                b = c
            else:
                a = c

        return pd.DataFrame({'Bisection Approximated root': approx_root,  'f(c) - [accuracy of root]': accuracy})


def sekant_method(f, x0, x1 = None, tol = 0.5e-4):

    approx_root = []
    accuracy = []

    while abs((x1 - (x1 - x0) / (f(x1) - f(x0)) * f(x1)) - x1) > tol:
        x_new = x1 - (x1 - x0) / (f(x1) - f(x0)) * f(x1)
        accuracy.append(abs(f(x_new)))
        approx_root.append(x_new)

        x0 = x1
        x1 = x_new

    return pd.DataFrame({'Sekant Approximated root': approx_root, 'f(c) - [accuracy of root]': accuracy})

def newtons_method(f, fprime, x0, tol = 0.5e-4):

    approx_root = []
    accuracy = []

    while abs((x0 - f(x0) / fprime(x0)) - x0) > tol:
        x_new = x0 - f(x0) / fprime(x0)
        accuracy.append(abs(f(x_new)))
        approx_root.append(x_new)
        x0 = x_new

    return pd.DataFrame({'Newtons Approximated root': approx_root, 'f(c) - [accuracy of root]': accuracy})


intervalls = [[-1.0,-0.5], [0.8,1.0],[1.6,2.0]]

#The roots for bisection method
biseect_root_1=bisection_method(f=f, a=intervalls[0][0],b=intervalls[0][1],tol=0.5e-10)
biseect_root_2=bisection_method(f=f, a=intervalls[1][0],b=intervalls[1][1],tol=0.5e-10)
biseect_root_3=bisection_method(f=f, a=intervalls[2][0],b=intervalls[2][1],tol=0.5e-10)

bisection_roots_consolidated_frame = pd.concat([biseect_root_1, biseect_root_2, biseect_root_3], axis=1)
print(bisection_roots_consolidated_frame)

Sekant_roots = [ sekant_method(f, n[0], n[1], tol=0.5e-10) for n in intervalls ]

Newton_roots = [newtons_method(f,fprime, n[0], tol=0.5e-10) for n in intervalls ]

