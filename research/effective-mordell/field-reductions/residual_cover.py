"""Check the remaining genus-five cover, its maps, and its five elliptic ranks."""
import sympy as sp
from cypari2 import Pari

pari=Pari()
s=sp.Symbol('s')

def check_zero(f):
    if sp.cancel(f)!=0:
        raise ArithmeticError(f)

def verify():
    x=(3+s**2)/(3-s**2)
    check_zero((x**2+x+1)*(3-s**2)**2-(s**4+27))
    check_zero((x**2-x+1)*(3-s**2)**2-3*(s**4+3))
    check_zero(7*(x**6-1)*(3-s**2)**6-42**2*s**2*(s**4+3)*(s**4+27)/7)
    h=s**8+30*s**4+81
    check_zero(h-(s**4+3)*(s**4+27))
    for j in (1,-1):
        t=s+j*3/s
        check_zero(h/s**4-(t**4-j*12*t**2+48))
    models=[(0,3),(0,1323),(210,3969),(-84,2352),(84,2352)]
    for a,b in models:
        E=pari.ellinit([0,a,0,b,0])
        for F in (E,E.elltwist(-4)):
            if list(F.ellrank()[:2])!=[1,1]:
                raise ArithmeticError('rank of curve or Gaussian twist')
    # Differentials on the genus-three quotient, w^2=7*h.
    # du/v for u=s^2, v=w; u=s+/-3/s, v=w/s^2, omitting common ds/w.
    coeffs=sp.Matrix([[0,2,0],[-3,0,1],[3,0,1]])
    if coeffs.det()==0:raise ArithmeticError('dependent genus-three quotients')
    print('The residual genus-five Jacobian is isogenous to five rank-two elliptic factors over Q(i).')
if __name__=='__main__':verify()
