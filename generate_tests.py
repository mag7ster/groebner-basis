from sympy import *
import random

x, y, z = symbols('x,y,z')
init_printing(use_unicode=False, wrap_line=False)


f = open("tests.txt", "w")

def deg():
    return random.randint(0, 6)

def coef():
    return random.randint(-10, 10)

def print_poly(poly):
    s = str(poly)
    f.write(s.replace('**', '^').replace(' ', '').replace('*', ''))
    f.write('\n')

def print_G(G):
    for poly in G:
        print_poly(poly)

for i in range(1000):

    f.write("test\n")

    G = []
    for _ in range(4):
        G.append(coef() * x**deg()*y**deg()*z**deg() + coef() * x**deg()*y**deg()*z**deg())

    print_G(G)

    f.write("Groebner\n")
    GR = groebner(G, x, y, z, domain=GF(998244353), order='grevlex')
    # print(GR)
    print_G(GR)

    f.write("calc\n")

f.close()
