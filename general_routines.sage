def mono2tuple(mono,gens=None):
    e = mono.exponents()[0]
    a = []
    if gens is not None:
        mapa = lambda i:gens[i]
    else:
        mapa = lambda i:i
    for i in range(len(e)):
        if e[i] != 0:
            a = a + ([mapa(i)]*e[i])
    a.reverse()
    return a

def poly2tuplelist(poly,gens):
    lista = []
    for m in poly.monomials():
        tupla = mono2tuple(m,gens)
        mult = abs(int(poly.coefficient(m)))
        
        lista += [tupla]*mult
    return lista

def _tuple2gen(i,j):
    return 'x'+str(i)+'_'+ str(j)

def getGenerators(N,d):
    gens = [(1,0)]
    r = 2
    for w in range(2,N+1):
        gens += [(w,k) for k in range(0,min(d,w-r)+1)]
    return gens
    
def degree_x1(mon):
    return mon.exponents()[0][0]
def remove_x1(poly):
    return sum([poly.coefficient(mon)*mon for mon in poly.monomials() if degree_x1(mon)==0 ])

def monFrom(weight):
    return R(_tuple2gen(weight,0))
    


