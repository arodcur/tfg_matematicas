
from __future__ import print_function
from ortools.sat.python import cp_model
import numpy as np
import math
import matplotlib.pyplot as plt
import time
                 
def solve_model(W,S,t,N,T,F,K,loc_x,loc_y):
    model = cp_model.CpModel()
    y = [model.NewIntVar(0, 1, '')for i in range(N)]
    z = [model.NewIntVar(0, T, '')for i in range(N)]
    x = [[[model.NewIntVar(0, 
                           0 if k not in F[i] or k not in F[j] 
                           else 1, '') 
                            for k in range(K)] 
                            for j in range(N)] 
                            for i in range(N)]
    for k in range(K):
        model.Add(sum(x[N-1][j][k] for j in range(N))==1)
        model.Add(sum(x[j][N-1][k] for j in range(N))==1)         
    for i in range(N-1):
        for f in F[i]:
            model.Add(sum(x[i][j][f] for j in range(N))==y[i])
            model.Add(sum(x[j][i][f] for j in range(N))==y[i])
    for j in range(N-1):
        for i in range(N-1):
            for f in F[i]:
                model.Add(z[j]>=z[i]+S[i]+t[i][j]
                *x[i][j][f]-T*(1-x[i][j][f]))
    model.Maximize(sum(y[i]*W[i] for i in range(N)))
    solver = cp_model.CpSolver()
    status = solver.Solve(model)
    if status ==cp_model.OPTIMAL: 
        for k in range(K):
            print("\n El recorrido del instrumento %i es:\n"%k)
            a=N-1
            b=0
            while(b!=N-1):
            
                for i in range(N):
                    if solver.Value(x[a][i][k])==1:
                        b=i
                if b<N-1:    
                    print("(" + str(loc_x[b]) + "," + str(loc_y[b]) + 
                            ") En el momento " + str(solver.Value(z[b])))
                if a <N-1 and b <N-1:
                    plt.plot([loc_x[a],loc_x[b]],
                             [loc_y[a],loc_y[b]],color="C%d"%k)
                a=b
            
            
    for i in range(N-1):
        if solver.Value(y[i])==1:
            plt.scatter(loc_x[i],loc_y[i],color="red")
    v_o=solver.ObjectiveValue()   
    return(v_o)






def main():

    rnd=np.random
    rnd.seed(2)
    e= 1000 #lado del grill
    K=5 #numero de instrumentos 
    T=70 #tiempo total
    n=25 #numero de estrellas
    v=100 #velocidad de movimiento de los intrumentos
    N=n+1
    w=[rnd.randint(1,10) for i in range(n)]
    s=[rnd.randint(5,50) for i in range(n)]
    loc_x=[rnd.randint(1,e) for i in range(n)]
    loc_y=[rnd.randint(1,e) for i in range(n)]
    k=e/K #altura de cada instrumento
    f=[math.ceil(loc_y[i]/k)-1 for i in range(n)]
    rand1=[rnd.randint(1,10000) for i in range(n)]
    g=[f[i]+1 if rand1[i]<=5000 else f[i]-1 for i in range(n)]
    for i in range(n):
        if g[i]==-1:
            g[i]=1
        if g[i]==K:
            g[i]=K-2
    rand2=[rnd.randint(1,10000) for i in range(n)]
    F=[[j for j in range(K)] if i==n else [f[i]] if rand2[i]<=6000 
        else sorted([f[i],g[i]])  for i in range(N)]
#    Las siguientes línes se usaron para crear la tabla
#    from prettytable import PrettyTable
#    tabla = PrettyTable(['Coordenadas', 'Beneficio','Tiempo de observación',
#                     'Instrumentos que deben observarlo'])
#
#    for i in range(n):
#        tabla.add_row([[loc_x[i],loc_y[i]], w[i], s[i], F[i]])
#    print(tabla)
    plt.figure()
    plt.scatter(loc_x,loc_y,color="green")
    for i in range(K+1):
        plt.hlines(k*i,0,1000)
    plt.vlines(0,0,1000)
    plt.vlines(1000,0,1000)
    d=[[math.sqrt((loc_x[j]-loc_x[i])**2+(loc_y[j]-loc_y[i])**2) *(1/v) 
    for i in range(n)] for j in range(n)]
    t = [[0 if n in (i,j) else round(d[i][j]) \
      for j in range(N)] for i in range(N)]
    W = [0 if i==n else w[i] for i in range(N)]
    S = [0 if i==n else s[i] for i in range(N)]

    
    Value=solve_model(W,S,t,N,T,F,K,loc_x,loc_y)

    return Value



if __name__ == '__main__':
    start_time= time.time()
    Value=main()
    print("\n El valor objetivo es ",Value)
    print("----%s seconds ----" %(time.time()-start_time))
    


