
from __future__ import print_function
from ortools.sat.python import cp_model
import numpy as np
import math
import matplotlib.pyplot as plt
import time

                 
def solve_model(W,D,E,M,T,K,loc_x,loc_y):
    model = cp_model.CpModel()
    xi=0
    xf=1000
    yi=0
    yf=550
    d=int(round((yf-yi)/5))
    y = [model.NewIntVar(0, 1, '')for i in range(M)]
    z = [model.NewIntVar(0, T, '')for i in range(M)]
    x = [[[model.NewIntVar(0, 1, '') 
                                for k in range(K)] 
                                for j in range(M)] 
                                for i in range(M)]
    
    for i in range(M-1):
        model.Add(y[i]*(loc_x[i]-xi)>=0)
        model.Add(y[i]*(xf-loc_x[i])>=0)
        model.Add(y[i]*(loc_y[i]-yi)>=0)
        model.Add(y[i]*(yf-loc_y[i])>=0)
    for i in range(M-1):
        for k in range(K):
            model.Add(sum(x[i][j][k] for j in range(M))*(loc_y[i]-(yi+k*d))>=0)
            model.Add(sum(x[i][j][k] for j in range(M))*((yi+(k+1)*d)-loc_y[i])>=0)
            model.Add(sum(x[j][i][k] for j in range(M))*(loc_y[i]-(yi+k*d))>=0)
            model.Add(sum(x[j][i][k] for j in range(M))*((yi+(k+1)*d)-loc_y[i])>=0)
    for k in range(K):
        model.Add(sum(x[M-1][j][k] for j in range(M))==1)
        model.Add(sum(x[j][M-1][k] for j in range(M))==1)         
    for i in range(M-1):
            model.Add(sum(x[i][j][k] for j in range(M) for k in range(K))==y[i])
            model.Add(sum(x[j][i][k] for j in range(M) for k in range(K))==y[i])
    for j in range(M-1):
        for i in range(M-1):
            for k in range(K):
                model.Add(z[j]>=z[i]+D[i]+E[i][j]
                *x[i][j][k]-T*(1-x[i][j][k]))
                
    model.Maximize(sum(y[i]*W[i] for i in range(M)))
    solver = cp_model.CpSolver()
    status = solver.Solve(model)
    if status ==cp_model.OPTIMAL: 
        for k in range(K):
            print("\n El recorrido del instrumento %i es:\n"%k)
            a=M-1
            b=0
            while(b!=M-1):
            
                for i in range(M):
                    if solver.Value(x[a][i][k])==1:
                        b=i
                if b<M-1:    
                    print(loc_x[b],loc_y[b])
                if a <M-1 and b <M-1:
                    plt.plot([loc_x[a],loc_x[b]],
                             [loc_y[a],loc_y[b]],color="C%d"%k)
                a=b
#            
    for i in range(M-1):
        if solver.Value(y[i])==1:
            plt.scatter(loc_x[i],loc_y[i],color="red")
#    for i in range(M):
#        for j in range(M):
#            for k in range(K):
#                if solver.Value(x[i][j][k])==1:
#                    a=-1111
#                    b=-1111
#                    c=-1111
#                    d=-1111
#                    if(i<M-1):
#                        a=loc_x[i]
#                        b=loc_y[i]
#                    if(j<M-1):
#                        c=loc_x[j]
#                        d=loc_y[j]                    
#                    print(a,b,c,d,k)
                    
    v_o=solver.ObjectiveValue()
    return(v_o)
    






def main():

    rnd=np.random
    rnd.seed(5)
    e1= 1000 #lado del grill
    e2=550
    K=5 #numero de instrumentos 
    T=70 #tiempo total
    m=50 #numero de estrellas
    v=100 #velocidad de movimiento de los intrumentos
    M=m+1
    w=[rnd.randint(1,10) for i in range(m)]
    d=[rnd.randint(5,50) for i in range(m)]
    loc_x=[rnd.randint(-e1/2,1.5*e1) for i in range(m)]
    loc_y=[rnd.randint(-e2/2,1.5*e2) for i in range(m)]
    k=e2/K
    plt.figure()
    plt.scatter(loc_x,loc_y,color="green")
    for i in range(K+1):
        plt.hlines(k*i,0,1000)
    plt.vlines(0,0,e2)
    plt.vlines(e1,0,e2)
    t=[[math.sqrt((loc_x[j]-loc_x[i])**2+(loc_y[j]-loc_y[i])**2) *(1/v) 
    for i in range(m)] for j in range(m)]
    E = [[0 if m in (i,j) else round(t[i][j]) \
      for j in range(M)] for i in range(M)]
    W = [0 if i==m else w[i] for i in range(M)]
    D = [0 if i==m else d[i] for i in range(M)]

    
    Value=solve_model(W,D,E,M,T,K,loc_x,loc_y)

    return Value



if __name__ == '__main__':
    start_time= time.time()
    Value=main()
    print("\n El valor objetivo es ",Value)
    print("----%s seconds ----" %(time.time()-start_time))

