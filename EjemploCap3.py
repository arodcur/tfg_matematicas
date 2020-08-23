
from __future__ import print_function
from ortools.sat.python import cp_model
import numpy as np
import math
import matplotlib.pyplot as plt
import time

                 
def solve_model(W,S,t,N,T,K,loc_x,loc_y,E,C_xi,C_xf,C_yi,C_yf):
    model = cp_model.CpModel()
    d=int(round((C_yf-C_yi)/(2*K)))
    y = [model.NewIntVar(0, 1, '')for i in range(N)]
    z = [model.NewIntVar(0, T, '')for i in range(N)]
    x = [[[model.NewIntVar(0, 1, '') 
                                for k in range(K)] 
                                for j in range(N)] 
                                for i in range(N)]
    
    for i in range(N-1):
        model.Add(y[i]*(loc_x[i]-C_xi)>=0)
        model.Add(y[i]*(C_xf-loc_x[i])>=0)
        model.Add(y[i]*(loc_y[i]-C_yi)>=0)
        model.Add(y[i]*(C_yf-loc_y[i])>=0)
    for i in range(N-1):
        for k in range(K):
            model.Add(sum(x[i][j][k] for j in range(N))
                      *(loc_y[i]-(C_yi+(2*k-E[i]+1)*d))>0)
            model.Add(sum(x[i][j][k] for j in range(N))
                      *((C_yi+(2*k+E[i]+1)*d)-loc_y[i])>0)
            model.Add(sum(x[j][i][k] for j in range(N))
                      *(loc_y[i]-(C_yi+(2*k-E[i]+1)*d))>0)
            model.Add(sum(x[j][i][k] for j in range(N))
                      *((C_yi+(2*k+E[i]+1)*d)-loc_y[i])>0)
    for k in range(K):
        model.Add(sum(x[N-1][j][k] for j in range(N))==1)
        model.Add(sum(x[j][N-1][k] for j in range(N))==1)
        
    for i in range(N-1):
        model.Add(sum(x[i][j][k] for j in range(N) 
                  for k in range(K))==y[i]*E[i])
        model.Add(sum(x[j][i][k] for j in range(N) 
                  for k in range(K))==y[i]*E[i])
        for k in range(K):
            model.Add(sum(x[i][j][k] for j in range(N))<=1)
            model.Add(sum(x[j][i][k] for j in range(N))<=1)
    for j in range(N-1):
        for i in range(N-1):
            for k in range(K):
                model.Add(z[j]>=z[i]+S[i]+t[i][j]
                *x[i][j][k]-T*(1-x[i][j][k]))
                
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
                    print(loc_x[b],loc_y[b])
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
    rnd.seed(5)
    K=5 #numero de instrumentos 
    T=70 #tiempo total
    n=50 #numero de estrellas
    v=100 #velocidad de movimiento de los intrumentos
    C_xi=0
    C_xf=1000
    C_yi=0
    C_yf=550
    e1=C_xf-C_xi
    e2=C_yf-C_yi
    N=n+1
    w=[rnd.randint(1,10) for i in range(n)]
    s=[rnd.randint(5,50) for i in range(n)]
    loc_x=[rnd.randint(-e1/2,1.5*e1) for i in range(n)]
    loc_y=[rnd.randint(-e2/2,1.5*e2) for i in range(n)]
    k=e2/K
    plt.figure()
    plt.scatter(loc_x,loc_y,color="green")
    for i in range(K+1):
        plt.hlines(k*i,0,1000)
    plt.vlines(0,0,e2)
    plt.vlines(e1,0,e2)
    d=[[math.sqrt((loc_x[j]-loc_x[i])**2+(loc_y[j]-loc_y[i])**2) *(1/v) 
    for i in range(n)] for j in range(n)]
    t = [[0 if n in (i,j) else round(d[i][j]) \
      for j in range(N)] for i in range(N)]
    W = [0 if i==n else w[i] for i in range(N)]
    S = [0 if i==n else s[i] for i in range(N)]
    E1=[rnd.randint(1,10000) for i in range(n)]
    E=[1 if E1[i]<=8000 else 2 for i in range(n)]
    
    Value=solve_model(W,S,t,N,T,K,loc_x,loc_y,E,C_xi,C_xf,C_yi,C_yf)

    return Value



if __name__ == '__main__':
    start_time= time.time()
    Value=main()
    print("\n El valor objetivo es ",Value)
    print("----%s seconds ----" %(time.time()-start_time))

