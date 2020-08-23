from __future__ import print_function
from ortools.sat.python import cp_model
import numpy as np
import math
import matplotlib.pyplot as plt
import time

def solve_model(W,S,t,N,T,K,loc_x,loc_y,E,Lx,Ly):
    model = cp_model.CpModel()
    e_xi=min(loc_x)
    e_xf=max(loc_x)
    e_yi=min(loc_y)
    e_yf=max(loc_y)
    d=int(round(Ly/(2*K)))
    y = [model.NewIntVar(0, 1, '')for i in range(N)]
    z = [model.NewIntVar(0, T, '')for i in range(N)]
    x = [[[model.NewIntVar(0, 1, '') 
                                for k in range(K)] 
                                for j in range(N)] 
                                for i in range(N)]
    co_x = model.NewIntVar(e_xi,e_xf, '')
    co_y = model.NewIntVar(e_yi,e_yf, '')

    
    for i in range(N-1):
        model.Add((loc_x[i]-(e_xf-e_xi)*(y[i]-1))>co_x)
        model.Add((loc_x[i]+(e_xf-e_xi)*(y[i]-1))<co_x+Lx)
        model.Add((loc_y[i]-(e_yf-e_yi)*(y[i]-1))>co_y)
        model.Add((loc_y[i]+(e_yf-e_yi)*(y[i]-1))<co_y+Ly)
    for i in range(N-1):
        for k in range(K):
            model.Add((loc_y[i]-(e_yf-e_yi)*(sum(x[i][j][k] for j in range(N))-1))>(co_y+(2*k-E[i]+1)*d))
            model.Add((loc_y[i]+(e_yf-e_yi)*(sum(x[i][j][k] for j in range(N))-1))<(co_y+(2*k+E[i]+1)*d))
            model.Add((loc_y[i]-(e_yf-e_yi)*(sum(x[j][i][k] for j in range(N))-1))>(co_y+(2*k-E[i]+1)*d))
            model.Add((loc_y[i]+(e_yf-e_yi)*(sum(x[j][i][k] for j in range(N))-1))<(co_y+(2*k+E[i]+1)*d))
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
    v_o=solver.ObjectiveValue()
    print(v_o)
    return(v_o)
    
def mostrarsol(W,S,t,N,T,K,loc_x,loc_y,E,Lx,Ly,loc_xo,loc_yo,ang):
    model = cp_model.CpModel()
    e_xi=min(loc_x)
    e_xf=max(loc_x)
    e_yi=min(loc_y)
    e_yf=max(loc_y)
    d=int(round(Ly/(2*K)))
    y = [model.NewIntVar(0, 1, '')for i in range(N)]
    z = [model.NewIntVar(0, T, '')for i in range(N)]
    x = [[[model.NewIntVar(0, 1, '') 
                                for k in range(K)] 
                                for j in range(N)] 
                                for i in range(N)]
    co_x = model.NewIntVar(e_xi,e_xf, '')
    co_y = model.NewIntVar(e_yi,e_yf, '')

    
    for i in range(N-1):
        model.Add((loc_x[i]-(e_xf-e_xi)*(y[i]-1))>co_x)
        model.Add((loc_x[i]+(e_xf-e_xi)*(y[i]-1))<co_x+Lx)
        model.Add((loc_y[i]-(e_yf-e_yi)*(y[i]-1))>co_y)
        model.Add((loc_y[i]+(e_yf-e_yi)*(y[i]-1))<co_y+Ly)
    for i in range(N-1):
        for k in range(K):
            model.Add((loc_y[i]-(e_yf-e_yi)*(sum(x[i][j][k] for j in range(N))-1))>(co_y+(2*k-E[i]+1)*d))
            model.Add((loc_y[i]+(e_yf-e_yi)*(sum(x[i][j][k] for j in range(N))-1))<(co_y+(2*k+E[i]+1)*d))
            model.Add((loc_y[i]-(e_yf-e_yi)*(sum(x[j][i][k] for j in range(N))-1))>(co_y+(2*k-E[i]+1)*d))
            model.Add((loc_y[i]+(e_yf-e_yi)*(sum(x[j][i][k] for j in range(N))-1))<(co_y+(2*k+E[i]+1)*d))
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
                    print(loc_xo[b],loc_yo[b])
                if a <N-1 and b <N-1:
                    plt.plot([loc_xo[a],loc_xo[b]],
                             [loc_yo[a],loc_yo[b]],color="C%d"%k)
                a=b
            
    for i in range(N-1):
        if solver.Value(y[i])==1:
            plt.scatter(loc_xo[i],loc_yo[i],color="red")
    xr=solver.Value(co_x)
    yr=solver.Value(co_y)
    x0=int(round(xr*np.cos(-ang)-yr*np.sin(-ang)))
    y0=int(round(xr*np.sin(-ang)+yr*np.cos(-ang)))
    x1=int(round((xr+Lx)*np.cos(-ang)-yr*np.sin(-ang)))
    y1=int(round((xr+Lx)*np.sin(-ang)+yr*np.cos(-ang)))
    
    plt.plot([x0,x1],[y0,y1])
    plt.plot([x0,Ly*np.cos(np.pi/2-ang)+x0],[y0,Ly*np.sin(np.pi/2-ang)+y0],color="black")
    plt.plot([x1,Ly*np.cos(np.pi/2-ang)+x1],[y1,Ly*np.sin(np.pi/2-ang)+y1],color="black")
    for i in range(K+1):
        plt.plot([2*d*i*np.cos(np.pi/2-ang)+x0,2*d*i*np.cos(np.pi/2-ang)+x1],[2*d*i*np.sin(np.pi/2-ang)+y0,2*d*i*np.sin(np.pi/2-ang)+y1],color="black")
                     
                   
    v_o=solver.ObjectiveValue()
    return(v_o)
    
    
def main():
    ang=8
    Ang=np.arange(0,np.pi,np.pi/ang)
    rnd=np.random
    rnd.seed(3)
    Lx= 1000 #lado del grill
    Ly=550
    K=5 #numero de instrumentos 
    T=100 #tiempo total
    n=50 #numero de estrellas
    v=100 #velocidad de movimiento de los intrumentos
    N=n+1
    w=[rnd.randint(1,10) for i in range(n)]
    s=[rnd.randint(5,50) for i in range(n)]
    loc_x=[rnd.randint(1,3*Lx) for i in range(n)]
    loc_y=[rnd.randint(1,3*Ly) for i in range(n)]
    plt.figure()
    plt.scatter(loc_x,loc_y,color="green")

    d=[[math.sqrt((loc_x[j]-loc_x[i])**2+(loc_y[j]-loc_y[i])**2) *(1/v) 
    for i in range(n)] for j in range(n)]
    t = [[0 if n in (i,j) else round(d[i][j]) \
      for j in range(N)] for i in range(N)]
    W = [0 if i==n else w[i] for i in range(N)]
    S = [0 if i==n else s[i] for i in range(N)]
    E1=[rnd.randint(1,10000) for i in range(n)]
    E=[1 if E1[i]<=8000 else 2 for i in range(n)]
    obj=0
    iobj=0
    for i in Ang:
        loc_xt=[]
        loc_yt=[]
        for j in range(n):
            loc_xt.append(int(round(loc_x[j]*np.cos(i)-loc_y[j]*np.sin(i))))
            loc_yt.append(int(round(loc_x[j]*np.sin(i)+loc_y[j]*np.cos(i))))
        Value=solve_model(W,S,t,N,T,K,loc_xt,loc_yt,E,Lx,Ly)
        if Value>obj:
            obj=Value
            iobj=i
    loc_xt=[]
    loc_yt=[]
    for j in range(n):
        loc_xt.append(int(round(loc_x[j]*np.cos(iobj)-loc_y[j]*np.sin(iobj))))
        loc_yt.append(int(round(loc_x[j]*np.sin(iobj)+loc_y[j]*np.cos(iobj))))
    Value=mostrarsol(W,S,t,N,T,K,loc_xt,loc_yt,E,Lx,Ly,loc_x,loc_y,iobj)
    return obj

    return Value



if __name__ == '__main__':
    start_time= time.time()
    Value=main()
    print("\n El valor objetivo es ",Value)
    print("----%s seconds ----" %(time.time()-start_time))

