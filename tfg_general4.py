
from __future__ import print_function
from ortools.sat.python import cp_model
import numpy as np
import math
import matplotlib.pyplot as plt
import time


def tr(x,y,a,l):
    xf=l*np.cos(a)+x
    yf=l*np.sin(a)+y
    f=[xf,yf]
    return f

def solve_model(W,D,E,M,T,K,loc_x,loc_y,V,e1,e2):
    model = cp_model.CpModel()
    xi=min(loc_x)
    xf=max(loc_x)
    yi=min(loc_y)
    yf=max(loc_y)
    d=int(round(e2/K))
    y = [model.NewIntVar(0, 1, '')for i in range(M)]
    z = [model.NewIntVar(0, T, '')for i in range(M)]
    x = [[[model.NewIntVar(0, 1, '') 
                                for k in range(K)] 
                                for j in range(M)] 
                                for i in range(M)]
    co_x = model.NewIntVar(xi,xf, '')
    co_y = model.NewIntVar(yi,yf, '')

    
    for i in range(M-1):
        model.Add((loc_x[i]-(xf-xi)*(y[i]-1))>=co_x)
        model.Add((loc_x[i]+(xf-xi)*(y[i]-1))<=co_x+e1)
        model.Add((loc_y[i]-(yf-yi)*(y[i]-1))>=co_y)
        model.Add((loc_y[i]+(yf-yi)*(y[i]-1))<=co_y+e2)
    for i in range(M-1):
        for k in range(K):
            model.Add((loc_y[i]-(yf-yi)*(sum(x[i][j][k] for j in range(M))-1))>=(co_y+(k-V[i]+1)*d))
            model.Add((loc_y[i]+(yf-yi)*(sum(x[i][j][k] for j in range(M))-1))<=(co_y+(k+V[i])*d))
            model.Add((loc_y[i]-(yf-yi)*(sum(x[j][i][k] for j in range(M))-1))>=(co_y+(k-V[i]+1)*d))
            model.Add((loc_y[i]+(yf-yi)*(sum(x[j][i][k] for j in range(M))-1))<=(co_y+(k+V[i])*d))
    for k in range(K):
        model.Add(sum(x[M-1][j][k] for j in range(M))==1)
        model.Add(sum(x[j][M-1][k] for j in range(M))==1)
        
    for i in range(M-1):
        model.Add(sum(x[i][j][k] for j in range(M) for k in range(K))==y[i]*V[i])
        model.Add(sum(x[j][i][k] for j in range(M) for k in range(K))==y[i]*V[i])
        for k in range(K):
            model.Add(sum(x[i][j][k] for j in range(M))<=1)
            model.Add(sum(x[j][i][k] for j in range(M))<=1)

    for j in range(M-1):
        for i in range(M-1):
            for k in range(K):
                model.Add(z[j]>=z[i]+D[i]+E[i][j]
                *x[i][j][k]-T*(1-x[i][j][k]))
                
    model.Maximize(sum(y[i]*W[i] for i in range(M)))
    solver = cp_model.CpSolver()
    status = solver.Solve(model)
#    if status ==cp_model.OPTIMAL: 
#        for k in range(K):
#            print("\n El recorrido del instrumento %i es:\n"%k)
#            a=M-1
#            b=0
#            while(b!=M-1):
#            
#                for i in range(M):
#                    if solver.Value(x[a][i][k])==1:
#                        b=i
#                if b<M-2:    
#                    print(loc_x[b],loc_y[b])
#                if a <M-2 and b <M-2:
#                    plt.plot([loc_x[a],loc_x[b]],
#                             [loc_y[a],loc_y[b]],color="C%d"%k)
#                a=b
#            
#    for i in range(M-1):
#        if solver.Value(y[i])==1:
#            plt.scatter(loc_x[i],loc_y[i],color="red")
#    x0=solver.Value(co_x)
#    y0=solver.Value(co_y)
#    plt.plot([x0,x0],[y0,y0+e2],color="black")
#    plt.plot([x0+e1,x0+e1],[y0,y0+e2],color="black")
#    for i in range(K+1):
#        plt.plot([x0,x0+e1],[y0+d*i,y0+d*i],color="black")
                     
                   
    v_o=solver.ObjectiveValue()
    print(v_o)
    return(v_o)
    
def mostrarsol(W,D,E,M,T,K,loc_x,loc_y,V,e1,e2,loc_xo,loc_yo,ang):
    model = cp_model.CpModel()
    xi=min(loc_x)
    xf=max(loc_x)
    yi=min(loc_y)
    yf=max(loc_y)
    d=int(round(e2/K))
    y = [model.NewIntVar(0, 1, '')for i in range(M)]
    z = [model.NewIntVar(0, T, '')for i in range(M)]
    x = [[[model.NewIntVar(0, 1, '') 
                                for k in range(K)] 
                                for j in range(M)] 
                                for i in range(M)]
    co_x = model.NewIntVar(xi,xf, '')
    co_y = model.NewIntVar(yi,yf, '')

    
    for i in range(M-1):
        model.Add((loc_x[i]-(xf-xi)*(y[i]-1))>=co_x)
        model.Add((loc_x[i]+(xf-xi)*(y[i]-1))<=co_x+e1)
        model.Add((loc_y[i]-(yf-yi)*(y[i]-1))>=co_y)
        model.Add((loc_y[i]+(yf-yi)*(y[i]-1))<=co_y+e2)
    for i in range(M-1):
        for k in range(K):
            model.Add((loc_y[i]-(yf-yi)*(sum(x[i][j][k] for j in range(M))-1))>=(co_y+(k-V[i]+1)*d))
            model.Add((loc_y[i]+(yf-yi)*(sum(x[i][j][k] for j in range(M))-1))<=(co_y+(k+V[i])*d))
            model.Add((loc_y[i]-(yf-yi)*(sum(x[j][i][k] for j in range(M))-1))>=(co_y+(k-V[i]+1)*d))
            model.Add((loc_y[i]+(yf-yi)*(sum(x[j][i][k] for j in range(M))-1))<=(co_y+(k+V[i])*d))
    for k in range(K):
        model.Add(sum(x[M-1][j][k] for j in range(M))==1)
        model.Add(sum(x[j][M-1][k] for j in range(M))==1)
        
    for i in range(M-1):
        model.Add(sum(x[i][j][k] for j in range(M) for k in range(K))==y[i]*V[i])
        model.Add(sum(x[j][i][k] for j in range(M) for k in range(K))==y[i]*V[i])
        for k in range(K):
            model.Add(sum(x[i][j][k] for j in range(M))<=1)
            model.Add(sum(x[j][i][k] for j in range(M))<=1)

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
                    print(loc_xo[b],loc_yo[b])
                if a <M-1 and b <M-1:
                    plt.plot([loc_xo[a],loc_xo[b]],
                             [loc_yo[a],loc_yo[b]],color="C%d"%k)
                a=b
            
    for i in range(M-1):
        if solver.Value(y[i])==1:
            plt.scatter(loc_xo[i],loc_yo[i],color="red")
    xr=solver.Value(co_x)
    yr=solver.Value(co_y)
    x0=int(round(xr*np.cos(-ang)-yr*np.sin(-ang)))
    y0=int(round(xr*np.sin(-ang)+yr*np.cos(-ang)))
    x1=int(round((xr+e1)*np.cos(-ang)-yr*np.sin(-ang)))
    y1=int(round((xr+e1)*np.sin(-ang)+yr*np.cos(-ang)))
    
    plt.plot([x0,x1],[y0,y1])
    plt.plot([x0,e2*np.cos(np.pi/2-ang)+x0],[y0,e2*np.sin(np.pi/2-ang)+y0],color="black")
    plt.plot([x1,e2*np.cos(np.pi/2-ang)+x1],[y1,e2*np.sin(np.pi/2-ang)+y1],color="black")
    for i in range(K+1):
        plt.plot([d*i*np.cos(np.pi/2-ang)+x0,d*i*np.cos(np.pi/2-ang)+x1],[d*i*np.sin(np.pi/2-ang)+y0,d*i*np.sin(np.pi/2-ang)+y1],color="black")
                     
                   
    v_o=solver.ObjectiveValue()
    return(v_o)
    
    
    
    
    
    



def main():
    ang=8
    Ang=np.arange(0,np.pi,np.pi/ang)
    rnd=np.random
    rnd.seed(2)
    e1= 1000 #lado del grill
    e2=550
    K=5 #numero de instrumentos 
    T=100 #tiempo total
    m=25 #numero de estrellas
    v=100 #velocidad de movimiento de los intrumentos
    M=m+1
    w=[rnd.randint(1,10) for i in range(m)]
    d=[rnd.randint(5,50) for i in range(m)]
    loc_x=[rnd.randint(1,3*e1) for i in range(m)]
    loc_y=[rnd.randint(1,3*e2) for i in range(m)]
    plt.figure()
    plt.xlim(max(max(loc_x),max(loc_y)))
    plt.ylim(max(max(loc_x),max(loc_y)))
    plt.scatter(loc_x,loc_y,color="green")
    

    t=[[math.sqrt((loc_x[j]-loc_x[i])**2+(loc_y[j]-loc_y[i])**2) *(1/v) 
    for i in range(m)] for j in range(m)]
    E = [[0 if m in (i,j) else round(t[i][j]) \
      for j in range(M)] for i in range(M)]
    W = [0 if i==m else w[i] for i in range(M)]
    D = [0 if i==m else d[i] for i in range(M)]
    V1=[rnd.randint(1,10000) for i in range(m)]
    V=[1 if V1[i]<=8000 else 2 for i in range(m)]
    
    obj=0
    iobj=0
    for i in Ang:
        loc_xt=[]
        loc_yt=[]
        for j in range(m):
            loc_xt.append(int(round(loc_x[j]*np.cos(i)-loc_y[j]*np.sin(i))))
            loc_yt.append(int(round(loc_x[j]*np.sin(i)+loc_y[j]*np.cos(i))))
        Value=solve_model(W,D,E,M,T,K,loc_xt,loc_yt,V,e1,e2)
        if Value>obj:
            obj=Value
            iobj=i
    loc_xt=[]
    loc_yt=[]
    for j in range(m):
        loc_xt.append(int(round(loc_x[j]*np.cos(iobj)-loc_y[j]*np.sin(iobj))))
        loc_yt.append(int(round(loc_x[j]*np.sin(iobj)+loc_y[j]*np.cos(iobj))))
    Value=mostrarsol(W,D,E,M,T,K,loc_xt,loc_yt,V,e1,e2,loc_x,loc_y,iobj)
    return obj




if __name__ == '__main__':
    start_time= time.time()
    Value=main()
    print("\n El valor objetivo es ",Value)
    print("----%s seconds ----" %(time.time()-start_time))
    

