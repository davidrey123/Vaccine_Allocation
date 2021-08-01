import time
import numpy as np

def read_data_node(filename):
    rty = np.sin(1)
    data = open(filename,'r')
    lines = data.readlines()
    lines = [line.rstrip('\n') for line in lines]
    data.close() 
    N = {}
    P = {}
    PM = {}
    node2dm = {}
    for line in lines:
        ls = line.split(',')
        i = int(ls[0])
        node2dm[i] = ls[1]
        P[i] = float(ls[2])
        PM[i] = int(ls[3])
        N[i] = []
        for j in range(4,len(ls)):
            N[i].append(int(ls[j]))    
    Ptot = sum(P[i] for i in N)
    print('Total population',Ptot)
    return N,P,PM,node2dm
    
def read_data_mobility(filename):
    data = open(filename,'r')
    lines = data.readlines()
    lines = [line.rstrip('\n') for line in lines]
    data.close()
    p = {}
    same = {}
    for line in lines:
        ls = line.split(',')    
        i = int(ls[0])    
        j = int(ls[1])    
        p[i,j] = float(ls[2])
        same[i,j] = int(ls[3])
    print('nb of arcs',len(p))
    return p,same
 
def read_data_dm(filename): 
    data = open('data_dm.txt','r')
    lines = data.readlines()
    lines = [line.rstrip('\n') for line in lines]
    data.close()  
    dm = {}
    popdm = {}
    S0 = {}
    I0 = {}
    bg = {}
    cap = {}
    for line in lines:
        ls = line.split(',')
        k = ls[0]
        popdm[k] = int(ls[1])
        S0[k] = int(ls[2])
        I0[k] = int(ls[3])
        bg[k] = (float(ls[4]),float(ls[5]))
        cap[k] = float(ls[6])
        dm[k] = []
        for i in range(7,len(ls)):
            dm[k].append(int(ls[i]))
    print('nb of dm',len(dm))
    return dm,popdm,S0,I0,bg,cap
    
def initialization(N,node2dm,popdm,S0,I0):
    S = {1:{i:0 for i in N}}
    I = {1:{i:0 for i in N}}
    R = {1:{i:0 for i in N}}
    D = {1:{i:0 for i in N}}    
    
    for i in N:
        k = node2dm[i]   
        S[1][i] = S0[k]/popdm[k]
        I[1][i] = I0[k]/popdm[k]    
        R[1][i] = 1.0 - S[1][i] - I[1][i]  
        D[1][i] = 0.0
        
    #---beta distribution (prior)
    a = {i:1.0 for i in N} 
    b = {i:1.0 for i in N} 
    
    #---node allocation upper bounds
    xub = {i:1.0 for i in N}
    
    #---total flow to population ratio
    ratio = 0.00927040328665764
    
    #---death rate
    lbd = 0.01
        
    return S,I,R,D,a,b,xub,ratio,lbd   

def initialization2(N,node2dm,popdm,S0,I0):
    S = {1:{i:0 for i in N}}
    I = {1:{i:0 for i in N}}
    R = {1:{i:0 for i in N}}
    D = {1:{i:0 for i in N}}    
    
    for i in N:
        k = node2dm[i]   
        S[1][i] = S0[k]/popdm[k]
        I[1][i] = I0[k]/popdm[k]    
        R[1][i] = 1.0 - S[1][i] - I[1][i]  
        D[1][i] = 0.0
        
    #---initial efficiency rate
    obs = {i:[0.5] for i in N}

    #---moving average size
    ma_size = 5
    
    #---node allocation upper bounds
    xub = {i:1.0 for i in N}
    
    #---total flow to population ratio
    ratio = 0.00927040328665764
    
    #---death rate
    lbd = 0.01   

    return S,I,R,D,obs,ma_size,xub,ratio,lbd
      
def simulation(t,x,te,N,S,I,R,D,node2dm,bg,lbd,p,ratio,dm):
    #---update local compartments using SIRD dynamics
    S[t+1] = {}
    I[t+1] = {}
    R[t+1] = {}
    D[t+1] = {}
    BS = {}
    pdmI = {k:{kp:0.0 for kp in dm} for k in dm}
    
    for k in dm:
        beta = bg[k][0]
        gamma = bg[k][1] 
        II = 0.0
        EI = 0.0
        for i in dm[k]:
            Stemp = 0.0
            Itemp = 0.0
            Rtemp = 0.0
            for j in N[i]:            
                Stemp += p[i,j]*(S[t][j]*(1-te[j]*x[j]) - S[t][i]*(1-te[i]*x[i]))
                Itemp += p[i,j]*(I[t][j] - I[t][i])
                Rtemp += p[i,j]*(R[t][j] + S[t][j]*te[j]*x[j] - R[t][i] - S[t][i]*te[i]*x[i])
                kp = node2dm[j]
                if kp == k:
                    II += ratio*p[i,j]*I[t][j]
                else:
                    pdmI[kp][k] += ratio*p[i,j]*I[t][j]
                    EI += ratio*p[i,j]*I[t][j]
            II += I[t][i] + beta*S[t][i]*I[t][i]*(1-te[i]*x[i]) - gamma*I[t][i]

            S[t+1][i] = S[t][i]*(1-te[i]*x[i]) - beta*S[t][i]*I[t][i]*(1-te[i]*x[i]) + ratio*Stemp
            I[t+1][i] = I[t][i] + beta*S[t][i]*I[t][i]*(1-te[i]*x[i]) - gamma*I[t][i] + ratio*Itemp
            R[t+1][i] = R[t][i] + S[t][i]*te[i]*x[i] + (1-lbd)*gamma*I[t][i] + ratio*Rtemp
            D[t+1][i] = 1.0 - S[t+1][i] - I[t+1][i] - R[t+1][i]
        if II < -1e-6:
            print(k,'negative II',II)        
        if EI+II < 1e-6:
            BS[k] = 0.0
        else:
            BS[k] = EI/(EI+II)
    return S,I,R,D,BS,pdmI  
    
def read_scenario(filename,dm,N):
    data = open(filename,'r')
    lines = data.readlines()
    lines = [line.rstrip('\n') for line in lines]
    data.close() 

    T = len(lines)-2

    #--dm-based data
    theta_true = {}

    #---node-based data
    theta_hat = {t:{} for t in range(1,T+1)}
    theta_env = {t:{} for t in range(1,T+1)}

    ls = lines[1].split(' ')
    cnt = 0
    for k in dm:
        theta_true[k] = float(ls[cnt])
        cnt += 1
    eps = float(ls[-1])

    for ind in range(2,T+2):
        ls = lines[ind].split(' ')
        cnt = 0
        for i in N:
            theta_env[ind-1][i] = float(ls[cnt])
            cnt += 1

    return theta_true,theta_env,theta_hat,T,eps
    
def set_budget(dm,P,cap,cB):
    #---per-period budget    
    B0 = {}
    for k in dm:
        Bmax = sum(P[i] for i in dm[k])
        B0[k] = cB*cap[k]*Bmax 
    return B0
    
def share_budget(t,dm,B0,BS,pdmI,cap,D,mode):
    B = {}
    
    if mode == 'none':
        B = B0
        
    elif mode == 'flow':
        den = {k:0.0 for k in dm}
        for k in dm:
            B[k] = B0[k]*(1 - BS[k])
            for kp in dm:
                den[k] += pdmI[kp][k]      
        for k in dm:
            if den[k] < 1e-6: continue            
            for kp in dm:
                if k == kp: continue
                B[kp] += B0[k]*BS[k]*pdmI[kp][k]/den[k]                
        
    elif mode == 'flow2cap':
        den = {k:0.0 for k in dm}
        for k in dm:
            B[k] = B0[k]*(1 - BS[k])
            for kp in dm:
                den[k] += pdmI[kp][k]/cap[kp]     
        for k in dm:
            if den[k] < 1e-6: continue
            for kp in dm:
                if k == kp: continue
                B[kp] += B0[k]*BS[k]*(pdmI[kp][k]/cap[kp])/den[k]
   
    elif mode == 'flow2capD':
        den = {k:0.0 for k in dm}
        for k in dm:
            Dkt = sum(D[t][i] for i in dm[k])
            Dktt = sum(D[t+1][i] for i in dm[k])
            if Dkt - Dktt <= 0:
                BS[k] = 0.0
            B[k] = B0[k]*(1 - BS[k])
            for kp in dm:
                den[k] += pdmI[kp][k]/cap[kp]     
        for k in dm:
            if den[k] < 1e-6: continue
            for kp in dm:
                if k == kp: continue
                B[kp] += B0[k]*BS[k]*(pdmI[kp][k]/cap[kp])/den[k]                
    
    else:
        print('incorrect mode')
        
    return B
    
def set_PM(PM,cPM):
    PM = {i:int(np.ceil(cPM*PM[i])) for i in PM}
    return PM
    
def optimize(k,t,xub,th,N,dm,S,I,p,same,bg,ratio,P,B):
    loss = {i:0.0 for i in dm[k]}
    for i in dm[k]:
        loss[i] = loss[i] - S[t][i]*th[i] + bg[k][0]*S[t][i]*I[t][i]*th[i]
        for j in N[i]:
            loss[i] = loss[i] + ratio*p[i,j]*S[t][i]*th[i]
            if same[i,j] == 1:
                loss[j] = loss[j] - ratio*p[i,j]*S[t][j]*th[j]

    v2w = {i:(-loss[i]/P[i]) for i in dm[k]}
    dm_N_sorted = sorted([i for i in dm[k]], key=lambda X: v2w[X], reverse=True)
    x = {i:0.0 for i in dm[k]}
    A = B[k]
    for i in dm_N_sorted:
        if P[i]*xub[i] <= A:
            x[i] = xub[i]
        else:
            x[i] = A/P[i]
        A = A - P[i]*x[i]
        if A < 1e-3: break
    return x  

def optimize_PB(k,t,xub,dm,P,B):
    dm_N_sorted = sorted([i for i in dm[k]], key=lambda X: P[X], reverse=True)
    x = {i:0.0 for i in dm[k]}
    A = B[k]
    for i in dm_N_sorted:
        if P[i]*xub[i] <= A:
            x[i] = xub[i]
        else:
            x[i] = A/P[i]
        A = A - P[i]*x[i]
        if A < 1e-3: break  
    return x   
    
def write_data_dm(filename,dm,P,S,I,R,D,T):
    Sopt = {k:[] for k in dm}
    Iopt = {k:[] for k in dm}
    Ropt = {k:[] for k in dm}
    Dopt = {k:[] for k in dm}

    for k in dm:
        for t in range(1,T+2):        
            Sopt[k].append(sum(P[i]*S[t][i] for i in dm[k]))
            Iopt[k].append(sum(P[i]*I[t][i] for i in dm[k]))
            Ropt[k].append(sum(P[i]*R[t][i] for i in dm[k]))
            Dopt[k].append(sum(P[i]*D[t][i] for i in dm[k]))

    file = open(filename,'w')
    for k in dm:
        file.write('%s' % k)
        for ind in range(T+1):
            ss = Sopt[k][ind]
            file.write(',%.1f' % ss)
        file.write('\n')
        file.write('%s' % k)
        for ind in range(T+1):
            ii = Iopt[k][ind]
            file.write(',%.1f' % ii)
        file.write('\n')
        file.write('%s' % k)
        for ind in range(T+1):
            rr = Ropt[k][ind]
            file.write(',%.1f' % rr)
        file.write('\n')
        file.write('%s' % k)
        for ind in range(T+1):
            dd = Dopt[k][ind]
            file.write(',%.1f' % dd)
        file.write('\n')    
    file.close()
    return
    
def write_data_node(filename,N,P,S,I,R,D,T):
    file = open(filename,'w')
    for i in N:
        file.write('%d' % i)
        for t in range(1,T+2):
            ss = P[i]*S[t][i]
            file.write(',%.1f' % ss)
        file.write('\n')
        file.write('%d' % i)
        for t in range(1,T+2):
            ii = P[i]*I[t][i]
            file.write(',%.1f' % ii)
        file.write('\n')
        file.write('%d' % i)
        for t in range(1,T+2):
            rr = P[i]*R[t][i]
            file.write(',%.1f' % rr)
        file.write('\n')
        file.write('%d' % i)
        for t in range(1,T+2):
            dd = P[i]*D[t][i]
            file.write(',%.1f' % dd)
        file.write('\n')    
    file.close()
    return
   
def optimize_scenario_TS(filename,N,P,PM,node2dm,p,same,dm,popdm,S0,I0,bg,cap,cB,cPM,mode):
    t0 = time.time()
    print(filename)
    theta_true,theta_env,theta_hat,T,eps = read_scenario(filename,dm,N)
    S,I,R,D,a,b,xub,ratio,lbd = initialization(N,node2dm,popdm,S0,I0)
    B0 = set_budget(dm,P,cap,cB)
    PM = set_PM(PM,cPM)
    xall = {}
    B = B0

    for t in range(1,T+1):
        tstart = time.time()        
        xdm = {k:{} for k in dm}
        xall[t] = {}
        
        #---sample model from prior  
        for i in N:
            theta_hat[t][i] = np.random.beta(a[i],b[i])    
        
        #---select action
        for k in dm:
            xdm[k] = optimize(k,t,xub,theta_hat[t],N,dm,S,I,p,same,bg,ratio,P,B)
            for i in dm[k]:
                xall[t][i] = xdm[k][i]
                
                #---update upper bound on x to account for re-vaccination / population mixing
                xub[i] = 1.0 - sum(xall[tt][i] for tt in range(max(t-PM[i]+1,1),t+1))

        #---apply action
        S,I,R,D,BS,pdmI = simulation(t,xall[t],theta_env[t],N,S,I,R,D,node2dm,bg,lbd,p,ratio,dm)
        
        #---budget sharing
        B = share_budget(t,dm,B0,BS,pdmI,cap,D,mode)
        
        #---update distributions    
        for i in N:
            if xall[t][i] > 0.0:
                for n in range(1):
                    draw = np.random.binomial(1,theta_env[t][i])
                    if draw == 1:
                        a[i] += 1
                    else:
                        b[i] += 1

    rt = time.time() - t0
    print('>>>\t%.1f\t%.1f\n' % (rt,rt/T))
    return S,I,R,D,xall,T
    
def optimize_scenario_GY(filename,N,P,PM,node2dm,p,same,dm,popdm,S0,I0,bg,cap,cB,cPM,mode):
    t0 = time.time()
    print(filename)
    theta_true,theta_env,theta_hat,T,eps = read_scenario(filename,dm,N)
    S,I,R,D,a,b,xub,ratio,lbd = initialization(N,node2dm,popdm,S0,I0)
    B0 = set_budget(dm,P,cap,cB)
    PM = set_PM(PM,cPM)
    xall = {}
    B = B0

    for t in range(1,T+1):
        tstart = time.time()        
        xdm = {k:{} for k in dm}
        xall[t] = {}
        
        #---take mean of prior
        for i in N:
            theta_hat[t][i] = a[i]/(a[i]+b[i])   
        
        #---select action
        for k in dm:
            xdm[k] = optimize(k,t,xub,theta_hat[t],N,dm,S,I,p,same,bg,ratio,P,B)
            for i in dm[k]:
                xall[t][i] = xdm[k][i]
                xub[i] = 1.0 - sum(xall[tt][i] for tt in range(max(t-PM[i]+1,1),t+1))

        #---apply action
        S,I,R,D,BS,pdmI = simulation(t,xall[t],theta_env[t],N,S,I,R,D,node2dm,bg,lbd,p,ratio,dm)
        
        #---budget sharing
        B = share_budget(t,dm,B0,BS,pdmI,cap,D,mode)
        
        #---update distributions    
        for i in N:
            #---only update prior distributions of nodes which were allocated vaccines
            if xall[t][i] > 0.0:
                for n in range(1):
                    draw = np.random.binomial(1,theta_env[t][i])
                    if draw == 1:
                        a[i] += 1
                    else:
                        b[i] += 1

    rt = time.time() - t0
    print('>>>\t%.1f\t%.1f\n' % (rt,rt/T))
    return S,I,R,D,xall,T
    
def optimize_scenario_MA(filename,N,P,PM,node2dm,p,same,dm,popdm,S0,I0,bg,cap,cB,cPM,mode):
    t0 = time.time()
    print(filename)
    theta_true,theta_env,theta_hat,T,eps = read_scenario(filename,dm,N)
    S,I,R,D,obs,ma_size,xub,ratio,lbd = initialization2(N,node2dm,popdm,S0,I0)
    B0 = set_budget(dm,P,cap,cB)
    PM = set_PM(PM,cPM)    
    xall = {}
    B = B0

    for t in range(1,T+1):
        tstart = time.time()        
        xdm = {k:{} for k in dm}
        xall[t] = {}
        
        #---determine effiency rates by moving average 
        for i in N:
            theta_hat[t][i] = sum(obs[i])/len(obs[i])  
        
        #---select action
        for k in dm:
            xdm[k] = optimize(k,t,xub,theta_hat[t],N,dm,S,I,p,same,bg,ratio,P,B)
            for i in dm[k]:
                xall[t][i] = xdm[k][i]
                xub[i] = 1.0 - sum(xall[tt][i] for tt in range(max(t-PM[i]+1,1),t+1))

        #---apply action
        S,I,R,D,BS,pdmI = simulation(t,xall[t],theta_env[t],N,S,I,R,D,node2dm,bg,lbd,p,ratio,dm)
        
        #---budget sharing
        B = share_budget(t,dm,B0,BS,pdmI,cap,D,mode)
        
        #---update observations    
        for i in N:
            if xall[t][i] > 0.0: 
                obs[i].append(theta_env[t][i])
                if len(obs[i]) > ma_size:
                    obs[i].pop(0)

    rt = time.time() - t0
    print('>>>\t%.1f\t%.1f\n' % (rt,rt/T))
    return S,I,R,D,xall,T  
    
def optimize_scenario_PB(filename,N,P,PM,node2dm,p,same,dm,popdm,S0,I0,bg,cap,cB,cPM,mode):

    file = open('debug.txt','w')

    t0 = time.time()
    print(filename)
    theta_true,theta_env,theta_hat,T,eps = read_scenario(filename,dm,N)
    S,I,R,D,obs,ma_size,xub,ratio,lbd = initialization2(N,node2dm,popdm,S0,I0)
    B0 = set_budget(dm,P,cap,cB)
    PM = set_PM(PM,cPM)    
    xall = {}
    B = B0

    for k in dm:
        file.write('%s\t%.3f\n' % (k,B[k]))    

    for t in range(1,T+1):    
        tstart = time.time()        
        xdm = {k:{} for k in dm}
        xall[t] = {}
        
        #---select action
        for k in dm:
            xdm[k] = optimize_PB(k,t,xub,dm,P,B)
            for i in dm[k]:
                xall[t][i] = xdm[k][i]
                
                #---update upper bound on x to account for re-vaccination / population mixing
                xub[i] = 1.0 - sum(xall[tt][i] for tt in range(max(t-PM[i]+1,1),t+1))

        #---apply action
        S,I,R,D,BS,pdmI = simulation(t,xall[t],theta_env[t],N,S,I,R,D,node2dm,bg,lbd,p,ratio,dm)
        
        #---budget sharing
        B = share_budget(t,dm,B0,BS,pdmI,cap,D,mode)
        
        for k in dm:
            file.write('%s\t%.3f\t%.3f\n' % (k,B[k],BS[k]))       

    rt = time.time() - t0
    print('>>>\t%.1f\t%.1f\n' % (rt,rt/T))
    
    file.close()    
    
    return S,I,R,D,xall,T
	
