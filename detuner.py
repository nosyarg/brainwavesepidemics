import matplotlib.pyplot as plt
from random import *
from scipy.optimize import minimize
timestep = .001
n = 10000
t = 0
lambd = 3/n
rho = .1
s = n-1
i = 1
r = 0
mode = 'lambda-'
def getdampingrat(lambdrho):
    dampingratiomax = 100
    maxtime = 1000
    if(min(lambdrho) < 0):
        return dampingratiomax
    lambd = lambdrho[0]
    rho = lambdrho[1]
    slist = []
    ilist = []
    rlist = []
    tlist = []
    timestep = .001
    n = 10000
    t = 0
    s = n-1
    i = 1
    r = 0
    scopy = s
    icopy = i
    rcopy = r
    while(t < maxtime):
        t += timestep
        s += (rho*r-lambd*s*i)*timestep
        i += (lambd*s*i-i)*timestep
        r += (i-rho*r)*timestep
        slist.append(s)
        ilist.append(i)
        rlist.append(r)
        tlist.append(t)

    plt.cla()
    plt.plot(tlist,slist, 'r')
    plt.plot(tlist,ilist, 'g')
    plt.plot(tlist,rlist, 'b')
    plt.savefig('deoutput.png')
    maxima = []
    minima = []
    for i in range(1,len(ilist)-1):
        if(ilist[i-1] < ilist[i] and ilist[i+1] < ilist[i]): 
            maxima.append(i)
        if(ilist[i-1] > ilist[i] and ilist[i+1] > ilist[i]): 
            minima.append(i)
    print('maxima')
    for i in maxima:
        print(ilist[i])
    print('minima')
    for i in minima:
        print(ilist[i])
    try:
        stablevel = (ilist[minima[-1]] + ilist[maxima[-1]])/2

        dampingratiomax = (ilist[maxima[1]] - stablevel)/(ilist[maxima[2]] - stablevel)
        dampingratiomin = (-ilist[minima[1]] + stablevel)/(-ilist[minima[2]] + stablevel)
    except:
        return dampingratiomax
    return dampingratiomax
    '''print("maxrat:" + str(dampingratiomax))
    print("minrat:" + str(dampingratiomin))
    print("lambda:" + str(lambd*n))
    print("rho:" + str(rho))
    print("mode:" + mode)
    if(mode == 'lambda-'):
        lambd /= dampingratiomax**(1/1000)
        if(dratold < dampingratiomax):
            mode = 'lambda+'
    elif(mode == 'lambda+'):
        lambd *= dampingratiomax**(1/1000)
        if(dratold < dampingratiomax):
            mode = 'rho+'
    elif(mode == 'rho+'):
        rho *= dampingratiomax**(1/100000)
        if(dratold < dampingratiomax):
            mode = 'rho-'
    elif(mode == 'rho-'):
        rho /= dampingratiomax**(1/100000)
        if(dratold < dampingratiomax):
            mode = 'lambda-'
    dratold = dampingratiomax'''
'''con = lambda x: min(x[0],x[1])
#nlc = NonlinearConstraint(con,0,2)
cons = ({'type': 'ineq', 'fun':con})
bnds = ((0,None))
result = minimize(getdampingrat,(3/n,.1),method='SLSQP',constraints=cons,bounds=bnds)'''
result = minimize(getdampingrat,(3/n,.1))

print(result.x)
print(result.fun)

'''
slist = []
ilist = []
rlist = []
tlist = []
t = 0
s = scopy
i = icopy
r = rcopy
posterior = 0
while(posterior == 0):
    while(t < 100):
        posterior = rho*r + lambd*s*i + i
        if(posterior ==0):
            print("early break")
            break
        t += 1/posterior
        rand = random()*posterior
        if(rand < rho*r):
            r -= 1
            s += 1
        elif(rand < rho*r+s*i*lambd):
            s -= 1
            i += 1
        elif(rand < rho*r + s*i*lambd + i):
            i -= 1
            r += 1
        else:
            print(s)
            print(i)
            print(r)
            print(rand)
            print(posterior)
            assert(1==0)
        slist.append(s)
        ilist.append(i)
        rlist.append(r)
        tlist.append(t)

plt.cla()
plt.plot(tlist,slist, 'r')
plt.plot(tlist,ilist, 'g')
plt.plot(tlist,rlist, 'b')
plt.savefig('stochdeoutput.png')
'''
