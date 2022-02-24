import matplotlib.pyplot as plt
from random import *
slist = []
ilist = []
rlist = []
tlist = []
timestep = .001
n = 10000
t = 0
lambd = 3/n
rho = .1
s = n-1
i = 1
r = 0
scopy = s
icopy = i
rcopy = r
while(t < 100):
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

slist = []
ilist = []
rlist = []
tlist = []
t = 0
s = scopy
i = icopy
r = rcopy

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
