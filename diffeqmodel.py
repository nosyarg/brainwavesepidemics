import matplotlib.pyplot as plt
slist = []
ilist = []
rlist = []
tlist = []
timestep = .1
t = 0
lambd = 1
rho = .05
n = 1000
s = (n-1)/n
i = 1
r = 0
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
#plt.xlim(tlist[-1]/3,2*tlist[-1]/3)
plt.plot(tlist,slist, 'r')
plt.plot(tlist,ilist, 'g')
plt.plot(tlist,rlist, 'b')
plt.savefig('deoutput.png')
