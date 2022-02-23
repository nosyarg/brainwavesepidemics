import networkx as nx
from random import *
from bisect import *
import matplotlib.pyplot as plt
def inlist(a, x):
    i = bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return 1
    return 0
r = 4
s = 5
n = s**2
mu = 8
u = 100000
p = mu/n
lambd = 2
rho = 1
def dist(n1,n2):
    return ((n1['loc'][0]-n2['loc'][0])**2 + (n1['loc'][1]-n2['loc'][1])**2)**(1/2)
def genlattice(s,mu,u):
    G = nx.Graph()
    for i in range(s):
        for j in range(s):
            G.add_node(i*s+j)
            G.nodes[i*s+j]['loc'] = (i,j)
            print(G.nodes[i*s+j]['loc'][0])
    priortotal = 0
    for n1 in G.nodes:
        for n2 in G.nodes:
            if(n1 > n2):
                priortotal += dist(G.nodes[n1],G.nodes[n2])**-u
    for n1 in G.nodes:
        for n2 in G.nodes:
            if(n1 > n2):
                edgeprob = mu*s**2*dist(G.nodes[n1],G.nodes[n2])**-u
                if random() < edgeprob:
                    G.add_edge(n1,n2)
    return G
while (r/n < .2):
        counter = 0
        infected = []
        susceptible = []
        si = []
        #G = nx.fast_gnp_random_graph(n,p)#set up the er graph and infect patient 0
        #G = genlattice(s,mu,u)
        G = nx.barabasi_albert_graph(n,mu)
        slist = []
        ilist = []
        rlist = []
        tlist = []
        t = 0
        r = 0
        i = 1
        s = n-1
        for j in range(len(G.nodes)):
                G.nodes[j]['state'] = 's'
        firstinfected = int(n*random())
        G.nodes[firstinfected]['state'] = 'i'
        infected = [firstinfected]
        susceptible = list(range(n))
        susceptible.remove(firstinfected)
        removed = []
        si = sorted(list(G.edges(0)))
        for edge in si:
            assert(G.nodes[edge[0]]['state'] == 's' or G.nodes[edge[1]]['state'] == 's')
        while(len(infected) > 0): #do the modified contact process algorithm for disease spread
                if(counter % 10 == 0):
                        mu = 0
                        for sus in susceptible:
                                mu += G.degree(sus)/len(susceptible)
                counter += 1
                if(random() < len(infected)/(len(infected)+(lambd)*len(si)+rho*r)):#node event
                        todie = choice(infected)
                        insort(removed,todie)
                        t += 1/(i+(lambd)*len(si) + rho*r)
                        i -= 1
                        r += 1
                        slist.append(s)
                        ilist.append(i)
                        rlist.append(r)
                        tlist.append(t)
                        #writefile.close()
                        G.nodes[todie]['state'] = 'r'
                        del infected[bisect_left(infected,todie)]
                        for node in G.neighbors(todie):
                                if(inlist(si,(todie,node))):
                                        del si[bisect_left(si,(todie,node))]
                                if(inlist(si,(node,todie))):
                                        del si[bisect_left(si,(node,todie))]
                        for edge in si:
                            assert(G.nodes[edge[0]]['state'] == 's' or G.nodes[edge[1]]['state'] == 's')
                elif(random()<lambd*len(si)/(rho*r+lambd*len(si))): #infection edge event 
                        t += 1/(i+(lambd)*len(si)+rho*r)
                        i += 1
                        s -= 1
                        slist.append(s)
                        ilist.append(i)
                        rlist.append(r)
                        tlist.append(t)
                        cross = choice(si)
                        if(G.nodes[cross[0]]['state']=='s'):
                                G.nodes[cross[0]]['state'] = 'i'
                                insort(infected,cross[0])
                                del susceptible[bisect_left(susceptible,cross[0])]
                                for node in G.neighbors(cross[0]):#change si list
                                        if(G.nodes[node]['state']=='i'):
                                                if(inlist(si,(cross[0],node))):
                                                        del si[bisect_left(si,(cross[0],node))]
                                                if(inlist(si,(node,cross[0]))):
                                                        del si[bisect_left(si,(node,cross[0]))]
                                        if(G.nodes[node]['state']=='s'):
                                                insort(si,(cross[0],node))
                        elif(G.nodes[cross[1]]['state']=='s'):
                                G.nodes[cross[1]]['state'] = 'i'
                                insort(infected,(cross[1]))
                                del susceptible[bisect_left(susceptible,cross[1])]
                                for node in G.neighbors(cross[1]):#change si list
                                        if(G.nodes[node]['state']=='i'):
                                                if(inlist(si,(cross[1],node))):
                                                        del si[bisect_left(si,(cross[1],node))]
                                                if(inlist(si,(node,cross[1]))):
                                                        del si[bisect_left(si,(node,cross[1]))]
                                        if(G.nodes[node]['state']=='s'):
                                                insort(si,(cross[1],node))
                        else:
                            print(G.nodes[cross[0]]['state'])
                            print(G.nodes[cross[1]]['state'])
                            for edge in si:
                                assert(G.nodes[edge[0]]['state'] == 's' or G.nodes[edge[1]]['state'] == 's')
                            assert(1==0)
                        for edge in si:
                            assert(G.nodes[edge[0]]['state'] == 's' or G.nodes[edge[1]]['state'] == 's')
                else: #RS event
                        assert(rho > 0)
                        todie = choice(removed)
                        t += 1/(i+(lambd)*len(si)+rho*r)
                        s += 1
                        r -= 1
                        slist.append(s)
                        ilist.append(i)
                        rlist.append(r)
                        tlist.append(t)
                        G.nodes[todie]['state'] = 's'
                        del removed[bisect_left(removed,todie)]
                        for node in G.neighbors(todie):
                                if((not inlist(si,(todie,node))) and (not inlist(si,(node,todie)))and (G.nodes[node]['state'] == 's')):
                                        insort(si,(node,todie))
                        for edge in si:
                            assert(G.nodes[edge[0]]['state'] == 's' or G.nodes[edge[1]]['state'] == 's')
                print(s)
                print(len(susceptible))
                assert(len(removed) == r)
                assert(len(susceptible) == s)
                for edge in si:
                    assert(G.nodes[edge[0]]['state'] == 's' or G.nodes[edge[1]]['state'] == 's')
        if(r/n %100 == 0): 
            plt.plot(tlist,slist)
            plt.plot(tlist,ilist)
            plt.plot(tlist,rlist)
            plt.show(block = False)
        survivors = 0
        print(r/n)

