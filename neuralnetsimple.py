import networkx as nx
from random import *
from bisect import *
from heapq import *
import matplotlib.pyplot as plt
'''def inlist(a, x):
    i = bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return 1
    return 0'''
#r = 4
def checkinvariants(G,si,infected,removed,susceptible):
    return 0
    for edge in si:
        assert(G.nodes[edge[0]]['state'] == 's' or G.nodes[edge[1]]['state'] == 's')
    for edge in si:
        assert(G.nodes[edge[0]]['state'] == 's' or G.nodes[edge[1]]['state'] == 's')
    for node in infected:
        assert(G.nodes[node]['state'] == 'i')
    for node in susceptible:
        assert(G.nodes[node]['state'] == 's')
    for node in removed:
        assert(G.nodes[node]['state'] == 'r')
    
removed = []
#s = 15
#n = s**2
n = 100
mu = 1
u = 100000
p = mu/n
lambd = 3/n
rho = .3
def dist(n1,n2):
    return ((n1['loc'][0]-n2['loc'][0])**2 + (n1['loc'][1]-n2['loc'][1])**2)**(1/2)
def genlattice(s,mu,u):
    G = nx.Graph()
    for i in range(s):
        for j in range(s):
            G.add_node(i*s+j)
            G.nodes[i*s+j]['loc'] = (i,j)
            #print(G.nodes[i*s+j]['loc'][0])
    priortotal = 0
    for n1 in G.nodes:
        for n2 in G.nodes:
            if(n1 > n2):
                priortotal += dist(G.nodes[n1],G.nodes[n2])**-u
    for n1 in G.nodes:
        for n2 in G.nodes:
            if(n1 > n2):
                edgeprob = mu*s**2*dist(G.nodes[n1],G.nodes[n2])**-u/priortotal
                if random() < edgeprob:
                    G.add_edge(n1,n2)
    return G

while (len(removed)/n < .2):
        counter = 0
        infected = []
        susceptible = []
        removed = []
        si = []
        #G = nx.fast_gnp_random_graph(n,p)#set up the er graph and infect patient 0
        #G = genlattice(s,mu,u)
        #G = nx.barabasi_albert_graph(n,mu)
        G = nx.complete_graph(n)
        plt.cla()
        nx.draw(G)
        plt.draw()
        plt.savefig('graph.png')
        slist = []
        ilist = []
        rlist = []
        tlist = []
        t = 0
        for j in range(len(G.nodes)):
                G.nodes[j]['state'] = 's'
        firstinfected = int(n*random())
        G.nodes[firstinfected]['state'] = 'i'
        infected = []
        heappush(infected,firstinfected)
        susceptiblelist = list(range(n))
        susceptible = []
        for num in susceptiblelist:
            heappush(susceptible,num)
        heappop(susceptible,firstinfected)
        silist = sorted(list(G.edges(firstinfected)))
        si = []
        for num in silist:
            heappush(si,num)
        while(len(infected) > 0): #do the modified contact process algorithm for disease spread
                checkinvariants(G,si,infected,removed,susceptible)
                if(counter % 10 == 0):
                        mu = 0
                        for sus in susceptible:
                                mu += G.degree(sus)/len(susceptible)
                counter += 1
                if(random() < len(infected)/(len(infected)+(lambd)*len(si)+rho*len(removed))):#node event
                        todie = choice(infected)
                        removed.append(todie)
                        t += 1/(len(infected)+(lambd)*len(si) + rho*len(removed))
                        slist.append(len(susceptible))
                        ilist.append(len(infected))
                        rlist.append(len(removed))
                        tlist.append(t)
                        G.nodes[todie]['state'] = 'r'
                        infected.remove(todie)
                        for node in G.neighbors(todie):
                                if(si.count((todie,node)) > 0):
                                        si.remove((todie,node))
                                if(si.count((node,todie)) > 0):
                                        si.remove((node,todie))
                        for edge in si:
                            assert(G.nodes[edge[0]]['state'] == 's' or G.nodes[edge[1]]['state'] == 's')
                elif(random()<lambd*len(si)/(rho*len(removed)+lambd*len(si))): #infection edge event 
                        t += 1/(len(infected)+(lambd)*len(si)+rho*len(removed))
                        slist.append(len(susceptible))
                        ilist.append(len(infected))
                        rlist.append(len(removed))
                        tlist.append(t)
                        cross = choice(si)
                        checkinvariants(G,si,infected,removed,susceptible)
                        if(G.nodes[cross[0]]['state']=='s'):
                                checkinvariants(G,si,infected,removed,susceptible)
                                G.nodes[cross[0]]['state'] = 'i'
                                infected.append(cross[0])
                                susceptible.remove(cross[0])
                                for node in G.neighbors(cross[0]):#change si list
                                        if(G.nodes[node]['state']=='i'):
                                                if(si.count((cross[0],node)) > 0):
                                                        si.remove((cross[0],node))
                                                if(si.count((node,cross[0])) > 0):
                                                        si.remove((node,cross[0]))
                                        if(G.nodes[node]['state']=='s'):
                                                si.append((cross[0],node))
                                checkinvariants(G,si,infected,removed,susceptible)
                        elif(G.nodes[cross[1]]['state']=='s'):
                                checkinvariants(G,si,infected,removed,susceptible)
                                G.nodes[cross[1]]['state'] = 'i'
                                infected.append(cross[1])
                                susceptible.remove(cross[1])
                                for node in G.neighbors(cross[1]):#change si list
                                        if(G.nodes[node]['state']=='i'):
                                                if(si.count((cross[1],node)) > 0):
                                                        si.remove((cross[1],node))
                                                if(si.count((node,cross[1])) > 0):
                                                        si.remove((node,cross[1]))
                                        if(G.nodes[node]['state']=='s'):
                                                si.append((cross[1],node))
                                checkinvariants(G,si,infected,removed,susceptible)
                        else:
                            #print(G.nodes[cross[0]]['state'])
                            #print(G.nodes[cross[1]]['state'])
                            checkinvariants(G,si,infected,removed,susceptible)
                            #for edge in si:
                            #    print(edge)
                            #    assert(G.nodes[edge[0]]['state'] == 's' or G.nodes[edge[1]]['state'] == 's')
                            assert(1==0)
                        #for edge in si:
                        #    assert(G.nodes[edge[0]]['state'] == 's' or G.nodes[edge[1]]['state'] == 's')
                        checkinvariants(G,si,infected,removed,susceptible)
                else: #RS event
                        assert(rho > 0)
                        todie = choice(removed)
                        t += 1/(len(infected)+(lambd)*len(si)+rho*len(removed))
                        slist.append(len(susceptible))
                        ilist.append(len(infected))
                        rlist.append(len(removed))
                        tlist.append(t)
                        G.nodes[todie]['state'] = 's'
                        removed.remove(todie)
                        susceptible.append(todie)
                        for node in G.neighbors(todie):
                                if((si.count((todie,node)) == 0) and (si.count((node,todie)) == 0) and (G.nodes[node]['state'] == 'i')):
                                        si.append((node,todie))
                        checkinvariants(G,si,infected,removed,susceptible)
                        #for edge in si:
                        #assert(G.nodes[edge[0]]['state'] == 's' or G.nodes[edge[1]]['state'] == 's')
                #print(s)
                #print(len(susceptible))
                checkinvariants(G,si,infected,removed,susceptible)
                #for edge in si:
                #    assert(G.nodes[edge[0]]['state'] == 's' or G.nodes[edge[1]]['state'] == 's')
                if(counter %400000 == 0): 
                    '''plt.plot(tlist,slist)
                    plt.plot(tlist,ilist)
                    plt.plot(tlist,rlist)
                    plt.savefig('output.png')'''
                    break
                    #plt.show(block = False)
        
        plt.cla()
        plt.xlim(tlist[-1]/3,2*tlist[-1]/3)
        plt.plot(tlist,slist, 'r')
        plt.plot(tlist,ilist, 'g')
        plt.plot(tlist,rlist, 'b')
        plt.savefig('output.png')
        plt.show(block = False)
        survivors = 0
        print(len(removed)/n)

