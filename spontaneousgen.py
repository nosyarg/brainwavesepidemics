import networkx as nx
from random import *
from bisect import *
import matplotlib.pyplot as plt
debug = 0
def inlist(a, x):
    i = bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return 1
    return 0
def checkinvariants(G,si,infected,removed,susceptible):
    if(debug):
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
n = 100000
mu = 2
u = 100000
maxtime = 400
p = mu/n
lambd = 200#infection rate
rho = .01#resusceptible rate
endo = .000001#endogenous infection rate
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
t = 0
while (t < maxtime):
        counter = 0
        infected = []
        susceptible = []
        removed = []
        si = []
        #G = nx.fast_gnp_random_graph(n,p)#set up the er graph and infect patient 0
        #G = genlattice(s,mu,u)
        G = nx.barabasi_albert_graph(n,mu)
        #G = nx.complete_graph(n)
        #plt.cla()
        #nx.draw(G)
        #plt.draw()
        #plt.savefig('graph.png')
        slist = []
        ilist = []
        rlist = []
        tlist = []
        t = 0
        for j in range(len(G.nodes)):
                G.nodes[j]['state'] = 's'
        firstinfected = int(n*random())
        G.nodes[firstinfected]['state'] = 'i'
        infected = [firstinfected]
        #heappush(infected,firstinfected)
        susceptible = sorted(list(range(n)))
        #heappop(susceptible,firstinfected)
        del susceptible[bisect_left(susceptible,firstinfected)]
        si = sorted(list(G.edges(firstinfected)))
        while(1): #do the modified contact process algorithm for disease spread
                checkinvariants(G,si,infected,removed,susceptible)
                counter += 1
                timestepsize = 1/(endo*len(susceptible) + len(infected)+(lambd)*len(si) + rho*len(removed))
                if(random() < len(infected)/(endo*len(susceptible) + len(infected)+(lambd)*len(si)+rho*len(removed))):#node event
                        todie = choice(infected)
                        insort(removed,todie)
                        t += timestepsize
                        slist.append(len(susceptible))
                        ilist.append(len(infected))
                        rlist.append(len(removed))
                        tlist.append(t)
                        G.nodes[todie]['state'] = 'r'
                        del infected[bisect_left(infected,todie)]
                        for node in G.neighbors(todie):
                                if(inlist(si,(todie,node))):
                                        del si[bisect_left(si,(todie,node))]
                                if(inlist(si,(node,todie))):
                                        del si[bisect_left(si,(node,todie))]
                        if(debug):
                            for edge in si:
                                assert(G.nodes[edge[0]]['state'] == 's' or G.nodes[edge[1]]['state'] == 's')
                elif(random()<lambd*len(si)/(endo*len(susceptible) + rho*len(removed)+lambd*len(si))): #infection edge event 
                        t += timestepsize
                        slist.append(len(susceptible))
                        ilist.append(len(infected))
                        rlist.append(len(removed))
                        tlist.append(t)
                        cross = choice(si)
                        checkinvariants(G,si,infected,removed,susceptible)
                        if(G.nodes[cross[0]]['state']=='s'):
                                checkinvariants(G,si,infected,removed,susceptible)
                                G.nodes[cross[0]]['state'] = 'i'
                                insort(infected,cross[0])
                                del susceptible[bisect_left(susceptible,cross[0])]
                                for node in G.neighbors(cross[0]):#change si list
                                        if(G.nodes[node]['state']=='i'):
                                                if(si.count((cross[0],node)) > 0):
                                                    del si[bisect_left(si,(cross[0],node))]
                                                if(si.count((node,cross[0])) > 0):
                                                    del si[bisect_left(si,(node,cross[0]))]
                                        if(G.nodes[node]['state']=='s'):
                                                insort(si,(cross[0],node))
                                checkinvariants(G,si,infected,removed,susceptible)
                        elif(G.nodes[cross[1]]['state']=='s'):
                                checkinvariants(G,si,infected,removed,susceptible)
                                G.nodes[cross[1]]['state'] = 'i'
                                insort(infected,cross[1])
                                del susceptible[bisect_left(susceptible,cross[1])]
                                for node in G.neighbors(cross[1]):#change si list
                                        if(G.nodes[node]['state']=='i'):
                                            if(inlist(si,(cross[1],node))):
                                                del si[bisect_left(si,(cross[1],node))]
                                            if(inlist(si,(node,cross[1]))):
                                                del si[bisect_left(si,(node,cross[1]))]
                                        if(G.nodes[node]['state']=='s'):
                                                insort(si,(cross[1],node))
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
                elif(random()<endo*len(susceptible)/(endo*len(susceptible) + rho*len(removed))): #Endogenous event
                        assert(endo > 0)
                        newinfected = choice(susceptible)
                        t += timestepsize
                        slist.append(len(susceptible))
                        ilist.append(len(infected))
                        rlist.append(len(removed))
                        tlist.append(t)
                        G.nodes[newinfected]['state'] = 'i'
                        del susceptible[bisect_left(susceptible,newinfected)]
                        insort(infected,newinfected)
                        for node in G.neighbors(newinfected):
                                if((not inlist(si, (newinfected,node))) and (not inlist(si,(node,newinfected))) and (G.nodes[node]['state'] == 's')):
                                        insort(si,(node,newinfected))
                                for node in G.neighbors(newinfected):#change si list
                                        if(G.nodes[node]['state']=='i'):
                                            if(inlist(si,(newinfected,node))):
                                                del si[bisect_left(si,(newinfected,node))]
                                            if(inlist(si,(node,newinfected))):
                                                del si[bisect_left(si,(node,newinfected))]
                        checkinvariants(G,si,infected,removed,susceptible)

                else: #RS event
                        assert(rho > 0)
                        todie = choice(removed)
                        t += timestepsize
                        slist.append(len(susceptible))
                        ilist.append(len(infected))
                        rlist.append(len(removed))
                        tlist.append(t)
                        G.nodes[todie]['state'] = 's'
                        del removed[bisect_left(removed,todie)]
                        insort(susceptible,todie)
                        for node in G.neighbors(todie):
                                if((not inlist(si, (todie,node))) and (not inlist(si,(node,todie))) and (G.nodes[node]['state'] == 'i')):
                                        insort(si,(node,todie))
                        checkinvariants(G,si,infected,removed,susceptible)
                        #for edge in si:
                        #assert(G.nodes[edge[0]]['state'] == 's' or G.nodes[edge[1]]['state'] == 's')
                #print(s)
                #print(len(susceptible))
                checkinvariants(G,si,infected,removed,susceptible)
                #for edge in si:
                #    assert(G.nodes[edge[0]]['state'] == 's' or G.nodes[edge[1]]['state'] == 's')
                if(t > maxtime): 
                    '''plt.plot(tlist,slist)
                    plt.plot(tlist,ilist)
                    plt.plot(tlist,rlist)
                    plt.savefig('output.png')'''
                    break
                    #plt.show(block = False)
        
        plt.cla()
        #plt.xlim(tlist[-1]/3,2*tlist[-1]/3)
        plt.plot(tlist,slist, 'r')
        plt.plot(tlist,ilist, 'g')
        plt.plot(tlist,rlist, 'b')
        plt.savefig('output.png')
        plt.show(block = False)
        survivors = 0
        print(len(removed)/n)
        maxima = []


