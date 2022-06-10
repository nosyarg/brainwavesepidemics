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
def runsim(
            genfunc = nx.barabasi_albert_graph,
            graphparams = [10000,2],
            lambd = 1,
            rho = 1,
            tmax = 200,
            maxnumtries = 10
            ):
    numtries = 0
    while ((len(removed)/n < .2 or t < maxtime) and numtries < maxnumtries):
            numtries += 1
            counter = 0
            infected = []
            susceptible = []
            removed = []
            si = []
            G = genfunc(*graphparams)
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
            susceptible = sorted(list(range(n)))
            del susceptible[bisect_left(susceptible,firstinfected)]
            si = sorted(list(G.edges(firstinfected)))
            while(len(infected) > 0): #do the modified contact process algorithm for disease spread
                    checkinvariants(G,si,infected,removed,susceptible)
                    counter += 1
                    if(random() < len(infected)/(len(infected)+(lambd)*len(si)+rho*len(removed))):#node event
                            todie = choice(infected)
                            insort(removed,todie)
                            t += 1/(len(infected)+(lambd)*len(si) + rho*len(removed))
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
                                checkinvariants(G,si,infected,removed,susceptible)
                                assert(1==0)
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
                            del removed[bisect_left(removed,todie)]
                            insort(susceptible,todie)
                            for node in G.neighbors(todie):
                                    if((not inlist(si, (todie,node))) and (not inlist(si,(node,todie))) and (G.nodes[node]['state'] == 'i')):
                                            insort(si,(node,todie))
                            checkinvariants(G,si,infected,removed,susceptible)
                    checkinvariants(G,si,infected,removed,susceptible)
                    if(t > maxtime): 
                        break
            if (numtries == maxnumtries):
                return 0
            else:
                return {'s':slist,'i':ilist,'r':rlist,'t':tlist}
                


