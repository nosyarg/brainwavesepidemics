import networkx as nx
from random import *
from bisect import *
def inlist(a, x):
    i = bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return 1
    return 0
r = 4
n = 1000
filepath = 'data/norewiresim.csv'
mu = 5
p = mu/n
lambd = 2
rho = 0
while (r/n < .2):
        counter = 0
        infected = []
        susceptible = []
        unexplored = []
        si = []
        G = nx.fast_gnp_random_graph(n,p)#set up the er graph and infect patient 0
        writefile = open(filepath,'w')
        writefile.write('n,p,lambda,rho,s,i,r,t,mu\n')
        writefile.close()
        t = 0
        r = 0
        i = 1
        s = n-1
        for j in range(len(G.nodes)):
                G.nodes[j]['state'] = 's'
        G.nodes[0]['state'] = 'i'
        infected = [0]
        susceptible = list(range(1,n))
        """
        for i in susceptible:
                unex = 1
                for j in G.neighbors(i):
                        if(inlist(infected,j)):
                                unex = 0
                if(unex):
                        unexplored.append(i)
        """
        si = sorted(list(G.edges(0)))
        while(len(infected) > 0): #do the modified contact process algorithm for disease spread
                if(counter % 10 == 0):
                        mu = 0
                        #print(len(susceptible))
                        for sus in susceptible:
                                mu += G.degree(sus)/len(susceptible)
                                """
                        for inf in infected:
                                mu += G.degree(inf)/len(infected)
                        """
                counter += 1
                if(random() < len(infected)/(len(infected)+(rho+lambd)*len(si))):#node event
                        todie = choice(infected)
                        t += 1/(i+(lambd+rho)*len(si))
                        i -= 1
                        r += 1
                        writefile = open(filepath,'a')
                        writefile.write(str(n) + ',' + str(p) + ',' + str(lambd) + ',' + str(rho) + ',' + str(s) + ',' + str(i) + ',' + str(r) + ',' + str(t) + ',' + str(mu) + '\n')
                        writefile.close()
                        G.nodes[todie]['state'] = 'r'
                        del infected[bisect_left(infected,todie)]
                        for node in G.neighbors(todie):
                                if(inlist(si,(todie,node))):
                                        del si[bisect_left(si,(todie,node))]
                                if(inlist(si,(node,todie))):
                                        del si[bisect_left(si,(node,todie))]
                elif(random()<lambd/(rho+lambd)): #infection edge event 
                        t += 1/(i+(lambd+rho)*len(si))
                        i += 1
                        s -= 1
                        writefile = open(filepath,'a')
                        writefile.write(str(n) + ',' + str(p) + ',' + str(lambd) + ',' + str(rho) + ',' + str(s) + ',' + str(i) + ',' + str(r) + ',' + str(t) + ',' + str(mu) + '\n')
                        writefile.close()
                        cross = choice(si)
                        if(G.nodes[cross[0]]['state']=='s'):
                                G.nodes[cross[0]]['state'] = 'i'
                                insort(infected,cross[0])
                                #print("WILL DELETE:")
                                #print(susceptible[bisect_left(susceptible,cross[0])])
                                del susceptible[bisect_left(susceptible,cross[0])]
                                for node in G.neighbors(cross[0]):#change si list
                                        if(G.nodes[node]['state']=='i'):
                                                if(inlist(si,(cross[0],node))):
                                                        del si[bisect_left(si,(cross[0],node))]
                                                if(inlist(si,(node,cross[0]))):
                                                        del si[bisect_left(si,(node,cross[0]))]
                                        if(G.nodes[node]['state']=='s'):
                                                insort(si,(cross[0],node))
                                                if(inlist(unexplored,node)):
                                                        del unexplored[bisect_left(unexplored,node)]
                        elif(G.nodes[cross[1]]['state']=='s'):
                                G.nodes[cross[1]]['state'] = 'i'
                                insort(infected,(cross[1]))
                                #print("WILL DELETE:")
                                #print(susceptible[bisect_left(susceptible,cross[1])])
                                del susceptible[bisect_left(susceptible,cross[1])]
                                for node in G.neighbors(cross[1]):#change si list
                                        if(G.nodes[node]['state']=='i'):
                                                if(inlist(si,(cross[1],node))):
                                                        del si[bisect_left(si,(cross[1],node))]
                                                if(inlist(si,(node,cross[1]))):
                                                        del si[bisect_left(si,(node,cross[1]))]
                                        if(G.nodes[node]['state']=='s'):
                                                insort(si,(cross[1],node))
                                                if(inlist(unexplored,node)):
                                                        del unexplored[bisect_left(unexplored,node)]
                        else:
                                 assert(1==0)
                else: #rewire edge event
                        t += 1/(i+(lambd+rho)*len(si))
                        writefile = open(filepath,'a')
                        writefile.write(str(n) + ',' + str(p) + ',' + str(lambd) + ',' + str(rho) + ',' + str(s) + ',' + str(i) + ',' + str(r) + ',' + str(t) + ',' + str(mu) + '\n')
                        writefile.close()
                        rewire = choice(si)
                        del si[bisect_left(si,rewire)]
                        G.remove_edge(rewire[0],rewire[1])
                        if(G.nodes[rewire[0]]['state'] == 's'):
                                svert = rewire[0]
                        else:
                                svert = rewire[1]
                        newvert = svert
                        while((newvert==svert) or (newvert in G.neighbors(svert))):#pick a valid vertex
                                newvert = int(n*random())
                        G.add_edge(newvert,svert)
                        if(G.nodes[newvert]['state'] == 'i'):
                                insort(si,(newvert,svert))
        survivors = 0
        print(r/n)
