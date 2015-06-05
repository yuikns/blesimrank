
import sys
import re
import networkx as nx
import time
import random
from collections import Counter
from collections import defaultdict
import math
import heapq



class ScalableSimRank:
    def __init__(self, dataset, T=11, P=11, Q=5, R1=100, R2=10000, K=20, theta=0.01):
        self.c = 0.6
        self.dataset = dataset
        self.networkFile = "data/"+self.dataset+".graph"
        self.outputFile = "result/"+self.dataset+".ssr"
        self.G = nx.Graph()
        self.newG = nx.Graph()
        self.T = T
        self.P = P
        self.Q = Q
        self.R1 = R1
        self.R2 = R2
        self.K = K
        self.theta = theta

        self.N = 0
        self.M = 0

        print "dataset=",self.dataset,", networkfile = ", self.networkFile,"outputfile = ", self.outputFile, "T = ", self.T, "P = ", self.P, "Q = ", self.Q

    def loadData(self):
        i = 0
        for line in open(self.networkFile, "r"):
            if i == 0:
                elements = re.compile('\t').split(line.strip())
                self.N = int(elements[0])
                self.M = int(elements[1])
                print "N=", self.N, ", M=", self.M
                for n in range(self.N):
                    self.G.add_node(n, weight=0)
                i+=1
                continue
            elements = re.compile('\t').split(line.strip())
            node1 = int(elements[0])
            node2 = int(elements[1])
            w = float(elements[2])
            self.G.add_edge(node1, node2, weight=w)
            self.G.node[node1]['weight'] += w
            self.G.add_edge(node2, node1, weight=w)
            self.G.node[node2]['weight'] += w

            self.newG = nx.Graph()

            i+=1
            if i % 10000 == 0:
                print i, "th edge"

        print "load edges done."
        print "#nodes = ", self.G.number_of_nodes(), ", #edges = ", self.G.number_of_edges()





    def calculateTransitionPros(self):

        self.transitionpros = [ [] for i in range(self.N)]

        for v in self.G.nodes():
            neighbors = self.G.neighbors(v)
            if len(neighbors) == 0:
                continue

            totalWeight = self.G.node[v]['weight']

            # print "totalWeight = ", totalWeight
            # for j in range(len(neighbors)):
            # 	print "to neighbor ", j, ", id = ", neighbors[j], ", weight = ", self.G.edge[v][neighbors[j]]['weight']

            self.transitionpros[v].append(self.G.edge[v][neighbors[0]]['weight']/totalWeight)
            for j in range(len(neighbors)-1):
                transitionProbability = self.G.edge[v][neighbors[j+1]]['weight']/totalWeight
                self.transitionpros[v].append( transitionProbability + self.transitionpros[v][j])

            # for j in range(len(neighbors)):
            # 	print "to neighbor ", j, ", id = ", neighbors[j], ", weight = ", self.transitionpros[v][j]


    def randomWalk(self,v):
        path = []
        currentNode = v
        t = 0
        while t < self.T:
            neighbors = self.G.neighbors(currentNode)
            if len(neighbors) == 0:
                break


            rand = random.random()
            for j in range(len(neighbors)):
                if self.transitionpros[v][j] >= rand:
                    break
            neighborIndex = j

            currentNode = neighbors[neighborIndex]
            path.append(currentNode)
            t+=1
        return path




    def generateBipartiteGraph(self):

        # The original nodes as the left nodes
        for n in range(self.N):
            self.newG.add_node(n)

        # Copy the original ndoes as the right nodes
        for n in range(self.N):
            self.newG.add_node(n+self.N)


        i = 0
        for v in self.G.nodes():
            neighbors = self.G.neighbors(v)
            if len(neighbors) == 0:
                continue
            # for j in range(len(neighbors)):
            # 	print "to neighbor ", j, ", id = ", neighbors[j], ", weight = ", self.G.edge[v][neighbors[j]]['weight'], ", accumulated weight =", self.transitionpros[v][j]
            for i in range(self.P):
                paths = []
                for q in range(self.Q):
                    paths.append(self.randomWalk(v))
                # print paths
                for t in range(self.T):
                    l = []
                    for q in range(self.Q):
                        l.append(paths[q][t])
                    # print "t=",t," l=",l
                    nodesMoreThanOneTime = set([x for x in l if l.count(x) > 1])
                    # print "common node = ", nodesMoreThanOneTime
                    for n in nodesMoreThanOneTime:
                        if n != v:
                            n += self.N
                            self.newG.add_edge(v, n)
            if i % 10000 == 0:
                print i, "th node"

    def calculateGammaBound(self):

        self.gamma = [[] for i in range(self.N)]
        for v in self.G.nodes():
            neighbors = self.G.neighbors(v)
            # print neighbors
            if len(neighbors) == 0:
                for t in range(self.T):
                    self.gamma[v].append(0)
            else:
                paths = []
                for r in range(self.R1):
                    paths.append(self.randomWalk(v))
                # print paths
                for t in range(self.T):
                    cnt = Counter()
                    for r in range(self.R1):
                        cnt[paths[r][t]]+=1
                    # print "t=",t,", cnt=",cnt
                    mu = 0
                    for path in cnt:
                        pathNumber = cnt[path]
                        mu +=(1-self.c)*(pathNumber * pathNumber) / (self.R1*self.R1)
                    self.gamma[v].append(math.sqrt(mu))

        # print self.gamma

    def calculateBetaBound(self):
        self.beta = [[] for i in range(self.N)]

        for v in self.G.nodes():
            neighbors = self.G.neighbors(v)
            if len(neighbors) == 0:
                for t in range(self.T):
                    self.beta[v].append(0)
            else:
                alpha = [0 for t in range(self.T)]
                paths = []
                for r in range(self.R2):
                    paths.append(self.randomWalk(v))

                for t in range(self.T):
                    cnt = Counter()
                    for r in range(self.R2):
                        cnt[paths[r][t]]+=1
                    for path in cnt:
                        pathNumber = cnt[path]
                        mu =(1-self.c)*pathNumber / self.R2
                        alpha[t]  = max(alpha[t], mu)
                # print alpha
                for d in range(self.T):
                    sum = 0
                    for t in range(d):
                        minDistance = d-t
                        maxDistance = d+t
                        maxAlpha = 0
                        for d_ in range(d-t, min(d+t+1, self.T)):
                            # print d, t, d_
                            maxAlpha = max(maxAlpha, alpha[d_])
                        sum += math.pow(self.c,t) * maxAlpha
                    self.beta[v].append(sum)

        # print self.beta



    def preprocess(self):
        self.loadData()
        print "Load data done."
        self.calculateTransitionPros()
        print "Calculate transition probabilities done."
        self.generateBipartiteGraph()
        print "Generate bipartite graph done."
        self.calculateGammaBound()
        print "Calcualte gamma bound done."
        # self.calculateBetaBound()



    def simrank(self, v, u, R):
        pathsOfv = []
        pathsOfu = []
        for r in range(R):
            pathsOfv.append(self.randomWalk(v))
            pathsOfu.append(self.randomWalk(u))

        eta = 0
        for t in range(self.T):
            cntOfv = Counter()
            for r in range(R):
                cntOfv[pathsOfv[r][t]]+=1
            # print cntOfv

            cntOfu = Counter()
            for r in range(R):
                cntOfu[pathsOfu[r][t]]+=1
            # print cntOfu

            cntInCommon = cntOfv & cntOfu
            # print cntInCommon
            for key in cntInCommon:
                alpha = cntOfv[key]
                beta = cntOfu[key]
                eta += (math.pow(self.c,t)*(1-self.c)*alpha*beta) / (R*R)

        return eta


    def queryOneNode(self, v):
        neighbors = self.newG.neighbors(v)
        if len(neighbors) == 0:
            return []
        else:
            # find candidates with common neighbors
            candidates = []
            for neighbor in neighbors:
                uIDs = self.newG.neighbors(neighbor)
                for u in uIDs:
                    if u != v:
                        candidates.append(u)
            candidates = set(candidates)
            # print candidates

            #prune by gamma bound and beta bound
            prunedCandidates = []
            for u in candidates:
                maxGammaBound = max(self.gamma[u])

                if maxGammaBound < self.theta:
                    continue
                else:
                    prunedCandidates.append(u)
            # print prunedCandidates

            # calculate simrank between v and the pruned candidates
            srList = []
            for candidate in prunedCandidates:
                sr = self.simrank(v, candidate, 100)
                if sr > self.theta:
                    sr = self.simrank(v, candidate, 100)
                    srList.append((sr,candidate))

            nLargestSR = heapq.nlargest(self.K, srList)
            return nLargestSR

    def queryOneNodeByName(self, name):
        queriedName = name.lower()
        print "query top simrank for ", queriedName
        nodemapfile = "data/"+self.dataset+".dict"

        name2id = {}
        id2name = {}
        for line in open(nodemapfile, "r"):
            elements = line.strip().split("\t")
            name = elements[0].lower()
            id = int(elements[1])
            name2id[name] = id
            id2name[id] = name


        queriedID = name2id[queriedName]
        print "userid = ", queriedID
        nLargestSR = self.queryOneNode(queriedID)
        for value, id in nLargestSR:
            name = id2name[id]
            print name,":", value, ", ",



    def queryAllNodes(self):

        print "query top simrank for all nodes"


        f = open(self.outputFile, "w")
        outputList = []
        for v in self.G.nodes():
            outputList.append((v,self.queryOneNode(v)))

            # output to file
            if (v % 100 == 0):
                print "==", v
                for v, output in outputList:
                    if len(output) > 0:
                        f.write(str(output[0][1])+":"+str(output[0][0]))
                        if len(output) > 1:
                            for i in range(1, len(output)):
                                f.write("\t"+str(output[i][1])+":"+str(output[i][0]))
                        f.write("\n")
                    else:
                        f.write("\n")
                outputList = []

        for v, output in outputList:
            if len(output) > 0:
                f.write(str(output[0][1])+":"+str(output[0][0]))
                if len(output) > 1:
                    for i in range(1, len(output)):
                        f.write("\t"+str(output[i][1])+":"+str(output[i][0]))
                f.write("\n")
            else:
                f.write("\n")


        f.close()


    def main(self):
        self.preprocess()
        self.query()




if len(sys.argv) < 9:
    print "usage: dataset T P Q R1 R2 K theta"
    exit()
#dataset
dataset=sys.argv[1]
#steps of random walks
T = int(sys.argv[2])
#times of iterations for each each node when generating bipartite graph
P = int(sys.argv[3])
#times of random walks for each node when generating bipartite graph
Q = int(sys.argv[4])
#times of random walks in Algorithm 1 and 3
R1 = int(sys.argv[5])
#times of random walks in algorithm 2
R2 = int(sys.argv[6])
#top-K nodes to be returned
K = int(sys.argv[7])
#the threshod of upper bound
theta = float(sys.argv[8])


print "dataset=", dataset, ",path length =", T, ", P=", P, ", Q=", Q, ", R1=", R1, ", R2=", R2, ", K=", K, ", theta=", theta
simrank = ScalableSimRank(dataset=dataset,T=T, P=P, Q=Q, R1=R1, R2=R2, K=K, theta=theta)
simrank.preprocess()
#simrank.queryOneNodeByName("jiawei han")
simrank.queryAllNodes()


