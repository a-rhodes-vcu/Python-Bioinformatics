from __future__ import division
from __future__ import print_function
from math import sqrt
import os

while True:

    var_linkage = input("Single or Complete Linkage?: ")
    if var_linkage !='Single':
        print("Please Try again")
        continue
    elif var_linkage !='Complete':
        print("Please Try again")
        continue
    else:
        break




class Tree:
    ################ these are global variables available to every method in the class  ###########
    nodes = {}

    distance = {}

    dummy = -9999999
    linkage = str

    bestDistance = nextBestDistance = dummy
    bestKeyLeft = bestKeyRight = nextBestKeyLeft = nextBestKeyRight = str

    metric = str

    def __init__(self, internal_name, gene_name, data, metric='Pearson',right=None,
                 left=None, nodes_represented=None,var_linkage='Complete'): #<----MESSAGE TO THE USER, PLEASE INDICATE SINGLE OR COMPLETE

        # PLEASE NOTE TO THE USER: THIS PROGRAM CAN ONLY TAKE EITHER SINGLE OR COMPLETE AS PARAMETERS
        #AS IT STANDS NOW, AVERAGE LINKING CLUSTERING IS NOT AN OPTION

        self.internal_name = internal_name
        self.name = gene_name
        self.right = right
        self.left = left

        self.taken = False

        self.data = data

        self.meh = var_linkage
        Tree.linkage = self.meh

        ################################    NODES REPRESENTED, CREATED A LIST OF EACH NODE AS IT'S MADE AND IT'S GENES
        ################################    WAS A DICT, DICT IS GONE
        if not nodes_represented:
            self.nodes_represented = [internal_name]
        else:
            self.nodes_represented = nodes_represented

############################################# CALLING THE DISTANCE FUNCTIONS, AND THE CLUSTERING FUNCTIONS  #########
        if metric == 'Euclidian':
            self.distance_function = self.__Euclidian_distance
            Tree.metric = 'Euclidian'
        elif metric == 'Pearson':
            self.distance_function = self.__Pearson_correlation
            Tree.metric = 'Pearson'

        for other_node in Tree.nodes.values():
            if other_node.taken is False: #and other_node.internal_name.startswith("N"):
                if metric == 'Single':
                    self.distance_function = self.__Single
                    Tree.metric = 'Single'
                elif metric == 'Complete':
                    self.distance_function = self.__Complete
                    Tree.metric = 'Complete'

        #print ("Instantiating node " + self.internal_name)

        for other_node in Tree.nodes.values():
            if other_node.taken is False:

                node_pair = (self.internal_name, other_node.internal_name)
                Tree.distance[node_pair] = self.distance_function(other_node)

                #print ("Distance between", self.internal_name, "and", other_node.internal_name, "is", Tree.distance[node_pair])
                Tree.distance[other_node.internal_name, self.internal_name] = Tree.distance[node_pair]

                if Tree.distance[node_pair] > Tree.bestDistance:

                        Tree.bestDistance = Tree.distance[node_pair]
                        Tree.bestKeyLeft = self.internal_name
                        Tree.bestKeyRight = other_node.internal_name
                        #print ("New champ is the pair", node_pair, "and the correlation is", Tree.bestDistance)

        Tree.nodes[self.internal_name] = self

#####################################   made the GTR file in class ###################################################
        if "GENE" not in self.internal_name:
            gtr_file = open("Bacillus_file.gtr","a")
            gtr_file.write("{}X\t{}X\t{}X\t{}\n".format(self.internal_name,Tree.bestKeyLeft,Tree.bestKeyRight,Tree.bestDistance))

#############################################   Distance    ###############################
    def __Euclidian_distance(self, other_node):
        x_vector = self.data
        y_vector = other_node.data

        deltas = 0

        for i in range(len(x_vector)):

            deltas += (x_vector[i] - y_vector[i])**2

        return -sqrt(deltas)


    def __Pearson_correlation(self, other_node):
        #did this in class
        x_data = self.data
        y_data = other_node.data
        xi_sum = 0
        yi_sum = 0
        xiyi_sum = 0
        x_2 = 0
        y_2 = 0
        for i in range(len(x_data)):
            xiyi_sum += x_data[i] * y_data[i]
            xi_sum += x_data[i]
            yi_sum += y_data[i]
            x_2 += x_data[i] ** 2
            y_2 += y_data[i] ** 2

        p = ((len(x_data)*xiyi_sum) - (xi_sum*yi_sum)) / ((sqrt((len(x_data)*x_2)-(xi_sum**2)))*(sqrt((len(y_data)*(y_2))-(yi_sum**2))))
        return p

################################################### Clustering  ##########################################
    def __Single(self, other):
        #this clustering method takes the best of the best
        dum= Tree.dummy
        left_child = self.nodes_represented
        right_child = other.nodes_represented
        for x in left_child:
            for y in right_child:
                value = Tree.distance[(x, y)]
                if value > dum:
                        dum = value
        return dum


    def __Complete(self, other):
        #this clustering method takes the best of the worst
        dum = Tree.dummy * -1
        left_child = self.nodes_represented
        right_child = other.nodes_represented
        for x in left_child:
            for y in right_child:
                value = Tree.distance[(x, y)]
                if value < dum:
                    dum = value
        return dum

################################################### Tree traversal, get a linear, untangled ordering of genes
    @staticmethod
    def DFS(root,new_file):
        #used code from a previous cource, also worked on this in class
        cdt_file = open("Bacillus_file.cdt", 'a')
        if root is not None:
            gene_vector = ""
            if "GENE" in root.internal_name:
                for k in root.data[0:]:
                    gene_vector = str(k) + "\t" + gene_vector
                new_file.write("{}X\t{}\t{}\t{}\n".format(root.internal_name,root.name,root.name,gene_vector))

            if root.left:
                Tree.DFS(Tree.nodes[root.left],new_file)

            if root.right:
                Tree.DFS(Tree.nodes[root.right],new_file)

    @staticmethod
    def __average(x, y):

        return [(x[i] + y[i]) / 2 for i in range(len(x))]

    @staticmethod
    def Cluster():
        for i in range(len(Tree.nodes) - 1):
            #print ("\nCreating node", i,"consisting of", Tree.bestKeyLeft, "and", Tree.bestKeyRight, "with distance", Tree.bestDistance)

            Tree.nodes[Tree.bestKeyLeft].taken = True
            Tree.nodes[Tree.bestKeyRight].taken = True

            leftData = Tree.nodes[Tree.bestKeyLeft].data
            rightData = Tree.nodes[Tree.bestKeyRight].data
            tempdata = Tree.__average(leftData, rightData)
            #print("New vector is", tempdata)
            tempname = "NODE" + str(i+1)

            Tree(tempname, tempname, tempdata, Tree.linkage, Tree.bestKeyLeft, Tree.bestKeyRight, nodes_represented=
                 Tree.nodes[Tree.bestKeyLeft].nodes_represented + Tree.nodes[Tree.bestKeyRight].nodes_represented)

            Tree.bestDistance = Tree.dummy
            for node_pair, cur_distance in Tree.distance.items():

                if cur_distance > Tree.bestDistance:
                    left_one, right_one = node_pair
                    if (Tree.nodes[left_one].taken is False) and (Tree.nodes[right_one].taken is False):
                        Tree.bestDistance = cur_distance
                        Tree.bestKeyLeft = left_one
                        Tree.bestKeyRight = right_one

        print ("Finished clustering: Please check the output files")
        return
    ########################################    Print out what nodes are being made with what genes. very helpful
    def node_dump(self):
        #surprising helpful method. I was able to know which nodes where being produced, and with each gene.
        print ("Node name is: " + self.internal_name)
        print ("Node to left is: " + str(self.left))
        print ("Node to right is: " + str(self.right))
        print ("Paired node?: " + str(self.taken))
        print ("Data associated with node: ", self.data)
        print ("Nodes this node represents: ", self.nodes_represented)
        return


def main():
    var_distance = input("Please enter Pearson or Euclidian: ")

    i = 0
    header = ''
    infile = open('BacillusData2.txt', 'r')
    for line in infile.readlines():
        line = line.strip()
        tempdata = line.split('\t')
        #print(tempdata)
        if tempdata[0] != 'GENE':

            Tree("GENE" + str(i+1), tempdata[0], [float(k) for k in tempdata[1:]], var_distance)
            i += 1
        elif tempdata[0] == 'GENE':
            header = ' '.join(tempdata[1:])
    #print(header)


############################################## FILE MAKING #################################

    if os.path.isfile("Bacillus_file.cdt"):
        os.remove("Bacillus_file.cdt")
    if os.path.isfile("Bacillus_file.gtr"):
        os.remove("Bacillus_file.gtr")

    bac_cdt = open("Bacillus_file.cdt", 'w')
    bac_cdt.write(
        "{}\t{}\t{}\t\n".format("GID", "Name", "Name2 " + header)) 

    Tree.Cluster()

    Tree.DFS(Tree.nodes["NODE"+str(i-2)],bac_cdt)

if __name__ == '__main__':
    main()