import sys
import time

class Binary_Search_Tree:
    """Base Class representing a Binary Search Tree Structure"""
    # Jet fuel can't melt binary trees #

    class _BST_Node:
        """Private class for storing linked nodes with values and references to their siblings"""
        def __init__(self, value):
            """Node Constructor with 3 attributes"""
            self._value = value     # value being added
            self._left = None       # left sibling node (if val less than parent)
            self._right = None      # right sibling node (if val greater than parent)

    def __init__(self):
        """Binary Tree Constructor (creates an empty binary tree)"""
        self._root = None     # Binary tree root

    def get_root(self):
        return self._root

    def insert_element(self, value):
        """Method to insert an element into the BST"""
        if self._root is None:                         # if there is a root in the tree
            self._root = Binary_Search_Tree._BST_Node(value)   # create a new root equal to value
        else:
            self._insert_element(value, self._root)  # insert an element, calls recursive function

    def _insert_element(self, value, node):
        """Private method to Insert elements recursively"""
        # if node._value == value:  # if we come across a value already in the list
        #    raise ValueError
        if value < node._value:
            if node._left is not None:  # if a left node child node exists
                return self._insert_element(value, node._left)
                # call the insert element function again to evaluate next node
            else:
                node._left = Binary_Search_Tree._BST_Node(value)
                # after recursively crawling through the tree add value to left leaf
        else:
            if node._right is not None:  # if a right node child node exists
                return self._insert_element(value, node._right)
                # call the insert element function again to evaluate next node
            else:
                node._right = Binary_Search_Tree._BST_Node(value)
                # after recursively crawling through the tree add value to right leaf

    def find(self, value):
        """Return the if value is in tree"""
        current_node = self._root
        while current_node is not None:
            if value == current_node._value:
                print( end ="True ")
                return "True"
            elif value < current_node._value:
                current_node = current_node._left
            else:
                current_node = current_node._right
        print(end ="False ")
        return "False"



if __name__ == '__main__':
    P1 = Binary_Search_Tree()

    size_of_data_array = sys.argv[1]
    size_of_query_array = sys.argv[2]

    a = int(size_of_data_array)  # cast to int
    b = int(size_of_query_array)  # cast to int
    Data_array = [0] * a  # inilitzing a new array with the length from the user input
    Query_array = [0] * b


    ####################################    adding values to arrays    ###################
    start = 0
    end = 0
    while True:
        for i in range(len(Data_array)):
            item = int(input())
            start = time.clock()
            P1.insert_element(item)
            end = time.clock()
        else:
            break

    new_query_array = []
    for i in range(len(Query_array)):
        item = int(input())
        new_query_array.append(item)

    print("Total Prep Time " + "%.2gs" % (end - start))
    while True:
        for i in range(len(new_query_array)):
            start = time.clock()
            P1.find(new_query_array[i])
            end = time.clock()
            print( "%.2gs" % (end - start))
        else:
            break

