import sys
import time
import random

    ############################ sort ##########################################
def msort(seq):
    if len(seq) == 1:
        return seq
    else:
        # recursion: break sequence down into chunks of 1
        mid = int(len(seq) / 2)
        left = msort(seq[:mid])
        right = msort(seq[mid:])

        i, j, k = 0, 0, 0  # i= left counter, j= right counter, k= master counter

        # run until left or right is out
        while i < len(left) and j < len(right):
            # if current left val is < current right val; assign to master list
            if left[i] < right[j]:
                seq[k] = left[i]
                i += 1
                k += 1
            # else assign right to master
            else:
                seq[k] = right[j]
                j += 1;
                k += 1

        remaining = left if i < j else right
        r = i if remaining == left else j

        while r < len(remaining):
            seq[k] = remaining[r]
            r += 1
            k += 1
        #print(seq)
        return seq


    ########################   Sequential search    #####################################

def seq_search(data_array, query_array):
    for i in range(len(query_array)):
        compare_arr(data_array, query_array[i])

def compare_arr(data_array, num):
    if num in data_array:
        start = time.clock()
        end = time.clock()
        print(end="True " + "%.2gs " % (end - start))
        return "True"
    else:
        start = time.clock()
        end = time.clock()
        print(end="False " + "%.2gs " % (end - start))
        return "False"


#########################      binary search    ########################################
def bin_search(data_array, query_array):
    for x in range(len(query_array)):
        compare_search(data_array, query_array[x])



def compare_search(data_array, num):
        left = 0
        right = len(data_array)
        #right1 = len(data_array)
        if left <= right:
            mid = int((left + right) / 2)
            if num < data_array[mid]:
                start = time.clock()
                left = mid - 1
                if num in data_array:
                    end = time.clock()
                    print(end="True " + "%.2gs " % (end - start))
                    return "True"
                elif num not in data_array:
                    start = time.clock()
                    end = time.clock()
                    print(end="False " + "%.2gs " % (end - start))
                    return "False"
            elif num > data_array[mid]:
                start = time.clock()
                left = mid + 1
                if num in data_array:
                    end = time.clock()
                    print(end="True " + "%.2gs " % (end - start))
                    return "True"
                elif num not in data_array:
                    start = time.clock()
                    end = time.clock()
                    print(end="False " + "%.2gs " % (end - start))
                    return "False"

            elif num == data_array[mid]:
                    start = time.clock()
                    end = time.clock()
                    print(end="True " + "%.2gs " % (end - start))
                    return "True"


###################################### user control #################################

size_of_data_array = sys.argv[1]
size_of_query_array = sys.argv[2]

a = int(size_of_data_array)  # cast to int
b = int(size_of_query_array)  # cast to int
Data_array = [0] * a  # inilitzing a new array with the length from the user input
Query_array = [0] * b

new_data_array = []
new_query_array = []
    ####################################    adding values to arrays    ###################
i = 0  # they say you dont need to initialize 'i' by zero, but i do it anyways
for i in range(len(Data_array)):
    responses = input()
    new_data_array.append(int(responses))  # cast to int
        # print(new_data_array)

for i in range(len(Query_array)):
    next_responses = input()
    new_query_array.append(int(next_responses))  # cast to int

#print(new_data_array)
#print(new_query_array)

"""
# test cases, 10,000 and 50,000
#10,000
r = []
for i in range (10000):
    r.append(random.randint(-10, 20))

s = []
for i in range (10000):
    s.append(random.randint(-10, 20))

bin_search(r,s)
seq_search(r,s)
"""


    #############################     call functions     ####################################


start_d = time.clock()
q_array = msort(new_query_array)
end_d = time.clock()
total_d = end_d - start_d

start_q = time.clock()
d_array = msort(new_data_array)
end_q = time.clock()
total_q = end_q - start_q

print("Total Prep Time " + "%.2gs" % (total_d + total_q))

bin_search(d_array, q_array)
print("\n")
seq_search(d_array, q_array)