import sys
import time

    ########################   Sequential search    #####################################
def seq_search(data_array, query_array):
    for i in range(len(query_array)):
        compare_arr(data_array, query_array[i])

def compare_arr(data_array, num):
    if num in data_array:
        start = time.clock()
        end = time.clock()
        print("True " + "%.2gs " % (end - start))
        return "True"
    else:
        start = time.clock()
        end = time.clock()
        print("False " + "%.2gs " % (end - start))
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
                    print("True " + "%.2gs " % (end - start))
                    return "True"
                elif num not in data_array:
                    start = time.clock()
                    end = time.clock()
                    print("False " + "%.2gs " % (end - start))
                    return "False"
            elif num > data_array[mid]:
                start = time.clock()
                left = mid + 1
                if num in data_array:
                    end = time.clock()
                    print("True " + "%.2gs " % (end - start))
                    return "True"
                elif num not in data_array:
                    start = time.clock()
                    end = time.clock()
                    print("False " + "%.2gs " % (end - start))
                    return "False"

            elif num == data_array[mid]:
                    start = time.clock()
                    end = time.clock()
                    print("True " + "%.2gs " % (end - start))
                    return "True"



#################################      user control ###########################


size_of_data_array = sys.argv[1]
size_of_query_array = sys.argv[2]

a = int(size_of_data_array)  # cast to int
b = int(size_of_query_array)  # cast to int
Data_array = [0] * a  # inilitzing a new array with the length from the user input
Query_array = [0] * b

new_data_array = []
new_query_array = []
    ####################################    adding values to arrays    ###################
i = 0  # they say you dont need to initialize 'i' by zero, but i do it anyways4
for i in range(len(Data_array)):
    responses = input()
    new_data_array.append(int(responses))  # cast to int
        # print(new_data_array)

for i in range(len(Query_array)):
    next_responses = input()
    new_query_array.append(int(next_responses))  # cast to int

#print(new_data_array)

#print(new_query_array)
    #############################     call functions     ####################################

seq_search(new_data_array, new_query_array)
print("\n")
bin_search(new_data_array,new_query_array)

