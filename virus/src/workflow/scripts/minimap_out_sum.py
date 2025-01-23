import sys
import csv

def write_dict(filename):
    query_dict = {}
    target_dict = {}
    with open(filename) as the_file:
        for line in the_file:

             query_contig = line.split('\t')[0]
             query_length = int(line.split('\t')[1])
             query_start = int(line.split('\t')[2])
             query_end = int(line.split('\t')[3])

             target_contig = line.split('\t')[5]
             target_length = int(line.split('\t')[6])
             target_start = int(line.split('\t')[7])
             target_end = int(line.split('\t')[8])

             if query_contig not in query_dict:
                 query_dict[query_contig] = [query_length, [[query_start, query_end]]]
             elif query_contig in query_dict:
                 query_dict[query_contig][-1].append([query_start, query_end])

             if target_contig not in target_dict:
                target_dict[target_contig] = [target_length, [[target_start, target_end]]]
             elif target_contig in target_dict:
                target_dict[target_contig][-1].append([target_start, target_end])

        dict_list = list([query_dict, target_dict])
    return dict_list

def mergeInterval(dict_list):

    query_merge_dict = {}
    target_merge_dict = {}

    query_dict = dict_list[0]
    target_dict = dict_list[1]

    for index in query_dict:
        contig = index
        length = query_dict[index][0]
        arr = query_dict[index][1]
        # print(index,": ", length, arr)
        arr.sort(key = lambda x: x[0])

        # array to hold the merged intervals
        m = []
        s = -10000
        max = -100000
        for i in range(len(arr)):
            a = arr[i]
            if a[0] > max:
                if i != 0:
                    m.append([s,max])
                max = a[1]
                s = a[0]
            else:
                if a[1] >= max:
                    max = a[1]

        #'max' value gives the last point of
        # that particular interval
        # 's' gives the starting point of that interval
        # 'm' array contains the list of all merged intervals

        if max != -100000 and [s, max] not in m:
            # interval = max - s
            m.append([s, max])
        # print("\n", index, ":", length, end = " ")
        mergeinterval = 0
        for i in range(len(m)):
            # print(m[i], end = " ")
            interval = m[i][1]-m[i][0]
            # print(interval, end = " ")
            mergeinterval += interval
            ratio = mergeinterval/length
        # print("mergeinterval", mergeinterval, ratio, end = " ")
        query_merge_dict[contig] = [length, mergeinterval, ratio]
    # print(query_merge_dict)
    with open(sys.argv[2], 'w') as csv_file:
        writer = csv.writer(csv_file)
        for key, value in query_merge_dict.items():
           writer.writerow([key] + value)

    for index in target_dict:
        contig = index
        length = target_dict[index][0]
        arr = target_dict[index][1]
        # print(index,": ", length, arr)
        arr.sort(key = lambda x: x[0])

        # array to hold the merged intervals
        m = []
        s = -10000
        max = -100000
        for i in range(len(arr)):
            a = arr[i]
            if a[0] > max:
                if i != 0:
                    m.append([s,max])
                max = a[1]
                s = a[0]
            else:
                if a[1] >= max:
                    max = a[1]

        #'max' value gives the last point of
        # that particular interval
        # 's' gives the starting point of that interval
        # 'm' array contains the list of all merged intervals

        if max != -100000 and [s, max] not in m:
            # interval = max - s
            m.append([s, max])
        # print("\n", index, ":", length, end = " ")
        mergeinterval = 0
        for i in range(len(m)):
            # print(m[i], end = " ")
            interval = m[i][1]-m[i][0]
            # print(interval, end = " ")
            mergeinterval += interval
            ratio = mergeinterval/length
        # print("mergeinterval", mergeinterval, ratio, end = " ")
        target_merge_dict[contig] = [length, mergeinterval, ratio]
    # print(target_merge_dict)
    with open(sys.argv[3], 'w') as csv_file:
        writer = csv.writer(csv_file)
        for key, value in target_merge_dict.items():
           writer.writerow([key] + value)

filename = sys.argv[1]
mergeInterval(write_dict(filename))
