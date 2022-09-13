import sys

import numpy
import numpy as np
import pandas as pd
import os
import spkmeansmodule

def k_means_pp(data_points, K):
    i = 1
    num = int(np.random.choice(range(len(data_points))))
    first_centroid = np.array(data_points[num])
    centroids = np.array([first_centroid])
    arr_of_index = [num]
    while i < K:
        arr_of_Dl = create_arr_of_Dl(centroids, data_points)
        i += 1
        num = int(np.random.choice(range(len(data_points)), p=probability(arr_of_Dl)))
        centroids = np.append(centroids, [data_points[num]], axis=0)
        arr_of_index.append(num)
    return arr_of_index


def find_the_Dl(centroids, vector):
    min_dis = float('inf')
    for i in range(len(centroids)):
        dis = np.sum(np.power((centroids[i] - vector), 2))
        if dis < min_dis:
            min_dis = dis
    return min_dis


def create_arr_of_Dl(centroids, data_points):
    arr_to_return = np.array([])
    for i in range(len(data_points)):
        arr_to_return = np.append(arr_to_return, find_the_Dl(centroids, data_points[i]))

    return arr_to_return


def probability(arr_of_Dl):
    arr_to_return = np.array([])
    for i in range(len(arr_of_Dl)):
        arr_to_return = np.append(arr_to_return, (arr_of_Dl[i] / np.sum(arr_of_Dl)))
    return arr_to_return


if __name__ == '__main__':
    arr = sys.argv
    if len(arr) != 4:
        print("Invalid Input!(1)")
        exit(1)

    else:
        try:
            k = int(arr[1])
            goal = arr[2]
            file_name = arr[3]
        except ValueError:
            print("Invalid Input!(2)")
            exit(1)

    np.random.seed(0)

    if goal == "spk":
        try:
            c_mat = spkmeansmodule.spk_module(file_name, k)

            data_points = numpy.array(c_mat)

            centroids = []

            index_of_centroid = k_means_pp(data_points, len(c_mat[0]))

            print(','.join(str(num) for num in index_of_centroid))

            for i in range(len(index_of_centroid)):
                centroids.append(c_mat[index_of_centroid[i]])

            list_of_centroids = spkmeansmodule.fit(len(data_points), len(data_points[0]), len(c_mat[0]),
                                                   300, 0, centroids, c_mat)

            for centroid in list_of_centroids:
                print(','.join(str('{:.4f}'.format(num)) for num in centroid))

        except IOError:
            print("Invalid Input!")
            exit(1)


    else:
        try:
            spkmeansmodule.wam_ddg_lnorm_jacobi(file_name, goal)

        except IOError:
            print("Invalid Input!")
            exit(1)
        except Exception:
            print("An Error Has Occurred")
            exit(1)