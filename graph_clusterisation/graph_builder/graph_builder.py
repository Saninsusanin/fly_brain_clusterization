import numpy as np

from math import ceil
from typing import List


def get_weight(distance, weights, step):
    if len(weights) == 1:
        return weights[0]

    left = 1
    right = len(weights)
    index = ceil(len(weights) / 2)

    while left <= right:
        if step * (index - 1) <= distance <= step * index:
            return weights[len(weights) - index]
        elif distance > step * index:
            left = index + 1
        elif distance < step * (index - 1):
            right = index - 1

        index = ceil((right - left) / 2) + left


def get_euclidean_distance(first_point, second_point):
    return (first_point[0] - second_point[0]) * (first_point[0] - second_point[0]) + \
           (first_point[1] - second_point[1]) * (first_point[1] - second_point[1]) + \
           (first_point[2] - second_point[2]) * (first_point[2] - second_point[2])


'''
    @brief this function updates graph
'''


def iterate_over_hemisphere(graph, hemisphere, mean_size_of_cluster, weights):
    first_vertex = 0
    first_point = np.zeros(3)
    second_point = np.zeros(3)
    step = mean_size_of_cluster / len(weights)

    while first_vertex < len(hemisphere):
        second_vertex = first_vertex + 1

        while second_vertex < len(hemisphere):
            first_point[0] = hemisphere.iloc[first_vertex, 1]
            first_point[1] = hemisphere.iloc[first_vertex, 2]
            first_point[2] = hemisphere.iloc[first_vertex, 3]
            second_point[0] = hemisphere.iloc[second_vertex, 1]
            second_point[1] = hemisphere.iloc[second_vertex, 2]
            second_point[2] = hemisphere.iloc[second_vertex, 3]
            distance = get_euclidean_distance(first_point, second_point)

            if distance <= mean_size_of_cluster:

                graph.append((int(hemisphere.iloc[first_vertex, 0]),
                              int(hemisphere.iloc[second_vertex, 0]),
                              get_weight(distance, weights, step)))

            second_vertex += 1

        print(first_vertex)
        first_vertex += 1


'''
    @param mean_size_of_cluster - radius of sphere which on average includes all vertices of cluster  
'''


def build_graph(data, mean_size_of_cluster, weights: List[int]):
    brain = data
    left_half = brain[brain.h == 'l']
    right_half = brain[brain.h == 'r']
    graph = []

    assert len(left_half) > 2, 'u have very small brain'
    assert len(right_half) > 2, 'u have very small brain'
    assert len(weights) != 0, 'not too many weights'
    assert 0 < mean_size_of_cluster, 'size of cluster is a positive number'

    iterate_over_hemisphere(graph, left_half, mean_size_of_cluster, weights)
    iterate_over_hemisphere(graph, right_half, mean_size_of_cluster, weights)

    return graph


def get_sliced_dataframe(data, percentage):
    sliced_data = None
    zones = data.zone.value_counts()

    for zone in list(zones.keys()):
        tmp = data[data.zone == zone]
        tmp = tmp.iloc[list(range(int(percentage * len(tmp))))]

        if sliced_data is None:
            sliced_data = tmp.copy(deep=True)
        else:
            sliced_data = sliced_data.append(tmp)

    return sliced_data
