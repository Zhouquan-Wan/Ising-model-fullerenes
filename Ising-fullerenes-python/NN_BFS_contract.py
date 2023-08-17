import numpy as np
from collections import deque
import os

def read_adjacent(path):
    with open(path, 'r') as f:
        lines = f.readlines()
    lines = lines[2:]
    adj_dist = {}
    for line in lines:
        columns = line.split()
        adj_dist[int(columns[-4])] = [int(num) for num in columns[-3:]]
    return adj_dist


def get_distance(adj_dist):
    visited = set()
    queue = deque([(1, 0)])
    dist_dict = {}
    while queue:
        node, distance = queue.popleft()
        if node not in visited:
            visited.add(node)
            if distance in dist_dict:
                dist_dict[distance].append(node)
            else:
                dist_dict[distance] = [node]
            for neighbor in adj_dist[node]:
                if neighbor not in visited:
                    queue.append((neighbor, distance + 1))
    return dist_dict


class Tensor:
    def __init__(self, ten, ten_idx_list, adj_dict):
        self.ten = ten.copy()
        self.ten_idx_list = ten_idx_list.copy()
        self.adj_dict = adj_dict.copy()

    def contract(self, ten_other, ten_idx_other):
        contract_axes_idx_list = []
        for ten_idx in self.adj_dict[ten_idx_other]:
            if ten_idx in self.ten_idx_list:
                contract_axes_idx_list.append(self.ten_idx_list.index(ten_idx))
        contract_axes_num = len(contract_axes_idx_list)
        contract_axes_idx_list_other = list(range(contract_axes_num))
        self.ten = np.tensordot(self.ten, ten_other, axes=(contract_axes_idx_list, contract_axes_idx_list_other))
        for axis_idx in reversed(sorted(contract_axes_idx_list)):
            self.ten_idx_list.pop(axis_idx)
        for axis_idx in range(contract_axes_num, ten_other.ndim):
            self.ten_idx_list.append(ten_idx_other)


def subdir_sort_key(dir_name):
    return int(dir_name[1:])


parent_directory_path = './data'
N_list = []
deg_list = []
subdir_list = sorted([subdir for subdir in os.listdir(parent_directory_path) if subdir.startswith('C')], key=subdir_sort_key)
for subdir in subdir_list:
    directory_path = os.path.join(parent_directory_path, subdir)
    file_name_list = [file for file in os.listdir(directory_path) if file.startswith('C')]
    for file_name in file_name_list:
        N = int(subdir[1:])
        if N < 180:
            N_list.append(N)
            file_path = os.path.join(directory_path, file_name)
            vertex = np.array([[[0, 1], [1, 1]], [[1, 1], [1, 0]]], dtype=object)
            adjacent_dict = read_adjacent(file_path)
            distance_dict = get_distance(adjacent_dict)
            distance_max = max(distance_dict)
            result = Tensor(vertex, [distance_dict[0][0], distance_dict[0][0], distance_dict[0][0]], adjacent_dict)
            for d in range(1, distance_max + 1):
                possible_ten_idx_list = distance_dict[d]
                for possible_ten_idx in possible_ten_idx_list:
                    result.contract(vertex, possible_ten_idx)
            deg_list.append(int(result.ten))
            print(file_path, ":", int(result.ten))

print("N_list =", N_list)
print("deg_list =", deg_list)




