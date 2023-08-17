import numpy as np
from itertools import product


def get_pentagon_energy(s_list, j1, j2):
    assert(len(s_list) == 5)
    s_list_raw = [2 * s - 1 for s in s_list]
    e = j1 * s_list_raw[0] * s_list_raw[1] + \
        j1 * s_list_raw[1] * s_list_raw[2] + \
        j1 * s_list_raw[2] * s_list_raw[3] + \
        j1 * s_list_raw[3] * s_list_raw[4] + \
        j1 * s_list_raw[4] * s_list_raw[0] + \
        j2 * s_list_raw[0] * s_list_raw[2] + \
        j2 * s_list_raw[1] * s_list_raw[3] + \
        j2 * s_list_raw[2] * s_list_raw[4] + \
        j2 * s_list_raw[3] * s_list_raw[0] + \
        j2 * s_list_raw[4] * s_list_raw[1]
    return e


def get_hexagon_energy(s_list, j1, j2):
    assert(len(s_list) == 6)
    s_list_raw = [2 * s - 1 for s in s_list]
    e = j1 * s_list_raw[0] * s_list_raw[1] + \
        j1 * s_list_raw[1] * s_list_raw[2] + \
        j1 * s_list_raw[2] * s_list_raw[3] + \
        j1 * s_list_raw[3] * s_list_raw[4] + \
        j1 * s_list_raw[4] * s_list_raw[5] + \
        j1 * s_list_raw[5] * s_list_raw[0] + \
        j2 * s_list_raw[0] * s_list_raw[2] + \
        j2 * s_list_raw[1] * s_list_raw[3] + \
        j2 * s_list_raw[2] * s_list_raw[4] + \
        j2 * s_list_raw[3] * s_list_raw[5] + \
        j2 * s_list_raw[4] * s_list_raw[0] + \
        j2 * s_list_raw[5] * s_list_raw[1]
    return e


j1 = 1
j2 = 1.5
print("j1 = ", j1)
print("j2 = ", j2)
ten_energy_pentagon = np.zeros([2, 2, 2, 2, 2])
ten_energy_hexagon = np.zeros([2, 2, 2, 2, 2, 2])

for s0, s1, s2, s3, s4 in product([0, 1], [0, 1], [0, 1], [0, 1], [0, 1]):
    ten_energy_pentagon[s0, s1, s2, s3, s4] = get_pentagon_energy([s0, s1, s2, s3, s4], j1, j2)

for s0, s1, s2, s3, s4, s5 in product([0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1]):
    ten_energy_hexagon[s0, s1, s2, s3, s4, s5] = get_hexagon_energy([s0, s1, s2, s3, s4, s5], j1, j2)

ten_energy_pentagon_min = ten_energy_pentagon.min()
ten_energy_hexagon_min = ten_energy_hexagon.min()

print("ground state energy pentagon = ", ten_energy_pentagon_min)
print("ground state energy hexagon = ", ten_energy_hexagon_min)
print("ground state energy lower bound = ", 12 * ten_energy_pentagon_min + 20 * ten_energy_hexagon_min)

tol = 1e-10
ten_pentagon = np.where(abs(ten_energy_pentagon - ten_energy_pentagon_min) < tol, 1, 0)
ten_hexagon = np.where(abs(ten_energy_hexagon - ten_energy_hexagon_min) < tol, 1, 0)

print("ground state pentagon:")
loc = np.where(ten_pentagon == 1)
print(list(zip(*loc)))

print("ground state hexagon:")
loc = np.where(ten_hexagon == 1)
print(list(zip(*loc)))

ten_temp1 = np.einsum('ABCDE,abcGAF,defHBG,ghiICH,jklJDI,mnoFEJ->abcdefghijklmno', ten_pentagon, ten_hexagon, ten_hexagon, ten_hexagon, ten_hexagon, ten_hexagon)

ten_temp2_1 = np.einsum('ABCDEFGHIJKLMNO,abcQBP,dRDCQ->APabcdREFGHIJKLMNO', ten_temp1, ten_hexagon, ten_pentagon)
ten_temp2_2 = np.einsum('APabcdREFGHIJKLMNO,efgSER,hTGFS->APabcdefghTHIJKLMNO', ten_temp2_1, ten_hexagon, ten_pentagon)
ten_temp2_3 = np.einsum('APabcdefghTHIJKLMNO,ijkUHT,lVJIU->APabcdefghijklVKLMNO', ten_temp2_2, ten_hexagon, ten_pentagon)
ten_temp2_4 = np.einsum('APabcdefghijklVKLMNO,mnoWKV,pXMLW->APabcdefghijklmnopXNO', ten_temp2_3, ten_hexagon, ten_pentagon)
ten_temp2 = np.einsum('APabcdefghijklmnopXNO,qrsYNX,tPAOY->abcdefghijklmnopqrst', ten_temp2_4, ten_hexagon, ten_pentagon)

result = np.einsum('abcdefghijklmnopqrst,stabcdefghijklmnopqr', ten_temp2, ten_temp2)
print("ground state energy lower bound degeneracy = ", result)
