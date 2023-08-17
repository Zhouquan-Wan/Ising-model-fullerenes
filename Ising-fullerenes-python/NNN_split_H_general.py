import numpy as np
from itertools import product


def get_pentagon_energy(s_list, j1, j2, ratio):
    assert(len(s_list) == 10)
    s_list_raw = [2 * s - 1 for s in s_list]
    e = j1 * s_list_raw[0] * s_list_raw[1] + \
        j1 * s_list_raw[1] * s_list_raw[2] + \
        j1 * s_list_raw[2] * s_list_raw[3] + \
        j1 * s_list_raw[3] * s_list_raw[4] + \
        j1 * s_list_raw[4] * s_list_raw[0] + \
        (1 - ratio) * j1 * s_list_raw[0] * s_list_raw[5] + \
        (1 - ratio) * j1 * s_list_raw[1] * s_list_raw[6] + \
        (1 - ratio) * j1 * s_list_raw[2] * s_list_raw[7] + \
        (1 - ratio) * j1 * s_list_raw[3] * s_list_raw[8] + \
        (1 - ratio) * j1 * s_list_raw[4] * s_list_raw[9] + \
        (1 - ratio) * j1 * s_list_raw[0] * s_list_raw[9] + \
        (1 - ratio) * j1 * s_list_raw[1] * s_list_raw[5] + \
        (1 - ratio) * j1 * s_list_raw[2] * s_list_raw[6] + \
        (1 - ratio) * j1 * s_list_raw[3] * s_list_raw[7] + \
        (1 - ratio) * j1 * s_list_raw[4] * s_list_raw[8] + \
        j2 * s_list_raw[0] * s_list_raw[2] + \
        j2 * s_list_raw[1] * s_list_raw[3] + \
        j2 * s_list_raw[2] * s_list_raw[4] + \
        j2 * s_list_raw[3] * s_list_raw[0] + \
        j2 * s_list_raw[4] * s_list_raw[1]
    return e


def get_hexagon_energy(s_list, j1, j2, ratio):
    assert (len(s_list) == 6)
    s_list_raw = [2 * s - 1 for s in s_list]
    e = ratio * j1 * s_list_raw[0] * s_list_raw[1] + \
        ratio * j1 * s_list_raw[1] * s_list_raw[2] + \
        ratio * j1 * s_list_raw[2] * s_list_raw[3] + \
        ratio * j1 * s_list_raw[3] * s_list_raw[4] + \
        ratio * j1 * s_list_raw[4] * s_list_raw[5] + \
        ratio * j1 * s_list_raw[5] * s_list_raw[0] + \
        j2 * s_list_raw[0] * s_list_raw[2] + \
        j2 * s_list_raw[1] * s_list_raw[3] + \
        j2 * s_list_raw[2] * s_list_raw[4] + \
        j2 * s_list_raw[3] * s_list_raw[5] + \
        j2 * s_list_raw[4] * s_list_raw[0] + \
        j2 * s_list_raw[5] * s_list_raw[1]
    return e


j1 = 1
j2_list = list(np.linspace(0, 1.5, 16)) + [2, 3, 4, 5, 10]
E_list = []
deg_list = []
for j2 in j2_list:
    if j2 / j1 > 0:
        ratio = min(j2 / j1, 1)
    else:
        ratio = 0
    print("j1 = ", j1)
    print("j2 = ", j2)

    ten_energy_pentagon = np.zeros([2] * 10)
    ten_energy_hexagon = np.zeros([2] * 6)
    for s_tup in product([0, 1], repeat=10):
        ten_energy_pentagon[s_tup] = get_pentagon_energy(s_tup, j1, j2, ratio)
    for s_tup in product([0, 1], repeat=6):
        ten_energy_hexagon[s_tup] = get_hexagon_energy(s_tup, j1, j2, ratio)

    ten_energy_pentagon_min = ten_energy_pentagon.min()
    ten_energy_hexagon_min = ten_energy_hexagon.min()
    print("ground state energy pentagon = ", ten_energy_pentagon_min)
    print("ground state energy hexagon = ", ten_energy_hexagon_min)
    print("ground state energy lower bound = ", 12 * ten_energy_pentagon_min + 20 * ten_energy_hexagon_min)
    E_list.append(12 * ten_energy_pentagon_min + 20 * ten_energy_hexagon_min)

    tol = 1e-10
    ten_pentagon = np.where(abs(ten_energy_pentagon - ten_energy_pentagon_min) < tol, 1, 0)
    ten_hexagon = np.where(abs(ten_energy_hexagon - ten_energy_hexagon_min) < tol, 1, 0)

    print("ground state pentagon:")
    loc = np.where(ten_pentagon == 1)
    print(list(zip(*loc)))
    print("num = ", len(list(zip(*loc))))

    print("ground state hexagon:")
    loc = np.where(ten_hexagon == 1)
    print(list(zip(*loc)))
    print("num = ", len(list(zip(*loc))))

    ten_pentagon = ten_pentagon.astype(float)
    ten_hexagon = ten_hexagon.astype(float)

    ten_temp1 = np.einsum('ABCDEFGHIJ,abFAJo,deGBFc,ghHCGf,jkIDHi,mnJEIl->abFcdeGfghHijkIlmnJo', ten_pentagon, ten_hexagon, ten_hexagon, ten_hexagon, ten_hexagon, ten_hexagon)

    ten_hex_pen = np.einsum('abcdef,ghijdklmec->abcgkhlimjef', ten_hexagon, ten_pentagon)

    ten_temp2_1 = np.einsum('ABCDEFGHIJKLMNOPQRST,abcdeUEDCBAY->abcdeUEFGHIJKLMNOPQRSTAY', ten_temp1, ten_hex_pen)
    ten_temp2_2 = np.einsum('abcdeUEFGHIJKLMNOPQRSTAY,efghiVIHGFEU->abcdefghiVIJKLMNOPQRSTAY', ten_temp2_1, ten_hex_pen)
    ten_temp2_3 = np.einsum('abcdefghiVIJKLMNOPQRSTAY,ijklmWMLKJIV->abcdefghijklmWMNOPQRSTAY', ten_temp2_2, ten_hex_pen)
    ten_temp2_4 = np.einsum('abcdefghijklmWMNOPQRSTAY,mnopqXQPONMW->abcdefghijklmnopqXQRSTAY', ten_temp2_3, ten_hex_pen)
    ten_half = np.einsum('abcdefghijklmnopqXQRSTAY,qrstaYATSRQX->abcdefghijklmnopqrst', ten_temp2_4, ten_hex_pen)

    result = np.einsum('abcdefghijklmnopqrst,stabcdefghijklmnopqr', ten_half, ten_half)
    deg_list.append(result)

    print("ground state energy lower bound degeneracy = ", result)
    print('\n')

print("E_list =", E_list)
print("degeneracy_list =", deg_list)