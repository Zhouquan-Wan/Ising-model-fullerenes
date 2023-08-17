import numpy as np


def read_adjacent(path):
    with open(path, 'r') as f:
        lines = f.readlines()
    lines = lines[2:]
    adj_dist = {}
    for line in lines:
        columns = line.split()
        adj_dist[int(columns[-4])] = [int(num) for num in columns[-3:]]
    return adj_dist


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

vertex = np.array([[[0, 1], [1, 1]], [[1, 1], [1, 0]]], dtype=object)
adjacent_dict = read_adjacent('./data/C500/C500-0.xyz')
contract_order_list = [1, 5, 4, 3, 2, 78, 69, 149, 166, 157, 135, 144, 127, 113, 122, 105, 91, 100, 83, 61, 77, 68, 150, 167, 158, 165, 156, 134, 143, 136, 145, 128, 112, 121, 114, 123, 106, 90, 99, 92, 101, 84, 62, 79, 70, 76, 67, 151, 169, 160, 168, 159, 164, 155, 133, 142, 137, 146, 138, 147, 129, 111, 120, 115, 124, 116, 125, 107, 89, 98, 93, 102, 94, 103, 85, 63, 81, 72, 80, 71, 82, 7, 8, 163, 154, 162, 153, 161, 152, 170, 27, 28, 148, 130, 139, 131, 140, 132, 141, 22, 23, 126, 108, 117, 109, 118, 110, 119, 17, 18, 104, 86, 95, 87, 96, 88, 97, 12, 13, 75, 66, 74, 65, 73, 64, 6, 188, 171, 189, 172, 191, 173, 185, 14, 15, 11, 254, 237, 255, 238, 257, 239, 251, 19, 20, 16, 298, 281, 299, 282, 301, 283, 295, 24, 25, 21, 342, 325, 343, 326, 345, 327, 339, 29, 30, 26, 214, 199, 208, 200, 209, 201, 210, 9, 10, 179, 187, 180, 190, 182, 184, 176, 267, 276, 259, 245, 253, 246, 256, 248, 250, 242, 311, 320, 303, 289, 297, 290, 300, 292, 294, 286, 347, 364, 355, 333, 341, 334, 344, 336, 338, 330, 377, 386, 369, 196, 205, 203, 212, 202, 211, 193, 223, 232, 215, 178, 186, 181, 183, 175, 266, 275, 268, 277, 260, 244, 252, 247, 249, 241, 310, 319, 312, 321, 304, 288, 296, 291, 293, 285, 348, 365, 356, 363, 354, 332, 340, 335, 337, 329, 376, 385, 378, 387, 370, 197, 206, 204, 213, 194, 222, 231, 224, 233, 216, 177, 192, 174, 265, 274, 269, 278, 270, 279, 261, 243, 258, 240, 309, 318, 313, 322, 314, 323, 305, 287, 302, 284, 349, 367, 358, 366, 357, 362, 353, 331, 346, 328, 375, 384, 379, 388, 380, 389, 371, 198, 207, 195, 221, 230, 225, 234, 226, 235, 217, 46, 50, 280, 262, 271, 263, 272, 264, 273, 42, 41, 45, 324, 306, 315, 307, 316, 308, 317, 37, 36, 40, 361, 352, 360, 351, 359, 350, 368, 57, 56, 60, 390, 372, 381, 373, 382, 374, 383, 52, 51, 55, 236, 218, 227, 219, 228, 220, 229, 47, 49, 456, 438, 447, 439, 448, 440, 449, 43, 44, 434, 416, 425, 417, 426, 418, 427, 38, 39, 405, 396, 404, 395, 403, 394, 412, 58, 59, 500, 482, 491, 483, 492, 484, 493, 53, 54, 478, 460, 469, 461, 470, 462, 471, 48, 441, 450, 445, 454, 446, 455, 437, 419, 428, 423, 432, 424, 433, 415, 393, 411, 402, 410, 401, 406, 397, 485, 494, 489, 498, 490, 499, 481, 463, 472, 467, 476, 468, 477, 459, 442, 451, 444, 453, 436, 420, 429, 422, 431, 414, 392, 409, 400, 407, 398, 486, 495, 488, 497, 480, 464, 473, 466, 475, 458, 443, 452, 435, 421, 430, 413, 391, 408, 399, 487, 496, 479, 465, 474, 457, 33, 32, 31, 35, 34]

result = Tensor(vertex, [contract_order_list[0], contract_order_list[0], contract_order_list[0]], adjacent_dict)
cnt = 0
for contract_ten in contract_order_list[1:]:
    result.contract(vertex, contract_ten)
    cnt = cnt + 1
    print(cnt, contract_ten, result.ten.ndim)

print(result.ten)