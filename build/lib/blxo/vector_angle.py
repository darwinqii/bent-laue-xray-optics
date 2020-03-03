import numpy as np


def angle_between_vectors(a, b):
    a = np.array(a)
    b = np.array(b)
    La = np.sqrt(a.dot(a))
    Lb = np.sqrt(b.dot(b))
    cos_alpha = a.dot(b) / (La * Lb)
    alpha_radians = np.arccos(cos_alpha)
    alpha_degrees = np.degrees(alpha_radians)
    return alpha_degrees


if __name__ == '__main__':
    print(angle_between_vectors((1, 2, 3), (2, 2, 4)))