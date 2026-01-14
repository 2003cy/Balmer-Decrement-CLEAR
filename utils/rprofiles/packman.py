import numpy as np

def calculate_angle(A, B, C):
    """
    Calculate the angle ABC (in degrees) given three points A, B, and C.
    B is the vertex of the angle.
    
    Parameters:
    A, B, C: Tuples representing (x, y) coordinates of the three points.
    
    Returns:
    Angle ABC in degrees.
    """
    BA = np.array(A) - np.array(B)
    BC = np.array(C) - np.array(B)
    
    cosine_angle = np.dot(BA, BC) / (np.linalg.norm(BA) * np.linalg.norm(BC))
    angle = np.arccos(np.clip(cosine_angle, -1.0, 1.0))
    
    return np.degrees(angle)


def double_packman(size, angle_left, angle_right):
    center = (size-2) / 2
    
    mask = np.zeros((size, size))
    for i in range(size):
        for j in range(size):
            angle_to_center = calculate_angle((i, j), (center, center), (center, 0))
            if (np.abs(angle_to_center) <= angle_left / 2 or np.abs(angle_to_center - 180) <= angle_right / 2):
                mask[i, j] = True
            else:
                mask[i, j] = False
    return mask