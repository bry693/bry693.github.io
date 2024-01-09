'''
STUDENT NAME: Bryanna Hernandez
Project 2
Fall 2022
COMP 313: Computer Graphics
Professor Schiffer
'''
# IMPORT LIBRARIES
import matplotlib.pyplot as plt
from math import sin, cos, radians
import sys
figure, axes = plt.subplots()
axes.set_aspect(1)
plt.axis('on')
plt.grid(True)



def objExport(file):
    file = open(file)

    list_of_lines = file.readlines()

    for line in list_of_lines:
        if line[0] == 'v':
            line_v = line.strip('v \n')
            pnt_vect = line_v.split(' ')

            for index, str in enumerate(pnt_vect):
                pnt_vect[index] = float(str)

            plane_x.append(pnt_vect[0])
            plane_y.append(pnt_vect[1])
            plane_z.append(pnt_vect[2])


        if line[0] == 'f':

            line_f = line.strip('f \n')

            face = line_f.split(' ')

            for index, str_num in enumerate(face):
                face[index] = int(str_num) - 1
            plane_f.append(face)

# PRESERVED ROTATED COORDINATES
rotd_x_coords = []
rotd_y_coords = []
rotd_z_coords = []
rotd_sequence_coords = []


# ROTATIONAL MATRIX FUNCTIONS

# Rotates a point vector from a given center a given number of degrees CC for + numbers on the x-axis
def rotx(center_vector, point_vector, angle):
    angle = radians(angle)

    rotx_mat = [
        [1, 0, 0],
        [0, cos(angle), -sin(angle)],
        [0, sin(angle), cos(angle)]
    ]

    xprod = rotx_mat[0][0] * point_vector[0] + rotx_mat[0][1] * point_vector[1] + rotx_mat[0][2] * point_vector[2]
    yprod = rotx_mat[1][0] * point_vector[0] + rotx_mat[1][1] * point_vector[1] + rotx_mat[1][2] * point_vector[2]
    zprod = rotx_mat[2][0] * point_vector[0] + rotx_mat[2][1] * point_vector[1] + rotx_mat[2][2] * point_vector[2]

    xg = xprod + center_vector[0]
    yg = yprod + center_vector[1]
    zg = zprod + center_vector[2]

    return [xg, yg, zg]


def roty(center_vector, point_vector, angle):
    angle = radians(angle)

    roty_mat = [
        [cos(angle), 0, sin(angle)],
        [0, 1, 0],
        [-sin(angle), 0, cos(angle)]
    ]

    xprod = roty_mat[0][0] * point_vector[0] + roty_mat[0][1] * point_vector[1] + roty_mat[0][2] * point_vector[2]
    yprod = roty_mat[1][0] * point_vector[0] + roty_mat[1][1] * point_vector[1] + roty_mat[1][2] * point_vector[2]
    zprod = roty_mat[2][0] * point_vector[0] + roty_mat[2][1] * point_vector[1] + roty_mat[2][2] * point_vector[2]

    xg = xprod + center_vector[0]
    yg = yprod + center_vector[1]
    zg = zprod + center_vector[2]

    return [xg, yg, zg]


def rotz(center_vector, point_vector, angle):
    angle = radians(angle)

    rotz_mat = [
        [cos(angle), -sin(angle), 0],
        [sin(angle), cos(angle), 0],
        [0, 0, 1]
    ]

    xprod = rotz_mat[0][0] * point_vector[0] + rotz_mat[0][1] * point_vector[1] + rotz_mat[0][2] * point_vector[2]
    yprod = rotz_mat[1][0] * point_vector[0] + rotz_mat[1][1] * point_vector[1] + rotz_mat[1][2] * point_vector[2]
    zprod = rotz_mat[2][0] * point_vector[0] + rotz_mat[2][1] * point_vector[1] + rotz_mat[2][2] * point_vector[2]

    xg = xprod + center_vector[0]
    yg = yprod + center_vector[1]
    zg = zprod + center_vector[2]

    return [xg, yg, zg]


# need to edit this
def plot_sequential_rotation(rotd_coords):
    count = 0
    for line in plane_f:
        if (len(line)//2 == 0):
            for i in range (0, len(line)):
                    plt.plot([rotd_coords[line[i]][0], rotd_coords[line[i+1]][0]], [rotd_coords[line[i]][1], rotd_coords[line[i+1]][1]], color='k')
        else:
            for i in range (0, len(line) - 1):
                    plt.plot([rotd_coords[line[i]][0], rotd_coords[line[i+1]][0]], [rotd_coords[line[i]][1], rotd_coords[line[i+1]][1]], color='k')
            plt.plot([rotd_coords[line[0]][0], rotd_coords[line[-1]][0]],
                     [rotd_coords[line[0]][1], rotd_coords[line[-1]][1]], color='k')

def plot_shape(vert_count, x_angle, y_angle, z_angle, center_vector, rot_seq):
    # GRID SETUP FOR SEQUENTIAL COMBINED ROTATIONS

    x_coords = plane_x
    y_coords = plane_y
    z_coords = plane_z

    # FIRST ROTATION
    # Takes values from the original lists of coord values (not point vectors!)
    # Passes them through each of the rotational functions and appends them into point vector matrices (list of lists)
    for vert_num in range(0, vert_count):
        '''
        # preserve x-rotated coords
        rotd_x_coords.append()
        # preserve y-rotated coords
        rotd_y_coords.append()
        # preserve z-rotated coords
        rotd_z_coords.append()
        '''
        rotd_x_coords.append(rotx(center_vector, (x_coords[vert_num], y_coords[vert_num], z_coords[vert_num]), x_angle))

        rotd_y_coords.append(roty(center_vector, (x_coords[vert_num], y_coords[vert_num], z_coords[vert_num]), y_angle))

        rotd_z_coords.append(rotz(center_vector, (x_coords[vert_num], y_coords[vert_num], z_coords[vert_num]), z_angle))

    # Rotate the selected one of six combination orders

    # SECOND AND THIRD ROTATIONS
    global rotd_sequence_coords
    if rot_seq == 'RxRyRz':
        rotd_sequence_coords = rotd_x_coords
        for vert_num in range(0, vert_count):
            rotd_sequence_coords[vert_num] = roty(center_vector, (
            rotd_sequence_coords[vert_num][0], rotd_sequence_coords[vert_num][1], rotd_sequence_coords[vert_num][2]),
                                                  y_angle)

        for vert_num in range(0, vert_count):
            rotd_sequence_coords[vert_num] = rotz(center_vector, (
            rotd_sequence_coords[vert_num][0], rotd_sequence_coords[vert_num][1], rotd_sequence_coords[vert_num][2]),
                                                  z_angle)

    elif rot_seq == 'RyRxRz':
        rotd_sequence_coords = rotd_y_coords
        for vert_num in range(0, vert_count):
            rotd_sequence_coords[vert_num] = rotx(center_vector, (
                rotd_sequence_coords[vert_num][0], rotd_sequence_coords[vert_num][1],
                rotd_sequence_coords[vert_num][2]),
                                                  x_angle)
        for vert_num in range(0, vert_count):
            rotd_sequence_coords[vert_num] = rotz(center_vector, (
                rotd_sequence_coords[vert_num][0], rotd_sequence_coords[vert_num][1],
                rotd_sequence_coords[vert_num][2]),
                                                  z_angle)

    elif rot_seq == 'RzRyRx':
        rotd_sequence_coords = rotd_z_coords
        for vert_num in range(0, vert_count):
            rotd_sequence_coords[vert_num] = roty(center_vector, (rotd_sequence_coords[vert_num][0],
                                                                  rotd_sequence_coords[vert_num][1],
                                                                  rotd_sequence_coords[vert_num][2]), y_angle)
        for vert_num in range(0, vert_count):
            rotd_sequence_coords[vert_num] = rotx(center_vector, (rotd_sequence_coords[vert_num][0],
                                                                  rotd_sequence_coords[vert_num][1],
                                                                  rotd_sequence_coords[vert_num][2]), x_angle)


    elif rot_seq == 'RzRxRy':
        rotd_sequence_coords = rotd_z_coords
        for vert_num in range(0, vert_count):
            rotd_sequence_coords[vert_num] = rotx(center_vector, (rotd_sequence_coords[vert_num][0],
                                                                  rotd_sequence_coords[vert_num][1],
                                                                  rotd_sequence_coords[vert_num][2]), x_angle)
        for vert_num in range(0, vert_count):
            rotd_sequence_coords[vert_num] = roty(center_vector, (rotd_sequence_coords[vert_num][0],
                                                                  rotd_sequence_coords[vert_num][1],
                                                                  rotd_sequence_coords[vert_num][2]), y_angle)


    elif rot_seq == 'RxRzRy':
        rotd_sequence_coords = rotd_x_coords
        for vert_num in range(0, vert_count):
            rotd_sequence_coords[vert_num] = rotz(center_vector, (rotd_sequence_coords[vert_num][0],
                                                                  rotd_sequence_coords[vert_num][1],
                                                                  rotd_sequence_coords[vert_num][2]), z_angle)
        for vert_num in range(0, vert_count):
            rotd_sequence_coords[vert_num] = roty(center_vector, (rotd_sequence_coords[vert_num][0],
                                                                  rotd_sequence_coords[vert_num][1],
                                                                  rotd_sequence_coords[vert_num][2]), y_angle)

    elif rot_seq == 'RyRzRx':
        rotd_sequence_coords = rotd_y_coords
        for vert_num in range(0, vert_count):
            rotd_sequence_coords[vert_num] = rotz(center_vector, (rotd_sequence_coords[vert_num][0],
                                                                  rotd_sequence_coords[vert_num][1],
                                                                  rotd_sequence_coords[vert_num][2]), z_angle)
        for vert_num in range(0, vert_count):
            rotd_sequence_coords[vert_num] = rotx(center_vector, (rotd_sequence_coords[vert_num][0],
                                                                  rotd_sequence_coords[vert_num][1],
                                                                  rotd_sequence_coords[vert_num][2]), x_angle)

    else:
        print('The rotation sequence is invalid.')
        sys.exit()

    plot_sequential_rotation(rotd_sequence_coords)


#change only these numbers to rotate the whole object
plane_center = [0, 0, 0]
rot_x = 20
rot_y = 60
rot_z = 90
seq_name = 'RyRzRx'

#here's where the rotation happens
plane_x, plane_y, plane_z, plane_f = [], [], [], []
objExport('plane2.obj')
plot_shape(1629, rot_x, rot_y, rot_z, plane_center, seq_name)
plane_x, plane_y, plane_z, plane_f = [], [], [], []
rotd_x_coords, rotd_y_coords, rotd_z_coords, rotd_y_coords = [], [], [], []
objExport('propeller.obj')
plot_shape(152,  rot_x, rot_y, rot_z, plane_center, seq_name)
plane_x, plane_y, plane_z, plane_f = [], [], [], []
rotd_x_coords, rotd_y_coords, rotd_z_coords, rotd_y_coords = [], [], [], []
objExport('wheel1.obj')
plot_shape(576, rot_x, rot_y, rot_z, plane_center, seq_name)
plane_x, plane_y, plane_z, plane_f = [], [], [], []
rotd_x_coords, rotd_y_coords, rotd_z_coords, rotd_y_coords = [], [], [], []
objExport('wheel2.obj')
plot_shape(576,  rot_x, rot_y, rot_z, plane_center, seq_name)
plane_x, plane_y, plane_z, plane_f = [], [], [], []
rotd_x_coords, rotd_y_coords, rotd_z_coords, rotd_y_coords = [], [], [], []
objExport('wheel3.obj')
plot_shape(576, rot_x, rot_y, rot_z, plane_center, seq_name)
plane_x, plane_y, plane_z, plane_f = [], [], [], []
rotd_x_coords, rotd_y_coords, rotd_z_coords, rotd_y_coords = [], [], [], []
objExport('wheel4.obj')
plot_shape(576, rot_x, rot_y, rot_z, plane_center, seq_name)
plt.show()
