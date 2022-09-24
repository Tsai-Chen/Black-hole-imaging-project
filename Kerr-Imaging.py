# -*- coding: utf-8 -*-
"""
Created on Tue May  3 19:49:02 2022

@author: 90394
"""

from math import sin, cos, sqrt, pi, tan
import numpy as np
import matplotlib.pyplot as plt

# Initial Parameters
M = 2
r_0 = 40
a = 1

# The main function
# Building Connection between initial angles and final angles
def init_angle_to_final_angle(theta_i, phi_i):
    t = [0]
    r = [r_0]
    theta = [pi / 2]
    phi = [0]

    # Conserved Quantities
    E = sqrt(1 - 2 * M / r_0)
    numerator1 = 2 * a * M * (2 * M - r_0) * sqrt(-(r_0 * (a ** 2 + r_0 * (-2 * M + r_0))) / (2 * M - r_0))
    numerator2 = sqrt(1 - 2 * M / r_0) * r_0 ** 2 * (a ** 2 + r_0 * (-2 * M + r_0)) * cos(theta_i) * sin(phi_i)
    denominator = sqrt(1 - 2 * M / r_0) * r_0 * (r_0 - 2 * M) * sqrt(r_0 * (r_0 + a ** 2 / (r_0 - 2 * M)))
    L = (numerator1 + numerator2) / denominator
    kappa = r_0 ** 2 * sin(theta_i) ** 2

    def Sigma(r, theta):
        return r ** 2 + a ** 2 * (cos(theta)) ** 2

    def Delta(r):
        return r ** 2 - 2 * M * r + a ** 2

    def T(r):
        T = -a * L + (a ** 2 + r ** 2) * E
        return T
# Four equations
    def dtdlambda(r, theta):
        RHS = -a ** 2 * E * (sin(theta)) ** 2 + L * a + ((a ** 2 + r ** 2) * (E * (r ** 2 + a ** 2) - L * a)) / Delta(r)
        return RHS / Sigma(r, theta)

    def dphidlambda(r, theta, flag):
        RHS = - flag * (a * E - L / (sin(theta)) ** 2) + a * T(r) / Delta(r)
        return RHS / Sigma(r, theta)

    def dthetadlambda(r, theta, flag):
        temp = kappa - (cos(theta)) ** 2 * ((a ** 2 * E ** 2) + L ** 2 / sin(theta) ** 2)
        try:
            RHS = flag * sqrt(temp)
        except:
            print(temp, "hi")
            print(sqrt(-1))
        return RHS / Sigma(r, theta)

    def drdlambda(r, theta, flag):
        temp = T(r) ** 2 - Delta(r) * ((E * a - L) ** 2 + kappa)
        try:
            RHS = flag * sqrt(temp)
        except:
            print(r)
            print(temp, "hi")
            print(sqrt(-1))
        return RHS / Sigma(r, theta)

    flag_r = -1
    flag_theta = -1
    flag_phi = 1
    if theta_i > 0:
        flag_theta = 1

    n = 0
    h = 0.04
    flag_stop = 0
    
    ## 4th order runge kutta method
    while True:

        if Delta(r[-1]) < 1e-5 or phi[-1] >= 4 * pi:
            flag_stop = 1
            break

        ## 0
        # r
        if abs(T(r[-1]) ** 2 - Delta(r[-1]) * ((E * a - L) ** 2 + kappa)) < 1e-10 or (
                T(r[-1]) ** 2 - Delta(r[-1]) * ((E * a - L) ** 2 + kappa) < 0 and T(r[-2]) ** 2 - Delta(r[-2]) * (
                (E * a - L) ** 2 + kappa) > 0):
            flag_r *= -1
            r.pop()
            t.pop()
            phi.pop()
            theta.pop()
            continue
        else:
            r0 = h * drdlambda(r[-1], theta[-1], flag_r)
        # theta
        if abs(kappa - (cos(theta[-1])) ** 2 * ((a ** 2 * E ** 2) + L ** 2 / sin(theta[-1]) ** 2)) < 1e-5 or kappa - (
        cos(theta[-1])) ** 2 * ((a ** 2 * E ** 2) + L ** 2 / sin(theta[-1]) ** 2) < 0:
            if abs(kappa - (cos(theta[-1])) ** 2 * ((a ** 2 * E ** 2) + L ** 2 / sin(theta[-1]) ** 2)) < 1e-28:
                theta0 = 0
            else:
                flag_theta *= -1
                r.pop()
                t.pop()
                phi.pop()
                theta.pop()
                continue
        else:
            theta0 = h * dthetadlambda(r[-1], theta[-1], flag_theta)
        t0 = h * dtdlambda(r[-1], theta[-1])
        phi0 = h * dphidlambda(r[-1], theta[-1], flag_phi)

        # r1 = h * drdlambda(r[-1]+r0/2, theta[-1]+theta0/2, flag_r)
        # theta1 = h * dthetadlambda(r[-1]+r0/2, theta[-1]+theta0/2,flag_theta)
        # t1 = h * dtdlambda(r[-1]+r0/2, theta[-1]+theta0/2)
        # phi1 = h * dphidlambda(r[-1]+r0/2, theta[-1]+theta0/2, flag_phi)
        # r2 = h * drdlambda(r[-1]+r1/2, theta[-1]+theta1/2, flag_r)
        # theta2 = h * dthetadlambda(r[-1]+r1/2, theta[-1]+theta1/2,flag_theta)
        # t2 = h * dtdlambda(r[-1]+r1/2, theta[-1]+theta1/2)
        # phi2 = h * dphidlambda(r[-1]+r1/2, theta[-1]+theta1/2, flag_phi)
        # r3 = h * drdlambda(r[-1]+r2, theta[-1]+theta2, flag_r)
        # theta3 = h * dthetadlambda(r[-1]+r2, theta[-1]+theta2,flag_theta)
        # t3 = h * dtdlambda(r[-1]+r2, theta[-1]+theta2)
        # phi3 = h * dphidlambda(r[-1]+r2, theta[-1]+theta2, flag_phi)
        ## 1
        # r
        if abs(T(r[-1] + r0 / 2) ** 2 - Delta(r[-1] + r0 / 2) * ((E * a - L) ** 2 + kappa)) < 1e-10 or (
                T(r[-1] + r0 / 2) ** 2 - Delta(r[-1] + r0 / 2) * ((E * a - L) ** 2 + kappa) < 0 and T(
                r[-2]) ** 2 - Delta(r[-2]) * ((E * a - L) ** 2 + kappa) > 0):
            flag_r *= -1
            r.pop()
            t.pop()
            phi.pop()
            theta.pop()
            continue
        else:
            r1 = h * drdlambda(r[-1] + r0 / 2, theta[-1] + theta0 / 2, flag_r)
        # theta
        if abs(kappa - (cos(theta[-1] + theta0 / 2)) ** 2 * (
                (a ** 2 * E ** 2) + L ** 2 / sin(theta[-1] + theta0 / 2) ** 2)) < 1e-5 or kappa - (
        cos(theta[-1] + theta0 / 2)) ** 2 * ((a ** 2 * E ** 2) + L ** 2 / sin(theta[-1] + theta0 / 2) ** 2) < 0:
            if abs(kappa - (cos(theta[-1] + theta0 / 2)) ** 2 * (
                    (a ** 2 * E ** 2) + L ** 2 / sin(theta[-1] + theta0 / 2) ** 2)) < 1e-28:
                theta1 = 0
            else:
                flag_theta *= -1
                r.pop()
                t.pop()
                phi.pop()
                theta.pop()
                continue
        else:
            theta1 = h * dthetadlambda(r[-1] + r0 / 2, theta[-1] + theta0 / 2, flag_theta)
        t1 = h * dtdlambda(r[-1] + r0 / 2, theta[-1] + theta0 / 2)
        phi1 = h * dphidlambda(r[-1] + r0 / 2, theta[-1] + theta0 / 2, flag_phi)

        ## 2
        # r
        if abs(T(r[-1] + r1 / 2) ** 2 - Delta(r[-1] + r1 / 2) * ((E * a - L) ** 2 + kappa)) < 1e-10 or (
                T(r[-1] + r1 / 2) ** 2 - Delta(r[-1] + r1 / 2) * ((E * a - L) ** 2 + kappa) < 0 and T(
                r[-2]) ** 2 - Delta(r[-2]) * ((E * a - L) ** 2 + kappa) > 0):
            flag_r *= -1
            r.pop()
            t.pop()
            phi.pop()
            theta.pop()
            continue
        else:
            r2 = h * drdlambda(r[-1] + r1 / 2, theta[-1] + theta1 / 2, flag_r)
        # theta
        if abs(kappa - (cos(theta[-1] + theta1 / 2)) ** 2 * (
                (a ** 2 * E ** 2) + L ** 2 / sin(theta[-1] + theta1 / 2) ** 2)) < 1e-5 or kappa - (
        cos(theta[-1] + theta1 / 2)) ** 2 * ((a ** 2 * E ** 2) + L ** 2 / sin(theta[-1] + theta1 / 2) ** 2) < 0:
            if abs(kappa - (cos(theta[-1] + theta1 / 2)) ** 2 * (
                    (a ** 2 * E ** 2) + L ** 2 / sin(theta[-1] + theta1 / 2) ** 2)) < 1e-28:
                theta2 = 0
            else:
                flag_theta *= -1
                r.pop()
                t.pop()
                phi.pop()
                theta.pop()
                continue
        else:
            theta2 = h * dthetadlambda(r[-1] + r1 / 2, theta[-1] + theta1 / 2, flag_theta)
        t2 = h * dtdlambda(r[-1] + r1 / 2, theta[-1] + theta1 / 2)
        phi2 = h * dphidlambda(r[-1] + r1 / 2, theta[-1] + theta1 / 2, flag_phi)

        ## 3
        # r
        if abs(T(r[-1] + r2) ** 2 - Delta(r[-1] + r2) * ((E * a - L) ** 2 + kappa)) < 1e-10 or (
                T(r[-1] + r2) ** 2 - Delta(r[-1] + r2) * ((E * a - L) ** 2 + kappa) < 0 and T(r[-2]) ** 2 - Delta(
                r[-2]) * ((E * a - L) ** 2 + kappa) > 0):
            flag_r *= -1
            r.pop()
            t.pop()
            phi.pop()
            theta.pop()
            continue
        else:
            r3 = h * drdlambda(r[-1] + r1 / 2, theta[-1] + theta1 / 2, flag_r)
        # theta
        if abs(kappa - (cos(theta[-1] + theta2)) ** 2 * (
                (a ** 2 * E ** 2) + L ** 2 / sin(theta[-1] + theta2) ** 2)) < 1e-5 or kappa - (
        cos(theta[-1] + theta2)) ** 2 * ((a ** 2 * E ** 2) + L ** 2 / sin(theta[-1] + theta2) ** 2) < 0:
            if abs(kappa - (cos(theta[-1] + theta2)) ** 2 * (
                    (a ** 2 * E ** 2) + L ** 2 / sin(theta[-1] + theta2) ** 2)) < 1e-28:
                theta3 = 0
            else:
                flag_theta *= -1
                r.pop()
                t.pop()
                phi.pop()
                theta.pop()
                continue
        else:
            theta3 = h * dthetadlambda(r[-1] + r2, theta[-1] + theta2, flag_theta)
        t3 = h * dtdlambda(r[-1] + r2, theta[-1] + theta2)
        phi3 = h * dphidlambda(r[-1] + r2, theta[-1] + theta2, flag_phi)

        r.append(r[-1] + (r0 + 2 * r1 + 2 * r2 + r3) / 6)
        t.append(t[-1] + (t0 + 2 * t1 + 2 * t2 + t3) / 6)
        theta.append(theta[-1] + (theta0 + theta1 + theta2 + theta3) / 6)
        phi.append(phi[-1] + (phi0 + phi1 * 2 + phi2 * 2 + phi3) / 6)

        n += 1

        if r[-1] >= 500:
            break
    if flag_stop == 1:
        return 0, 0, 1
        # 1 = fallen into BH

    if flag_stop == 0:
        return theta[-1] - pi / 2, pi - phi[-1], 0

# Parameters for image processing
Rs = 40
grid = 256

# Processing Image
img = plt.imread('berkeley.png')
h = len(img)
w = len(img[1])
print(h,w)
def angle_to_pixel(theta, phi, h, w):
    return Rs * tan(phi) + w / 2, h / 2 - (Rs / cos(phi) * tan(theta))

# 'Grid'
angular_grid_phi = list(np.linspace(-pi/3, pi/3, grid))
angular_grid_theta = list(np.linspace(-pi/3, pi/3, grid))
pixel_grid_to_angular_grid = []
for i in range(grid):
    row = []
    for j in range(grid):
        row.append((float(angular_grid_theta[i]),float(angular_grid_phi[j])))
    pixel_grid_to_angular_grid.append(row)

processed_picture = []
count = 0
for i in range(grid):
    colomn = []
    for j in range(grid):
        theta_i, phi_i = pixel_grid_to_angular_grid[i][j]
        theta_f, phi_f, flag_BH = init_angle_to_final_angle(theta_i, phi_i)
#        theta_f, phi_f, flag_BH =  theta_i, phi_i,0 
        pixel_x_pre, pixel_y_pre = angle_to_pixel(theta_f, phi_f, h, w)
        pixel_x = int(pixel_x_pre)
        pixel_y = int(pixel_y_pre)
        if pixel_x >= 0 and pixel_x < w and pixel_y >= 0 and pixel_y < h and flag_BH == 0:
            color = [img[pixel_y][pixel_x][0],img[pixel_y][pixel_x][1],img[pixel_y][pixel_x][2]]
        else:
            if flag_BH == 1:
                color=[0,0,0]
            else:
                color = [0,100,0]
        colomn.append(color)

        count += 1
        if count % grid ==0:
            print('count = {} / {}'.format(count, grid**2))
    processed_picture.insert(0,colomn)

plt.imshow(processed_picture)
plt.show()