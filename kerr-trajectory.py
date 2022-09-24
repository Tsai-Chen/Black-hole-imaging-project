# -*- coding: utf-8 -*-
"""
Created on Fri May  6 09:34:45 2022

@author: 90394
"""


from math import sin, cos, sqrt, pi
from numpy import linspace, arctan
import numpy as np
import sympy
import matplotlib.pyplot as plt
M = 2
r_0 = 10
a = 1
process = 0

theta_i_list = np.linspace(-pi/6, pi/6, 7)
theta_i_list = [pi/6, pi/12, 0, -pi/12, -pi/6]
phi_i_list = np.linspace(-pi/6, pi/6, 7)
phi_i_list = [pi/3, -pi/3]

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

for theta_i in theta_i_list:
    for phi_i in phi_i_list:

        E = sqrt(1 - 2*M/r_0)
        numerator1 = 2*a*M*(2*M-r_0)*sqrt(-(r_0*(a**2+r_0*(-2*M+r_0)))/(2*M-r_0))
        numerator2 = sqrt(1-2*M/r_0)*r_0**2*(a**2+r_0*(-2*M+r_0))*cos(theta_i)*sin(phi_i)
        denominator = sqrt(1-2*M/r_0)*r_0*(r_0-2*M)*sqrt(r_0*(r_0+a**2/(r_0-2*M)))
        L = (numerator1+numerator2)/denominator
        kappa = r_0**2*sin(theta_i)**2

        def Sigma(r, theta):
            return r**2 + a**2*(cos(theta))**2
        def Delta(r):
            return r**2 - 2*M*r + a**2
        def T(r):
            T = -a*L + (a**2 + r**2)*E
            return T
        def dtdlambda(r, theta):
            RHS = -a**2*E*(sin(theta))**2 + L*a + ((a**2 + r**2)*(E*(r**2 + a**2) - L*a))/Delta(r)
            return RHS/Sigma(r, theta)
        def dphidlambda(r, theta, flag):
            RHS = - flag*(a*E - L/(sin(theta))**2) + a*T(r)/Delta(r)
            return RHS/Sigma(r, theta)

        def dthetadlambda(r, theta, flag):
            temp = kappa - (cos(theta))**2*((a**2*E**2) + L**2/sin(theta)**2)
            try:
                RHS = flag*sqrt(temp)
            except:
                print(temp,"hi")
                print(sqrt(-1))
            return RHS/Sigma(r, theta)
        def drdlambda(r, theta, flag):
            temp = T(r)**2 - Delta(r)*((E*a - L)**2 + kappa)
            try:
                RHS = flag* sqrt(temp)
            except:
                print(r)
                print(temp,"hi")
                print(sqrt(-1))
            return RHS/Sigma(r, theta)

        t = [0]
        r = [r_0]
        theta = [pi/2]
        phi = [0]

        flag_r = -1
        flag_theta = -1
        flag_phi = 1
        if theta_i>0:
            flag_theta = 1
            
        n = 0
        h = 0.001
        ## 4th order runge kutta method
        while n < 400000:
        
            if Delta(r[-1]) < 1e-8:
                print('BH')
                break
            
            ## 0
            #r
            if abs(T(r[-1])**2 - Delta(r[-1])*((E*a - L)**2 + kappa)) < 1e-10 or (T(r[-1])**2 - Delta(r[-1])*((E*a - L)**2 + kappa) <0 and T(r[-2])**2 - Delta(r[-2])*((E*a - L)**2 + kappa) >0 ):
                flag_r*=-1
                r.pop()
                t.pop()
                phi.pop()
                theta.pop()
                continue
            else:
                r0 = h * drdlambda(r[-1], theta[-1], flag_r)
            #theta
            if abs(kappa - (cos(theta[-1]))**2*((a**2*E**2)+L**2/sin(theta[-1])**2)) < 1e-5 or kappa - (cos(theta[-1]))**2*((a**2*E**2)+L**2/sin(theta[-1])**2)<0:
                    if abs(kappa - (cos(theta[-1]))**2*((a**2*E**2)+L**2/sin(theta[-1])**2)) < 1e-28:
                        theta0=0
                    else:
                        flag_theta*=-1
                        r.pop()
                        t.pop()
                        phi.pop()
                        theta.pop()
                        continue
            else:
                theta0 = h * dthetadlambda(r[-1], theta[-1],flag_theta)
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
            #r
            if abs(T(r[-1]+r0/2)**2 - Delta(r[-1]+r0/2)*((E*a - L)**2 + kappa)) < 1e-10 or (T(r[-1]+r0/2)**2 - Delta(r[-1]+r0/2)*((E*a - L)**2 + kappa) <0 and T(r[-2])**2 - Delta(r[-2])*((E*a - L)**2 + kappa) >0 ):
                flag_r*=-1
                r.pop()
                t.pop()
                phi.pop()
                theta.pop()
                continue
            else:
                r1 = h * drdlambda(r[-1]+r0/2, theta[-1]+theta0/2, flag_r)
            #theta
            if abs(kappa - (cos(theta[-1]+theta0/2))**2*((a**2*E**2)+L**2/sin(theta[-1]+theta0/2)**2)) < 1e-5 or kappa - (cos(theta[-1]+theta0/2))**2*((a**2*E**2)+L**2/sin(theta[-1]+theta0/2)**2)<0:
                    if abs(kappa - (cos(theta[-1]+theta0/2))**2*((a**2*E**2)+L**2/sin(theta[-1]+theta0/2)**2)) < 1e-28:
                        theta1=0
                    else:
                        flag_theta*=-1
                        r.pop()
                        t.pop()
                        phi.pop()
                        theta.pop()
                        continue
            else:
                theta1 = h * dthetadlambda(r[-1]+r0/2, theta[-1]+theta0/2, flag_theta)
            t1 = h * dtdlambda(r[-1]+r0/2, theta[-1]+theta0/2)
            phi1 = h * dphidlambda(r[-1]+r0/2, theta[-1]+theta0/2, flag_phi)

            ## 2
            #r
            if abs(T(r[-1]+r1/2)**2 - Delta(r[-1]+r1/2)*((E*a - L)**2 + kappa)) < 1e-10 or (T(r[-1]+r1/2)**2 - Delta(r[-1]+r1/2)*((E*a - L)**2 + kappa) <0 and T(r[-2])**2 - Delta(r[-2])*((E*a - L)**2 + kappa) >0 ):
                flag_r*=-1
                r.pop()
                t.pop()
                phi.pop()
                theta.pop()
                continue
            else:
                r2 = h * drdlambda(r[-1]+r1/2, theta[-1]+theta1/2, flag_r)
            #theta
            if abs(kappa - (cos(theta[-1]+theta1/2))**2*((a**2*E**2)+L**2/sin(theta[-1]+theta1/2)**2)) < 1e-5 or kappa - (cos(theta[-1]+theta1/2))**2*((a**2*E**2)+L**2/sin(theta[-1]+theta1/2)**2)<0:
                    if abs(kappa - (cos(theta[-1]+theta1/2))**2*((a**2*E**2)+L**2/sin(theta[-1]+theta1/2)**2)) < 1e-28:
                        theta2=0
                    else:
                        flag_theta*=-1
                        r.pop()
                        t.pop()
                        phi.pop()
                        theta.pop()
                        continue
            else:
                theta2 = h * dthetadlambda(r[-1]+r1/2, theta[-1]+theta1/2, flag_theta)
            t2 = h * dtdlambda(r[-1]+r1/2, theta[-1]+theta1/2)
            phi2 = h * dphidlambda(r[-1]+r1/2, theta[-1]+theta1/2, flag_phi)

            ## 3
            #r
            if abs(T(r[-1]+r2)**2 - Delta(r[-1]+r2)*((E*a - L)**2 + kappa)) < 1e-10 or (T(r[-1]+r2)**2 - Delta(r[-1]+r2)*((E*a - L)**2 + kappa) <0 and T(r[-2])**2 - Delta(r[-2])*((E*a - L)**2 + kappa) >0 ):
                flag_r*=-1
                r.pop()
                t.pop()
                phi.pop()
                theta.pop()
                continue
            else:
                r3 = h * drdlambda(r[-1]+r1/2, theta[-1]+theta1/2, flag_r)
            #theta
            if abs(kappa - (cos(theta[-1]+theta2))**2*((a**2*E**2)+L**2/sin(theta[-1]+theta2)**2)) < 1e-5 or kappa - (cos(theta[-1]+theta2))**2*((a**2*E**2)+L**2/sin(theta[-1]+theta2)**2)<0:
                    if abs(kappa - (cos(theta[-1]+theta2))**2*((a**2*E**2)+L**2/sin(theta[-1]+theta2)**2)) < 1e-28:
                        theta3=0
                    else:
                        flag_theta*=-1
                        r.pop()
                        t.pop()
                        phi.pop()
                        theta.pop()
                        continue
            else:
                theta3 = h * dthetadlambda(r[-1]+r2, theta[-1]+theta2, flag_theta)
            t3 = h * dtdlambda(r[-1]+r2, theta[-1]+theta2)
            phi3 = h * dphidlambda(r[-1]+r2, theta[-1]+theta2, flag_phi)

            r.append(r[-1] + (r0+2*r1+2*r2+r3)/6)
            t.append(t[-1] + (t0+2*t1+2*t2+t3)/6)
            theta.append(theta[-1] + (theta0+theta1+theta2+theta3)/6)
            phi.append(phi[-1] + (phi0+phi1*2+phi2*2+phi3)/6)

            n += 1

        ## plot the trajectory
        x = []
        y = []
        z = []
        for j in range(len(t)-2):
            x.append(r[j]*sin(theta[j])*cos(phi[j]))
            y.append(r[j]*sin(theta[j])*sin(phi[j]))
            z.append(r[j]*cos(theta[j]))
        ax.plot(x, y, z, color = 'b', linewidth = 1) # light
        process += 1
        print('process: ' + str(process) + '/' + str(len(theta_i_list)*len(phi_i_list)))

ax.scatter(x[0], y[0], z[0], color = 'r') # observer       

radius = M + sqrt(M**2-a**2)

# plot the black hole
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = radius * np.outer(np.cos(u), np.sin(v))
y = radius * np.outer(np.sin(u), np.sin(v))
z = radius * np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_wireframe(x, y, z, rstride=10, cstride=10, color='red', linewidth = 1.5)

ax.set_title("Light Trajectory")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_xlim((-5,10))
ax.set_ylim((-5,10))
ax.set_zlim((-5,5))

plt.show()