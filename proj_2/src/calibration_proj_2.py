#!/usr/bin/env python3
# ENPM Project 2 Fall AY2021
# Script for calibration of Kuka Robot

import rospy
from std_msgs.msg import Float64
from numpy.linalg import inv
from sympy import *
from sympy import sin, cos, atan2, Matrix, Subs, N
from sympy import diff #to take partial derivative
import sympy as sp
from sympy import nsolve
import math
import numpy as np
from sympy.physics.mechanics import dynamicsymbols, Point, ReferenceFrame
from sympy.physics.vector import init_vprinting
init_vprinting(use_latex='mathjax', pretty_print=False)
from sympy.matrices.common import a2idx
init_printing(use_unicode=True)
import matplotlib as mpl
import matplotlib.pyplot as plt

d2r=np.deg2rad

#---------------------PART 1-----------------#

#---Forward Kinematics (HW#3)------#
theta1=symbols('theta1')
theta2=symbols('theta2')
theta3=symbols('theta3')
theta4=symbols('theta4')
theta5=symbols('theta5')
theta6=symbols('theta6') 
theta7=symbols('theta7')

#--Measurements from Figure 2 (in meters)--
d02=0.36 #joint 0 to 2
d01=d02/2 #180mm #joint 0 to 1
d12=d02/2 #joint 1 to 2
d24=0.42 #joint 2 to 4
d23=d24/2 #joint 2 to 3
d34=d24/2 #joint 3 to 4
d04=0.780 #joint 0 to 4
d47=0.505 #joint 4 to 7
d45=0.201 #joint 4 to 5
d56=0.1985 #joint 5 to 6
d46=d45+d56 #joint 4 to 6 #399.5 mm
d67=0.1055 #joint 6 to 7
d26=d24+d46
tool=0.1 #10 cm calibration tool

#---Homogenous Transformation Matrices with Joint 5 locked----
A1=Matrix([[cos(theta1),0,-sin(theta1),0],[sin(theta1),0,cos(theta1),0],[0,-1,0,d01],[0,0,0,1]])
A2=Matrix([[cos(theta2),0,sin(theta2),0],[sin(theta2),0,-cos(theta2),0],[0,1,0,d12],[0,0,0,1]])
A3=Matrix([[cos(theta3),0,-sin(theta3),0],[sin(theta3),0,cos(theta3),0],[0,-1,0,d23],[0,0,0,1]])
A4=Matrix([[cos(theta4),0,sin(theta4),0],[sin(theta4),0,-cos(theta4),0],[0,1,0,d34],[0,0,0,1]])
A5=Matrix([[cos(theta5),0,-sin(theta5),0],[sin(theta5),0,cos(theta5),0],[0,-1,0,d45],[0,0,0,1]])
A6=Matrix([[cos(theta6),0,sin(theta6),0],[sin(theta6),0,-cos(theta6),0],[0,1,0,d56],[0,0,0,1]])
A7=Matrix([[cos(theta6),-sin(theta7),0,0],[sin(theta7), cos(theta7),0,0],[0,0,1,(d67)],[0,0,0,1]])
Atool=Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,tool],[0,0,0,1]])

#---Final Transformation np.array (Base to Frame 7)---
H=A1*A2*A3*A4*A5*A6*A7

#---Final Transformation np.array (Base to Tool tip)---
H_tool=A1*A2*A3*A4*A5*A6*A7*Atool
# print(H.shape) #4x4
# print("\nForward Kinematics Transformation is: \n",np.array(H.evalf(2)))

#--End Effector (calibration tool tip) position---
#X is first 3 elements of 4th column (d)
# X = Matrix( H[[0,1,2],:][:,3] ) #no tool attached
X = Matrix( H_tool[[0,1,2],:][:,3] ) #d_x, d_y, and d_z
# d_x=H[0,3]
# d_y=H[1,3]
# d_z=H[2,3]
# print(X)

#----Solve for joint angles based on EE position---
# -------Using Inv Kinematics-----
#-------------starting EE at goal 1
goal_1=[-0.2, 0.7, 0.0]
goal_2=[0.2, 0.7, 0.0]
goal_3=[0.2, 0.4, 0.0]
goal_4=[-0.2, 0.4, 0.0]
d_x=goal_1[0]
d_y=goal_1[1]
d_z=goal_1[2]

t1=acos(d_x/d46)
b=tool+d67-d02
d=sqrt(d_y*d_y+b*b)
alpha=acos((d24*d24+d*d-d46*d46)/(2*d24*d))
beta2=acos((d46*d46+d*d-d24*d24)/(2*d46*d))
gamma2=acos((d24*d24+d46*d46-d*d)/(2*d24*d46))
phi=atan2(b,d_y)
t2=-(np.pi/2-(alpha+phi))
t3=0 #lock out joint 3
t4=-(np.pi-beta2)
t5=0 #lock out joint 5
t6=gamma2-phi-(np.pi/2)
t7=0 #no need to rotate joint 7; angle is irrelevant to completion of task

#-----Plot EE Position--------
plt.figure(800)
plt.axis([-0.5, 0.5, 0.3, 0.8])
plt.title("ENPM 662 Project#2\nJerry Pittman (117707120)\nMaitreya Kulkarni (UMDID)")
plt.xlabel('X-axis of Workspace (m)')
plt.ylabel('Y-axis of Workspace(m)')
plt.text(-0.2,0.7,"Goal 1")
plt.plot(-0.2,0.7,'bo')
plt.text(0.2,0.7,"Goal 2")
plt.plot(0.2,0.7,'bo')
plt.text(0.2,0.4,"Goal 3")
plt.plot(0.2,0.4,'bo')
plt.text(-0.2,0.4,"Goal 4")
plt.plot(-0.2,0.4,'bo')
plt.grid()

plt.figure(1)
plt.title("Joint 1 Angles")
plt.ylabel('Joint Angle (rad)')
plt.xlabel('time (sec)')
plt.grid()

plt.figure(2)
plt.title("Joint 2 Angles")
plt.ylabel('Joint Angle (rad)')
plt.xlabel('time (sec)')
plt.grid()

plt.figure(3)
plt.title("Joint 3 Angles")
plt.ylabel('Joint Angle (rad)')
plt.xlabel('time (sec)')
plt.grid()

plt.figure(4)
plt.title("Joint 4 Angles")
plt.ylabel('Joint Angle (rad)')
plt.xlabel('time (sec)')
plt.grid()

plt.figure(5)  #Joint 5
plt.title("Joint 5 Angles")
plt.ylabel('Joint Angle (rad)')
plt.xlabel('time (sec)')
plt.grid()

plt.figure(6)
plt.title("Joint 6 Angles")
plt.ylabel('Joint Angle (rad)')
plt.xlabel('time (sec)')
plt.grid()

plt.figure(7)
plt.title("Joint 7 Angles")
plt.ylabel('Joint Angle (rad)')
plt.xlabel('time (sec)')
plt.grid()

msg = """
Moving Kuka for Calibration of Tool!
---------------------------
(Positions are in mm)
Position pt is (x, y, z)
Position 1 is (-200, 600, 0)
Position 2 is (200, 600, 0)
Position 3 is (200, 400, 0)
Position 4 is (-200, 400, 0)
Allowed error is 10 mm

CTRL-C to quit
"""

rospy.init_node('kuka_calibrate', anonymous=True)
  
# PositionJointInterface_J1_controller
# type: effort_controllers/JointPositionController
# joint: iiwa_joint_1

#Create Publishers to Joint Controllers
pub_J1 = rospy.Publisher('/proj_2/PositionJointInterface_J1_controller/command', Float64, queue_size=10)     
pub_J2 = rospy.Publisher('/proj_2/PositionJointInterface_J2_controller/command', Float64, queue_size=10) 
pub_J3 = rospy.Publisher('/proj_2/PositionJointInterface_J3_controller/command', Float64, queue_size=10) 
pub_J4 = rospy.Publisher('/proj_2/PositionJointInterface_J4_controller/command', Float64, queue_size=10) 
pub_J5 = rospy.Publisher('/proj_2/PositionJointInterface_J5_controller/command', Float64, queue_size=10) 
pub_J6 = rospy.Publisher('/proj_2/PositionJointInterface_J6_controller/command', Float64, queue_size=10) 
pub_J7 = rospy.Publisher('/proj_2/PositionJointInterface_J7_controller/command', Float64, queue_size=10) 
    
time=0.0
# delta_t=0.5
error=0.01 #10 cm
cal=0.01 #calibrate every 10cm
i=2 #iterate calibration rows
    
def kuka_calibrate():
    goal_1=[-0.2, 0.7, 0.0]
    goal_2=[0.2, 0.7, 0.0]
    goal_3=[0.2, 0.4, 0.0]
    goal_4=[-0.2, 0.4, 0.0]
    d_x=goal_1[0]
    d_y=goal_1[1]
    d_z=goal_1[2]

    #while position of EE is not goal_3 or 4   
    while (d_x!=goal_4[0] and d_y!=goal_4[1]) or (d_x!=goal_3[0] and d_y!=goal_3[1]):
    #calibrating between goal 1 and 2
        if d_y==goal_1[1]: #700
            if d_x==goal_2[0]:
                d_y=goal_2[1]-cal #move down a row
                print("\nRow complete!\nNow calibrating to the right!")
                
            d_x=d_x+cal
        elif d_x!=goal_1[0] and d_y!=goal_2[1]-error: #calibrating between goal 2 and sub1 goal 1
            if d_x==goal_1[0] and d_y==goal_1[1]-error:
                d_y=goal_1[1]-cal #move down a row
                print("\nRow complete!\nNow calibrating to the left!")
                
            d_x=d_x-cal
        elif d_y==goal_1[1]-(i*cal): #calibrating between goal 1-(i*cal) and goal 2-(i*cal)
            if d_x==goal_2[0]:
                d_y=d_y-cal #move down a row (-y)
                print("\nRow complete!\nNow calibrating to the right!")
                i=i+1
                
            d_x=d_x+cal #move toward right (+x)
        elif d_x!=goal_1[0] and d_y!=goal_2[1]-(i*cal):#calibrating between goal 2-(i*cal) and  goal 1-i*cal
            if d_x==goal_1[0] and d_y==goal_1[1]-(i*cal):
                d_y=d_y-cal #move down a row (-y)
                print("\nRow complete!\nNow calibrating to the left!")
                i=i+1
                
            d_x=d_x-cal #move towards left (-x)
            
        #Graph EE Position
        plt.figure(800)
        plt.scatter(d_x, d_y)
    
        #Joint angles Should update automatically based on
        # change in d_x and d_y

        #Get Time from Ros
        seconds = rospy.get_time()
        
        #plot Joint Angles
        plt.figure(1) #Joint 1
        plt.scatter(seconds, t1)
        plt.figure(2) #Joint 2
        plt.scatter(seconds, t2)
        plt.figure(3)  #Joint 3
        plt.scatter(seconds, t3)
        plt.figure(4)  #Joint 4
        plt.scatter(seconds, t4)
        plt.figure(5)  #Joint 5
        plt.scatter(seconds, t5)
        plt.figure(6)  #Joint 6
        plt.scatter(seconds, t6)
        plt.figure(7)  #Joint 7
        plt.scatter(seconds, t7)

        #publish joint angle commands
        pub_J1.publish(t1)
        pub_J2.publish(t2) 
        pub_J3.publish(t3) 
        pub_J4.publish(t4) 
        pub_J5.publish(t5)
        pub_J6.publish(t6) 
        pub_J7.publish(t7) 
        
        if (d_x==goal_4[0] and d_y==goal_4[1]) or (d_x==goal_3[0] and d_y==goal_3[1]):
            print("Calibration Complete!")

if __name__=="__main__":
    try:
        print (msg)
        kuka_calibrate()

    except rospy.ROSInterruptException:
        pass

    finally:
        #publish joint angle commands
        pub_J1.publish(t1) 
        pub_J2.publish(t2)
        pub_J3.publish(t3) 
        pub_J4.publish(t4) 
        pub_J5.publish(t5) 
        pub_J6.publish(t6) 
        pub_J7.publish(t7) 
