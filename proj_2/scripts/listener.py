#!/usr/bin/env python3
#ENPM662 Fall AY2021
##Project 2 Script Subscriber

import rospy
from std_msgs.msg import Float64

def callback(data):
    rospy.loginfo("I heard %s", data.data)

def listener():

    rospy.init_node('listener', anonymous=True)
    
    rospy.Subscriber("/iiwa/PositionJointInterface_J1_controller/command", Float64, callback)
    rospy.Subscriber("/iiwa/PositionJointInterface_J2_controller/command", Float64, callback)
    rospy.Subscriber("/iiwa/PositionJointInterface_J3_controller/command", Float64, callback)
    rospy.Subscriber("/iiwa/PositionJointInterface_J4_controller/command", Float64, callback)
    rospy.Subscriber("/iiwa/PositionJointInterface_J5_controller/command", Float64, callback)
    rospy.Subscriber("/iiwa/PositionJointInterface_J6_controller/command", Float64, callback)
    rospy.Subscriber("/iiwa/PositionJointInterface_J7_controller/command", Float64, callback)

    # spin() simply keeps python from exiting until this node is stopped
    rospy.spin()

if __name__ == '__main__':
    listener()
