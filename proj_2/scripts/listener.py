#!/usr/bin/env python3

##Project 2 Script Subscriber

import rospy
from std_msgs.msg import Float64

def callback(data):
    rospy.loginfo("I heard %s", data.data)

def listener():

    rospy.init_node('listener', anonymous=True)
    # pub_J1 = rospy.Publisher('/proj_2/PositionJointInterface_J1_controller/command', Float64, queue_size=10) 
    rospy.Subscriber("/proj_2/PositionJointInterface_J1_controller/command", Float64, callback)
    rospy.Subscriber("/proj_2/PositionJointInterface_J2_controller/command", Float64, callback)
    rospy.Subscriber("/proj_2/PositionJointInterface_J3_controller/command", Float64, callback)
    rospy.Subscriber("/proj_2/PositionJointInterface_J4_controller/command", Float64, callback)
    rospy.Subscriber("/proj_2/PositionJointInterface_J6_controller/command", Float64, callback)
    rospy.Subscriber("/proj_2/PositionJointInterface_J7_controller/command", Float64, callback)

    # spin() simply keeps python from exiting until this node is stopped
    rospy.spin()

if __name__ == '__main__':
    listener()
