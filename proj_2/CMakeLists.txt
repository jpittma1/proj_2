cmake_minimum_required(VERSION 2.8.3)

project(proj_2)

find_package(catkin REQUIRED COMPONENTS
  cmake_modules
  std_msgs
  geometry_msgs
  roscpp
  rospy
  roslaunch
)

catkin_package(
	INCLUDE_DIRS include
	LIBRARIES	 ${PROJECT_NAME}
	CATKIN_DEPENDS std_msgs geometry_msgs roscpp rospy
)

find_package(roslaunch)

foreach(dir config launch meshes urdf)
	install(DIRECTORY ${dir}/
		DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}/${dir})
endforeach(dir)

install(
  PROGRAMS
    src/calibration_proj2.py
    scripts/listener.py
  DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)
