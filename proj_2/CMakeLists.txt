cmake_minimum_required(VERSION 2.8.3)

project(proj_2)

find_package(catkin REQUIRED COMPONENTS
  iiwa_ros
  iiwa_msgs
  cmake_modules
  roscpp
  rospy
  roslaunch
)

catkin_package(
	INCLUDE_DIRS 
	LIBRARIES	 ${PROJECT_NAME}
	CATKIN_DEPENDS iiwa_msgs iiwa_ros roscpp rospy
)

find_package(roslaunch)

foreach(dir config launch meshes urdf)
	install(DIRECTORY ${dir}/
		DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}/${dir})
endforeach(dir)

