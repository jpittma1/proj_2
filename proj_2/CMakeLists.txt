cmake_minimum_required(VERSION 3.0.2)

project(proj_2)

find_package(catkin REQUIRED)

catkin_package()

find_package(roslaunch)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  std_msgs
)

# catkin_python_setup()

# catkin_install_python(PROGRAMS src/calibration.py
#   DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
# )

catkin_install_python(PROGRAMS  scripts/listener.py src/calibration_proj2.py
  DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)

foreach(dir config launch meshes urdf)
	install(DIRECTORY ${dir}/
		DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}/${dir})
endforeach(dir)
