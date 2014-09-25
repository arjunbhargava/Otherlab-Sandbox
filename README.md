Otherlab-Sandbox
================

This folder contains path planning code for Kuka robot as they weld bicycle frames. There are several hidden files which contain earlier iterations of the code. 

The folders are organized as follows: 
	
Otherlab/ 
	sandbox/
		kinematic_chain/
			This folder contains the c++ code with python wrappers for the simulation code. 
			It is possible to plan the path for a single robot or both of the robots. 
			effector_robot_1.cpp contains the most recent version of the single planner. 
			robot_2_planner.cpp contains the most recent version of the planner for the second robot (it's dependent on the first robot being planned, so don't try to run this code on its own.) 
			The other important files in this folder are RobotSystem.cpp and RobotSystem.h. This is a class that contains the data for each of the robot systems, such a the mesh data, locations, and other useful metadata. 
			Fully detailed descriptions of the code are in the folders themselves. 

		Visualization/ 
			The contents of this folder can be sorted roughly into 2 kinds of files - general rendering/handling files, and the ik files. 
			Important files: 
				effectorik.py - The inverse kinematics which treat a point on the bike frame as the effector. 
				plan_1.py - Plans and renders the first robot
				timed_plan.py - Calls plan_1.py and then computes a new path for the second robot based on the path of the first robot. 			 
		
