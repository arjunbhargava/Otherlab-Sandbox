#!/usr/bin/env python
from __future__ import division
import sys
import math
from gui import *
from geode.value import parser, PropManager
from geode import *
from geode.openmesh import *
from geode.random import *
from OpenGL import GL, GLUT
from numpy import *
from props import props
from effectorik import *
from math import *
from scatter_points import *
from kinematic_chain import *
from plan_1 import *

#Data structures for paths

#These are the paths returned by the C++ path planning function. The organization is a vector<vector<vector>>>,
#with the outermost vector containing the entire path, the middle vector containing each submotion, and the 
#innermost vector having dimension of the number of DOF of the robots. 
angle_path1 = [] 
angle_path2 = []

#Unstructued paths are the solutions to the path planning problem, but organized in a way such that there's only 2 layers,
# a vector of containing 6 dimensional vectors. 
unstructured_path1 = []
unstructured_path2 = []


#Initialization of parameters
target_axis = array([0, 0, -1]) #for the second robot. 
initial_rotation = pi #So that they are facing each other. 
location = array([1200, -600, 0]) #Initial offset of the origin of the second robot. 

#step1 is the index of the first point we want to touch, step2 is the index of the last point. 
step1 = 0
step2 = 1

#Initialization of robots. Here we're initializing the second robot, since the first one gets done inside plan_1.py. 
arm2_frame = Frames(location,Rotation.from_angle_axis(initial_rotation,array([0,0,1])))
armset2 = System([Arm('a3',6, bf =arm2_frame)], generate_targets(target_point, array([0, 0, -1]), 0,), 2)
armset_dummy2 = System([Arm('a4',6, bf=arm2_frame)], generate_targets(target_point, array([0, 0, -1]), 0),2)


#The planning method returns a series of steps, where a full path is composed of several motions. Each of these motions 
#is converted to one very long motion in this method, such that the rendering can be done more easily.  
def adjustPaths(): 
  global unstructured_path1
  global unstructured_path2
  print "we doin it"
  
  for i, motion1 in enumerate(angle_path1):
    for j, motion2 in enumerate(angle_path2):
      if i == j:
       if len(motion1) > len(motion2):
          for k in range(0, len(motion1) - len(motion2)):
            motion2.append(motion2[-1])
          angle_path2[j] = motion2
       else:
          for k in range(0, len(motion2) - len(motion1)):
            motion1.append(motion1[-1])
          angle_path1[i] = motion1

  unstructured_path1 = convertUnstructured(angle_path1)
  unstructured_path2 = convertUnstructured(angle_path2)
  print len(unstructured_path1), len(unstructured_path2)


#The paths returned by the RRT are not guaranteed to be hte sam elegnth for each of the motions. This routine 
#picks the shorter of the 2 motions for each step and extends it until they are the same length. This is basically just 
#to get the timing of the motions correct. 
def timePaths():    
  global angle_path1 
  global angle_path2
  for i, motion1 in enumerate(angle_path1):
    for j, motion2 in enumerate(angle_path2):
      if i == j:
        temp_motion1 = []
        temp_motion2 = []
        pos1 = 0
        pos2 = 0
        flag = 0
        counter = 0
        while flag == 0:
          counter = counter + 1
          temp_motion1.append(motion1[pos1])
          temp_motion2.append(motion2[pos2])
          if pos1 < len(motion1)-1:
            pos1 = pos1 + 1
          if pos2 < len(motion2)-1:
            pos2 = pos2 + 1
          if temp_motion1[-1][6]==1 and temp_motion2[-1][6]==1:
            flag = 1
          if pos1 <= len(motion1) - 1 and pos2 <= len(motion2) -1:
            if motion1[pos1][6] > motion2[pos2][6]:
              pos1 = pos1 - 1
            elif motion2[pos2][6] > motion1[pos1][6]: 
              pos2 = pos2 - 1
        angle_path1[i] = temp_motion1
        angle_path2[i] = temp_motion2


#Render update here, uses the unstructured paths. 
def update():
  total_time = len(unstructured_path1)
  frame = props.get("frame")()
  state1 = unstructured_path1[int(frame)]
  state2 = unstructured_path2[int(frame)]
  print state1[6], state2[6];
  for arm in armset1.arms:
    for i, node in enumerate(arm.nodes):
        node.rotation.set(float(state1[i]) * 180/pi)
  for arm in armset2.arms:
    for j , node in enumerate(arm.nodes):
        node.rotation.set(float(state2[j]) * 180/pi)


#First we set up the problem, then we solve it. 


#Here I control the timesteps that I'm interested in solving. 
def setup_2(targets): 
  global step
  angle_goals = [];
  points = []
  goal_list = [];

  for t in range(step1, step2): 
    print t
    counter = 0
    #The targets are acquired from plan_1.py, so the second robot knows where to go...
    target_point = targets[t - step1] 
    armset_dummy2.set_bike_point(target_point)
    goal_list.append(target_point)
    coarse_samples = []

    #This is the point where you could add more possible target points, as in the commented out loop before. 
    #For this implementation I have the second robot only solving for a single solution. 
   
    #for i in linspace(-pi, pi, 11):
#      axis_r = Rotation.from_angle_axis(i, target_axes[t])
#      sample_axis = axis_r * initial_axis
      #if normals[t][2] > 0:
    #target_axis = -asarray(normals[t])
    #else: 
    #target_axis = array([0, 0, -1])

    target_axis = array([0, 0, -1]) #We set the target axis to point down (the robot flange will point along this vector)
    tp = target_point - (200 * target_axis) #Calculate the offset from the target point along the axis of interest. 
    triangle = generate_targets(tp, target_axis, 0) #Generate a target based on the point and axis. 
    points.append(triangle[0])
    points.append(triangle[1]) #Two points is sufficient to constrain the solution (a line..)
    armset_dummy2.set_target(triangle[0:2]) #Set the targets for the IK solver. 
    goal, e = armset_dummy2.solve_ik(counter) #e is the error from the IK. 
    #if magnitude(e) < .01: #Filter the errors here... Sometimes it doesn't converge to a solution, not sure why. 
    coarse_samples.append(append(goal, 1.))
    angle_goals.append(coarse_samples)

  #Visualize where we're doing the initial sampling
  points.append(target_point)
  armset2.set_target(points)
  armset2.set_bike_point(target_point)

 # scatter_points(points) #If you want to see where the goal points are.
  return goal_list, angle_goals


def solve_2(first_path, targets):
  print "Solving 2"
  global angle_path2
  global location
  global initial_rotation
  global unstructured_path2

  #Calculate all of the target points
  goal_list, angle_goals = setup_2(targets) 
  #Get the initial offsets for the robot meshes
  robot_origins = asarray(get_origins())

  #Get the robot meshes into an array
  for arm in armset2.arms:
    robot_meshes = arm.submeshes()

  #Parameter path refers to the path from the first robot, but organized in a way that the C++ code accepts it. 
  parameter_path = []
  for i, motion in enumerate(first_path):
    print motion
    solution = []
    for j, submotion in enumerate(motion):
      submotion.append(1/(len(motion)-1) * j)
      solution.append(asarray(submotion))
    parameter_path.append(asarray(solution))


#Solve the second robot here. 
  angle_path2 = plan_2(6, asarray(goal_list), robot_origins, robot_meshes, [bike_mesh()], .1, 150, 100., initial_rotation, location, .01, angle_goals,
   0, 100, asarray(parameter_path))


#This method sets up the environment and also calls plan_1.py through getAnglePath. Not too much going on here. 
def view():
  global step1
  global step2
  global angle_path1
  global angle_path2
  global unstructured_path1 
  global unstructured_path2
  # Set up main window, with both robots
  app = QEApp(sys.argv,True)
  main = MainWindow(props)
  ks = KukaScene(armset1)
  ks2 = KukaScene(armset2)
  main.view.add_scene("kuka",ks)
  main.view.add_scene("kuk2", ks2)
  main.view.show_all(True)
  #We compute the solution for the first robot
  angle_path1, targets = getAnglePath(step1, step2)
  #print angle_path1
  #Then given the axes and previous path, solve the second robot. 
  solve_2(angle_path1, targets)
  for i, motion in enumerate(angle_path2):
    print motion

  timePaths()
  adjustPaths()
  #More rendering stuff
  print angle_path2
  props.get("last_frame").set(len(unstructured_path2))
  main.resize_timeline(100);
  l = listen(props.get("frame"),update)
  main.init()
  app.run()


def main():
  #Pulls the properties for props
  parser.parse(props,'visualizer')
  view()
 
if __name__=='__main__':
  main()
