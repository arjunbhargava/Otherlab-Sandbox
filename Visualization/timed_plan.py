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
from follow import *
from math import *
from scatter_points import *
from kinematic_chain import *
from plan_1 import *

#Data structures for paths
angle_path1 = []
angle_path2 = [];
unstructured_path1 = []
unstructured_path2 = []
frame1 = 0
frame2 = 0
#Initialization of parameters
target_axis = [0, 0, -1]
initial_rotation = pi
location = array([1600, -600, 0])
step1 = 0
step2 = 1

#Initialization of robots
arm2_frame = Frames(location,Rotation.from_angle_axis(initial_rotation,array([0,0,1])))
armset2 = System([Arm('a3',6, effector_functions, bf =arm2_frame)], generate_targets(target_point, target_axis, 0, False), 1)
armset_dummy2 = System([Arm('a4',6,effector_functions, bf=arm2_frame)], generate_targets(target_point, target_axis, 0, False),1)



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

def setup_2(target_axes, target_point): 
  global step
  angle_goals = [];
  points = []
  goal_list = [];

  for t in range(step1, step2): #Check the logic here...
    print t
    counter = 0
 #   target_point = target_point + 100 * target_axes[t]
    target_point =  array([700, 0, 1500])
    goal_list.append(target_point)
    coarse_samples = []
    for i in linspace(-pi, pi, 13):
#      axis_r = Rotation.from_angle_axis(i, target_axes[t])
#      sample_axis = axis_r * initial_axis
      triangle = generate_targets(target_point, target_axis, i, False)
      points.append(triangle[0])
      points.append(triangle[1])
      points.append(triangle[2])
      armset_dummy2.set_target(triangle)
      goal, timeout = armset_dummy2.solve_ik(counter)
      coarse_samples.append(append(goal, 1.))
    angle_goals.append(coarse_samples)

  #Visualize where we're doing the initial sampling
  pp = target_point - array([0,0,100])
  points.append(pp)
 # scatter_points(points)
  return goal_list, angle_goals

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

def solve_2(first_path, target_axes):
  print "Solving 2"
  global angle_path2
  global location
  global initial_rotation
  global unstructured_path2

  #Calculate all of the target points
  goal_list, angle_goals = setup_2(target_axes, target_point) 
  #Get the initial offsets for the robot meshes
  robot_origins = asarray(get_origins())

  #Get the robot meshes into an array
  for arm in armset2.arms:
    robot_meshes = arm.submeshes()

  parameter_path = []
  for i, motion in enumerate(first_path):
    print motion
    solution = []
    for j, submotion in enumerate(motion):
      submotion.append(1/(len(motion)-1) * j)
      solution.append(asarray(submotion))
    parameter_path.append(asarray(solution))

  print parameter_path[0]
  angle_path2 = plan_2(6, asarray(goal_list), robot_origins, robot_meshes, [bike_mesh()], .1, 150, 100., initial_rotation, location, .01, angle_goals,
   0, 100, asarray(parameter_path))

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
  angle_path1, target_axes, target_point = getAnglePath(step1, step2)
  #print angle_path1
  #Then given the axes and previous path, solve the second robot. 
  solve_2(angle_path1, target_axes)
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
