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

angle_path2 = [];
target_point = [1200,0, 1400]
# target_axis = [0, 0, 1]
initial_rotation = pi
armset2 = System([Arm('a3',6, effector_functions)], generate_targets(target_point, target_axis, 0), 1)
armset_dummy2 = System([Arm('a4',6,effector_functions)], generate_targets(target_point, target_axis, 0),1)

# def update():
#   frame = props.get("frame")()
#   state1 = angle_path1[int(frame)]
#   state2 = angle_path2[int(frame)]
#   for arm in armset1.arms:
#     for i, node in enumerate(arm.nodes):
#       node.rotation.set(float(state[i]) * 180/pi)
#   for arm in obstacle_arm.arms:
#     for j , node in enumerate(arm.nodes):
#       node.rotation.set(float(state[j+6] * 180/pi))

def setup_2(target_axes): 
  global angle_path
  global location
  angle_goals = [];
  points = []
  goal_list = [];

  for t in range(0, 1): #Check the logic here...
    print t
    for cone_angle in linspace(0, pi/24, 3):
      initial_axis = [target_axes[t][2] * sin(cone_angle), 0, cos(cone_angle)]
      counter = 0
      goal_list.append(target_point)
      coarse_samples = []
      for i in linspace(-pi, pi, 10):
        axis_r = Rotation.from_angle_axis(i, target_axis)
        sample_axis = axis_r * initial_axis
        triangle = generate_targets(target_point, sample_axis, i)
        print triangle
        points.append(triangle[0])
        points.append(triangle[1])
        points.append(triangle[2])
        armset_dummy2.set_target(triangle)
        goal, timeout = armset_dummy2.solve_ik(counter)
        coarse_samples.append(goal)
      angle_goals.append(coarse_samples)

  #Visualize where we're doing the initial sampling
  pp = target_point + array([0,0,100])
  points.append(pp)
#  scatter_points(points)
  return goal_list, angle_goals


def solve_2(first_path, target_axes):
  global angle_path2
  goal_list, angle_goals = setup_2(target_axes) 
  #Get the initial offsets for the robot meshes
  robot_origins = asarray(get_origins())

  #Get the robot meshes into an array
  for arm in armset1.arms:
    robot_meshes = arm.submeshes()
 
   angle_path2 = plan_2(6, goal_list, robot_origins, robot_meshes, 
     [bike_mesh()], .1, 50 , 10., initial_rotation, array([0, 0, 0]), .01, angle_goals, 1, 100, first_path)


def view():
  # Set up main window
  app = QEApp(sys.argv,True)
  main = MainWindow(props)
  ks = KukaScene(armset1)
  ks2 = KukaScene(armset2)
  main.view.add_scene("kuka",ks)
  main.view.add_scene("kuk2", ks2)
  main.view.show_all(True)
  first_path, target_axes = getAnglePath()
  solve_2(first_path, target_axes)
  #Set up the actual targets for the RRT
 # solve()
  #props.get("last_frame").set(len(angle_path)) 
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
