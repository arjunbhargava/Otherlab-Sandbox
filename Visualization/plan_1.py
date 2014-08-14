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

angle_path = [];
target_point = [1200,0, 1400]
target_axis = [0, 0, 1]
initial_rotation = 0 
unstructured_path = []

def generate_targets(point, axis, rotation):
  r = Rotation.from_angle_axis(rotation, axis)
  axis = normalized(asarray(axis))
  axis2 = unit_orthogonal_vector(axis);
  axis3 = normalized(cross(axis2, axis))
  p1 = [];
  p1.append(point + r*(axis2 * sqrt(25**2 + 25**2)))
  p1.append(point + r*(axis2 * -25 + axis3 * -25))
  p1.append(point + r*(axis2 * -25 + axis3 * 25))
  return asarray(p1);

armset1 = System([Arm('a1',6, effector_functions)], generate_targets(target_point, target_axis, 0), 0)
armset_dummy = System([Arm('a2',6,effector_functions)], generate_targets(target_point, target_axis, 0),0)

def update_1():
  frame = props.get("frame")()
  state = unstructured_path[int(frame)]
  for arm in armset1.arms:
    for i, node in enumerate(arm.nodes):
      node.rotation.set(float(state[i]) * 180/pi)

def setup_1(): 
  global angle_path
  global location
  angle_goals = [];
  points = []
  goal_list = [];
  
  all_goals, rotation_axes, rotation_angles, target_axes, normals = armset1.pointNormals()
  target_axes = []
  all_goals = getContours()
  for t in range(0, 1): #Check the logic here...
    print t
#    if rand() < 1:
    sample_axis = array([0, 0, 1])
#    else: sample_axis = array([0, 0, -1])

    target_axes.append(sample_axis)
    counter = 0
    touch_point = all_goals[t]
    print touch_point
    goal_list.append(target_point)  
    coarse_samples = []
    counter = 0
    for theta in linspace(-pi, pi, 13): 
      r = Rotation.from_angle_axis(theta, sample_axis)
      p1 = target_point - r * touch_point
      triangle = generate_targets(p1, sample_axis, theta)
      # points.append(p1)
      points.append(triangle[0])
      points.append(triangle[1])
      points.append(triangle[2])
      armset_dummy.set_target(triangle)
      goal, timeout = armset_dummy.solve_ik(counter)
      coarse_samples.append(goal)
    angle_goals.append(coarse_samples)
  pp = target_point + array([0,0,100])
  points.append(pp)
  armset1.set_target(points)
  #Visualize where we're doing the initial sampling
#  scatter_points(points)

  return goal_list, angle_goals, target_axes


def solve_1():
  global angle_path
  global unstructured_path
  print "Solving 1"
  goal_list, angle_goals, target_axes = setup_1() 
  #Get the initial offsets for the robot meshes
  robot_origins = asarray(get_origins())

  #Get the robot meshes into an array
  for arm in armset1.arms:
    robot_meshes = arm.submeshes()
 
  angle_path = plan_1(6, goal_list, robot_origins, robot_meshes, 
    [bike_mesh()], .1, 50 , 10., initial_rotation, array([0, 0, 0]), .01, angle_goals, 0, 100)

  for i in range(0, len(angle_path)):
    for j in range(0, len(angle_path[i])):
      unstructured_path.append(asarray(angle_path[i][j]))
  return angle_path, target_axes

def getAnglePath():
  return solve_1()


def view_1():
  # Set up main window
  solve_1()
  app = QEApp(sys.argv,True)
  main = MainWindow(props)
  ks = KukaScene(armset1)
  main.view.add_scene("kuka",ks)
  main.view.show_all(True)

  #Set up the actual targets for the RRT
  props.get("last_frame").set(len(unstructured_path)) 
  main.resize_timeline(100);
  l = listen(props.get("frame"),update_1)
  main.init()
  app.run()


def main():
  #Pulls the properties for props
  parser.parse(props,'visualizer')
  view_1()
 
if __name__=='__main__':
  main()
