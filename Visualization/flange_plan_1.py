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
target_point = [600,0, 1100]
target_axis = [0, 0, 1]
initial_rotation = 0 
unstructured_path1 = []
plane = True
smoothing = 0

def convertUnstructured(angle_path):
  unstructured_path = [] 
  for i in range(0, len(angle_path)):
    for j in range(0, len(angle_path[i])):
      unstructured_path.append(asarray(angle_path[i][j]))
  return unstructured_path

def generate_targets(point, axis, rotation, sample_axis, tp):
  sample_axis = normalized(sample_axis)
  r = Rotation.from_angle_axis(rotation, axis)
  axis = normalized(asarray(axis))
  axis2 = unit_orthogonal_vector(axis);
  axis3 = normalized(cross(axis2, axis))
  p1 = [];
  p1.append(point + r*(axis2 * sqrt(25**2 + 25**2)))
  p1.append(point + r*(axis2 * -25 + axis3 * -25))
  p1.append(point + r*(axis2 * -25 + axis3 * 25))
 # print sample_axis
  angle2 = acos(clamp(dot(axis,sample_axis), -1, 1))
  #print angle2
  rotation_axis2 = cross(axis, sample_axis)
  #print rotation_axis2
  r2 = Rotation.from_angle_axis(angle2, rotation_axis2)

  #if flip:
 #   for i in range(0, len(p1)):
  orientation_vector = asarray(tp) - asarray(point)#p1[i]
  #  # angle = acos(clamp(orientation_vector, array([0, 0, 1]), -1, 1))  
  # r2 = Rotation.from_angle_axis(a2, normalized(cross(orientation_vector, axis))) 
  p1 = r2 * (asarray(p1) - asarray(tp)) + asarray(tp)
  return asarray(p1);

armset1 = System([Arm('a1',6, effector_functions)], generate_targets(target_point, target_axis, 0, array([0, 0, 1]), target_point), 0)
armset_dummy = System([Arm('a2',6,effector_functions)], generate_targets(target_point, target_axis, 0, array([0, 0, 1]), target_point),0)

def update_1():
  frame = props.get("frame")()
  state = unstructured_path1[int(frame)]
  for arm in armset1.arms:
    for i, node in enumerate(arm.nodes):
      node.rotation.set(float(state[i]) * 180/pi)

def setup_1(step1, step2): 
  global angle_path
  global location
  angle_goals = [];
  points = []
  goal_list = [];
  
  flip = False
  all_goals, rotation_axes, rotation_angles, target_axes, normals = armset1.pointNormals()
  this_target = []
  #.all_goals = getContours()
  targets = []
  for t in range(step1, step2):
    touch_point = all_goals[t]
    sample_axis = target_axes[t]
    counter = 0
    print sample_axis
    # if sample_axis[2] == -1:
    #   flip = True
    # # else: flip = False
    # sample_axis = array([0, -1, 0])

    #print touch_point
    coarse_samples = []
    counter = 0
    #sample_axis = array([0, 0, 1])
    # if sample_axis[2] == -1:
    #   target_point = array([700, 0, 800])
    # else:
    #   target_point = array([600, 0, 1300])


    targets.append(target_point - 100*sample_axis)
    this_target = target_point - 200*sample_axis
    goal_list.append(this_target)  

    for theta in linspace(0, 2*pi,11): 
      r = Rotation.from_angle_axis(theta, array([0, 0, 1]))
      p1 = this_target - r * touch_point
     # p1 = [0, 0, 0]
      triangle = generate_targets(p1, array([0, 0, 1]), theta, sample_axis, this_target)
      armset_dummy.set_target(triangle)
      points.append(triangle[0])
      points.append(triangle[1])
      points.append(triangle[2])
      goal, e = armset_dummy.solve_ik(counter)
      print magnitude(e)
      if magnitude(e) < .01:
        coarse_samples.append(goal)
    angle_goals.append(coarse_samples)
  #pp = this_target + sample_axis*100
  points.append(target_point)
  armset1.set_target(points)
  #Visualize where we're doing the initial sampling
  #scatter_points(points)
  return goal_list, angle_goals, normals, targets


def solve_1(step1, step2):
  global angle_path
  global unstructured_path1
  print "Solving 1"
  goal_list, angle_goals, normals, targets = setup_1(step1, step2) 

  #Get the initial offsets for the robot meshes
  robot_origins = asarray(get_origins())

  #Get the robot meshes into an array
  for arm in armset1.arms:
    robot_meshes = arm.submeshes()
 
  angle_path = plan_1(6, goal_list, robot_origins, robot_meshes, 
    [bike_mesh()], .1, 150 , 150., initial_rotation, array([0, 0, 0]), .01, angle_goals, smoothing, 100)

  unstructured_path1 = convertUnstructured(angle_path)

  return angle_path, normals, targets

def getAnglePath(step1, step2):
  return solve_1(step1, step2)


def view_1():
  # Set up main window
  solve_1(0,5)
  app = QEApp(sys.argv,True)
  main = MainWindow(props)
  ks = KukaScene(armset1)
  main.view.add_scene("kuka",ks)
  main.view.show_all(True)

  #Set up the actual targets for the RRT
  props.get("last_frame").set(len(unstructured_path1)) 
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
