#!/usr/bin/env python
from __future__ import division
import sys
import math
import random
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



#Initializing parameters for planning the path of the first robot. 
angle_path = [];
target_point = [700,0, 1100] #Point in space where we want the middle of the effector to be. 
target_axis = [0, 0, 1] #Axis laong which we want the normal of the ffector. 
initial_rotation = 0 
unstructured_path1 = [] 
#Flags for whether or not we want smoothing of the paths. For the algorithm to run fast you don't want smoothing on
smoothing = 0. 
#Specify which robot we're working with for the IK solver. 
robot_number = 1


#Converts the path returned from the path planning (a vector<vector<vector>>>) into a single level array whose 
#elements are a 6D array containing the angles for each joint. 
def convertUnstructured(angle_path):
  unstructured_path = [] 
  for i in range(0, len(angle_path)):
    for j in range(0, len(angle_path[i])):
      unstructured_path.append(asarray(angle_path[i][j]))
  return unstructured_path


#Generates 3 points around the point in space where want the effector to be. These targets are the point itself, 
#a point along the specified axis, and then a third point in the plane created by that normal. This would need to be
#constrained more depending on the location of the other robot, and there's a place to do that later in the code. 
def generate_targets(point, axis, rotation):
  r = Rotation.from_angle_axis(rotation, axis)
  axis = normalized(asarray(axis))
  axis2 = unit_orthogonal_vector(axis);
  p1 = [];
  p1.append(point)
  p1.append(point + axis*25)
  p1.append(point + r*(axis2*25))
  return asarray(p1);


#Create the first robot. We generate an armset as a system that's used for rendering,  and a dummy armset which is used for
#computing the IK. 

armset1 = System([Arm('a1',6)], generate_targets(target_point, target_axis, 0), 1)
armset_dummy = System([Arm('a2',6)], generate_targets(target_point, target_axis, 0),1)


#This updates the rendering for the single robot. 
def update_1():
  frame = props.get("frame")()
  state = unstructured_path1[int(frame)]
  for arm in armset1.arms:
    for i, node in enumerate(arm.nodes):
      node.rotation.set(float(state[i]) * 180/pi)


#Based on the soliution path returned by the path planner, 
#we recompute the locations of the final targets, and pass these to the planner for the second robot. 
def getTargetsFromPath():
  targets =[]
  for i, step in enumerate(angle_path):
    state = step[-1]
    for arm in armset_dummy.arms:
      for j, node in enumerate(arm.nodes):
        node.rotation.set(float(state[j]) * 180/pi)
    targets.append(armset_dummy.effectors()[0])
  return targets

#The general procedure is to set up the targets/acquire the IK solutions, and then in a separate routine
#send the relevant data to the C++ function that creates an RRT and finds paths. 

def setup_1(step1, step2):  #The parameters here are the indices of the points we want to find paths between. 
  global angle_path
  global location

  angle_goals = [];
  points = []
  targets = []
  all_goals = getContours()

  for t in range(step1, step2): #Time control here.
    armset_dummy.set_bike_point(all_goals[t]) #Acquire the appropriate point on the bike that we want to treat as an effector
    sample_axis = target_axis #Normal axis from the effector, we will sample in space maintaining this normal. 
    coarse_samples = []
    counter = 0
    for x in linspace(-200, 200, 4): #Generate a grid of points
      for y in linspace(-700, 700, 4):
        counter = counter+1
        this_target = target_point + array([x, y, 0])  #Can add a z here, we're just sampling the entire space basicaly. 
       # for theta in linspace(2*pi/3, 4*pi/3, 3):  #3 points become the target - the point, a point above this first point, and
       #then a third point in the plane created by the normal. The location of this third point is specified by the an angle in the plane. 
        triangle = generate_targets(this_target, sample_axis, pi) #Generate the targets, here choosing the angle to be pi(facing the other robot)
        armset_dummy.set_target(triangle) #Sets the target so the IK knows what to solve.
        points.append(triangle[0]) #"points" is a data structure used for plotting the sampling if you want. 
        points.append(triangle[1])
        points.append(triangle[2])
        goal, e = armset_dummy.solve_ik(0)
        print magnitude(e) #there are some convergence issues here, so filter bad solutions. 
        if magnitude(e) < .1:
          coarse_samples.append(goal)
    angle_goals.append(random.shuffle(coarse_samples))
  points.append(target_point)
  armset1.set_target(points) #Set the target for the robot that's actually rendered. 
  armset1.set_bike_point(all_goals[t]) 

  return angle_goals, targets

#Gathers all the relevant parameters and then passes them to the plan_1 function, which is a C++ wrapper. 
def solve_1(step1, step2):
  global angle_path
  global unstructured_path1
  print "Solving 1"
  angle_goals, targets = setup_1(step1, step2) 

  #Get the initial offsets for the robot meshes
  robot_origins = asarray(get_origins())

  #Get the robot meshes into an array
  for arm in armset1.arms:
    robot_meshes = arm.submeshes()
 
 #Computes the path. Returns it in a path/motion/step structure, where each layer is a vector. 
  angle_path = plan_1(6, step2-step1, robot_origins, robot_meshes, 
    [bike_mesh()], .1, 150 , 150., initial_rotation, array([0, 0, 0]), .01, angle_goals, smoothing, 100)

  targets = getTargetsFromPath()
  unstructured_path1 = convertUnstructured(angle_path)
  print targets
  return angle_path,targets

def getAnglePath(step1, step2):
  return solve_1(step1, step2)


def view_1():
  # Set up main window
  solve_1(0,1)
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
