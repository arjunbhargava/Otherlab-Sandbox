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
from scatter3d_demo import *
from kinematic_chain import *

angle_path = [];



def generate_targets(point, axis):
  axis = asarray(axis)
  axis = normalized(axis)
  axis2 = unit_orthogonal_vector(axis);
  axis3 = normalized(cross(axis2, axis))
  p1 = [];
  p1.append(point + axis2 * sqrt(1**2 + 1**2))
  p1.append(point + axis2 * -1+ axis3 * -1)
  p1.append(point + axis2 * -1 + axis3 * 1)
  return asarray(p1);


target_point = [1000, 0, 1000];
target_axis = [0, 0,1]
location = array([2400, 0,0])

initial_rotation = pi

armset1 = System([Arm('a1',6, effector_functions)], generate_targets(target_point, target_axis), 0)
armset_dummy = System([Arm('a2',6,effector_functions)], generate_targets(target_point, target_axis),0)
   
obstacle_arm_frame = Frames(location,Rotation.from_angle_axis(initial_rotation,array([0,0,1])))
obstacle_arm = System([Arm('obstacle',6,effector_functions2,bf=obstacle_arm_frame)], generate_targets(target_point, target_axis),1)
armset_dummy2 = System([Arm('obstacle2',6,effector_functions2,bf=obstacle_arm_frame)], generate_targets(target_point, target_axis),1)



def view():
  # Set up main window
  global angle_path
  global location

  all_goals = getContour()
  goal_list = [];
  angle_goals = [];
  for t in range(0, 1):#len(goal_list)):
    touch_point = all_goals[t];
    print touch_point
    p1 = array([1250, 0, 1400]) + touch_point #500 + 400*cos(pi/10 * t), 1425 + 300*sin(pi/10 *t)])
    goal_list.append(p1)
    coarse_samples = []
    target_list = []
    for cone_angle in linspace(pi/12, pi/24, 3):
      initial_axis = [sin(cone_angle), 0, cos(cone_angle)]
      counter = 0
      for i in linspace(0, 2* pi, 10):
        axis_r = Rotation.from_angle_axis(i, target_axis)
        sample_axis = axis_r * initial_axis
        armset_dummy.set_target(generate_targets(p1, sample_axis))
        armset_dummy2.set_target(generate_targets(p1, sample_axis))
        goal = armset_dummy.solve_ik(counter)
        goal2 = armset_dummy2.solve_ik(counter)
        final_configuration = append(goal, goal2)
        coarse_samples.append(final_configuration)
        angle_goals.append(coarse_samples)
        counter = counter + 1
  print goal_list
  print angle_goals

  robot_origins = asarray(get_origins())
  robot_origins2 = asarray(get_origins());
  both_origins = array([robot_origins, robot_origins2])

  both_meshes = []
  for arm in armset1.arms:
    both_meshes.append(arm.submeshes())
  for arm in obstacle_arm.arms:
    both_meshes.append(arm.submeshes())


#static Nested<real> plan(unsigned int links, double robot_number, vector<Array<real>> goalState, 
#Array<Vector<real,3>,2> parsed_offsets, vector<vector<Ref<TriMesh>>> robotMeshes, vector<Ref<TriMesh>> obstacleMeshes, 
#double resolution, double range, double solve_time, double initial_angle, initial locations, tolerance) 

  angle_path = sample_path(6,2, goal_list, both_origins, both_meshes, 
    [bike_mesh()], .1, 50 , 30., initial_rotation, [array([0, 0, 0]), location], .01, angle_goals, 0, 100)
#  props.get("last_frame").set(len(angle_path))
    # main.resize_timeline(100);
  plot_points(angle_path, goal_list)
  #kukaview = main.add_view("kuka")
 #Turn on the listener for props.get("frame"), then calls update (timeline updating)
  #l = listen(props.get("frame"),update)

  # Launch
 # main.init()
  #app.run()

def main():
  #Pulls the properties for props
  parser.parse(props,'visualizer')
  view()
 
if __name__=='__main__':
  main()



