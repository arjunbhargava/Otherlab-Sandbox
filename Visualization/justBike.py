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
from bikefollow import *
from math import *
from scatter_points import *
from kinematic_chain import *

angle_path = [];
target_point = [1000,0, 1400]
target_axis = [0, 0, 1]
initial_rotation = 0 


def generate_targets(point, axis, rotation):
  r = Rotation.from_angle_axis(rotation, axis)
  axis = normalized(asarray(axis))
  axis2 = unit_orthogonal_vector(axis);
  axis3 = normalized(cross(axis2, axis))
  p1 = [];
  p1.append(point + r*(axis2 * sqrt(25**2 + 25**2)))
  p1.append(point + r*(axis2 * -25+ axis3 * -25))
  p1.append(point + r*(axis2 * -25 + axis3 * 25))
  return asarray(p1);

armset1 = System([Arm('a1',6, effector_functions)], generate_targets(target_point, target_axis, 0), 0)
armset_dummy = System([Arm('a2',6,effector_functions)], generate_targets(target_point, target_axis, 0),0)

def view():
  global angle_path
  global location

  app = QEApp(sys.argv,True)
  main = MainWindow(props)
  ks = KukaScene(armset1)
  main.view.add_scene("kuka",ks)
  main.view.show_all(True)

  points = []
  all_goals, rotation_axes, rotation_angles, target_axes, normals = armset1.pointNormals() 
  points = getContours();
  armset1.set_target(points)

  main.resize_timeline(100);
  main.init()
  app.run()


def main():
  parser.parse(props,'visualizer')
  view()

if __name__=='__main__':
  main()
