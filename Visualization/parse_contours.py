from __future__ import division
import sys
from gui import *
import geode
import scipy
from geode import *
from numpy import *
from numpy.linalg import *
from math import *
from geode.openmesh import *
from geode.vector import *
from scipy import *
from OpenGL import GL, GLUT
from scipy.optimize import fmin_tnc
from props import props
from follow import *
from scatter3d_demo import *
from kinematic_chain import *

def getContour():
	file = 'pts.txt'
	origins = [ [float(i) for i in l.strip().split(' ')] for l in open(file)]
	return origins

	
