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



#Reads the origins of each of the meshes, returns them in an array. 
def get_origins():
  base = './origins_kr16'
  return [ [float(i) for i in l.strip().split(' ')] for l in open(base+'.txt')]


#This is the input function for the touch off points for the path planner. 
def getContours():
  touch = 'pts.txt' #Input would be changed here...
  contours = [ [float(i) for i in l.strip().split(' ')] for l in open(touch)]  
  for i, val in enumerate(contours):
 #   contours[i] = contours[i] + array([-100, 0,-300]) If the bike is translated relative to the flange, needs to be changed. 
    contours[i] = Rotation.from_angle_axis(pi/2, [0, 0 , 1]) * contours[i]
  return contours
  

#Reads in the bike mesh, and rotates it to align it with the initial configuration of the robots
#All manipulations are relative to this initial position
def bike_mesh():
  bike = TriMesh()
  bike.read("bike_mesh.om")
 # bike.translate(array([-100, 0, -300])) #If you want to move the bike relative to the flange do that here
  bike.set_X(Rotation.from_angle_axis(pi/2, [0, 0, 1]) * bike.X())
  bike.request_face_normals()
  bike.request_vertex_normals()
  bike.update_normals()
  return bike


#Rendering function, fairly self explanatory
def render_mesh_at(i,tm):
  c = wheel_color(pi*i)
  GL.glColor3f(c[0],c[1],c[2])
  render_mesh(tm,ctFace,False)

#No major changes from the initial implementation of the node object. 
class Node(object):
  def __init__(self,m,offset,r,p,bf):
    assert offset is not None
    assert m is not None
    self.base_frame = bf
    self.offset = Prop("offset",offset)
    self.mesh = Prop("mesh",m)
    self.axis = Prop("axis",None)
    self.negate_rotation = Prop("negate",False)
    self.rotation = r

    m.request_face_normals()
    m.request_vertex_normals()
    m.update_normals()

    self.parent = Prop("parent_node",p)

  @cache_method
  def world_axis(self):
    wa = self.frame().r*self.axis()
    wa *= 1./norm(wa)
    return wa

  @cache_method
  def frame(self):
    assert self.axis() is not None
    f = self.parent().frame() if self.parent() is not None else self.base_frame
    factor = -1. if self.negate_rotation() else 1.
    r = Rotation.from_angle_axis(factor*self.rotation()*pi/180.,f.r*self.axis())
    return Frames(f.t+f.r*self.offset(),r*f.r)


#Arm object which contains all of the metadata for the robot meshes. 
class Arm(object):
  def __init__(self,n,dof,bf=Frames(array([0.,0.,0.]),Rotation.from_angle_axis(0.,array([1.,0.,0.])))):
    self.name = n
    self.dof = dof
    self.nodes = []
    self.base_frame = bf
    self.axes = []
    for d in range(self.dof):
      self.axes.append(props.add(self.name+"axis%s"%d,0.))
    this_origins = get_origins()
    t = array([0,0,0])
    for i in range(6):
      o = array(this_origins[i]) - t
      t = this_origins[i]
      m = TriMesh()
      bf ='kr16/'
      m.read(bf+'axis%s.stl'%(i+1))
      r = self.axes[i]
      n = Node(m,o,r,self.nodes[i-1] if i>0 else None, self.base_frame)
      self.nodes.append(n)

    self.nodes[0].axis.set((0,0,1))
    self.nodes[0].rotation.set_min(-185.).set_max(185.).set_step(1)
    self.nodes[0].negate_rotation.set(True)

    self.nodes[1].axis.set((0,1,0))
    self.nodes[1].rotation.set_min(-65.).set_max(125.).set_step(1)

    self.nodes[2].axis.set((0,1,0))
    self.nodes[2].rotation.set_min(-210.).set_max(64.).set_step(1)

    self.nodes[3].axis.set((1,0,0))
    self.nodes[3].rotation.set_min(-165.).set_max(165.).set_step(1)
    self.nodes[3].negate_rotation.set(True)

    self.nodes[4].axis.set((0,1,0))
    self.nodes[4].rotation.set_min(-130.).set_max(130.).set_step(1)

    self.nodes[5].axis.set((1,0,0))
    self.nodes[5].rotation.set_min(-350.).set_max(350.)
    self.nodes[5].rotation.set_step(.5)
    self.nodes[5].negate_rotation.set(True)

    self.base_mesh = TriMesh()
    base_fn = 'kr16/base.stl'
    self.base_mesh.read(base_fn)
    self.base_mesh.transform(self.base_frame.matrix())
    self.base_mesh.request_face_normals()
    self.base_mesh.request_vertex_normals()
    self.base_mesh.update_normals()

  @cache_method
  def effectors(self):
    return self.nodes[-1]

  @cache_method
  def mesh(self):
    m = TriMesh()
    m.request_face_normals()
    m.request_vertex_normals()
    m.add_mesh(self.base_mesh)
    for node in self.nodes:
      tm=node.mesh().copy()
      tm.transform(node.frame().matrix())
      m.add_mesh(tm)
    m.update_normals()
    return m

  #Return each piece of the robot mesh in an array, used for collision checking in the planning code. 
  def submeshes(self):
    meshes = []
    for node in self.nodes:
      tm = node.mesh().copy()
      meshes.append(tm)
    return meshes

  @cache_method
  def lists(self):
    ls = []
    ls.append(cache_render(curry(render_mesh_at,0,self.base_mesh)))
    for i,n in enumerate(self.nodes):
      tm = n.mesh()
      ls.append(cache_render(curry(render_mesh_at,i+1,tm)))
    return ls


#Class that has all the methods for computing the inverse kinematics. 

class System():

  #Requires the correct arm set, an appropriate target (3 points here), and the appropriate robot number (1 or 2). 
  def __init__(self,arms,target, number):
    self.arms = arms
    self.target_points = target
    self.number = number
    self.bike_point = []
    self.bm = bike_mesh()
    self.bft = self.bm.face_tree()

  #Usually 3 points, can be easily changed. 
  def set_target(self, new_target): 
    self.target_points = new_target

  #Returns the targets
  def targets(self):
    return self.target_points

  #This is the point on the bike we are interested in computing IK on. 
  def set_bike_point(self, bike_point):
    self.bike_point = bike_point

  #Computes the locations of the effectors, around a given "bike_point"
  def effectors(self):
    node = self.arms[0].effectors()
    this_effectors = [] 
    if self.number == 1: # For the first robot we need 3 points. Note the offset from the flange. 
      this_effectors.append(node.frame().t + node.frame().r*(array([1,0,0])*100 + self.bike_point)) #Point 1
      bike_normal = normalized(self.pointNormal(self.bike_point))
      this_effectors.append(node.frame().t + node.frame().r*(array([1.,0,0])*100  + self.bike_point + bike_normal*25)) #Point 2
      flange_normal = node.frame().r * array([1, 0, 0])
      angle_between = acos(clamp(dot(bike_normal, flange_normal), -1, 1))
      orthogonal_projection = normalized(flange_normal) - (magnitude(flange_normal) * cos(angle_between) * bike_normal)
      this_effectors.append(node.frame().t + node.frame().r*(array([1.,0,0])*100  + self.bike_point + orthogonal_projection*25)) #Point 3
    else: #for the second robot we just align the normals. 
      this_effectors.append(node.frame().t)
      this_effectors.append(node.frame().t + node.frame().r*(array([1.,0,0])*25))
    return this_effectors

  def Jacobian(self):
    eff = self.effectors()
    nn = self.n_nodes()
    J = zeros(shape=(len(eff),3,nn))
    j_offset = 0
    i_offset = 0
    for arm in self.arms:
      eff = self.effectors()
      for i,node in enumerate(arm.nodes):
        for j,e in enumerate(eff):
          delta = e - node.frame().t
          J[j+j_offset,:,i+i_offset] = cross(node.world_axis(),delta)
          if node.negate_rotation():
            J[j+j_offset,:,i+i_offset]*=-1
      i_offset+=len(arm.nodes)
      j_offset+=len(eff)

    return J

  def n_nodes(self):
    return sum([len(arm.nodes) for arm in self.arms])

  def axis_props(self):
    l = []
    for arm in self.arms:
      l.extend(arm.axes)
    return l

#Clamp function, just renamed here. 
  def this_clamp(self, w,d,norm_type=None):
    n = norm(w,ord=norm_type)
    if n < d:
      return w
    else:
      return w*d/n

#IK solver. Pretty much unmodified except that it calls the relevant functions for getting effectors/targets
  def solve_ik(self, initialize):
    gamma_max = pi*.5
    if initialize < 1:
      for i in range(0, self.n_nodes()):
        p = self.axis_props()[i]
        p.set(0)

    e = self.targets()- self.effectors()
    for i,v in enumerate(e):
      e[i] = self.this_clamp(v, 500.)

    iters = 0
    while norm(e) > .1 and iters < 3e2:
      iters+=1
      e = self.targets()-self.effectors()
      for i,v in enumerate(e):
        e[i] = self.this_clamp(v, 500.)
      e = ndarray.flatten(e)

      if norm(e) == 0:
        return

      J = self.Jacobian()
      rho = geode.vector.magnitudes(J,axis=1)
      nn = self.n_nodes()
      U,s,V = svd(J.reshape(3*len(self.effectors()),nn))

      U = asarray(U)
      V = asarray(V)
      Umags = magnitudes(U)
      phisum = zeros(nn)
      S = ndarray(shape=(U.shape[1],V.shape[0]))
      S[:len(s),:len(s)] = diag(s)

      for i in range(U.shape[1]):
        Ui = U[:,i]
        alpha = dot(Ui,e)
        Ni = sum(array([Umags])[:,i])
        if i >= S.shape[1] or s[i] == 0:
          continue
        #ith singular value
        sigi_inv = 1./s[i]
        #ith row of J; note: different from SLDS paper, which is wrong.
        vv = V[i,:]
        abv = abs(vv)
        Mi = 0.
        #sum Mil for all effectors
        for l in range(len(self.effectors())):
          Mil = sigi_inv*dot(rho[l],abv)
          Mi += Mil
        if Mi == 0:
          continue

        #p = axis_props()[i]
        #gamma_max=min(min(abs(p.min-p()),abs(p.max-p())),pi*.25)
        gamma_i = min(1.,Ni/Mi)*gamma_max
        phi_i = self.this_clamp(sigi_inv*alpha*vv,gamma_i,inf)
        phisum+=phi_i

      dt = self.this_clamp(phisum,gamma_max,inf)
      for i,val in enumerate(dt):
        p = self.axis_props()[i]
        p.set(p() + val*180/pi)
      
    initial_theta = array([p()*pi/180 for p in self.axis_props()])
    return initial_theta, e

#Renders the bike frame, offsets it from the flange by 100 units. 
  def bike(self):
    frame = self.arms[0].effectors().frame()
    bike = bike_mesh()    
    bike.set_X(frame * bike.X())
    bike.translate(frame.r * (array([1,0,0]) * 100))
    bike.request_face_normals()
    bike.request_vertex_normals()
    bike.update_normals()
    return bike

#Returns the normal of the closest point from the bike mesh. 
  def pointNormal(self, point):
    close_point = self.bft.closest_point(point, 200)
    return self.bm.smooth_normal(close_point[1], close_point[2])

#Just a bunch of rendering. 
class KukaScene(Scene):
  def __init__(self,system):
      self.system = system
      
  @cache_method
  def lists(self):
    return [arm.lists() for arm in self.system.arms]

  def render(self,*args):
    for i, e in enumerate(self.system.effectors()):
      GL.glPointSize(20)
      GL.glBegin(GL.GL_POINTS)
      v = [0., 0., 0.];
      v[i%3]= 1
      GL.glColor3f(v[0], v[1], v[2]) 
      GL.glVertex3f(e[0],e[1],e[2])
      GL.glEnd()
      GL.glEnable(GL.GL_DEPTH_TEST)

    for i, e in enumerate(self.system.targets()):
      # if i == len(self.system.targets())-1:
      #   GL.glPointSize(20)
      # else:
      GL.glPointSize(5)
      GL.glBegin(GL.GL_POINTS)
      v = [0., 0., 0.];
      v[i%3]= 1
      GL.glColor3f(v[0], v[1], v[2]) 
      GL.glVertex3f(e[0],e[1],e[2])
      GL.glEnd()
      GL.glEnable(GL.GL_DEPTH_TEST)

    if self.system.number == 1:
      render_mesh_at(0, self.system.bike())

    for armid,armlists in enumerate(self.lists()):
      for id,list in enumerate(armlists):
        with gl_scope():
          if id!=0:
            c = wheel_color(pi*id+.05)*.8
            GL.glColor4f(c[0]*1.2,c[1]*1.2,c[2]*1.2,.2)
            GL.glLineWidth(4)
            GL.glBegin(GL.GL_LINES)
            p = self.system.arms[armid].nodes[id-1].frame().t
            GL.glVertex3f(p[0],p[1],p[2])
            v = self.system.arms[armid].nodes[id-1].world_axis()
            scl = 400
            GL.glColor4f(c[0],c[1],c[2],.7)
            GL.glVertex3f(p[0]+scl*v[0],p[1]+scl*v[1],p[2]+scl*v[2])
            GL.glEnd()
            gl_mult_mx(self.system.arms[armid].nodes[id-1].frame().matrix())
          list().call()

  def bounding_box(self):
    b = Box3d((0,0,0),(0,0,0))
    for arm in self.system.arms:
      for n in arm.nodes:
        tm = n.mesh().copy()
        tm.transform(n.frame().matrix())
        b = Box3d(minimum(b.min,tm.bounding_box().min),maximum(b.max,tm.bounding_box().max))
    return b

