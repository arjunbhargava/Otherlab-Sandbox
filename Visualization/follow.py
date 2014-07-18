from __future__ import division
import sys
from gui import *
import geode
from geode import *

from numpy import *
from numpy.linalg import *

from geode.openmesh import *
from geode.vector import *

from OpenGL import GL, GLUT

from props import props

base = './origins_kr16'
origins = [ [float(i) for i in l.strip().split(' ')] for l in open(base+'.txt')]

def get_origins():
  return origins

def render_mesh_at(i,tm):
  c = wheel_color(pi*i)
  GL.glColor3f(c[0],c[1],c[2])
  render_mesh(tm,ctFace,False)

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


def end_effector(node):
  return node.frame().t + node.frame().r*array([1.,0,0])*100.

def normal(node):
  return node.frame().t + node.frame().r*array([1.,0,0])*50.

def up(node):
  return end_effector(node) + node.frame().r*array([0,0,1])*50.

class Arm(object):
  def __init__(self,n,dof,eff_functions,bf=Frames(array([0.,0.,0.]),Rotation.from_angle_axis(0.,array([1.,0.,0.])))):
    self.name = n
    self.dof = dof
    self.nodes = []
    self.base_frame = bf
    self.axes = []
    self.effector_functions = eff_functions
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
    return [  f(self.nodes[-1]) for f in self.effector_functions]

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

class System():
  def __init__(self,arms):
    self.arms=arms

  def effectors(self):
    e = []
    for arm in self.arms:
      e.extend(arm.effectors())
    return e

  #@cache
  def Jacobian(self):
    eff = self.effectors()
    nn = self.n_nodes()
    J = zeros(shape=(len(eff),3,nn))
    j_offset = 0
    i_offset = 0
    for arm in self.arms:
      eff = arm.effectors()
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

  @cache
  def targets():
    t = 0#props.get("frame")()*.01
    p1 = array([817,400*cos(t),1425 + 300*sin(t)])
    return array([p1,p1+array([-50.,0,0.]),p1+array([0.,0.,50.])])
    #return array([p1,p1+array([0,0,50.]),p1,p1-array([0,0,50.])])

  def clamp(self, w,d,norm_type=None):
    n = norm(w,ord=norm_type)
    if n < d:
      return w
    else:
      return w*d/n

  def update_effector(self):
    gamma_max = pi*.25

    e = self.targets()- self.effectors()
    for i,v in enumerate(e):
      e[i] = self.clamp(v, 350.)

    iters = 0
    while norm(e) > .05 and iters < 1e2:
      iters+=1
      e = self.targets()-self.effectors()
      for i,v in enumerate(e):
        e[i] = self.clamp(v, 350.)
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
      #print 'U : ',U.shape
      #print 's: ',s
      #print 'V : ',V.shape
      S = ndarray(shape=(U.shape[1],V.shape[0]))
      #print 'S: ',S.shape
      S[:len(s),:len(s)] = diag(s)
      #assert allclose(Jm,dot(U,dot(S,V)))

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
        phi_i = self.clamp(sigi_inv*alpha*vv,gamma_i,inf)
        phisum+=phi_i

      dt = self.clamp(phisum,gamma_max,inf)
      #print 'dt: ',dt
      for i,val in enumerate(dt):
        p = self.axis_props()[i]
        p.set(p() + val*180/pi)
      
      #return array(angles * pi/180)
    print self.arms[0]
    return array([p()*pi/180 for p in self.axis_props()])

    #f = self.arms[0].nodes[-1].frame() 
    #ang = f.r.euler_angles()
    #return ang

#frame2 = Frames(array([2043.,0.,0.]),Rotation.from_angle_axis(pi,array([0.,0.,1.])))

effector_functions = [end_effector,normal,up]

#arms = [Arm('welding',6,effector_functions),Arm('fixturing',6,effector_functions,frame2) ]

class KukaScene(Scene):
  def __init__(self,system):
      self.system = system
      
  @cache_method
  def lists(self):
    return [arm.lists() for arm in self.system.arms]

  def render(self,*args):
    for t in self.system.targets():
      GL.glPointSize(10)
      GL.glBegin(GL.GL_POINTS)
      GL.glColor3f(1,0,0)
      GL.glVertex3f(t[0],t[1],t[2])
      GL.glEnd()
    for e in self.system.effectors():
      GL.glDisable(GL.GL_DEPTH_TEST)
      GL.glPointSize(5)
      GL.glBegin(GL.GL_POINTS)
      GL.glColor3f(0,1,0)
      GL.glVertex3f(e[0],e[1],e[2])
      GL.glEnd()
      GL.glEnable(GL.GL_DEPTH_TEST)

    GL.glDisable(GL.GL_DEPTH_TEST)
    GL.glPointSize(50)
    GL.glBegin(GL.GL_POINTS)
    GL.glColor3f(0,0,1)
    GL.glVertex3f(716, 397, 1427)
    GL.glVertex3f(1455, 0, 1320)
    GL.glEnd()
    GL.glEnable(GL.GL_DEPTH_TEST)

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

