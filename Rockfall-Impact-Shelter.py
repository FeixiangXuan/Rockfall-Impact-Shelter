####################################Numerical Study of Rockfall Impact on a Shelter#######################################
#                                                                                                                        #
#                                              <Editor:Feixiang Xuan_China>                                              #
#                                                                                                                        #
#                                       <Supervisor:Prof.Francesco Calvetti_Italy>                                       #
#                                                                                                                        #
#                                        <Special_thanks:Dr.Václav Šmilauer_Czech>                                       #
#                                                                                                                        #
#                                   <Special_thanks:Dr.Emmanouil Vrairaktaris_Greece>                                    #
#                                                                                                                        #
##########################################################################################################################

# 感谢玄承涛从小到大在学习上的帮助，在人生道路上的指引
# 此致，敬礼 :)
                                  

print('Study the impact force onto the RC shelter')

# Import the using modules
from woo.core import *
from woo.dem import *
from woo.fem import *
import woo
import woo.gl
import math
from math import pi
from minieigen import *
from past.builtins import execfile
import woo.pack, woo.utils, numpy

# Set Api
woo.master.usesApi=10103

# Set the graphic options
S=woo.master.scene=Scene(fields=[DemField(gravity=(0,0,-10))],dtSafety=0.8)
S.gl.membrane.node=True
S.gl=woo.gl.GlSetup(
    woo.gl.Gl1_Membrane(uScale=0,relPhi=0,refConf=False),
    woo.gl.Gl1_DemField(shape=woo.gl.Gl1_DemField.shapeNonSpheres,colorBy=woo.gl.Gl1_DemField.colorDisplacement,vecAxis='norm',colorBy2=woo.gl.Gl1_DemField.colorDisplacement),
)
#S.gl.demField.colorRange2.mnmx=(0,2.)

# Define the soil grains (num=20000,mechanical property as below)
if 1:
 sp=woo.pack.SpherePack()
 sp.makeCloud((-1.4,-1.4,0.1),(1.4,1.4,0.65),rMean=.025,rRelFuzz=.3,periodic=False,num=20000,)
 # add to the simulation
 sp.toSimulation(S,mat=FrictMat(young=200e6,density=2650,ktDivKn=0.25,tanPhi=0.3))

# Define the RC plate (Dimension=2.8m*2.8m,thickness=0.2m,mechanical property as below)
plate_thickness=0.2
ht=plate_thickness/2.
xmax,ymax=2.8,2.8
xdiv,ydiv=11,11
mat=FrictMat(young=30e9,tanPhi=.3,ktDivKn=.25,density=2500)
ff=woo.pack.gtsSurface2Facets(woo.pack.sweptPolylines2gtsSurface([[(x,y,0) for x in numpy.linspace(-1.4,xmax/2,num=xdiv)] for y in numpy.linspace(-1.4,ymax/2,num=ydiv)]),flex=True,halfThick=ht,mat=mat)
S.dem.par.add(ff,nodes=True)
 
# Define the box which hold the soil grains
wall=woo.dem.Wall.makeBox(
  ((-1.4,-1.4,0.1),(1.4,1.4,1.1)),
  (1,1,0,1,1,0)
)
S.dem.par.add(wall)

# Set the boundary condition of the RC plate (Fix the translational D.o.F at the four corner of it)
for n in S.dem.nodes:
    n.dem.blocked=''
    if (n.pos[0]==-1.4) and (n.pos[1]==-1.4): n.dem.blocked='xyz'
    if (n.pos[0]==1.4) and (n.pos[1]==-1.4): n.dem.blocked='xyz'
    if (n.pos[0]==1.4) and (n.pos[1]==1.4): n.dem.blocked='xyz'
    if (n.pos[0]==-1.4) and (n.pos[1]==1.4): n.dem.blocked='xyz'

# Set the fake inertia and mass to the lumped node of the CST
for n in S.dem.nodes: DemData.setOriMassInertia(n)

# block the rotational D.o.F of the spheres
for i in range(len(sp)):
    S.dem.par[i].blocked=''
    if i in range(len(sp)): S.dem.par[i].blocked='XYZ'

# radius of the boulder
r=0.3
# weight of the boulder
m=100
# volume of the boudler
v=4/3*pi*r**3
# density of the boulder
d=m/v
# create material of boulder
mat_boulder=woo.dem.FrictMat(
                               young=200e6, 
                               density=d, 
                               tanPhi=0.3, 
                               ktDivKn=0.25 
)
# create boulder shape
sphere=woo.dem.Sphere.make(
                          (0,0,(0.53+r)),
                          radius=r, 
                          wire=False, 
                          color=2,                
                          mat=mat_boulder, 
                          fixed=False 
)
# according to the law of conservation of energy
# MgH=1/2*M*V^2
g=9.8
H=50
v=math.sqrt(2*g*H)

# Find the supports nodes
i=0
for i in range(len(S.dem.nodes)):
    # find the corner nodes (reaction)
    if (S.dem.nodes[i].pos[0]==-1.4) and (S.dem.nodes[i].pos[1]==-1.4): 
       print('c_1',i)
       c_1=i
    elif (S.dem.nodes[i].pos[0]==1.4) and (S.dem.nodes[i].pos[1]==-1.4): 
       print('c_2',i)
       c_2=i 
    elif (S.dem.nodes[i].pos[0]==1.4) and (S.dem.nodes[i].pos[1]==1.4): 
       print('c_3',i)
       c_3=i 
    elif (S.dem.nodes[i].pos[0]==-1.4) and (S.dem.nodes[i].pos[1]==1.4): 
       print('c_4',i)
       c_4=i 
    else: i=i+1

# find the diagonal nodes (check the displ)
i=0
for i in range(len(S.dem.nodes)):
    if (S.dem.nodes[i].pos[0]>-1.4 and S.dem.nodes[i].pos[0]<-0.84) and (S.dem.nodes[i].pos[1]>-1.4 and S.dem.nodes[i].pos[1]<-0.84) and (S.dem.nodes[i].pos[2]==0): 
       print('d_5',i)
       d_5=i
    elif (S.dem.nodes[i].pos[0]==-0.84) and (S.dem.nodes[i].pos[1]==-0.84) and (S.dem.nodes[i].pos[2]==0): 
       print('d_4',i)
       d_4=i
    elif (S.dem.nodes[i].pos[0]>-0.84 and S.dem.nodes[i].pos[0]<-0.28) and (S.dem.nodes[i].pos[1]>-0.84 and S.dem.nodes[i].pos[1]<-0.28) and (S.dem.nodes[i].pos[2]==0): 
       print('d_3',i)
       d_3=i
    elif (S.dem.nodes[i].pos[0]>-0.56 and S.dem.nodes[i].pos[0]<-0.) and (S.dem.nodes[i].pos[1]>-0.56 and S.dem.nodes[i].pos[1]<-0.) and (S.dem.nodes[i].pos[2]==0): 
       print('d_2',i)
       d_2=i
    elif (S.dem.nodes[i].pos[0]==0) and (S.dem.nodes[i].pos[1]==0): 
       print('d_1',i)
       d_1=i
    else: i=i+1
    
# find the symmetry diagonal nodes (compare the displ)
i=0
for i in range(len(S.dem.nodes)):
    if (S.dem.nodes[i].pos[0]>-1.4 and S.dem.nodes[i].pos[0]<-0.84) and (S.dem.nodes[i].pos[1]>0.84 and S.dem.nodes[i].pos[1]<1.4) and (S.dem.nodes[i].pos[2]==0): 
       print('d_5_s',i)
       d_5_s=i
    elif (S.dem.nodes[i].pos[0]==-0.84) and (S.dem.nodes[i].pos[1]>0.56 and S.dem.nodes[i].pos[1]<0.84) and (S.dem.nodes[i].pos[2]==0): 
       print('d_4_s',i)
       d_4_s=i
    elif (S.dem.nodes[i].pos[0]>-0.84 and S.dem.nodes[i].pos[0]<-0.28) and (S.dem.nodes[i].pos[1]>0.28 and S.dem.nodes[i].pos[1]<0.56) and (S.dem.nodes[i].pos[2]==0): 
       print('d_3_s',i)
       d_3_s=i
    elif (S.dem.nodes[i].pos[0]>-0.56 and S.dem.nodes[i].pos[0]<-0.) and (S.dem.nodes[i].pos[1]>0. and S.dem.nodes[i].pos[1]<0.28) and (S.dem.nodes[i].pos[2]==0): 
       print('d_2_s',i)
       d_2_s=i
    elif (S.dem.nodes[i].pos[0]==0) and (S.dem.nodes[i].pos[1]==0) and (S.dem.nodes[i].pos[2]==0): 
       print('d_1_s',i)
       d_1_s=i
    else: i=i+1

# Find the overall force and calculate the stress increment at these 8 areas (8 specific area)
# radius_1
r_1=0.15
# radius_2
r_2i=0.25
r_2o=0.35
# radius_3
r_3i=0.55
r_3o=0.65
# radius_4
r_4i=0.85
r_4o=0.95
# radius_5
r_5i=1.15
r_5o=1.25
# raidus_new
r_6=0.45
# area_1
a_1=pi*r_1**2
# area_2
a_2=pi*r_2o**2-pi*r_2i**2
# area_3
a_3=pi*r_3o**2-pi*r_3i**2
# area_4
a_4=pi*r_4o**2-pi*r_4i**2
# area_5
a_5=pi*r_5o**2-pi*r_5i**2
# area_new
a_6=pi*r_2i**2-pi*r_1**2
a_7=pi*r_6**2-pi*r_2o**2
a_8=pi*r_3i**2-pi*r_6**2
             
# determine the top layer height
h=[]
for i in range(len(sp)):
    a=S.dem.par[i].pos[2]
    h=h+[a]
        
# find the exact particle
for i in range(len(sp)):
    if S.dem.par[i].pos[2]==max(h):
       print('the heightest particle is',i)
       h_n=i

print('the top layer height is ',max(h))

# Create the engine
S.engines=DemField.minimalEngines(damping=.25,verletDist=-0.01)+[IntraForce([In2_Membrane_ElastMat(bending=True)])]+[PyRunner(50,'check(S)')]

# From the Check(S), plotting the results and adding the rockfall
def check(S):
  if S.step==100:
      S.dem.par.add(sphere)
      p1=S.dem.par[-1]
      p1.vel=(0,0,-v)
  elif S.step>100:
      # total force
      total_force=0
      # stress
      force_1=0
      force_2=0
      force_3=0
      force_4=0
      force_5=0
      force_6=0
      force_7=0
      force_8=0
      stress_1=0
      stress_2=0
      stress_3=0
      stress_4=0
      stress_5=0
      stress_6=0
      stress_7=0
      stress_8=0
      for c in S.dem.con:
          # filter out the contacts with (par+par)+(par+membrane)
          if not isinstance (c.pA.shape,woo.dem.Wall) and not isinstance (c.pB.shape,woo.dem.Wall):
              # filter out the contacts with (par+membrane)
              if isinstance (c.pA.shape,woo.fem.Membrane) or isinstance (c.pB.shape,woo.fem.Membrane):
                  total_force=total_force+abs(c.phys.force[0])
                  # filter the contacts in the specific area
                  # area_1    
                  if (c.pB.pos[0]**2+c.pB.pos[1]**2)**0.5<r_1:     
                      force_1=force_1+abs(c.phys.force[0])
                      stress_1=force_1/a_1
              
                  # area_2    
                  if (c.pB.pos[0]**2+c.pB.pos[1]**2)**0.5<r_2o and (c.pB.pos[0]**2+c.pB.pos[1]**2)**0.5>r_2i:     
                      force_2=force_2+abs(c.phys.force[0])
                      stress_2=force_2/a_2
              
                  # area_3    
                  if (c.pB.pos[0]**2+c.pB.pos[1]**2)**0.5<r_3o and (c.pB.pos[0]**2+c.pB.pos[1]**2)**0.5>r_3i:     
                      force_3=force_3+abs(c.phys.force[0])
                      stress_3=force_3/a_3
              
                  # area_4    
                  if (c.pB.pos[0]**2+c.pB.pos[1]**2)**0.5<r_4o and (c.pB.pos[0]**2+c.pB.pos[1]**2)**0.5>r_4i:     
                      force_4=force_4+abs(c.phys.force[0])
                      stress_4=force_4/a_4
              
                  # area_5    
                  if (c.pB.pos[0]**2+c.pB.pos[1]**2)**0.5<r_5o and (c.pB.pos[0]**2+c.pB.pos[1]**2)**0.5>r_5i:
                      force_5=force_5+abs(c.phys.force[0])
                      stress_5=force_5/a_5
              
                  # area_6
                  if (c.pB.pos[0]**2+c.pB.pos[1]**2)**0.5<r_2i and (c.pB.pos[0]**2+c.pB.pos[1]**2)**0.5>r_1:
                      force_6=force_6+abs(c.phys.force[0])
                      stress_6=force_6/a_6
    
                  # area_7
                  if (c.pB.pos[0]**2+c.pB.pos[1]**2)**0.5<r_6 and (c.pB.pos[0]**2+c.pB.pos[1]**2)**0.5>r_2o:
                      force_7=force_7+abs(c.phys.force[0])
                      stress_7=force_7/a_7
              
                  # area_8
                  if (c.pB.pos[0]**2+c.pB.pos[1]**2)**0.5<r_3i and (c.pB.pos[0]**2+c.pB.pos[1]**2)**0.5>r_6:
                      force_8=force_8+abs(c.phys.force[0])
                      stress_8=force_8/a_8
      
      S.plot.addData(t=S.time,i=S.step,
      # rock boulder impact force [kN]
      uf=S.dem.par[-1].f[2]/1000,
      # reaction force [kN]
      r1=S.dem.nodes[c_1].dem.force[2]/1000,r2=S.dem.nodes[c_2].dem.force[2]/1000,r3=S.dem.nodes[c_3].dem.force[2]/1000,r4=S.dem.nodes[c_4].dem.force[2]/1000,
      # displacement [mm]
      d1=S.dem.nodes[d_1].pos[2]*1000,d2=S.dem.nodes[d_2].pos[2]*1000,d3=S.dem.nodes[d_3].pos[2]*1000,d4=S.dem.nodes[d_4].pos[2]*1000,d5=S.dem.nodes[d_5].pos[2]*1000,
      # displacement symmetry [mm]
      ds1=S.dem.nodes[d_1_s].pos[2]*1000,ds2=S.dem.nodes[d_2_s].pos[2]*1000,ds3=S.dem.nodes[d_3_s].pos[2]*1000,ds4=S.dem.nodes[d_4_s].pos[2]*1000,ds5=S.dem.nodes[d_5_s].pos[2]*1000, 
      # the top soil layer [m]
      h=S.dem.par[h_n].pos[2],
      # overall force on the plate [kN]
      t_f=total_force/1000,
      # the stress [kPa]
      s_1=stress_1/1000,
      s_2=stress_6/1000,
      s_3=stress_2/1000,
      s_4=stress_7/1000,
      s_5=stress_8/1000,
      s_6=stress_3/1000,
      s_7=stress_4/1000,
      s_8=stress_5/1000)
      # Save the txt
      # reaction force
      S.plot.saveDataTxt('support reaction',vars=('t','r1','r2','r3','r4'))
      # displ
      S.plot.saveDataTxt('displ',vars=('t','d1','d2','d3','d4','d5'))
      # displ_symmetry
      S.plot.saveDataTxt('displ_s',vars=('t','ds1','ds2','ds3','ds4','ds5'))
      # rock boulder height
      S.plot.saveDataTxt('rock boulder height',vars=('t','h'))
      # rockboulder impact force
      S.plot.saveDataTxt('impact force',vars=('t','uf'))
      # overall force on the plate
      S.plot.saveDataTxt('overall force',vars=('t','t_f'))
      # the stress increment on the plate 
      S.plot.saveDataTxt('stress increment',vars=('t','s_1','s_2','s_3','s_4'))
     


S.saveTmp()
import woo.qt
woo.qt.Controller()
#woo.qt.View() 

