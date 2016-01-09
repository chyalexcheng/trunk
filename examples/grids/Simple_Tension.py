""" Author: Hongyang Cheng <chyalexcheng@gmail>
    Test #1: prepare membrane strucutre using chained cylinders
"""

# import modules
from yade import qt, plot
import numpy as np

#####################
##  Key controlls  ##
#####################

# type of membrane material
GSType = 'PP'
color = [245./255,245./255,245./255]
# global damping
damp = 0.2
width = 0.1
# assumed radius
rGrid = 5.e-3
# discretization per cylinder
nL = 5
# factor for greater GridCo-GridCo stiffness
stif = 1e1
# density
rho = 443.96

# tensile modulus of membrane material
if GSType == 'PE':
	maxStrainStress = [0.17035, 37.971e6]
	thick = .251e-3
if GSType == 'PP':
	maxStrainStress = [0.1425, 80.8e6]
	thick = .393e-3
young = maxStrainStress[1]/maxStrainStress[0]
sigTmax = maxStrainStress[1]

## material parameters for external behavior
# consider 10 cm width to compute 2D modulus and strength
young *= width
sigTmax *= width
# m2i: membrane-interface
E_m2i = stif*young; v_m2i = 0.33; phi_m2i = radians(34)

## material parameters for external behavior
# rescale 2D modulus and strength to fit with assumed radius, !!! obmit pi here (code bug)
young *= thick/(pi*rGrid**2)
sigTmax *= thick/(pi*rGrid**2)
# m: membrane, i:interface
E_m = young; v_m = 0.; phi_m = 0.; sigTmax_m = sigTmax; sigSmax_m = sigTmax
E_i = 0.   ; v_i = 0.; phi_i = 0.; sigTmax_i = sigTmax; sigSmax_i = sigTmax

#################
##  Functions  ##
#################

def addPlotData():
	T = []
	e = O.bodies[4].state.displ().norm()/O.bodies[4].state.refPos.norm()
	for i in range(4):
		inter = O.interactions[i,i+1]
		if inter.isReal:
			T.append(inter.phys.normalForce.norm()/thick/width)
		else:
			O.pause()
			print 'membrane broke...'
	plot.addData(e = e, T1 = T[0], T2 = T[1], T3 = T[2], T4 = T[3])

def gridNode(center,radius,dynamic=None,fixed=False,wire=False,color=None,highlight=False,material=-1):
	b=Body()
	b.shape=GridNode(radius=radius,color=color if color else utils.randomColor(),wire=wire,highlight=highlight)
	#V=(4./3)*math.pi*radius**3	# will be overriden by the connection
	V=0.
	geomInert=(2./5.)*V*radius**2	# will be overriden by the connection
	utils._commonBodySetup(b,V,Vector3(geomInert,geomInert,geomInert),material,pos=center,dynamic=dynamic,fixed=fixed)
	b.aspherical=False
	b.bounded=False
	b.mask=0	#avoid contact detection with the nodes. Manual interaction will be set for them in "gridConnection" below.
	return b

def gridConnection(id1,id2,radius,wire=False,color=None,highlight=False,material=-1,mask=1,cellDist=None):
	b=Body()
	b.shape=GridConnection(radius=radius,color=color if color else utils.randomColor(),wire=wire,highlight=highlight)
	sph1=O.bodies[id1] ; sph2=O.bodies[id2]
	i=createInteraction(id1,id2)
	nodeMat=sph1.material
	b.shape.node1=sph1 ; b.shape.node2=sph2
	sph1.shape.addConnection(b) ; sph2.shape.addConnection(b)
	if(O.periodic):
		if(cellDist!=None):
			i.cellDist=cellDist
		segt=sph2.state.pos + O.cell.hSize*i.cellDist - sph1.state.pos
	else: segt=sph2.state.pos - sph1.state.pos
	L=segt.norm()
	V=0.5*L*math.pi*radius**2
	geomInert=(2./5.)*V*radius**2
	utils._commonBodySetup(b,V,Vector3(geomInert,geomInert,geomInert),material,pos=sph1.state.pos,dynamic=False,fixed=True)
	sph1.state.mass = sph1.state.mass + V*nodeMat.density
	sph2.state.mass = sph2.state.mass + V*nodeMat.density
	for k in [0,1,2]:
		sph1.state.inertia[k] = sph1.state.inertia[k] + geomInert*nodeMat.density
		sph2.state.inertia[k] = sph2.state.inertia[k] + geomInert*nodeMat.density
	b.aspherical=False
	if O.periodic:
		i.phys.unp= -(sph2.state.pos + O.cell.hSize*i.cellDist - sph1.state.pos).norm() + sph1.shape.radius + sph2.shape.radius
		b.shape.periodic=True
		b.shape.cellDist=i.cellDist
	else:
		i.phys.unp= -(sph2.state.pos - sph1.state.pos).norm() + sph1.shape.radius + sph2.shape.radius	
	i.geom.connectionBody=b
	I=math.pi*(2.*radius)**4/64.
	E=nodeMat.young
	print i.phys.kn, i.phys.ks, i.phys.kr, i.phys.ktw
	i.phys.kn=E*math.pi*(radius**2)/L
	# use the correct normalAdhesion and shear Adhesion
	i.phys.normalAdhesion = nodeMat.normalCohesion*math.pi*(radius**2)
	i.phys.shearAdhesion = nodeMat.shearCohesion*math.pi*(radius**2)
	# do not use shear and rolling/bending stiffness
	#~ i.phys.kr=E*I/L
	#~ i.phys.ks=12.*E*I/(L**3)
	#~ G=E/(2.*(1+nodeMat.poisson))
	#~ i.phys.ktw=2.*I*G/L
	b.mask=mask
	print i.phys.kn, i.phys.ks, i.phys.kr, i.phys.ktw
	return b
	
###############
##  Engines  ##
###############

## Engines need to be defined first since the function gridConnection creates the interaction
O.engines=[
	ForceResetter(),
	InsertionSortCollider([
		Bo1_GridConnection_Aabb(),
	]),
	InteractionLoop(
		[Ig2_GridNode_GridNode_GridNodeGeom6D(),
		 Ig2_GridConnection_GridConnection_GridCoGridCoGeom(),
		 ],
		[Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(
			setCohesionNow=False,
			setCohesionOnNewContacts=False),
		 Ip2_FrictMat_FrictMat_FrictPhys(),
		 ],
		[Law2_ScGeom6D_CohFrictPhys_CohesionMoment(),
       Law2_GridCoGridCoGeom_FrictPhys_CundallStrack(),
       ]
	),
	NewtonIntegrator(gravity=(0,0,0),damping=damp,label='newton'),
	PyRunner(command='addPlotData()',iterPeriod=5000,label='plotter'),
]

#################
##  Materials  ##
#################

# interface
O.materials.append(CohFrictMat(young=E_i,poisson=v_i,density=rho,frictionAngle=phi_i \
	,normalCohesion=sigTmax_i,shearCohesion=sigSmax_i,momentRotationLaw=False,label='iMat'))
# membrane
O.materials.append(CohFrictMat(young=E_m,poisson=v_m,density=rho,frictionAngle=phi_m \
	,normalCohesion=sigTmax_m,shearCohesion=sigSmax_m,momentRotationLaw=False,label='mMat'))
# membrane to interface
O.materials.append(FrictMat(young=E_m2i,poisson=v_m2i,density=rho,frictionAngle=phi_m2i,label='m2iMat'))	

###################
##  Build model  ##
###################
# create nodes
dL = 0.06
start = Vector3(0,0,0)
ids = []
for i in range(5):
	ids.append(O.bodies.append(
		gridNode(start+Vector3(dL*i,0,0),rGrid,wire=False,fixed=False,material='mMat',color=[1,1,0])))
# create connection between membrane nodes
for i,j in zip( ids[:-1], ids[1:]):
	O.bodies.append(gridConnection(i,j,rGrid,wire=True,material='m2iMat',color=[1,1,0]))
## set boundary condition
# fix both ends
for i in [ids[0],ids[-1]]:
	O.bodies[i].state.blockedDOFs = 'xyzXYZ'
# apply constant velocity on the plus side
O.bodies[ids[-1]].state.vel = Vector3(0.01*dL,0,0)

# run with plot
plot.plots = {'e':('T1','T2','T3','T4')}
O.dt = 0.5*PWaveTimeStep()
O.saveTmp()
O.run()

