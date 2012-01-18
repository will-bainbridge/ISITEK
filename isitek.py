#!/usr/bin/python

################################################################################

import numpy
import os
import cPickle as pickle
import scipy.misc
import scipy.sparse
import scipy.sparse.linalg
import scipy.special
import sys
import time

class Struct:
	def __init__(self, **keywords):
		self.__dict__.update(keywords)

class Timer(object):
	def __init__(self, name=None, multiline=False):
		self.name = name
		self.multiline = multiline
	def __enter__(self):
		self.start = time.time()
		if self.name:
			print '%s ...' % self.name ,
			if self.multiline:
				print
		sys.stdout.flush()
	def __exit__(self, type, value, traceback):
		if self.multiline:
			print ' ...' ,
		print 'done in %.3f s' % (time.time() - self.start)
		
################################################################################

def nodegrid(a,b):
	return [ x.T for x in numpy.meshgrid(a,b) ]

def dot_sequence(*args):
	if len(args) == 1: return args[0]
	else: return numpy.dot( args[0] , dot_sequence(*args[1:]) )

def string_multiple_replace(string,dict):
	for s,r in dict.iteritems():
		string = string.replace(s,r)
	return string

################################################################################

def read_data_file(data_filename):

	file = open(data_filename,'rb')
	data = pickle.load(file)
	file.close()

	node = data['node']
	face = data['face']
	element = data['element']
	boundary = data['boundary']
	u = data['u']
	order = data['order']

	return node,face,element,boundary,u,order

#------------------------------------------------------------------------------#

def read_input_file(input_filename):

	geometry_filename = []
	order = []
	boundary = []
	constants = []
	flux = []
	mesh_size = []
	wind = []
	initial = []

	file = open(input_filename,'r')

	for line in file.readlines():

		lineparts = line.split()

		if len(lineparts) >= 2 and lineparts[0] == 'geometry_filename':
			geometry_filename = lineparts[1]

		if len(lineparts) >= 2 and lineparts[0] == 'order':
			order = numpy.array([ int(x) for x in lineparts[1:] ])

		if len(lineparts) >= 4 and lineparts[0] == 'boundary':
			boundary.append(Struct(
				face = sum([ list(z) if len(z) == 1 else range(*z) for z in [ tuple( int(y) for y in x.split(':') ) for x in lineparts[1].split(',') ] ],[]) ,
				variable = int(lineparts[2]) ,
				condition = tuple(sum([ x == y for x in lineparts[3] ]) for y in 'xy') ,
				value = float(lineparts[4]) if len(lineparts) >= 5 else 0.0 ))

		if len(lineparts) >= 2 and lineparts[0] == 'initial':
			initial = lineparts[1:]

		if len(lineparts) >= 2 and lineparts[0] == 'constants':
			constants = lineparts[1]

		if len(lineparts) >= 6 and lineparts[0] == 'flux':
			flux.append(Struct(
				equation = int(lineparts[1]) ,
				variable = [ int(x) for x in lineparts[2].split(',') ] ,
				direction = 0*(lineparts[3] == 'x') + 1*(lineparts[3] == 'y') ,
				differential = [ tuple( sum([ x == y for x in z ]) for y in 'xy' ) for z in lineparts[4].split(',') ] ,
				power = [ int(x) for x in lineparts[5].split(',') ] ,
				constant = lineparts[6] ,
				method = lineparts[7] ))

		if len(lineparts) >= 2 and lineparts[0] == 'wind':
			wind = eval( 'lambda n,u,v:' + lineparts[1] , {'numpy':numpy} , {} )

		if len(lineparts) >= 2 and lineparts[0] == 'mesh_size':
			mesh_size = int(lineparts[1])

	file.close()

	if len(constants):
		constants = dict([ (y[0],float(y[1])) for y in [ x.split('=') for x in constants.split(';') ] ])
	else:
		constants = {}

	if len(flux):
		for i in range(0,len(flux)):
			flux[i].constant = eval(flux[i].constant,{},constants)

	if len(initial):
		replace = {'pi':'numpy.pi','cos(':'numpy.cos(','sin(':'numpy.sin('}
		for i in range(0,len(initial)):
			initial[i] = eval( 'lambda x,y:' + string_multiple_replace(initial[i],replace) , {'numpy':numpy} , constants )

	return geometry_filename,order,boundary,initial,flux,wind,mesh_size

#------------------------------------------------------------------------------#
		
def element_sequential_indices(e,element,face):
	
	n = len(element[e].face)

	polyline = numpy.array([ list(face[f].node) for f in element[e].face ])
	polynode = numpy.unique(polyline)
	ones = numpy.ones((n,1))

	connect = 1*(ones*polyline[:,0] == (ones*polynode).T) + 2*(ones*polyline[:,1] == (ones*polynode).T)

	side = [0]*n
	vertex = [0]*n
	for i in range(1,n):
		temp = connect[connect[:,side[i-1]] == (int(not vertex[i-1])+1),:].flatten() * (numpy.arange(0,n) != side[i-1])
		side[i] = temp.nonzero()[0][0]
		vertex[i] = temp[side[i]]-1

	return [side,vertex]

#------------------------------------------------------------------------------#

def read_geometry(geometry_filename):
		
	# read the geometry file
	file = open(geometry_filename,'r')
	data = file.readlines()
	file.close()

	# generate the mesh structures
	i = 0
	while i < len(data):
		if data[i].strip().split()[0] == 'NODES':
			nn = int(data[i].strip().split()[1])
			node = [ Struct(x=(0.0,0.0)) for _ in range(0,nn) ]
			for n in range(0,nn):
				node[n].x = tuple( [ float(x) for x in data[i+1+n].strip().split() ] )
			i += nn
		elif data[i].strip().split()[0] == 'FACES':
			nf = int(data[i].strip().split()[1])
			face = [ Struct(node=(0,0),border=[],size=1.0,normal=(0.0,0.0),centre=(0.0,0.0),boundary=[],Q=[]) for temp in range(0,nf) ]
			for f in range(0,nf):
				face[f].node = tuple( [ int(x) for x in data[i+1+f].strip().split() ] )
			i += nf
		elif data[i].strip().split()[0] == 'CELLS' or data[i].strip().split()[0] == 'ELEMENTS':
			ne = int(data[i].strip().split()[1])
			element = [ Struct(face=[],orientation=[],size=1.0,area=0.0,centre=(0.0,0.0),unknown=[],V=[],P=[],Q=[],W=[],X=[]) for temp in range(0,ne) ]
			for e in range(0,ne):
				element[e].face = [ int(x) for x in data[i+1+e].strip().split() ]
			i += ne
		else:
			i += 1

	# generate borders
	for e in range(0,ne):
		for f in element[e].face:
			face[f].border.append(e)

	# additional element geometry
	for e in range(0,ne):

		s,t = element_sequential_indices(e,element,face)

		index = [ face[element[e].face[i]].node[j] for i,j in zip(s,t) ]
		cross = [ node[index[i-1]].x[0]*node[index[i]].x[1]-node[index[i]].x[0]*node[index[i-1]].x[1] for i in range(0,len(element[e].face)) ]
		
		element[e].area = 0.5*sum(cross)
		element[e].centre = tuple([ sum([ (node[index[i-1]].x[j]+node[index[i]].x[j])*cross[i] for i in range(0,len(element[e].face)) ])/(6*element[e].area) for j in range(0,2) ])

		element[e].orientation = [ 2*t[i]-1 for i in s ]

		if element[e].area < 0.0:
			element[e].area = -element[e].area
			element[e].orientation = [ -x for x in element[e].orientation ]

		element[e].size = numpy.sqrt(element[e].area)

	# additional face geometry
	for f in range(0,nf):
		face[f].normal = ( -node[face[f].node[1]].x[1]+node[face[f].node[0]].x[1] , +node[face[f].node[1]].x[0]-node[face[f].node[0]].x[0] )
		face[f].size = 0.5*numpy.sqrt(numpy.dot(face[f].normal,face[f].normal))
		face[f].centre = tuple([ 0.5*(node[face[f].node[1]].x[i]+node[face[f].node[0]].x[i]) for i in range(0,2) ])

	# return
	return node,face,element

#------------------------------------------------------------------------------#

def assign_boundaries():
	nv = len(order)
	for f in range(0,len(face)):
		face[f].boundary = [ [] for v in range(0,nv) ]
	for b in range(0,len(boundary)):
		for f in boundary[b].face:
			face[f].boundary[boundary[b].variable].append(b)

#------------------------------------------------------------------------------#

def generate_unknowns():
	nv = len(order)
	np = order*(order+1)/2
	nu = 0
	for e in range(0,len(element)):
		element[e].unknown = [[]]*nv
		for v in range(0,nv):
			element[e].unknown[v] = range(nu,nu+np[v])
			nu += np[v]
	return numpy.zeros(nu)

#------------------------------------------------------------------------------#

def generate_constants(order):

	max_order = max(order)

	ng = 2*max_order-1
	gauss_locations,gauss_weights = [ x.real for x in scipy.special.orthogonal.p_roots(ng) ]

	#nh = 7
	#hammer_locations = numpy.array([
	#	[0.101286507323456,0.101286507323456],[0.797426958353087,0.101286507323456],[0.101286507323456,0.797426958353087],
	#	[0.470142064105115,0.470142064105115],[0.059715871789770,0.470142064105115],[0.470142064105115,0.059715871789770],
	#	[0.333333333333333,0.333333333333333]])
	#hammer_weights = 0.5 * numpy.array([
	#	0.125939180544827,0.125939180544827,0.125939180544827,0.132394152788506,0.132394152788506,0.132394152788506,
	#	0.225000000000000])

	#nh = 9
	#hammer_locations = numpy.array([
	#	[0.437525248383384,0.437525248383384],[0.124949503233232,0.437525248383384],[0.437525248383384,0.124949503233232],
	#	[0.165409927389841,0.037477420750088],[0.037477420750088,0.165409927389841],[0.797112651860071,0.165409927389841],
	#	[0.165409927389841,0.797112651860071],[0.037477420750088,0.797112651860071],[0.797112651860071,0.037477420750088]])
	#hammer_weights = 0.5 * numpy.array([
	#	0.205950504760887,0.205950504760887,0.205950504760887,0.063691414286223,0.063691414286223,0.063691414286223,
	#	0.063691414286223,0.063691414286223,0.063691414286223])

	nh = 12
	hammer_locations = numpy.array([
		[0.063089014491502,0.063089014491502],[0.873821971016996,0.063089014491502],[0.063089014491502,0.873821971016996],
		[0.249286745170910,0.249286745170910],[0.501426509658179,0.249286745170910],[0.249286745170910,0.501426509658179],
		[0.310352451033785,0.053145049844816],[0.053145049844816,0.310352451033785],[0.636502499121399,0.310352451033785],
		[0.310352451033785,0.636502499121399],[0.053145049844816,0.636502499121399],[0.636502499121399,0.053145049844816]])
	hammer_weights = 0.5 * numpy.array([
		0.050844906370207,0.050844906370207,0.050844906370207,0.116786275726379,0.116786275726379,0.116786275726379,
		0.082851075618374,0.082851075618374,0.082851075618374,0.082851075618374,0.082851075618374,0.082851075618374])

	taylor_coefficients = numpy.array([])
	taylor_powers = numpy.zeros((0,2),dtype=int)
	for i in range(0,2*max_order):
		taylor_coefficients = numpy.append(taylor_coefficients,scipy.misc.comb(i*numpy.ones(i+1),numpy.arange(0,i+1))/scipy.misc.factorial(i))
		taylor_powers = numpy.append(taylor_powers,numpy.array([range(i,-1,-1),range(0,i+1)],dtype=int).T,axis=0)

	powers_taylor = numpy.zeros((2*max_order,2*max_order),dtype=int)
	for i in range(0,taylor_powers.shape[0]): powers_taylor[taylor_powers[i][0]][taylor_powers[i][1]] = i
		
	factorial = scipy.misc.factorial(numpy.arange(0,2*max_order))

	return gauss_locations,gauss_weights,hammer_locations,hammer_weights,taylor_coefficients,taylor_powers,powers_taylor,factorial

#------------------------------------------------------------------------------#

def basis(x,y,element,n,differential):

	if taylor_powers[n,0] < differential[0] or taylor_powers[n,1] < differential[1]:
		return numpy.zeros(x.shape)
	p = taylor_powers[n]
	q = taylor_powers[n]-differential

	constant = taylor_coefficients[n] / numpy.power( element.size , sum(p) )
	constant = constant * factorial[p[0]] * factorial[p[1]] / ( factorial[q[0]] * factorial[q[1]] )

	return constant * numpy.power(x-element.centre[0],q[0]) * numpy.power(y-element.centre[1],q[1])

#------------------------------------------------------------------------------#

def derivative_transform_matrix(A,order):

	n = order*(order+1)/2

	D = numpy.zeros((n,n))
	D[0,0] = 1.0

	for i in range(0,order-1):

		old = numpy.nonzero(numpy.sum(taylor_powers,axis=1) == i)[0]
		temp = numpy.append( taylor_powers[old,:] + [1,0] , taylor_powers[old[taylor_powers[old,0] == 0],:] + [0,1] , axis=0 )
		new = powers_taylor[temp[:,0],temp[:,1]]
		index = nodegrid(old,old)

		D[nodegrid(new,new)] = numpy.append(
				A[0,0] * numpy.append( D[index] , numpy.zeros((i+1,1)) , axis=1 ) +
				A[0,1] * numpy.append( numpy.zeros((i+1,1)) , D[index] , axis=1 ) ,
				A[1,0] * numpy.append( D[old[-1],[old]] , [[0]] , axis=1 ) +
				A[1,1] * numpy.append( [[0]] , D[old[-1],[old]] , axis=1 ) , axis=0 )

	return D

#------------------------------------------------------------------------------#

def calculate_element_matrices():

	nf = len(face)
	ne = len(element)
	nv = len(order)

	max_order = max(order)

	np = numpy.array([ len(x) for x in element[0].unknown ])
	max_np = max(np)

	ng = len(gauss_weights)
	nh = len(hammer_weights)

	# initialise
	if do.pre:
		for e in range(0,ne):
			element[e].V = numpy.zeros((max_np,max_np))
			element[e].P = numpy.zeros((max_np,(len(element[e].face)-2)*nh,max_np))
			element[e].Q = [ numpy.zeros((ng,max_np)) for i in range(0,len(element[e].face)) ]
			element[e].W = numpy.zeros((len(element[e].face)-2)*nh)
			element[e].X = numpy.zeros(((len(element[e].face)-2)*nh,2))
		for f in range(0,nf):
			face[f].Q = [ [] for v in range(0,nv) ]

	# element vandermonde matrices
	if do.pre:
		for e in range(0,ne):
			for i in range(0,max_np):
				for j in range(0,max_np):
					element[e].V[i,j] = basis(numpy.array(element[e].centre[0]),numpy.array(element[e].centre[1]),element[e],i,taylor_powers[j])

	# triangulation and element area quadrature
	for e in range(0,ne):

		# triangulate
		nt = len(element[e].face)-2
		v = numpy.zeros((nt,3),dtype=int)
		v[:,0] = face[element[e].face[0]].node[0]
		j = 0
		for i in range(0,len(element[e].face)):
			f = element[e].face[i]
			o = int(element[e].orientation[i] < 0)
			v[j][1:] = numpy.array(face[f].node)[[1-o,o]]
			j += not any(v[j][1:] == v[j][0])
			if j >= nt: break

		# integration locations in and area of the triangles
		element[e].X = numpy.zeros(((len(element[e].face)-2)*nh,2))
		area = numpy.zeros(nt)
		for i in range(0,nt):
			d = numpy.array([ [ node[v[i][j]].x[k] - node[v[i][0]].x[k] for k in range(0,2) ] for j in range(1,3) ])
			element[e].X[i*nh:(i+1)*nh] = ( numpy.ones((nh,1))*node[v[i][0]].x +
					hammer_locations[:,0][numpy.newaxis].T*d[0] +
					hammer_locations[:,1][numpy.newaxis].T*d[1] )
			area[i] = numpy.cross(d[0,:],d[1,:])

		# integration weights
		element[e].W = (numpy.array([area]).T*hammer_weights).flatten()

	# element FEM numerics matrices
	if do.pre:
		for e in range(0,ne):

			# basis function values at the integration points
			for i in range(0,max_np):
				for j in range(0,max_np):
					element[e].P[i][:,j] = basis(element[e].X[:,0],element[e].X[:,1],element[e],j,taylor_powers[i])

	# element DG-FEM numerics matrices
	if do.pre:
		for e in range(0,ne):
			for i in range(0,len(element[e].face)):

				f = element[e].face[i]
			
				# integration locations along the face
				temp = gauss_locations[numpy.newaxis].T
				x = 0.5*(1.0-temp)*node[face[f].node[0]].x + 0.5*(1.0+temp)*node[face[f].node[1]].x

				# basis function values at the integration points
				for j in range(0,max_np):
					element[e].Q[i][:,j] = basis(x[:,0],x[:,1],element[e],j,[0,0])

	# face IDG-FEM numerics matrices
	for f in range(0,nf):

		# adjacent element and boundaries
		a = numpy.array(face[f].border)
		na = len(a)
		b = numpy.array(face[f].boundary,dtype=object)
		nb = [ len(i) for i in b ]

		if do.pre or (do.re and any(b)):

			# rotation to face coordinates
			R = numpy.array([[face[f].normal[0],face[f].normal[1]],[-face[f].normal[1],face[f].normal[0]]])
			R /= numpy.sqrt(numpy.dot(face[f].normal,face[f].normal))

			# face locations
			x = 0.5*(1.0-gauss_locations[numpy.newaxis].T)*node[face[f].node[0]].x + 0.5*(1.0+gauss_locations[numpy.newaxis].T)*node[face[f].node[1]].x
			y = face[f].centre + numpy.dot( x - face[f].centre , R.T )
			w = gauss_weights

			# adjacent integration locations
			xa = [ element[a[i]].X for i in range(0,na) ]
			ya = [ face[f].centre + numpy.dot( xa[i] - face[f].centre , R.T ) for i in range(0,na) ]
			wa = numpy.append(element[a[0]].W,element[a[1]].W) if na == 2 else element[a[0]].W

			for v in range(0,nv):

				# face basis indices
				temp = nodegrid(range(0,2*order[v]),range(0,2*order[v])) # NOTE # not sufficient for boundary faces with 2 bordering elements
				face_taylor = powers_taylor[ numpy.logical_and( temp[0] + na*temp[1] < na*order[v] + nb[v] , temp[1] < order[v] ) ]

				# number of interpolation unknowns
				ni = len(face_taylor)

				# matrices
				P = numpy.zeros((na*nh,na*np[v]))
				for j in range(0,np[v]):
					for k in range(0,na):
						P[k*nh:(1+k)*nh,j+k*np[v]] = basis(xa[k][:,0],xa[k][:,1],element[a[k]],j,[0,0])

				Q = numpy.zeros((na*nh,ni))
				for j in range(0,ni):
					for k in range(0,na):
						Q[k*nh:(k+1)*nh,j] = basis(ya[k][:,0],ya[k][:,1],face[f],face_taylor[j],[0,0])

				A = dot_sequence( P.T , numpy.diag(wa) , Q )
				B = dot_sequence( P.T , numpy.diag(wa) , P )

				# boundary parts
				if nb[v]:

					dA = numpy.zeros((nb[v]*order[v],ni))
					for i in range(0,nb[v]):
						for j in range(0,ni):
							for k in range(0,order[v]):
								dA[k+i*order[v],j] = basis(
										numpy.array(face[f].centre[0]),
										numpy.array(face[f].centre[1]),
										face[f],face_taylor[j],[0,k])
								# TODO # change [0,k] here for different BCs

					dB = numpy.zeros((nb[v]*order[v],nb[v]))
					for i in range(0,nb[v]): dB[i*order[v],i] = 1.0

					A = numpy.append( A , dA , axis=0 )
					B = numpy.append( numpy.append( B , numpy.zeros((B.shape[0],nb[v])) , axis=1 ) ,
							numpy.append( numpy.zeros((nb[v]*order[v],B.shape[1])) , dB , axis=1 ) ,
							axis=0 )

				# solve interpolation problem
				D = numpy.linalg.solve(A,B)

				# interpolated values
				F = numpy.zeros((ng,ni))
				face[f].Q[v] = numpy.zeros((np[v],ng,D.shape[1]))
				for j in range(0,np[v]):
					for k in range(0,ni):
						F[:,k] = basis(y[:,0],y[:,1],face[f],face_taylor[k],taylor_powers[j])
					face[f].Q[v][j] = numpy.dot( F , D )

				# transform differentials to x and y
				T = derivative_transform_matrix(numpy.linalg.inv(R),order[v])
				for j in range(0,ng): face[f].Q[v][:,j] = numpy.dot( T , face[f].Q[v][:,j] )

#------------------------------------------------------------------------------#

def initialise_unknowns():

	ne = len(element)
	np = [ len(x) for x in element[0].unknown ]
	nv = len(order)

	max_order = max(order)
	max_order_sq = max_order*max_order
	max_np = max(np)

	for e in range(0,ne):

		x = element[e].centre

		delta = numpy.linspace(-0.1*element[e].size/2,0.1*element[e].size/2,max_order)

		dx = [ temp.flatten() for temp in nodegrid(delta,delta) ]

		p = [ taylor_powers[0:max_np,i] for i in range(0,2) ]

		M = ((numpy.ones((max_np,1)) * dx[0]).T ** (numpy.ones((max_order_sq,1)) * p[0]) *
				(numpy.ones((max_np,1)) * dx[1]).T ** (numpy.ones((max_order_sq,1)) * p[1]) *
				(numpy.ones((max_order_sq,1)) * (scipy.misc.comb(p[0]+p[1],p[0])/scipy.misc.factorial(p[0]+p[1]))))

		inv_M = numpy.linalg.pinv(M)
		inv_V = numpy.linalg.inv(element[e].V)

		for v in range(0,nv):
			u[element[e].unknown[v]] = dot_sequence( inv_V , inv_M , initial[v](dx[0]+x[0],dx[1]+x[1]) )[0:np[v]]

#------------------------------------------------------------------------------#

def generate_system():

	ne = len(element)
	nf = len(face)
	ng = len(gauss_weights)
	nh = len(hammer_weights)
	np = [ len(x) for x in element[0].unknown ]
	nv = len(order)

	max_np = max(np)
	sum_np = sum(np)
	sum_np_sq = sum_np*sum_np

	# local dense jacobian
	L = Struct(i=[],x=[])

	# csr system jacobian
	J = Struct(np=0,p=[],ni=0,i=[],x=[])
	for e in range(0,ne): J.ni += sum_np_sq
	for f in range(0,nf): J.ni += (len(face[f].border) == 2)*2*sum_np_sq
	J.p = numpy.zeros(u.shape[0]+1,dtype=int)
	J.i = numpy.zeros(J.ni,dtype=int)
	J.x = numpy.zeros(J.ni)
	J.ni = 0

	# function vector
	F = numpy.zeros(u.shape)

	# differentiating matrices
	D = [ numpy.zeros((max_np,max_np)) for i in range(0,max_np) ]
	for i in range(0,max_np):
		index = ( (taylor_powers[0:max_np,0] - taylor_powers[i,0] >= 0) & (taylor_powers[0:max_np,1] - taylor_powers[i,1] >= 0) )
		D[i][0:sum(index),index] = numpy.eye(sum(index))
	
	for e in range(0,ne):

		# number of sides
		ns = len(element[e].face)

		# adjacent elements
		adj = - numpy.ones(ns,dtype=int)
		for i in range(0,ns):
			temp = numpy.array(face[element[e].face[i]].border)
			temp = temp[temp != e]
			if len(temp): adj[i] = temp[0]
		n_adj = sum(adj >= 0)
		i_adj = numpy.arange(0,ns)[adj >= 0]

		# local matrices to add to the system
		L.i = numpy.zeros((sum_np,(1+n_adj)*sum_np),dtype=int)
		L.i[:,0:sum_np] = numpy.tile( sum(element[e].unknown,[]) , (sum_np,1) )
		for i in range(0,n_adj): L.i[:,(i+1)*sum_np:(i+2)*sum_np] = numpy.tile( sum(element[adj[i_adj[i]]].unknown,[]) , (sum_np,1) )
		L.x = numpy.zeros(L.i.shape)

		# indices into the local matrices
		index_e = [ numpy.arange(sum(np[:v]),sum(np[:v+1]))[numpy.newaxis] for v in range(0,nv) ]
		index_a = [ [] for i in range(0,ns) ]
		for i in range(0,n_adj):
			index_a[i_adj[i]] = [ numpy.array([
				range(sum(np[:v]),sum(np[:v+1])) +
				range((i+1)*sum_np+sum(np[:v]),(i+1)*sum_np+sum(np[:v+1])) ])
				for v in range(0,nv) ]

		# mass matrices
		M = dot_sequence( element[e].P[powers_taylor[0,0]].T , numpy.diag(element[e].W) , element[e].P[powers_taylor[0,0]] )
		inv_M = [ numpy.linalg.inv(M[0:np[v],0:np[v]]) for v in range(0,nv) ]

		# loop over fluxes
		for term in flux:

			nt = len(term.variable)

			# numbers of equations
			npe = np[term.equation]

			# direction index
			direction = powers_taylor[1-term.direction,term.direction]

			# powers
			P = numpy.array(term.power)[numpy.newaxis].T

			# equation matrix
			A = - term.constant * dot_sequence( inv_M[term.equation] , element[e].P[direction][:,0:npe].T , numpy.diag(element[e].W) )

			# calculate the coefficients and values
			B = [ [] for t in range(0,nt) ]
			X = numpy.zeros((nt,nh))
			for t,v,differential in zip(range(0,nt),term.variable,term.differential):
				B[t] = element[e].P[powers_taylor[differential]][:,0:np[v]]
				X[t,:] = numpy.dot( B[t] , u[element[e].unknown[v]] )

			# add to the local jacobian
			Y = X ** P
			for t,v in zip(range(0,nt),term.variable):
				temp = numpy.copy(Y)
				temp[t,:] = P[t] * X[t,:] ** (P[t]-1)
				L.x[index_e[term.equation].T,index_e[v]] += dot_sequence( A , numpy.diag(numpy.prod(temp,axis=0)) , B[t] )

			# add to the function vector
			F[element[e].unknown[term.equation]] += numpy.dot( A , numpy.prod(Y,axis=0) )

			# face components
			for i in range(0,ns):

				f = element[e].face[i]
				a = adj[i]
				b = numpy.array(face[f].boundary,dtype=object)

				# face normal
				normal = element[e].orientation[i] * numpy.array(face[f].normal)

				# corresponding face index
				if a >= 0: j = numpy.arange(0,len(element[a].face))[numpy.array(element[a].face) == f]

				# wind
				if a >= 0 and ('u' in term.method):
					ui = [ dot_sequence( gauss_weights , element[e].Q[i][:,0:np[v]] , u[element[e].unknown[v]] ) for v in range(0,nv) ]
					uo = [ dot_sequence( gauss_weights , element[a].Q[j][:,0:np[v]] , u[element[a].unknown[v]] ) for v in range(0,nv) ]
					w = wind( normal , ui , uo )
				else:
					w = True

				# equation matrix
				A = normal[term.direction] * term.constant * dot_sequence( inv_M[term.equation] ,
						element[e].Q[i][:,0:npe].T , numpy.diag(0.5*gauss_weights) )

				# calculate the coefficients and values
				B = [ [] for t in range(0,nt) ]
				X = numpy.zeros((nt,ng))
				for t,v,differential,method in zip(range(0,nt),term.variable,term.differential,term.method):

					# where there is an adjacent element
					if a >= 0:

						# interpolated flux
						if method == 'i' or len(b[v]):
							if face[f].border[0] == e: temp = numpy.array(range(0,2*np[v]))
							else: temp = numpy.array(range(np[v],2*np[v])+range(0,np[v]))
							B[t] = face[f].Q[v][powers_taylor[differential]][:,temp]

						# averaged flux
						elif method == 'a':
							B[t] = 0.5*numpy.append(element[e].Q[i][:,0:np[v]],element[a].Q[j][:,0:np[v]],axis=1)

						# upwind flux
						elif method == 'u':

							B[t] = numpy.zeros((ng,2*np[v]))
							if w: B[t][:,0:np[v]] += element[e].Q[i][:,0:np[v]]
							else: B[t][:,np[v]:2*np[v]] += element[a].Q[j][:,0:np[v]]

						# values
						X[t,:] = numpy.dot( B[t] , numpy.append(u[element[e].unknown[v]],u[element[a].unknown[v]]) )

					# interpolated flux where there is no adjacent element
					else:
						B[t] = face[f].Q[v][powers_taylor[differential]][:,0:np[v]]
						X[t,:] = numpy.dot( B[t] , u[element[e].unknown[v]] )

					# interpolated flux at boundaries
					if len(b[v]):
						for k in range(0,len(b[v])):
							X[t,:] += boundary[b[v][k]].value * face[f].Q[v][powers_taylor[differential]][:,(1+(a>=0))*np[v]+k]
				
				# add to the local jacobian
				Y = X ** P
				for t,v in zip(range(0,nt),term.variable):
					temp = numpy.copy(Y)
					temp[t,:] = P[t] * X[t,:] ** (P[t]-1)
					L.x[index_e[term.equation].T,index_a[i][v] if a >= 0 else index_e[v]] += dot_sequence(
							A , numpy.diag(numpy.prod(temp,axis=0)) , B[t] )

				# add to the function vector
				F[element[e].unknown[term.equation]] += numpy.dot( A , numpy.prod(Y,axis=0) )

		# add dense local jacobian to csr global jacobian
		J.p[J.np:J.np+L.i.shape[0]] = numpy.arange(J.ni,J.ni+L.i.size,L.i.shape[1])
		J.np += L.i.shape[0]
		J.i[J.ni:J.ni+L.i.size] = L.i.flatten()
		J.x[J.ni:J.ni+L.i.size] = L.x.flatten()
		J.ni += L.i.size

	# return the global system
	J.p[J.np] = J.ni
	return [ scipy.sparse.csr_matrix((J.x,J.i,J.p)) , F ]

#------------------------------------------------------------------------------#

def write_display_file(display_filename,n):

	nv = len(order)
	np = numpy.array([ len(x) for x in element[0].unknown ])

	Q = numpy.linalg.inv(numpy.array([[1,-1,-1,1],[1,1,-1,-1],[1,1,1,1],[1,-1,1,-1]]))
		
	file = open(display_filename,'w')
	
	for e in range(0,len(element)):
		
		s,t = element_sequential_indices(e,element,face)
		
		for i in range(0,len(element[e].face)):
			
			quad = numpy.array( [ element[e].centre ,
				face[element[e].face[s[i-1]]].centre ,
				node[face[element[e].face[s[i]]].node[t[i]]].x ,
				face[element[e].face[s[i]]].centre ] )

			a = numpy.dot(Q,quad)
			
			mesh = numpy.append( numpy.mgrid[0:n+1,0:n+1]*(2.0/n)-1.0 , numpy.zeros((nv,n+1,n+1)) , axis=0 )
			mesh[0:2] = [ a[0,j] + a[1,j]*mesh[0] + a[2,j]*mesh[1] + a[3,j]*mesh[0]*mesh[1] for j in range(0,2) ]

			for j in range(0,max(np)):
				phi = basis(mesh[0],mesh[1],element[e],j,[0,0])
				for v in numpy.arange(0,nv)[j < np]:
					mesh[2+v] += u[element[e].unknown[v][j]]*phi

			file.write( '\n\n'.join([ '\n'.join([ ' '.join(['%e']*(2+nv)) % tuple(mesh[:,i,j]) for j in range(0,n+1) ]) for i in range(0,n+1) ]) + '\n\n\n' )
			
	file.close()

#------------------------------------------------------------------------------#

def write_data_file(data_filename):
	file = open(data_filename,'wb')
	pickle.dump({'node':node,'face':face,'element':element,'boundary':boundary,'order':order,'u':u},file,protocol=pickle.HIGHEST_PROTOCOL)
	file.close()
	
################################################################################
	
path = sys.argv[1]
action = sys.argv[2].lower()

directory = os.path.dirname(path)
name = os.path.basename(path)

input_filename = directory + os.sep + name + '.input'
data_filename = directory + os.sep + name + '.data'
display_filename = directory + os.sep + name + '.display'

do = Struct(pre = 'p' in action , re = 'r' in action , init = 'i' in action , solve = 's' in action , display = 'd' in action )

#------------------------------------------------------------------------------#

if not do.pre:
	with Timer('reading data from "%s"' % data_filename):
		node,face,element,boundary,u,order = read_data_file(data_filename)

with Timer('reading input from "%s"' % input_filename):
	input_data = read_input_file(input_filename)
	if do.pre:
		geometry_filename = directory + os.sep + input_data[0]
		order = input_data[1]
	if do.pre or do.re:
		boundary = input_data[2]
	if do.init:
		initial = input_data[3]
	if do.solve:
		for i in range(0,len(boundary)): boundary[i].value = input_data[2][i].value
		flux = input_data[4]
		wind = input_data[5]
	if do.display:
		mesh_size = input_data[6]

with Timer('generating constants'):
	(gauss_locations,gauss_weights,
			hammer_locations,hammer_weights,
			taylor_coefficients,taylor_powers,powers_taylor,
			factorial) = generate_constants(order)

if do.pre:
	with Timer('reading and processing geometry from "%s"' % geometry_filename):
		node,face,element = read_geometry(geometry_filename)
	with Timer('generating unknowns'):
		u = generate_unknowns()

if do.pre or do.re:
	with Timer('assigning boundaries to faces'):
		assign_boundaries()
	with Timer('calculating element matrices'):
		calculate_element_matrices()

if do.init:
	with Timer('initialising the unknowns'):
		initialise_unknowns()

if do.solve:

	with Timer('iterating',True):

		ne = len(element)
		nv = len(order)

		index = [ numpy.zeros(u.shape,dtype=bool) for v in range(0,nv) ]
		for e in range(0,ne):
			for v in range(0,nv):
				index[v][element[e].unknown[v]] = True

		unit = numpy.zeros(u.shape)
		#for e in range(0,ne): unit[element[e].unknown[0][0]] = 1.0

		for iteration in range(0,1):
			J,f = generate_system()
			print '    ' + ' '.join([ '%.4e' % numpy.max(numpy.abs(f[i])) for i in index ])
			u += scipy.sparse.linalg.spsolve(J,-f-unit)

if do.display:
	with Timer('saving display to "%s"' % display_filename):
		write_display_file(display_filename,mesh_size)

if do.pre or do.re or do.init or do.solve:
	with Timer('saving data to "%s"' % data_filename):
		write_data_file(data_filename)

################################################################################
