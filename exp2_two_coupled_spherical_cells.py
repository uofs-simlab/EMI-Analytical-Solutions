from firedrake import *
from irksome import *
import fractional_step as fs
import time

start_time = time.time()

# Change to the target mesh
mesh = Mesh("meshes/two_coupled_cells_0.1.msh")

"""
------------------
Parameters
------------------
"""
# Cell model
v_rest = 5
Rm = 1
Cm = 1
C12 =1
Rgap =1
# Tissue model
sigmae = 1
sigmai = 1


"""
------------------
Function spaces
------------------
"""
Ue = FunctionSpace(mesh, 'CG', 1)
Ui = VectorFunctionSpace(mesh, 'CG', 1, dim=2)
I = VectorFunctionSpace(mesh, 'CG', 1, dim=3)
V = VectorFunctionSpace(mesh, 'CG', 1, dim=2)
Wi = FunctionSpace(mesh, 'CG', 1)
W = Ue * Ui * I * V * Wi
w = Function(W, name="w")
ue, ui, i, vn, wn = split(w)
t_ue, t_ui, t_i, t_v, t_w = TestFunctions(W)
ui1, ui2 = split(ui)
t_ui1, t_ui2 = split(t_ui)
v1, v2 = split(vn)
t_v1, t_v2 = split(t_v)
i1, i2, i12 = split(i)
t_i1, t_i2, t_i12 = split(t_i)

"""
------------------
Mesh definition
------------------
"""
boundary = 2
ec_space = 1
cell1 = 3
cell2 = 4
gap = 5
gap_boundary = 8
ic_boundary_cell1 = 6
ic_boundary_cell2 = 7

ti = 1/4
t = Constant(ti)
x, y, z = SpatialCoordinate(mesh)
radius = sqrt(x**2 + y**2 + z**2)

"""
------------------
EMI equations
------------------
"""
F_emi = inner(sigmae * grad(ue), grad(t_ue)) * dx(ec_space) - inner(i1, t_ue)('-') * dS(cell1) - inner(i2, t_ue)('-') * dS(cell2) + \
	inner(sigmai * grad(ui1), grad(t_ui1)) * dx(cell1) + inner(i1, t_ui1)('-')  * dS(cell1) + inner(i12, t_ui1)('-') * dS(gap) + \
	inner(sigmai * grad(ui2), grad(t_ui2)) * dx(cell2) + inner(i2, t_ui2)('-')  * dS(cell2) - inner(i12, t_ui2)('-') * dS(gap) + \
	((ui1 - ue - v1) * t_v1)('-') * dS(cell1) + \
	((ui2 - ue - v2) * t_v2)('-') * dS(cell2) + \
	((ui1 - ui2 - wn) * t_w)('-') * dS(gap) + \
	(Dt(v1) * t_i1)('-') *dS(cell1) - (i1 / Cm * t_i1)('-') * dS(cell1) + \
	(Dt(v2) * t_i2)('-') * dS(cell2) - (i2 / Cm * t_i2)('-') * dS(cell2) + \
	(Dt(wn) * t_i12)('-') * dS(gap) - (i12 / C12 * t_i12)('-') * dS(gap)

"""
------------------
Initial conditions
------------------
"""
# W
w.sub(4).assign(-20*exp(-ti))

# V1,V2
w.sub(3).interpolate(as_vector([
	(5/181) * (181 + 145 * exp(-ti) + 4 * exp(-ti/10) * (9 * cos(ti) + 10 * sin(ti))),
	(5/181) * (181 + 869 * exp(-ti) + 4 * exp(-ti/10) * (9 * cos(ti) + 10 * sin(ti)))
]))

# Im1, Im2, Im12
w.sub(2).interpolate(as_vector([2*(exp(-ti/10))*cos(ti), 2*(exp(-ti/10))*cos(ti), 0]))

# Ui1, Ui2
w.sub(1).interpolate(as_vector([
	(5/181) * (181 + 145 * exp(-ti) + (2 * exp(-ti/10) * ((905 + 18 * radius) * cos(ti) + 20 * radius * sin(ti))) / radius),
	(5/181) * (181 + 869 * exp(-ti) + (exp(-ti/10) * (2 * (905 + 18 * radius) * cos(ti) + 40 * radius * sin(ti))) / radius)
]))

#Ue
w.sub(0).interpolate((50*exp(-ti/10)*cos(ti))/radius)

"""
------------------
Boundary conditions
------------------
"""
#Uapp
r2_bc = (50*exp(-t/10)*cos(t))/6
bc = DirichletBC(W.sub(0), r2_bc, [boundary])

# Ui1, Ui2
r_inner = 3
r01_bc = as_vector([(5/181) * (181 + 145 * exp(-t) + (2 * exp(-t/10) * ((905 + 18 * r_inner) * cos(t) + 20 * r_inner * sin(t))) / r_inner), 0])
r02_bc = as_vector([0, (5/181) * (181 + 869 * exp(-t) + (exp(-t/10) * (2 * (905 + 18 * r_inner) * cos(t) + 40 * r_inner * sin(t))) / r_inner)])
bc1 = DirichletBC(W.sub(1), r01_bc, [ic_boundary_cell1])
bc2 = DirichletBC(W.sub(1), r02_bc, [ic_boundary_cell2])

#W 12
r05_bc = -20*exp(-t)
bc3 = DirichletBC(W.sub(4), r05_bc, [ic_boundary_cell1])
bc4 = DirichletBC(W.sub(4), r05_bc, [ic_boundary_cell2])

# J12
r03_bc = as_vector([0, 0, 0])
bc5 = DirichletBC(W.sub(2), r03_bc, [ic_boundary_cell1])
bc6 = DirichletBC(W.sub(2), r03_bc, [ic_boundary_cell2])


"""
------------------
Cell model
------------------
"""
def f_explicit(t, w):
	dw = Function(w)
	dw.assign(0)
	dw.sub(3).interpolate((-1/(Cm *Rm))*(w.sub(3)- as_vector([v_rest,v_rest])))
	dw.sub(4).assign((-1/(C12*Rgap))*w.sub(4))
	return dw

"""
------------------
OS implementation
------------------
"""

OS_method="Godunov"
integrator_1 = 'FE'
integrator_2 = BackwardEuler()

N = 160
tf = 8
dt = Constant((tf-ti)/N)
delta_t = (tf-ti)/N #For print purposes

path_save="results/data"
fname = path_save+"/exp2_"+OS_method+"_tf_"+str(tf)+"_dt_"+str(delta_t)+"_N_"+str(N)+".h5"

fs.fractional_step(
	[f_explicit, F_emi],
	dt, w, t, tf, OS_method, {(2,): integrator_2, (1,): integrator_1}, 
	solver_parameters={2: {'solver_parameters': {
		'snes_type': 'ksponly', 
		'ksp_rtol': 1e-7, 
		'pc_type': 'ilu', 
		'pc_factor_shift_type': 'nonzero', }, 
		'bcs': [bc, bc1, bc2, bc3, bc4, bc5, bc6]}},
	fname=fname)

end_time = time.time()
print(f"Time taken: {end_time - start_time:.6f} seconds")
