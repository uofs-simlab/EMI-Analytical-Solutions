from firedrake import *
from irksome import *
import fractional_step as fs
import time

start_time = time.time()

# Change to the target mesh
mesh = Mesh("meshes/single_circular_cell_0.10.msh")

"""
------------------
Parameters
------------------
"""
# Cell model
v_rest = 5
Rm = 1
Cm = 1
# Tissue model
sigmae = 1
sigmai = 1

"""
------------------
Function spaces
------------------
"""
Ue = FunctionSpace(mesh, 'CG', 1)
Ui = FunctionSpace(mesh, 'CG', 1)
I = FunctionSpace(mesh, 'CG', 1)
V = FunctionSpace(mesh, 'CG', 1)
W = Ue * Ui * I * V
w = Function(W, name="w")
ue, ui, i, vn = split(w)
t_ue, t_ui, t_i, t_v = TestFunctions(W)

"""
------------------
Mesh definition
------------------
"""
boundary = 2
ec_space = 1
cell1 = 3
ic_boundary = 4

ti = 1/4
t = Constant(ti)
x, y = SpatialCoordinate(mesh)
radius = sqrt(x**2 + y**2)

"""
------------------
EMI equations
------------------
"""
F_emi = inner(sigmae * grad(ue), grad(t_ue)) * dx(ec_space) - inner(i, t_ue)('-') * dS(cell1) + \
	inner(sigmai * grad(ui), grad(t_ui)) * dx(cell1) + inner(i, t_ui)('-') * dS(cell1) + \
	((ui - ue - vn) * t_v)('-') * dS(cell1) + \
	(Dt(vn) * t_i)('-') * dS(cell1) - ((i / Cm) * t_i)('-') * dS(cell1)

"""
------------------
Initial conditions
------------------
"""
w.sub(3).assign(5 + 14*exp(-ti) + cos(ti) - sin(ti)) #V
w.sub(2).assign(-2*sin(ti)) #i
w.sub(1).interpolate(10+14*exp(-ti)+cos(ti)-sin(ti)+10*ln(radius)*sin(ti)) #Ui
w.sub(0).interpolate(5+10*ln(radius)*sin(ti)) #Ue

"""
------------------
Boundary conditions
------------------
"""
# Ue
r2_bc = 5 + 10*ln(6)*sin(t)
bc = DirichletBC(W.sub(0), r2_bc, [boundary])

# Ui1, Ui2
r_inner = 3
r0_bc = 10+14*exp(-t)+cos(t)-sin(t)+10*ln(r_inner)*sin(t)
bc1 = DirichletBC(W.sub(1), r0_bc, [ic_boundary])

"""
------------------
Cell model
------------------
"""
def f_explicit(t, w):
	dw = Function(w)
	dw.assign(0)
	dw.sub(3).assign((-1/(Cm *Rm))*(w.sub(3)- v_rest))
	return dw

"""
------------------
OS implementation
------------------
"""

OS_method="Godunov"
integrator_1 = 'FE'
integrator_2 = BackwardEuler()

N = 112
tf = 7
dt = Constant((tf-ti)/N)
delta_t = (tf-ti)/N #For print purposes

path_save="results/data"
fname = path_save+"/exp1_"+OS_method+"_tf_"+str(tf)+"_dt_"+str(delta_t)+"_N_"+str(N)+".h5"

fs.fractional_step(
	[f_explicit, F_emi],
	dt, w, t, tf, OS_method, {(2,): integrator_2, (1,): integrator_1}, 
	solver_parameters={2: {'solver_parameters': {
		'snes_type': 'ksponly',
		'ksp_rtol': 1e-7,
		'pc_type': 'ilu',
		'pc_factor_shift_type': 'nonzero', }, 
		'bcs': [bc, bc1]}},
	fname=fname)

end_time = time.time()
print(f"Time taken: {end_time - start_time:.6f} seconds")
