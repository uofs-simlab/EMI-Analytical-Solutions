from firedrake import *
from irksome import *
import fractional_step as fs
import time

start_time = time.time()

# Change to the target mesh
mesh = Mesh("meshes/manufactured_solution_0.035.msh")
mesh.init()

"""
------------------
Parameters
------------------
"""
#Cell model
Rm = 1
Cm = 1
C12 =1
Rgap =1
#Tissue model
sigmae = 4
sigmai = 1

"""
------------------
Mesh definition
------------------
"""
N_cells = 4
cell1 = 3
dirichlet_bc = 1
dirichlet_bc2 = 2

ti = 0
t = Constant(ti)
x, y, z = SpatialCoordinate(mesh)

## Determine gap junctions
cells = [cell1 + i for i in range(N_cells)]
extracellular = 1
boundary_map = {}
cell_map = {}
f = mesh.cell_to_facets

for i in cells:
    s = mesh.cell_subset(i)
    boundaries = set(f.data[s.indices][:,:,1].flatten())
    boundaries.discard(-1)
    for b in boundaries:
        if b not in boundary_map:
            boundary_map[b] = []
        boundary_map[b].append(i)
    cell_map[i] = boundaries

gap_map = {}
lowest_idx = 0
for bound in boundary_map:
    c = boundary_map[bound]
    if len(c) == 1:
        continue
    gap_map[bound] = (lowest_idx, min(c), max(c))
    lowest_idx += 1

n_cells = N_cells
n_gaps = len([x for x in boundary_map if len(boundary_map[x]) > 1])


"""
------------------
Function spaces
------------------
"""
Ue = FunctionSpace(mesh, 'CG', 1)
Ui = VectorFunctionSpace(mesh, 'CG', 1, dim=n_cells)
# Place all I1, I2, I3... in a single function since all of them are 
# equal to zero and reduce the computational demand.
I = FunctionSpace(mesh, 'CG', 1)
V = VectorFunctionSpace(mesh, 'CG', 1, dim=n_cells)
Wi = VectorFunctionSpace(mesh, 'CG', 1, dim=max(2, n_gaps))
W = MixedFunctionSpace([Ue, Ui, I, V, Wi])

w = Function(W, name="w")

ue, ui, i , vn, wn= split(w)

pe, p_i, t_i, t_v, t_w = TestFunctions(W)

ui = split(ui)
vi = split(vn)
wi = split(wn)

t_v = split(t_v)
p_i = split(p_i)
t_w = split(t_w)

"""
------------------
EMI equations (2)
------------------
"""
Fpe = inner(sigmae * grad(ue), grad(pe)) * dx(extracellular) - 192*(pi**2)*cos(4*pi*x)*cos(4*pi*y)*cos(4*pi*z) * pe * dx(extracellular)
for k in cells:
    Fpe -= inner(i, pe)('+') * dS(k)

Fp_i = 0
for k in cells:
    Fp_i += inner(sigmai * grad(ui[k - cell1]), grad(p_i[k - cell1])) * dx(k)
    Fp_i += inner(i, p_i[k-cell1])('+') * dS(k)
    Fp_i -= 48*(1 + (k - cell1+1)*exp(-t))*(pi**2)*cos(4*pi*x)*cos(4*pi*y)*cos(4*pi*z) * p_i[k - cell1] * dx(k) #Forcing terms

for bound in gap_map:
    idx, celli, cell2 = gap_map[bound]

    Fp_i += inner(i, p_i[celli-cell1])('+') * dS(bound) - inner(i, p_i[cell2 - cell1])('+') * dS(bound)

Fmem = 0
for k in cells:
    Fmem += ((ui[k-cell1] - ue - vi[k-cell1]) * t_v[k-cell1])('+') * dS(k)
    Fmem += (Dt(vi[k-cell1]) * t_i)('+') * dS(k) - (i / Cm * t_i)('+') * dS(k)

Fgap = 0
for bound in gap_map:
    idx, celli, cell2 = gap_map[bound]
    Fgap += ((ui[celli - cell1] - ui[cell2 - cell1] - wi[idx]) * t_w[idx])('+') * dS(bound)
    Fgap += (Dt(wi[idx]) * t_i)('+') * dS(bound) - (i / C12 * t_i)('+') * dS(bound)


"""
------------------
Boundary conditions
------------------
"""
bc = DirichletBC(W.sub(0), 0, [dirichlet_bc])
bc2 = DirichletBC(W.sub(0), 0, [dirichlet_bc2])

"""
------------------
Initial conditions
------------------
"""
w.sub(4).interpolate(as_vector([
    -exp(-ti)*cos(4*pi*x)*cos(4*pi*y)*cos(4*pi*z),
    -2*exp(-ti)*cos(4*pi*x)*cos(4*pi*y)*cos(4*pi*z),
    -2*exp(-ti)*cos(4*pi*x)*cos(4*pi*y)*cos(4*pi*z),
    -exp(-ti)*cos(4*pi*x)*cos(4*pi*y)*cos(4*pi*z)
]))  # [W12, ...]
w.sub(3).interpolate(as_vector([
    1*exp(-ti)*cos(4*pi*x)*cos(4*pi*y)*cos(4*pi*z),
    2*exp(-ti)*cos(4*pi*x)*cos(4*pi*y)*cos(4*pi*z),
    3*exp(-ti)*cos(4*pi*x)*cos(4*pi*y)*cos(4*pi*z),
    4*exp(-ti)*cos(4*pi*x)*cos(4*pi*y)*cos(4*pi*z)
]))  # [v1, ...]
w.sub(2).assign(0) #[I]
w.sub(1).interpolate(as_vector([
    (1 + 1*exp(-ti))*cos(4*pi*x)*cos(4*pi*y)*cos(4*pi*z),
    (1 + 2*exp(-ti))*cos(4*pi*x)*cos(4*pi*y)*cos(4*pi*z),
    (1 + 3*exp(-ti))*cos(4*pi*x)*cos(4*pi*y)*cos(4*pi*z),
    (1 + 4*exp(-ti))*cos(4*pi*x)*cos(4*pi*y)*cos(4*pi*z)
])) #[Ui1, ...]
w.sub(0).interpolate(cos(4*pi*x)*cos(4*pi*y)*cos(4*pi*z)) #Ue


"""
------------------
Cell model
------------------
"""
def f_explicit(t, w):
    dw = Function(w)
    dw.assign(0)
    dw.sub(3).assign((-1/(Cm *Rm))*w.sub(3))
    dw.sub(4).assign((-1/(C12*Rgap))*w.sub(4))
    return dw

"""
------------------
OS implementation
------------------
"""

N = 160
tf = 1
dt = Constant((tf-ti)/N)
delta_t = (tf-ti)/N #For print purposes

OS_method="Godunov"
integrator_1 = 'FE'
integrator_2 = BackwardEuler()

path_save="results/data"
fname = path_save+"/exp3_"+OS_method+"_tf_"+str(tf)+"_dt_"+str(delta_t)+"_N_"+str(N)+".h5"

fs.fractional_step(
    [f_explicit, Fpe + Fp_i + Fmem + Fgap], 
    dt, w, t, tf, OS_method, {(2,): integrator_2, (1,): integrator_1}, 
    solver_parameters={2: {'solver_parameters': {
        'snes_type': 'ksponly', 
        'ksp_rtol': 1e-7, 
        'pc_type': 'ilu', 
        'pc_factor_shift_type': 'nonzero'}, 
        'bcs': [bc, bc2]}},
        fname=fname)

end_time = time.time()
print(f"Time taken: {end_time - start_time:.6f} seconds")
