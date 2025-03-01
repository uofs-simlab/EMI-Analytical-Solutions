from firedrake import *
import math
import csv

def format_scientific(value):
	return "{:.3e}".format(value)

def get_error(fname, n):
	path_save = "data/"
	error_list=[fname, n]
	fname_path= path_save+fname
	with CheckpointFile(fname_path, 'r') as afile:
		mesh = afile.load_mesh()
		timestepping_history = afile.get_timestepping_history(mesh, "w")
		num_times = len(timestepping_history['index'])
		w = afile.load_function(mesh, "w", idx=num_times-1)
		t = afile.get_attr("/times", str(num_times-1))

	#Mesh definitions
	N_cells = 4
	boundary = 2
	ec_space = 1
	cell1 = 3
	cells = [cell1 + i for i in range(N_cells)]
	boundary_map = {}
	cell_map = {}
	x, y, z = SpatialCoordinate(mesh)
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
	w_coef = [-1, -2, -2, -1]

	#For ueAprox
	ueAprox = w.sub(0)
	ueExact_expr = cos(4*pi*x)*cos(4*pi*y)*cos(4*pi*z)
	ueExact = Function(ueAprox.function_space()).interpolate(ueExact_expr)
	ue_error = ueAprox - ueExact
	L2_ue_error = sqrt(assemble(inner(ue_error, ue_error) * dx(ec_space)))
	error_list.append(format_scientific(L2_ue_error))
	
	for k in cells:
		#For uiAprox
		uiAprox = w.sub(1)[k-cell1]
		uiExact_expr = (1 + (k-cell1+1)*exp(-t))*cos(4*pi*x)*cos(4*pi*y)*cos(4*pi*z)
		uiExact = Function(ueAprox.function_space()).interpolate(uiExact_expr)
		ui_error = uiAprox - uiExact
		L2_ui_error = sqrt(assemble(inner(ui_error, ui_error) * dx(k)))
		error_list.append(format_scientific(L2_ui_error))

		# For vAprox
		vAprox = w.sub(3)[k-cell1]
		vExact_expr = (k-cell1+1)*exp(-t)*cos(4*pi*x)*cos(4*pi*y)*cos(4*pi*z)
		vExact = Function(ueAprox.function_space()).interpolate(vExact_expr)
		v_error = vAprox - vExact
		L2_v_error = sqrt(assemble(inner(v_error, v_error)* dS(k)))
		error_list.append(format_scientific(L2_v_error))

	# For wAprox
	for bound in gap_map:
		idx, celli, cell2 = gap_map[bound]
		wAprox = w.sub(4)[idx]
		wExact_expr = w_coef[idx]*exp(-t)*cos(4*pi*x)*cos(4*pi*y)*cos(4*pi*z)
		wExact = Function(ueAprox.function_space()).interpolate(wExact_expr)
		w_error = wAprox - wExact
		L2_w_error = sqrt(assemble(inner(w_error, w_error)* dS(bound)))
		error_list.append(format_scientific(L2_w_error))

	return error_list

h5names = {
	10:"exp3_Godunov_tf_1_dt_0.1_N_10.h5",
	20:"exp3_Godunov_tf_1_dt_0.05_N_20.h5",
	40:"exp3_Godunov_tf_1_dt_0.025_N_40.h5",
	80:"exp3_Godunov_tf_1_dt_0.0125_N_80.h5",
	160:"exp3_Godunov_tf_1_dt_0.00625_N_160.h5"
	}

csv_file = 'errors_exp3.csv'

data = [["File_name", "n", "e_ue", 
		"e_ui1", "e_ui2","e_ui3", "e_ui4", 
		"e_v1", "e_v2", "e_v3", "e_v4", 
		"e_w12", "e_w13", "e_w24", "e_w34"]]

for n, fname in h5names.items():
	data_i= get_error(fname, n)
	data.append(data_i)

with open(csv_file, mode='w', newline='') as file:
	writer = csv.writer(file)
	writer.writerows(data)

print("Error CSV file created successfully")
