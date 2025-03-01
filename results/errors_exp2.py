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
	boundary = 2
	ec_space = 1
	cell1 = 3
	cell2 = 4
	gap = 5
	gap_boundary = 8
	ic_boundary_cell1 = 6
	ic_boundary_cell2 = 7

	# Define the radius
	x, y, z = SpatialCoordinate(mesh)
	radius = (x**2 + y**2 + z**2)**0.5
	annular_region = conditional(And(radius >= 4, radius <= 5.0), 1.0, 0)

	ueExact_expr = (50*exp(-t/10)*cos(t))/radius
	ui1Exact_expr = (5/181) * (181 + 145 * exp(-t) + (2 * exp(-t/10) * ((905 + 18 * radius) * cos(t) + 20 * radius * sin(t))) / radius)
	ui2Exact_expr = (5/181) * (181 + 869 * exp(-t) + (exp(-t/10) * (2 * (905 + 18 * radius) * cos(t) + 40 * radius * sin(t))) / radius)
	v1Exact_expr = (5/181) * (181 + 145 * exp(-t) + 4 * exp(-t/10) * (9 * cos(t) + 10 * sin(t)))
	v2Exact_expr = (5/181) * (181 + 869 * exp(-t) + 4 * exp(-t/10) * (9 * cos(t) + 10 * sin(t)))
	w12Exact_expr = -20*exp(-t)

	#For ueAprox
	ueAprox = w.sub(0)
	ueExact = Function(ueAprox.function_space()).interpolate(ueExact_expr)
	ue_error = ueAprox - ueExact
	L2_ue_error = sqrt(assemble(inner(ue_error, ue_error) * dx(ec_space)))
	max_ue_error=Function(w.sub(0).function_space())
	annular_region1 = conditional(And(radius >= 5, radius <= 6.0), 1.0, 0)
	max_ue_error.interpolate(annular_region1*abs(ue_error))
	Linf_ue_error = max(abs(max_ue_error.dat.data))
	error_list.append(format_scientific(L2_ue_error))
	#error_list.append(format_scientific(Linf_ue_error))

	#For ui1Aprox
	ui1Aprox = w.sub(1)[0]
	ui1Exact = Function(ueAprox.function_space()).interpolate(ui1Exact_expr)
	ui1_error = ui1Aprox - ui1Exact
	L2_ui1_error = sqrt(assemble(annular_region * inner(ui1_error, ui1_error) * dx(cell1)))
	error_list.append(format_scientific(L2_ui1_error))

	#For ui2Aprox
	ui2Aprox = w.sub(1)[1]
	ui2Exact = Function(ueAprox.function_space()).interpolate(ui2Exact_expr)
	ui2_error = ui2Aprox - ui2Exact
	L2_ui2_error = sqrt(assemble(annular_region * inner(ui2_error, ui2_error) * dx(cell2)))
	error_list.append(format_scientific(L2_ui2_error))

	# For v1Aprox
	v1Aprox = w.sub(3)[0]
	v1Exact = Function(ueAprox.function_space()).interpolate(v1Exact_expr)
	v1_error = v1Aprox - v1Exact
	L2_v1_error = sqrt(assemble(inner(v1_error, v1_error)* dS(cell1)))
	error_list.append(format_scientific(L2_v1_error))

	# For v2Aprox
	v2Aprox = w.sub(3)[1]
	v1Exact = Function(ueAprox.function_space()).interpolate(v2Exact_expr)
	v2_error = v2Aprox - v1Exact
	L2_v2_error = sqrt(assemble(inner(v2_error, v2_error)* dS(cell2)))
	error_list.append(format_scientific(L2_v2_error))

	# For w12Aprox
	w12Aprox = w.sub(4)
	w12Exact = Function(w12Aprox.function_space()).interpolate(w12Exact_expr)
	w12_error = w12Aprox - w12Exact
	L2_w12_error = sqrt(assemble(annular_region * inner(w12_error, w12_error)* dS(gap)))
	max_w12_error=Function(w.sub(4).function_space())
	max_w12_error.interpolate(annular_region*abs(w12_error))
	Linf_w12_error = max(abs(max_w12_error.dat.data))
	error_list.append(format_scientific(L2_w12_error))
	#error_list.append(format_scientific(Linf_w12_error))

	return error_list

h5names = {
	10:"exp2_Godunov_tf_8_dt_0.775_N_10.h5",
	20:"exp2_Godunov_tf_8_dt_0.3875_N_20.h5",
	40:"exp2_Godunov_tf_8_dt_0.19375_N_40.h5",
	80:"exp2_Godunov_tf_8_dt_0.096875_N_80.h5",
	160:"exp2_Godunov_tf_8_dt_0.0484375_N_160.h5"
	}

csv_file = 'errors_exp2.csv'

data = [["File_name", "n", "e_ue", "e_u1", "e_u1", "e_v1", "e_v2", "e_w12"]]

for n, fname in h5names.items():
	data_i= get_error(fname, n)
	data.append(data_i)

with open(csv_file, mode='w', newline='') as file:
	writer = csv.writer(file)
	writer.writerows(data)

print("Error CSV file created successfully")
