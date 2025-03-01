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
	internal_cell_boundary = 4

	# Define the radius
	x, y = SpatialCoordinate(mesh)
	radius = (x**2 + y**2)**0.5
	
	ueExact_expr = 5 + 10*ln(radius)*sin(t)
	uiExact_expr = 10+14*exp(-t)+cos(t)-sin(t)+10*ln(radius)*sin(t)
	vExact_expr = 5 + 14*exp(-t) + cos(t) - sin(t)
	iExact_expr = -2*sin(t)

	#For ueAprox
	ueAprox = w.sub(0)
	ueExact = Function(ueAprox.function_space()).interpolate(ueExact_expr)
	ue_error = ueAprox - ueExact
	L2_ue_error = sqrt(assemble(inner(ue_error, ue_error) * dx(ec_space)))
	error_list.append(format_scientific(L2_ue_error))

	#For uiAprox
	uiAprox = w.sub(1)
	uiExact = Function(uiAprox.function_space()).interpolate(uiExact_expr)
	ui_error = uiAprox - uiExact
	L2_ui_error = sqrt(assemble(inner(ui_error, ui_error) * dx(cell1)))
	error_list.append(format_scientific(L2_ui_error))

	# For vAprox
	vAprox = w.sub(3)
	vExact = Function(vAprox.function_space()).interpolate(vExact_expr)
	v_error = vAprox - vExact
	L2_v_error = sqrt(assemble(inner(v_error, v_error)* dS(cell1)))
	error_list.append(format_scientific(L2_v_error))

	return error_list


h5names = {
	7:"exp1_Godunov_tf_7_dt_0.9642857142857143_N_7.h5",
	14:"exp1_Godunov_tf_7_dt_0.48214285714285715_N_14.h5",
	28:"exp1_Godunov_tf_7_dt_0.24107142857142858_N_28.h5",
	56:"exp1_Godunov_tf_7_dt_0.12053571428571429_N_56.h5",
	112:"exp1_Godunov_tf_7_dt_0.060267857142857144_N_112.h5"
	}

csv_file = 'errors_exp1.csv'

data = [["File_name", "n", "e_ue", "e_u1", "e_v"]]

for n, fname in h5names.items():
	data_i= get_error(fname, n)
	data.append(data_i)

with open(csv_file, mode='w', newline='') as file:
	writer = csv.writer(file)
	writer.writerows(data)

print("Error CSV file created successfully")
