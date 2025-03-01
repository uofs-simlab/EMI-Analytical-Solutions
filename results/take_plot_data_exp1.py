import csv
import numpy as np
from firedrake import *

path_save = "data/"
point = (0, 5)
h5names = {
	7:"exp1_Godunov_tf_7_dt_0.9642857142857143_N_7.h5",
	14:"exp1_Godunov_tf_7_dt_0.48214285714285715_N_14.h5",
	28:"exp1_Godunov_tf_7_dt_0.24107142857142858_N_28.h5",
	56:"exp1_Godunov_tf_7_dt_0.12053571428571429_N_56.h5",
	112:"exp1_Godunov_tf_7_dt_0.060267857142857144_N_112.h5"
}

for n, fname in h5names.items():
	fname_path = path_save + fname
	with CheckpointFile(fname_path, 'r') as afile:
		mesh = afile.load_mesh()
		timestepping_history = afile.get_timestepping_history(mesh, "w")
		num_times = len(timestepping_history['index'])

		times_aprox = [afile.get_attr("/times", str(i)) for i in range(num_times)]
		ue_aprox_values = []
		ui_aprox_values = []
		v_aprox_values = []

		for i in range(num_times):
			w = afile.load_function(mesh, "w", idx=i)
			ue_aprox_values.append(w.sub(0)(point))
			ui_aprox_values.append(w.sub(1)(point))
			v_aprox_values.append(w.sub(3)(point))

		csv_filename = f"plot_data_exp1/data_n_{n}.csv"
		with open(csv_filename, mode='w', newline='') as csvfile:
			writer = csv.writer(csvfile)
			writer.writerow(["Time", "Ue", "Ui", "V"])  # Writing header
			for t, ue, ui, v in zip(times_aprox, ue_aprox_values, ui_aprox_values, v_aprox_values):
				writer.writerow([t, ue, ui, v])

		print(f"Data for n={n} saved to {csv_filename}")
