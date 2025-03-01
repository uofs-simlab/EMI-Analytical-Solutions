import csv
import numpy as np
from firedrake import *

path_save = "data/"

h5names = {
	10:"exp2_Godunov_tf_8_dt_0.775_N_10.h5",
	20:"exp2_Godunov_tf_8_dt_0.3875_N_20.h5",
	40:"exp2_Godunov_tf_8_dt_0.19375_N_40.h5",
	80:"exp2_Godunov_tf_8_dt_0.096875_N_80.h5",
	160:"exp2_Godunov_tf_8_dt_0.0484375_N_160.h5"
}

for n, fname in h5names.items():
	fname_path = path_save + fname
	with CheckpointFile(fname_path, 'r') as afile:
		mesh = afile.load_mesh()
		timestepping_history = afile.get_timestepping_history(mesh, "w")
		num_times = len(timestepping_history['index'])

		times_aprox = [afile.get_attr("/times", str(i)) for i in range(num_times)]
		ue_aprox_values = []
		ui1_aprox_values = []
		ui2_aprox_values = []
		v1_aprox_values = []
		v2_aprox_values = []
		w12_aprox_values = []

		for i in range(num_times):
			w = afile.load_function(mesh, "w", idx=i)
			ue_aprox_values.append(w.sub(0).at([5, 0, 0]))
			ui1_aprox_values.append(w.sub(1).at([-5, 0, 0])[0])
			ui2_aprox_values.append(w.sub(1).at([5, 0, 0])[1])
			v1_aprox_values.append(w.sub(3).at([-5, 0, 0])[0])
			v2_aprox_values.append(w.sub(3).at([5, 0, 0])[1])
			w12_aprox_values.append(w.sub(4).at([0, 4.5,0]))

		csv_filename = f"plot_data_exp2/data_n_{n}.csv"
		with open(csv_filename, mode='w', newline='') as csvfile:
			writer = csv.writer(csvfile)
			writer.writerow(["Time", "Ue", "Ui1", "Ui2", "V1", "V2", "W12"])  # Writing header
			for t, ue, ui1, ui2, v1, v2, w12 in zip(times_aprox, ue_aprox_values, ui1_aprox_values, ui2_aprox_values, v1_aprox_values, v2_aprox_values, w12_aprox_values):
				writer.writerow([t, ue, ui1, ui2, v1, v2, w12])

		print(f"Data for n={n} saved to {csv_filename}")
