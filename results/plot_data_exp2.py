import matplotlib.pyplot as plt
import csv
import numpy as np
import matplotlib.lines as mlines

plt.rcParams.update({
	"text.usetex": True,
	"font.family": "serif",
	"font.serif": ["Times New Roman"],
	"axes.labelsize": 16,
	"xtick.labelsize": 12,
	"ytick.labelsize": 12,
	"legend.fontsize": 12,
})

line_styles = {20: "--", 40: "-." , 80: (0, (5, 2, 1, 2, 1, 2)) , 160: ":"}
colors = {
	"Ue": "cornflowerblue",
	"Ui1": "darkorchid",
	"Ui2": "deepskyblue",
	"V1": "tomato",
	"V2": "darkorange",
	"W12": "firebrick"
}
cl_mapping = {20: 0.28, 40: 0.20, 80: 0.14, 160: 0.10}

t_smooth = np.linspace(0, 7, 500)
radius = 5
ue_exact_smooth = (50*np.exp(-t_smooth/10)*np.cos(t_smooth))/radius
ui1_exact_smooth = (5/181) * (181 + 145 * np.exp(-t_smooth) + (2 * np.exp(-t_smooth/10) * ((905 + 18 * radius) * np.cos(t_smooth) + 20 * radius * np.sin(t_smooth))) / radius)
ui2_exact_smooth = (5/181) * (181 + 869 * np.exp(-t_smooth) + (np.exp(-t_smooth/10) * (2 * (905 + 18 * radius) * np.cos(t_smooth) + 40 * radius * np.sin(t_smooth))) / radius)
v1_exact_smooth = (5/181) * (181 + 145 * np.exp(-t_smooth) + 4 * np.exp(-t_smooth/10) * (9 * np.cos(t_smooth) + 10 * np.sin(t_smooth)))
v2_exact_smooth = (5/181) * (181 + 869 * np.exp(-t_smooth) + 4 * np.exp(-t_smooth/10) * (9 * np.cos(t_smooth) + 10 * np.sin(t_smooth)))
w12_exact_smooth = -20*np.exp(-t_smooth)

"""
-----------
Plot 1
-----------
"""
t_zoom = 2.5
delta_t = 0.2
fig, axs = plt.subplots(1, 2, figsize=(12, 6), dpi=300, gridspec_kw={'width_ratios': [2, 1]})

# Main plot
ax_main = axs[0]

# Store handles for the exact solutions
ue_line, = ax_main.plot(t_smooth, ue_exact_smooth, color=colors["Ue"], linewidth=2, label=r"$u_e$")
ui1_line, = ax_main.plot(t_smooth, ui1_exact_smooth, color=colors["Ui1"], linewidth=2, label=r"$u_{i}^{(1)}$")
ui2_line, = ax_main.plot(t_smooth, ui2_exact_smooth, color=colors["Ui2"], linewidth=2, label=r"$u_{i}^{(2)}$")

# Store handles for numerical solutions
numerical_handles = []

ax_zoom = axs[1]
ax_zoom.set_xlim(t_zoom - delta_t, t_zoom + delta_t)
ax_zoom.set_ylim(-8, 4)

for csv_file in ["data_n_20.csv", "data_n_40.csv", "data_n_80.csv", "data_n_160.csv"]:
	times = []
	ue_values = []
	ui1_values = []
	ui2_values = []

	with open("plot_data_exp2/"+csv_file, mode='r') as file:
		reader = csv.reader(file)
		next(reader)  # Skip header
		for row in reader:
			times.append(float(row[0]))
			ue_values.append(float(row[1]))
			ui1_values.append(float(row[2]))
			ui2_values.append(float(row[3]))

	n_steps = int(csv_file.split('_')[2].split('.')[0])
	cl = cl_mapping[n_steps]
	ue_num, = ax_main.plot(times, ue_values, linestyle=line_styles[n_steps], color=colors["Ue"], label=rf"$h={cl}, n_\mathrm{{f}}={n_steps}$")
	ui1_num, = ax_main.plot(times, ui1_values, linestyle=line_styles[n_steps], color=colors["Ui1"], label=rf"$h={cl}, n_\mathrm{{f}}={n_steps}$")
	ui2_num, = ax_main.plot(times, ui2_values, linestyle=line_styles[n_steps], color=colors["Ui2"], label=rf"$h={cl}, n_\mathrm{{f}}={n_steps}$")
	numerical_handles.extend([ue_num, ui1_num, ui2_num])
	
	ue_zoom_exact = (50*np.exp(-np.array([t_zoom - delta_t, t_zoom + delta_t])/10)*np.cos([t_zoom - delta_t, t_zoom + delta_t]))/radius
	ui1_zoom_exact = (5/181) * (181 + 145 * np.exp(-np.array([t_zoom - delta_t, t_zoom + delta_t])) + (2 * np.exp(-np.array([t_zoom - delta_t, t_zoom + delta_t])/10) * ((905 + 18 * radius) * np.cos([t_zoom - delta_t, t_zoom + delta_t]) + 20 * radius * np.sin([t_zoom - delta_t, t_zoom + delta_t]))) / radius)
	ui2_zoom_exact = (5/181) * (181 + 869 * np.exp(-np.array([t_zoom - delta_t, t_zoom + delta_t])) + (np.exp(-np.array([t_zoom - delta_t, t_zoom + delta_t])/10) * (2 * (905 + 18 * radius) * np.cos([t_zoom - delta_t, t_zoom + delta_t]) + 40 * radius * np.sin([t_zoom - delta_t, t_zoom + delta_t]))) / radius)

	ax_zoom.plot([t_zoom - delta_t, t_zoom + delta_t], ue_zoom_exact, color=colors["Ue"], label=r"Exact $u_e$")
	ax_zoom.plot([t_zoom - delta_t, t_zoom + delta_t], ui1_zoom_exact, color=colors["Ui1"], label=r"Exact $u_{i}^{(1)}$")
	ax_zoom.plot([t_zoom - delta_t, t_zoom + delta_t], ui2_zoom_exact, color=colors["Ui2"], label=r"Exact $u_{i}^{(2)}$")

	for func_values, color in zip([ue_values, ui1_values, ui2_values], [colors["Ue"], colors["Ui1"], colors["Ui2"]]):
		before = [(t, v) for t, v in zip(times, func_values) if t < t_zoom - delta_t]
		after = [(t, v) for t, v in zip(times, func_values) if t > t_zoom + delta_t]
		if before and after:
			(t_before, v_before) = before[-1]
			(t_after, v_after) = after[0]
			value_left = v_before + (v_after - v_before) * (t_zoom - delta_t - t_before) / (t_after - t_before)
			value_right = v_before + (v_after - v_before) * (t_zoom + delta_t - t_before) / (t_after - t_before)
			ax_zoom.plot(times, func_values, linestyle=line_styles[n_steps], color=color)

ax_main.set_xlabel(r"$t$")
ax_main.set_ylabel(r"$u_e, u_{i}^{(1)}, u_{i}^{(2)}$")


dummy_handle = mlines.Line2D([], [], color='none', marker='', linestyle='')
first_col_handles = [ue_line, ui1_line, ui2_line , dummy_handle, dummy_handle, dummy_handle]
remaining_handles = numerical_handles
handles = first_col_handles + remaining_handles
legend = ax_main.legend(handles=handles, loc="upper right", ncol=3)
ax_main.add_artist(legend)

ax_zoom.set_xlabel(r"$t$")
ax_zoom.set_title("Zoom arround $t=2.5$")
ax_zoom.set_xticks([t_zoom - delta_t, t_zoom, t_zoom + delta_t])
plt.tight_layout()
plt.savefig('plot1_exp2.pdf', format='pdf')

"""
-----------
Plot 2
-----------
"""
t_zoom = 2
delta_t = 0.2
fig, axs = plt.subplots(1, 2, figsize=(12, 6), dpi=300, gridspec_kw={'width_ratios': [2, 1]})

# Main plot
ax_main = axs[0]

# Store handles for the exact solutions
v1_line, = ax_main.plot(t_smooth, v2_exact_smooth, color=colors["V1"], linewidth=2, label=r"$v^{(1)}$")
v2_line, = ax_main.plot(t_smooth, v2_exact_smooth, color=colors["V2"], linewidth=2, label=r"$v^{(2)}$")
w12_line, = ax_main.plot(t_smooth, w12_exact_smooth, color=colors["W12"], linewidth=2, label=r"$w^{(1,2)}$")

# Store handles for numerical solutions
numerical_handles = []

ax_zoom = axs[1]
ax_zoom.set_xlim(t_zoom - delta_t, t_zoom + delta_t)
ax_zoom.set_ylim(-6, 15)

for csv_file in ["data_n_20.csv", "data_n_40.csv", "data_n_80.csv", "data_n_160.csv"]:
	times = []
	v1_values = []
	v2_values = []
	w12_values = []

	with open("plot_data_exp2/"+csv_file, mode='r') as file:
		reader = csv.reader(file)
		next(reader)  # Skip header
		for row in reader:
			times.append(float(row[0]))
			v1_values.append(float(row[4]))
			v2_values.append(float(row[5]))
			w12_values.append(float(row[6]))

	n_steps = int(csv_file.split('_')[2].split('.')[0])
	cl = cl_mapping[n_steps]
	v1_num, = ax_main.plot(times, v1_values, linestyle=line_styles[n_steps], color=colors["V1"], label=rf"$h={cl}, n_\mathrm{{f}}={n_steps}$")
	v2_num, = ax_main.plot(times, v2_values, linestyle=line_styles[n_steps], color=colors["V2"], label=rf"$h={cl}, n_\mathrm{{f}}={n_steps}$")
	w12_num, = ax_main.plot(times, w12_values, linestyle=line_styles[n_steps], color=colors["W12"], label=rf"$h={cl}, n_\mathrm{{f}}={n_steps}$")
	numerical_handles.extend([v1_num, v2_num, w12_num])

	v1_zoom_exact = (5/181) * (181 + 145 * np.exp(-np.array([t_zoom - delta_t, t_zoom + delta_t])) + 4 * np.exp(-np.array([t_zoom - delta_t, t_zoom + delta_t])/10) * (9 * np.cos([t_zoom - delta_t, t_zoom + delta_t]) + 10 * np.sin([t_zoom - delta_t, t_zoom + delta_t])))
	v2_zoom_exact = (5/181) * (181 + 869 * np.exp(-np.array([t_zoom - delta_t, t_zoom + delta_t])) + 4 * np.exp(-np.array([t_zoom - delta_t, t_zoom + delta_t])/10) * (9 * np.cos([t_zoom - delta_t, t_zoom + delta_t]) + 10 * np.sin([t_zoom - delta_t, t_zoom + delta_t])))
	w12_zoom_exact = -20*np.exp(-np.array([t_zoom - delta_t, t_zoom + delta_t]))

	ax_zoom.plot([t_zoom - delta_t, t_zoom + delta_t], v1_zoom_exact, color=colors["V1"], label=r"Exact $v^{(1)}$")
	ax_zoom.plot([t_zoom - delta_t, t_zoom + delta_t], v2_zoom_exact, color=colors["V2"], label=r"Exact $v^{(2)}$")
	ax_zoom.plot([t_zoom - delta_t, t_zoom + delta_t], w12_zoom_exact, color=colors["W12"], label=r"Exact $w^{(1,2)}$")

	# Extract and interpolate zoomed region
	for func_values, color in zip([v1_values, v2_values, w12_values], [colors["V1"], colors["V2"], colors["W12"]]):
		before = [(t, v) for t, v in zip(times, func_values) if t < t_zoom - delta_t]
		after = [(t, v) for t, v in zip(times, func_values) if t > t_zoom + delta_t]

		if before and after:
			(t_before, v_before) = before[-1]
			(t_after, v_after) = after[0]

			value_left = v_before + (v_after - v_before) * (t_zoom - delta_t - t_before) / (t_after - t_before)
			value_right = v_before + (v_after - v_before) * (t_zoom + delta_t - t_before) / (t_after - t_before)

			ax_zoom.plot(times, func_values, linestyle=line_styles[n_steps], color=color)

ax_main.set_xlabel(r"$t$")
ax_main.set_ylabel(r"$v^{(1)}, v^{(2)}, w^{(1,2)}$")

legend1_handles = [v1_line, v2_line, w12_line] + numerical_handles[:6]
legend1 = ax_main.legend(handles=legend1_handles, loc="upper right", ncol=3)
legend2_handles = numerical_handles[6:]
legend2 = ax_main.legend(handles=legend2_handles, loc="lower right", ncol=2)
ax_main.add_artist(legend1)

ax_zoom.set_xlabel(r"$t$")
ax_zoom.set_title("Zoom arround $t=2$")
ax_zoom.set_xticks([t_zoom - delta_t, t_zoom, t_zoom + delta_t])
plt.tight_layout()
plt.savefig('plot2_exp2.pdf', format='pdf')
