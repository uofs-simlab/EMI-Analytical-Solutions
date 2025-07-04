import matplotlib.pyplot as plt
import csv
import numpy as np

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
    "axes.labelsize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
})

line_styles = {14:"--" , 28: "-.", 56: (0, (5, 2, 1, 2, 1, 2)), 112: ":"}
colors = {"Ue": "cornflowerblue", "Ui": "darkorchid", "V": "tomato"}
cl_mapping = {14: 0.28, 28: 0.20, 56: 0.14, 112: 0.10}

t_smooth = np.linspace(0, 7, 500)
radius = 5

ue_exact_smooth = 5 + 10 * np.log(radius) * np.sin(t_smooth)
ui_exact_smooth = 10 + 14 * np.exp(-t_smooth) + np.cos(t_smooth) - np.sin(t_smooth) + 10 * np.log(radius) * np.sin(t_smooth)
v_exact_smooth = 5 + 14 * np.exp(-t_smooth) + np.cos(t_smooth) - np.sin(t_smooth)

t_zoom = 1.0
delta_t = 0.2
fig, axs = plt.subplots(1, 2, figsize=(12, 6), dpi=300, gridspec_kw={'width_ratios': [2, 1]})

# Main plot
ax_main = axs[0]

ue_line, = ax_main.plot(t_smooth, ue_exact_smooth, color=colors["Ue"], linewidth=2, label=r"$u_e$")
ui_line, = ax_main.plot(t_smooth, ui_exact_smooth, color=colors["Ui"], linewidth=2, label=r"$u_{i}^{(1)}$")
v_line, = ax_main.plot(t_smooth, v_exact_smooth, color=colors["V"], linewidth=2, label=r"$v^{(1)}$")

numerical_handles = []

ax_zoom = axs[1]
ax_zoom.set_xlim(t_zoom - delta_t, t_zoom + delta_t)
ax_zoom.set_ylim(5, 30)

for csv_file in ["data_n_14.csv", "data_n_28.csv", "data_n_56.csv", "data_n_112.csv"]:
	times = []
	ue_values = []
	ui_values = []
	v_values = []

	with open("plot_data_exp1/"+csv_file, mode='r') as file:
		reader = csv.reader(file)
		next(reader)  # Skip header
		for row in reader:
			times.append(float(row[0]))
			ue_values.append(float(row[1]))
			ui_values.append(float(row[2]))
			v_values.append(float(row[3]))

	n_steps = int(csv_file.split('_')[2].split('.')[0])
	cl = cl_mapping[n_steps]

	ue_num, = ax_main.plot(times, ue_values, linestyle=line_styles[n_steps], color=colors["Ue"], label=rf"$h={cl}, n_\mathrm{{f}}={n_steps}$")
	ui_num, = ax_main.plot(times, ui_values, linestyle=line_styles[n_steps], color=colors["Ui"], label=rf"$h={cl}, n_\mathrm{{f}}={n_steps}$")
	v_num, = ax_main.plot(times, v_values, linestyle=line_styles[n_steps], color=colors["V"], label=rf"$h={cl}, n_\mathrm{{f}}={n_steps}$")
	numerical_handles.extend([ue_num, ui_num, v_num])

	ue_zoom_exact = 5 + 10 * np.log(radius) * np.sin([t_zoom - delta_t, t_zoom + delta_t])
	ui_zoom_exact = 10 + 14 * np.exp(-np.array([t_zoom - delta_t, t_zoom + delta_t])) + np.cos([t_zoom - delta_t, t_zoom + delta_t]) - np.sin([t_zoom - delta_t, t_zoom + delta_t]) + 10 * np.log(radius) * np.sin([t_zoom - delta_t, t_zoom + delta_t])
	v_zoom_exact = 5 + 14 * np.exp(-np.array([t_zoom - delta_t, t_zoom + delta_t])) + np.cos([t_zoom - delta_t, t_zoom + delta_t]) - np.sin([t_zoom - delta_t, t_zoom + delta_t])

	ax_zoom.plot([t_zoom - delta_t, t_zoom + delta_t], ue_zoom_exact, color=colors["Ue"], label=r"Exact $u_e$")
	ax_zoom.plot([t_zoom - delta_t, t_zoom + delta_t], ui_zoom_exact, color=colors["Ui"], label=r"Exact $u_{i}^{(1)}$")
	ax_zoom.plot([t_zoom - delta_t, t_zoom + delta_t], v_zoom_exact, color=colors["V"], label=r"Exact $v^{(1)}$")

	for func_values, color in zip([ue_values, ui_values, v_values], [colors["Ue"], colors["Ui"], colors["V"]]):
		before = [(t, v) for t, v in zip(times, func_values) if t < t_zoom - delta_t]
		after = [(t, v) for t, v in zip(times, func_values) if t > t_zoom + delta_t]

		if before and after:
			(t_before, v_before) = before[-1]
			(t_after, v_after) = after[0]

			value_left = v_before + (v_after - v_before) * (t_zoom - delta_t - t_before) / (t_after - t_before)
			value_right = v_before + (v_after - v_before) * (t_zoom + delta_t - t_before) / (t_after - t_before)

			ax_zoom.plot(times, func_values, linestyle=line_styles[n_steps], color=color)

ax_main.set_xlabel(r"$t$")
ax_main.set_ylabel(r"$u_e, u_{i}^{(1)}, v^{(1)}$")

legend1_handles = [ue_line, ui_line, v_line] + numerical_handles[:6]
legend1 = ax_main.legend(handles=legend1_handles, loc="upper right", ncol=3, bbox_to_anchor=(1, 1))

legend2_handles = numerical_handles[6:]
legend2 = ax_main.legend(handles=legend2_handles, loc="lower left", ncol=2, bbox_to_anchor=(0, 0))

ax_main.add_artist(legend1)

ax_zoom.set_xlabel(r"$t$")
ax_zoom.set_title("Zoom arround $t=1$")
ax_zoom.set_xticks([t_zoom - delta_t, t_zoom, t_zoom + delta_t])

plt.tight_layout()
plt.savefig('plot_exp1.pdf', format='pdf')
