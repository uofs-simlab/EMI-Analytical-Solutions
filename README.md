# EMI Analytical Solutions

This repository contains the code to generate and handle results for the paper **"Analytical Solutions for the Extracellular-Membrane-Intracellular (EMI) Model"**.

## 📂 Repository Structure

The repository includes the following key components:

### 🔹 Experiment Scripts

| File                                  | Description                                        |
| ------------------------------------- | -------------------------------------------------- |
| `exp1_single_circular_cell.py`        | Simulates a single circular cell.                  |
| `exp2_two_coupled_spherical_cells.py` | Simulates two coupled spherical cells.             |
| `exp3_manufactured_solution.py`       | Implements a manufactured solution for validation. |

### 🔹 Meshes

The `meshes/` folder contains `.geo` files for generating meshes required for the experiments. To generate a mesh, use the following command:

```sh
gmsh -<dim> <file_name>.geo -o <output_file_name>.msh
```

where `<dim>` is the dimension of the mesh (e.g., 2 for 2D, 3 for 3D).

### 🔹 Results

The `results/` folder contains scripts to process the results. Since the data files are large, you can download them from [Zenodo](https://doi.org/10.5281/zenodo.14948002) and place them in a `data/` folder inside the repository.

Once the data is available, you can use the following scripts to process it:

| File                       | Description                                                                      |
| -------------------------- | -------------------------------------------------------------------------------- |
| `errors_exp<n>.py`         | Generates an Excel dataset with the $L^2$-norm error values presented in the paper. |
| `take_plot_data_exp<n>.py` | Extracts data for plotting.                                                      |
| `plot_data_exp<n>.py`      | Generates the plots shown in the paper.                                          |

---

## 🛠 Requirements

The following dependencies are required to run the simulations:

- `gmsh`
- `firedrake`
- `irksome`
- `pythOS`

---

## 🏗 Installation

You can install the dependencies manually, but for convenience, the repository includes `container.def`, a definition file for building a container using **Apptainer** (formerly Singularity) with Ubuntu Linux.

### 🔹 Building the Container

Run the following command to build the container:

```sh
singularity build container.sif container.def
```

Once built, the `container.sif` file will be available for execution.

### 🔹 Running the Container

To launch the container, use the following command:

```sh
apptainer shell --writable-tmpfs --bind ./:/container --bind <path-to-repo>/EMI-Analytical-Solutions:<path-to-repo>/EMI-Analytical-Solutions container.sif
```

Replace `<path-to-repo>` with the actual path to your repository.

### 🔹 Inside the Container

Once inside, navigate to the container directory:

```sh
cd /container
```

Then, activate the **Firedrake** virtual environment:

```sh
source /opt/firedrake/bin/activate
```

You're now ready to run simulations and process results!

---



