

# Steps

Here we use the Double Mach reflection problem as test case.

### 1. Generating mesh
```bash
cd pressio-demoapps/examples
python ../meshing_scripts/create_full_mesh_for.py --problem doublemach2d_s5 -n 200 50 --outDir .
```

### 2. Running the FOM simulation
```bash
python generate_data.py -s 5 --dt 0.001 -T 0.25
# -s  sampling-frequency
# -dt time-step-size
# -T  simulation-time
```
This will generate two files:
```
rhs_snapshots.txt
state_snapshots.txt
```
To make sure things are right, you can want to make a plot:
```bash
# this script will read state_snapshots.txt and plot density at last step
python make_plot.py
```

### 3. Generate sample mesh points
```bash
python psample.py
```
This will generate a file `/samplemesh/sample_mesh_gids.dat` which contain the
indices of the sample mesh cells picked by the psample method.

### 4. Create the stencil mesh and visualize the sample mesh points
```bash
cd pressio-demoapps/meshing_scripts/
python create_sample_mesh.py --outdir ../examples/samplemesh --fullMeshDir ../examples/
python plot_mesh.py --wdir ../examples/samplemesh -p show
```
This will show a plot with the stencil mesh.
