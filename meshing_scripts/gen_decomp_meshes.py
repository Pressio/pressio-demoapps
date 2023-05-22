import os
import subprocess

mesh_script = "./create_full_mesh.py"

# -----

# double shock reflection
xbounds = [0.0, 4.0]
ybounds = [0.0, 1.0]
zbounds = [None, None]

Nx = 200
Ny = 50
Nz = 0

ndom_x = 2
ndom_y = 2
ndom_z = 1

stencil = 3

overlap = 2

outdir_base = "/home/crwentl/research/code/pressio-demoapps/build-debug-noprint/tests_cpp/eigen_2d_euler_double_mach_reflection_schwarz/firstorder/meshes/mesh_2x2"

# -----

ndim = 1

assert int(overlap/2) == (overlap/2), "Overlap must be even for now (FIX)"
assert Nx > 0
assert Ny >= 0
assert Nz >= 0
assert ndom_x > 0
assert all([bound is not None for bound in xbounds])
assert xbounds[1] > xbounds[0]
if Ny > 0:
    assert ndom_y > 0
    assert all([bound is not None for bound in ybounds])
    assert ybounds[1] > ybounds[0]
    ndim += 1
    if Nz > 0:
        assert ndom_z > 0
        assert all([bound is not None for bound in zbounds])
        assert zbounds[1] > zbounds[0]
        ndim += 1
else:
    assert Nz == 0

def prep_dim(N, ndom, bounds):
    d = (bounds[1] - bounds[0]) / N
    N_dom = [int(N / ndom)] * ndom
    fill = N - int(N / ndom) * ndom
    for idx in range(fill):
        N_dom[idx] += 1

    return d, N_dom 

def prep_dom_dim(dom_idx, ndom, N_dom, offset, bounds, d):
    
    # cells
    n = N_dom[dom_idx]
    if ndom > 1:
        if (dom_idx == 0) or (dom_idx == (ndom - 1)):
            n += offset
        else:
            n += 2 * offset
    
    # bounds
    bound = [None, None]
    if dom_idx == 0:
        bound[0] = bounds[0]
    else:
        bound[0] = (sum(N_dom[:dom_idx]) - offset) * d
    if (dom_idx == ndom - 1):
        bound[1] = bounds[1]
    else:
        bound[1] = (sum(N_dom[:dom_idx+1]) + offset) * d

    return n, bound

def main():
    
    offset = int(overlap / 2)
 
    dx, Nx_dom = prep_dim(Nx, ndom_x, xbounds)
    if (Ny > 0):
        dy, Ny_dom = prep_dim(Ny, ndom_y, ybounds)
        if (Nz > 0):
            dz, Nz_dom = prep_dim(Nz, ndom_z, zbounds)

    dom_count = 0
    xbound = [None, None]
    ybound = [None, None]
    zbound = [None, None]
    for dom_z in range(ndom_z):
        for dom_y in range(ndom_y):
            for dom_x in range(ndom_x):

                print("Domain " + str(dom_count))
                print("Cells: ", end='') 

                # subdomain dimensions
                nx, xbound = prep_dom_dim(
                    dom_x, ndom_x, Nx_dom, offset, xbounds, dx
                )
                print(str(nx), end='')
                if Ny > 0:
                    ny, ybound = prep_dom_dim(
                        dom_y, ndom_y, Ny_dom, offset, ybounds, dy
                    )
                    print(" x " + str(ny), end='')
                    if Nz > 0:
                        nz, zbound = prep_dom_dim(
                            dom_z, ndom_z, Nz_dom, offset, zbounds, dz
                        )
                        print(" x " + str(nz), end='')
                print("")

                print("x-bounds: (" + str(xbound[0]) + ", " + str(xbound[1]) + ")")
                if Ny > 0:
                    print("y-bounds: (" + str(ybound[0]) + ", " + str(ybound[1]) + ")")
                    if Nz > 0:
                        print("z-bounds: (" + str(zbound[0]) + ", " + str(zbound[1]) + ")")

                # subdomain mesh subdirectory
                outdir = os.path.join(outdir_base, "domain_" + str(dom_count))
                if not os.path.isdir(outdir):
                    os.makedirs(outdir)
                print("Output directory: " + outdir)

                if Ny == 0:
                    args = ("python3", mesh_script,
                    "-n", str(nx),
                    "--outDir", outdir,
                    "--stencilSize", str(stencil),
                    "--bounds", str(xbound[0]), str(xbound[1]))
                elif Nz == 0:
                    args = ("python3", mesh_script,
                    "-n", str(nx), str(ny),
                    "--outDir", outdir,
                    "--stencilSize", str(stencil),
                    "--bounds", str(xbound[0]), str(xbound[1]),
                                str(ybound[0]), str(ybound[1]),
                    )
                else:
                    args = ("python3", mesh_script,
                    "-n", str(nx), str(ny), str(nz),
                    "--outDir", outdir,
                    "--stencilSize", str(stencil),
                    "--bounds", str(xbound[0]), str(xbound[1]),
                                str(ybound[0]), str(ybound[1]),
                                str(zbound[0]), str(zbound[1]),
                    )

                popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

                dom_count += 1

    with open(os.path.join(outdir_base, "info_domain.dat"), "w") as f:
        f.write("dim %8d\n" % ndim)
        f.write("ndomX %8d\n" % ndom_x)
        if Ny > 0:
            f.write("ndomY %8d\n" % ndom_y)
            if Nz > 0:
                f.write("ndomZ %8d\n" % ndom_y)

        f.write("overlap %8d\n" % overlap)

if __name__ == "__main__":
    main()
