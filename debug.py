"""Script for debugging the algo."""
import pygalfunc as pgf
import pygalview as pgv

mesh = pgf.loadPolyMesh(
    pgf.absPath(
        pgf.var_string("/home/rnjth94/buffer/parametrization/bimba.obj")))
pgv.show("original", mesh)

paired = pgf.loadPolyMesh(
    pgf.absPath(
        pgf.var_string(
            "/home/rnjth94/buffer/parametrization/bimba_paired.obj")))
pgv.show("paired", paired)

subd = pgf.loadPolyMesh(
    pgf.absPath(
        pgf.var_string(
            "/home/rnjth94/buffer/parametrization/bimba_subdivided.obj")))
pgv.show("subdivided", subd)

with open("temp/chord_err.txt", "r") as f:
    err = [float(line) for line in f.readlines()]

indices = sorted(list(range(len(err))), key=lambda i: err[i])

for ci in indices[:10]:
    with open(f"temp/chord{ci}.txt", "r") as f:
        findices = [int(line) for line in f.readlines()]
        faces = pgf.var_int(findices)
        chord = pgf.subMesh(subd, faces)
        colors = [(1., 0., 0.) for _ in range(len(findices))]
        chord = pgf.meshWithVertexColors(chord, pgf.var_vec3(colors))
        pgv.show(f"chord{ci}", chord)

pgv.runCommands("""
zoomextents
meshedges on
""")
