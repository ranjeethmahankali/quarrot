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

for ci in range(880, 900):
    with open(f"/home/rnjth94/buffer/parametrization/chord{ci}.txt", "r") as f:
        faces = pgf.var_int([int(line) for line in f.readlines()])
    chord = pgf.subMesh(subd, faces)
    pgv.show(f"chord{ci}", chord)

pgv.runCommands("""
zoomextents
meshedges on
""")
