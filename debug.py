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

pgv.runCommands("""
zoomextents
meshedges on
""")
