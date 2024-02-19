"""Script for debugging the algo."""
import pygalfunc as pgf
import pygalview as pgv

relpath = pgf.var_string("/home/rnjth94/buffer/parametrization/bimba.obj")
mesh = pgf.loadPolyMesh(pgf.absPath(relpath))
pgv.show("original", mesh)
