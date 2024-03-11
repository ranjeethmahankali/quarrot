"""Script for debugging the algo."""
import pygalfunc as pgf
import pygalview as pgv

subd = pgf.loadPolyMesh(
    pgf.absPath(
        pgf.var_string(
            "/home/rnjth94/buffer/parametrization/bimba_subdivided.obj")))
pgv.show("subdivided", subd)

paths = pgf.var_string([
    f"/home/rnjth94/buffer/parametrization/bimba_collapsed{ci}.obj"
    for ci in range(6)
])
collapsedList = pgf.loadPolyMesh(pgf.absPath(paths))
pgv.print("Number of meshes", pgf.listLength(collapsedList))
index = pgv.slideri32("Index", 0, 5, 0)
collapsed = pgf.listItem(collapsedList, index)
pgv.print("Path", pgf.listItem(paths, index))
pgv.show("collapsed", collapsed)

pgv.runCommands("""
zoomextents
meshedges on
""")
