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
collapsed = pgf.listItem(collapsedList, pgv.slideri32("Index", 0, 10, 0))
pgv.show("collapsed", collapsed)

pgv.runCommands("""
zoomextents
meshedges on
""")
