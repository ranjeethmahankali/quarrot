"""Script for debugging the algo."""
import pygalfunc as pgf
import pygalview as pgv

nfiles = 20
testpath = "/home/rnjth94/dev/quarrot/temp/"

subd = pgf.loadPolyMesh(
    pgf.absPath(pgf.var_string(testpath + "subdivided.obj")))
pgv.show("subdivided", subd)

paths = pgf.var_string(
    [testpath + f"collapsed{ci}.obj" for ci in range(nfiles)])

collapsedList = pgf.loadPolyMesh(pgf.absPath(paths))
pgv.print("Number of meshes", pgf.listLength(collapsedList))
index = pgv.slideri32("Index", 0, nfiles - 1, 0)
collapsed = pgf.listItem(collapsedList, index)
pgv.print("Path", pgf.listItem(paths, index))
pgv.show("collapsed", collapsed)

pgv.runCommands("""
zoomextents
meshedges on
""")
