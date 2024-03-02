#include <quarrot.h>

namespace quarrot {

void subdivide(Mesh& mesh)
{
  if (std::all_of(mesh.faces_begin(), mesh.faces_end(), [&mesh](FaceH fh) {
        return std::distance(mesh.cfv_begin(fh), mesh.cfv_end(fh)) == 4;
      })) {
    // The mesh is already a quad-only mesh.
    return;
  }
  OpenMesh::EPropHandleT<VertH> everts;
  OpenMesh::FPropHandleT<VertH> fverts;
  mesh.add_property(everts);
  mesh.add_property(fverts);
  // Compute and store edge midpoints.
  for (EdgeH eh : mesh.edges()) {
    mesh.property(everts, eh) = mesh.add_vertex(
      0.5 * (mesh.point(mesh.to_vertex_handle(mesh.halfedge_handle(eh, 0))) +
             mesh.point(mesh.to_vertex_handle(mesh.halfedge_handle(eh, 1)))));
  }
  // Compute and store face centers.
  for (FaceH fh : mesh.faces()) {
    mesh.property(fverts, fh) = mesh.add_vertex(
      std::accumulate(mesh.cfv_begin(fh),
                      mesh.cfv_end(fh),
                      glm::dvec3(0., 0., 0.),
                      [&](glm::dvec3 sum, VertH vh) { return sum + mesh.point(vh); }) /
      double(std::distance(mesh.cfv_begin(fh), mesh.cfv_end(fh))));
  }
  // Replace faces with subdivided quads.
  std::vector<std::array<VertH, 4>> quads;  // Temporary buffer to store subdivided faces.
  for (FaceH fh : mesh.faces()) {
    quads.clear();
    std::transform(
      mesh.cfh_begin(fh), mesh.cfh_end(fh), std::back_inserter(quads), [&](HalfH he) {
        return std::array<VertH, 4> {
          {mesh.from_vertex_handle(he),
           mesh.property(everts, mesh.edge_handle(he)),
           mesh.property(fverts, mesh.face_handle(he)),
           mesh.property(everts, mesh.edge_handle(mesh.prev_halfedge_handle(he)))}};
      });
    mesh.delete_face(fh, false);
    for (const auto& quad : quads) {
      mesh.add_face(quad.data(), quad.size());
    }
  }
  mesh.remove_property(everts);
  mesh.remove_property(fverts);
  mesh.delete_isolated_vertices();
  mesh.garbage_collection();
}

}  // namespace quarrot
