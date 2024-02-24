#include <quarrot.h>
#include <OpenMesh/Core/Utils/Property.hh>
#include <algorithm>
#include <cassert>
#include <glm/gtx/norm.hpp>
#include <queue>
#include <stdexcept>

namespace quarrot {

void simplify(Mesh& mesh)
{
  throw std::logic_error("Not Implemented");
}

static bool isTriangle(const Mesh& mesh, FaceH fh)
{
  if (!fh.is_valid()) {
    return false;
  }
  return std::distance(mesh.cfv_begin(fh), mesh.cfv_end(fh)) == 3;
}

static void mark_fixed_edges(Mesh& mesh)
{
  for (EdgeH eh : mesh.edges()) {
    auto he  = mesh.halfedge_handle(eh, 0);
    auto ohe = mesh.opposite_halfedge_handle(he);
    mesh.status(eh).set_locked(!isTriangle(mesh, mesh.face_handle(he)) ||
                               !isTriangle(mesh, mesh.face_handle(ohe)));
  }
}

static std::array<VertH, 4> quad_from_diagonal(const Mesh& mesh, EdgeH eh)
{
  auto he  = mesh.halfedge_handle(eh, 0);
  auto ohe = mesh.opposite_halfedge_handle(he);
  return {{mesh.from_vertex_handle(he),
           mesh.to_vertex_handle(mesh.next_halfedge_handle(ohe)),
           mesh.to_vertex_handle(he),
           mesh.to_vertex_handle(mesh.next_halfedge_handle(he))}};
}

static void calc_beta(Mesh& mesh, const OpenMesh::EPropHandleT<double>& beta)
{
  static constexpr double SQRT_3 = 1.73205080757;
  /*
               C
              /|\
             / | \
            /  |  \
           D   |   B
            \  |  /
             \ | /
              \|/
               A
   */
  for (EdgeH eh : mesh.edges()) {
    if (mesh.status(eh).locked()) {
      continue;
    }
    std::array<VertH, 4>      verts = quad_from_diagonal(mesh, eh);
    std::array<glm::dvec3, 4> points;
    for (size_t i = 0; i < 4; ++i) {
      points[i] = mesh.point(verts[i]);
    }
    std::array<double, 4> alphas;
    for (size_t i = 0; i < 4; ++i) {
      const auto& a     = points[i];
      const auto& b     = points[(i + 1) % 4];
      const auto& c     = points[(i + 2) % 4];
      auto        ca    = a - c;
      auto        cb    = b - c;
      glm::dvec3  norm  = glm::cross(ca, cb);
      double      dot   = glm::length(norm);
      glm::dvec3  enorm = mesh.normal(mesh.face_handle(mesh.halfedge_handle(eh, 0))) +
                         mesh.normal(mesh.face_handle(mesh.halfedge_handle(eh, 1)));
      if (glm::dot(enorm, norm) < 0.) {
        dot = -dot;
      }
      alphas[i] =
        2. * SQRT_3 * dot / (glm::length2(ca) + glm::length2(b - a) + glm::length2(cb));
    }
    std::sort(alphas.begin(), alphas.end());
    mesh.property(beta, eh) = (alphas[0] * alphas[1]) / (alphas[2] * alphas[3]);
  }
}

static double calc_edge_beta_star(Mesh&                                 mesh,
                                  EdgeH                                 eh,
                                  const OpenMesh::EPropHandleT<double>& beta,
                                  const OpenMesh::EPropHandleT<double>& beta_star)
{
  if (mesh.status(eh).locked()) {
    return 0.;
  }
  auto                 he        = mesh.halfedge_handle(eh, 0);
  auto                 ohe       = mesh.opposite_halfedge_handle(he);
  std::array<EdgeH, 4> neighbors = {{mesh.edge_handle(mesh.next_halfedge_handle(ohe)),
                                     mesh.edge_handle(mesh.prev_halfedge_handle(ohe)),
                                     mesh.edge_handle(mesh.next_halfedge_handle(he)),
                                     mesh.edge_handle(mesh.prev_halfedge_handle(he))}};
  return mesh.property(beta, eh) -
         std::accumulate(
           neighbors.begin(), neighbors.end(), 0., [&](double total, EdgeH ne) {
             return total + ((mesh.status(ne).locked() || mesh.status(ne).deleted())
                               ? 0.
                               : mesh.property(beta, ne));
           });
}

static void calc_beta_star(Mesh&                                 mesh,
                           const OpenMesh::EPropHandleT<double>& beta,
                           const OpenMesh::EPropHandleT<double>& beta_star)
{
  for (EdgeH eh : mesh.edges()) {
    mesh.property(beta_star, eh) = calc_edge_beta_star(mesh, eh, beta, beta_star);
  }
}

void pair_triangles(Mesh& mesh, double gamma)
{
  /*
   Lo, S. H.. “Generating quadrilateral elements on plane and over curved surfaces.”
   Computers & Structures 31 (1989): 421-426.
  */
  mesh.update_face_normals();
  mark_fixed_edges(mesh);
  OpenMesh::EPropHandleT<double> beta, beta_star;
  mesh.add_property(beta);
  mesh.add_property(beta_star);
  // Compute beta values and initial beta-star values of edges.
  calc_beta(mesh, beta);
  calc_beta_star(mesh, beta, beta_star);
  // Priority queue of diagonals to be removed.
  using Pair = std::pair<double, EdgeH>;
  std::priority_queue<Pair> queue;
  for (EdgeH eh : mesh.edges()) {
    const auto& status = mesh.status(eh);
    if (status.locked() || status.deleted()) {
      continue;
    }
    double b = mesh.property(beta, eh);
    if (b < gamma) {
      continue;
    }
    queue.push(std::make_pair(mesh.property(beta_star, eh), eh));
  }
  // Remove diagonals
  while (!queue.empty()) {
    auto [bstar, eh] = queue.top();
    queue.pop();
    if (mesh.status(eh).deleted()) {
      continue;
    }
    else if (mesh.property(beta, eh) < gamma) {
      mesh.status(eh).set_locked(true);
      continue;
    }
    else if (mesh.status(eh).locked()) {
      continue;
    }
    // Delete edge and replace with quad.
    auto verts = quad_from_diagonal(mesh, eh);
    mesh.delete_edge(eh, false);
    mesh.add_face(verts.data(), verts.size());
    // Update neighbors.
    std::array<HalfH, 4> neighbors;
    {
      auto he   = mesh.halfedge_handle(eh, 0);
      auto ohe  = mesh.opposite_halfedge_handle(he);
      neighbors = {{mesh.next_halfedge_handle(ohe),
                    mesh.prev_halfedge_handle(ohe),
                    mesh.next_halfedge_handle(he),
                    mesh.prev_halfedge_handle(he)}};
    }
    for (HalfH nhe : neighbors) {
      mesh.status(mesh.edge_handle(nhe)).set_locked(true);
      auto ohe = mesh.opposite_halfedge_handle(nhe);
      if (!mesh.face_handle(ohe).is_valid()) {
        continue;
      }
      auto h1                      = mesh.next_halfedge_handle(ohe);
      auto e1                      = mesh.edge_handle(h1);
      auto e2                      = mesh.edge_handle(mesh.next_halfedge_handle(h1));
      mesh.property(beta_star, e1) = calc_edge_beta_star(mesh, e1, beta, beta_star);
      mesh.property(beta_star, e2) = calc_edge_beta_star(mesh, e2, beta, beta_star);
    }
  }
  mesh.remove_property(beta);
  mesh.remove_property(beta_star);
  mesh.delete_isolated_vertices();
  mesh.garbage_collection();
}

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
