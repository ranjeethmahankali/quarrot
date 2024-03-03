#include <dset.h>
#include <algorithm>
#include <cassert>
#include <cstdint>

#include <OpenMesh/Core/Utils/Property.hh>

#include <debug.h>
#include <open_mesh_adaptor.h>
#include <probabilistic_quadircs.h>
#include <quarrot.h>

namespace quarrot {

using Quadric = pq::quadric<pq::math<double, glm::dvec3, glm::dvec3, glm::dmat3>>;

void calc_face_centroids(Mesh& mesh, OpenMesh::FPropHandleT<glm::dvec3> centerProp)
{
  for (FaceH fh : mesh.faces()) {
    if (mesh.status(fh).deleted()) {
      continue;
    }
    mesh.property(centerProp, fh) = mesh.calc_face_centroid(fh);
  }
}

void calc_vert_quadrics(Mesh&                              mesh,
                        OpenMesh::VPropHandleT<Quadric>    qprop,
                        OpenMesh::FPropHandleT<glm::dvec3> fcenterprop)
{
  for (VertH vh : mesh.vertices()) {
    double stddev =
      0.25 * std::accumulate(
               mesh.cve_begin(vh), mesh.cve_end(vh), 0., [&](double total, EdgeH eh) {
                 return mesh.calc_edge_sqr_length(eh) + total;
               });
    Quadric q;
    for (HalfH he : mesh.voh_range(vh)) {  // This only covers half of the triangles.
      FaceH fh = mesh.face_handle(he);
      if (!fh.is_valid()) {
        continue;
      }
      q += Quadric::probabilistic_triangle_quadric(mesh.point(vh),
                                                   mesh.point(mesh.to_vertex_handle(he)),
                                                   mesh.property(fcenterprop, fh),
                                                   stddev);
    }
    for (HalfH he : mesh.vih_range(vh)) {  // The other half of of the triangles.
      FaceH fh = mesh.face_handle(he);
      if (!fh.is_valid()) {
        continue;
      }
      q +=
        Quadric::probabilistic_triangle_quadric(mesh.point(vh),
                                                mesh.property(fcenterprop, fh),
                                                mesh.point(mesh.from_vertex_handle(he)),
                                                stddev);
    }
    mesh.property(qprop, vh) = q;
  }
}

struct Props
{
  OpenMesh::FPropHandleT<std::array<int, 2>> m_chord_indices;
  OpenMesh::FPropHandleT<glm::dvec3>         m_face_centers;
  OpenMesh::VPropHandleT<Quadric>            m_vert_quadrics;
  Mesh&                                      m_mesh;

  explicit Props(Mesh& mesh)
      : m_mesh(mesh)
  {
    mesh.add_property(m_chord_indices);
    for (FaceH fh : mesh.faces()) {  // Initialize all indices to -1.
      mesh.property(m_chord_indices, fh).fill(-1);
    }
    mesh.add_property(m_face_centers);
    calc_face_centroids(mesh, m_face_centers);
    mesh.add_property(m_vert_quadrics);
    calc_vert_quadrics(mesh, m_vert_quadrics, m_face_centers);
  }

  ~Props()
  {
    m_mesh.remove_property(m_chord_indices);
    m_mesh.remove_property(m_face_centers);
    m_mesh.remove_property(m_vert_quadrics);
  }
};

struct PolychordTraverse
{
  std::vector<FaceH>   m_faces;
  DisjointSets         m_vsets;
  std::vector<Quadric> m_quadrics;
  double               m_err_q  = 0.;
  double               m_err_d2 = 0.;

  explicit PolychordTraverse(size_t numVertices)
      : m_vsets(numVertices)
  {}

  bool valid() const { return !m_faces.empty(); }

  void copy_quadrics(Mesh& mesh, OpenMesh::VPropHandleT<Quadric> qprop)
  {
    m_quadrics.clear();
    m_quadrics.reserve(mesh.n_vertices());
    const std::vector<Quadric>& src = mesh.property(qprop).data_vector();
    std::copy(src.begin(), src.end(), std::back_inserter(m_quadrics));
  }

  void record_collapse(Mesh& mesh, HalfH he)
  {
    VertH from = mesh.from_vertex_handle(he);
    VertH to   = mesh.to_vertex_handle(he);
    // Unite the sets.
    m_vsets.unite(uint32_t(from.idx()), uint32_t(to.idx()));
    // Sum the quadrics.
    Quadric q              = m_quadrics[from.idx()] + m_quadrics[to.idx()];
    m_quadrics[from.idx()] = q;
    m_quadrics[to.idx()]   = q;
    // Update the max errors;
    m_err_q  = std::max(m_err_q, q(q.minimizer()));
    m_err_d2 = std::max(m_err_d2, mesh.calc_edge_sqr_length(he));
  }

  void reset(Mesh& mesh, const Props& props)
  {
    m_faces.clear();
    m_vsets.reset();
    copy_quadrics(mesh, props.m_vert_quadrics);
    m_err_q  = 0.;
    m_err_d2 = 0.;
  }

  double valence_error(Mesh& mesh) { throw std::logic_error("Not Implemented"); }
};

struct PolychordInfo
{
  HalfH  m_start_he;
  double m_err_q = 0.;
  double m_err_d = 0.;
  double m_err_v = 0.;

  PolychordInfo(HalfH he, double eq, double ed, double ev)
      : m_start_he(he)
      , m_err_q(eq)
      , m_err_d(ed)
      , m_err_v(ev)
  {}
};

void walk_polychord(Mesh&              mesh,
                    HalfH              he,
                    const Props&       props,
                    int                index,
                    PolychordTraverse& ptraverse)
{
  HalfH start   = he;
  FaceH curface = mesh.face_handle(he);
  ptraverse.reset(mesh, props);
  do {
    auto& indices = mesh.property(props.m_chord_indices, curface);
    ptraverse.m_faces.push_back(curface);
    // Unite the vertices as if the halfedge is being collapsed.
    ptraverse.record_collapse(mesh, he);
    if (indices[0] != -1 && indices[1] != -1) {
      ptraverse.m_faces.clear();
      break;
    }
    he = mesh.opposite_halfedge_handle(
      mesh.next_halfedge_handle(mesh.next_halfedge_handle(he)));
    curface = mesh.face_handle(he);
  } while (curface.is_valid() && he != start);
  // Mark all the faces as traversed.
  for (FaceH fh : ptraverse.m_faces) {
    auto& indices = mesh.property(props.m_chord_indices, fh);
    if (indices[0] == -1) {
      indices[0] = index;
    }
    else if (indices[1] == -1) {
      indices[1] = index;
    }
  }
}

void polychord_collapse(Mesh& mesh)
{
  Props props(mesh);
  // Find all polychords.
  std::vector<PolychordInfo> chords;
  debug::clear("chord_err");
  PolychordTraverse ptraverse(mesh.n_vertices());  // temporary storage.
  for (FaceH fh : mesh.faces()) {
    const std::array<int, 2>& indices = mesh.property(props.m_chord_indices, fh);
    if (indices[0] == -1) {
      HalfH he = *mesh.cfh_begin(fh);
      walk_polychord(mesh, he, props, int(chords.size()), ptraverse);
      if (ptraverse.valid()) {
        // debug
        debug::write_faces(ptraverse.m_faces, "chord" + std::to_string(chords.size()));
        debug::append(ptraverse.m_err_q, "chord_err");
        // debug
        chords.emplace_back(he, ptraverse.m_err_q, std::sqrt(ptraverse.m_err_d2), 0.);
      }
    }
    if (indices[1] == -1) {
      HalfH he = mesh.next_halfedge_handle(*mesh.cfh_begin(fh));
      walk_polychord(mesh, he, props, int(chords.size()), ptraverse);
      if (ptraverse.valid()) {
        // debug
        debug::write_faces(ptraverse.m_faces, "chord" + std::to_string(chords.size()));
        debug::append(ptraverse.m_err_q, "chord_err");
        // debug
        chords.emplace_back(he, ptraverse.m_err_q, std::sqrt(ptraverse.m_err_d2), 0.);
      }
    }
  }
  // Debug-start
  std::cout << "Number of chords found: " << chords.size() << std::endl;
  // Debug-end
  throw std::logic_error("Not Implemented");
}

void simplify(Mesh& mesh)
{
  polychord_collapse(mesh);
  throw std::logic_error("Not Implemented");
}

}  // namespace quarrot
