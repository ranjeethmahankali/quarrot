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
  OpenMesh::FPropHandleT<std::array<int, 2>> chordIndices;
  OpenMesh::FPropHandleT<glm::dvec3>         faceCenters;
  OpenMesh::VPropHandleT<Quadric>            vertQuadric;
  Mesh&                                      mMesh;

  explicit Props(Mesh& mesh)
      : mMesh(mesh)
  {
    mesh.add_property(chordIndices);
    for (FaceH fh : mesh.faces()) {  // Initialize all indices to -1.
      mesh.property(chordIndices, fh).fill(-1);
    }
    mesh.add_property(faceCenters);
    calc_face_centroids(mesh, faceCenters);
    mesh.add_property(vertQuadric);
    calc_vert_quadrics(mesh, vertQuadric, faceCenters);
  }

  ~Props()
  {
    mMesh.remove_property(chordIndices);
    mMesh.remove_property(faceCenters);
    mMesh.remove_property(vertQuadric);
  }
};

struct PolychordTraverse
{
  std::vector<FaceH>   faces;
  DisjointSets         vsets;
  std::vector<Quadric> quadrics;
  double               max_err = 0.;

  explicit PolychordTraverse(size_t numVertices)
      : vsets(numVertices)
  {}

  bool valid() const { return !faces.empty(); }

  void copy_quadrics(Mesh& mesh, OpenMesh::VPropHandleT<Quadric> qprop)
  {
    quadrics.clear();
    quadrics.reserve(mesh.n_vertices());
    const std::vector<Quadric>& src = mesh.property(qprop).data_vector();
    std::copy(src.begin(), src.end(), std::back_inserter(quadrics));
  }

  void unite_quadrics(int a, int b)
  {
    Quadric q   = quadrics[a] + quadrics[b];
    quadrics[a] = q;
    quadrics[b] = q;
    max_err     = std::max(max_err, q(q.minimizer()));
  }

  void reset(Mesh& mesh, const Props& props)
  {
    faces.clear();
    vsets.reset();
    copy_quadrics(mesh, props.vertQuadric);
    max_err = 0.;
  }
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
    auto& indices = mesh.property(props.chordIndices, curface);
    ptraverse.faces.push_back(curface);
    // Unite the vertices as if the halfedge is being collapsed.
    VertH v0 = mesh.from_vertex_handle(he);
    VertH v1 = mesh.to_vertex_handle(he);
    ptraverse.vsets.unite(uint32_t(v0.idx()), uint32_t(v1.idx()));
    ptraverse.unite_quadrics(v0.idx(), v1.idx());
    if (indices[0] != -1 && indices[1] != -1) {
      ptraverse.faces.clear();
      break;
    }
    he = mesh.opposite_halfedge_handle(
      mesh.next_halfedge_handle(mesh.next_halfedge_handle(he)));
    curface = mesh.face_handle(he);
  } while (curface.is_valid() && he != start);
  // Mark all the faces as traversed.
  for (FaceH fh : ptraverse.faces) {
    auto& indices = mesh.property(props.chordIndices, fh);
    if (indices[0] == -1) {
      indices[0] = index;
    }
    else if (indices[1] == -1) {
      indices[1] = index;
    }
  }
}

struct PolychordInfo
{
  HalfH  start_he;
  double err_q = 0.;
  double err_d = 0.;
  double err_v = 0.;

  PolychordInfo(HalfH he, double eq, double ed, double ev)
      : start_he(he)
      , err_q(eq)
      , err_d(ed)
      , err_v(ev)
  {}
};

void polychord_collapse(Mesh& mesh)
{
  Props props(mesh);
  // Find all polychords.
  std::vector<PolychordInfo> chords;
  PolychordTraverse          ptraverse(mesh.n_vertices());  // temporary storage.
  for (FaceH fh : mesh.faces()) {
    const std::array<int, 2>& indices = mesh.property(props.chordIndices, fh);
    if (indices[0] == -1) {
      HalfH he = *mesh.cfh_begin(fh);
      walk_polychord(mesh, he, props, int(chords.size()), ptraverse);
      if (ptraverse.valid()) {
        // debug
        // debug
        chords.emplace_back(he, 0., ptraverse.max_err, 0.);
      }
    }
    if (indices[1] == -1) {
      HalfH he = mesh.next_halfedge_handle(*mesh.cfh_begin(fh));
      walk_polychord(mesh, he, props, int(chords.size()), ptraverse);
      if (ptraverse.valid()) {
        chords.emplace_back(he, 0., ptraverse.max_err, 0.);
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
