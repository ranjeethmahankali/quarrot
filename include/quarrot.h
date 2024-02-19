#pragma once

#include <open_mesh_adaptor.h>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Traits.hh>

namespace quarrot {
struct MeshTraits : public OpenMesh::DefaultTraits
{
  typedef glm::dvec3 Point;
  typedef glm::dvec3 Normal;

  VertexAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal);
  HalfedgeAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::PrevHalfedge);
  EdgeAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal);
  FaceAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal);
};

using Mesh  = OpenMesh::PolyMesh_ArrayKernelT<MeshTraits>;
using VertH = Mesh::VertexHandle;
using HalfH = Mesh::HalfedgeHandle;
using EdgeH = Mesh::EdgeHandle;
using FaceH = Mesh::FaceHandle;

void simplify(Mesh& mesh);

void pair_triangles(Mesh& mesh, double gamma = 0.2);

void subdivide(Mesh& mesh);

}  // namespace quarrot
