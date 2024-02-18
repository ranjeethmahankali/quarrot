#pragma once

#include <OpenMeshAdaptor.h>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Traits.hh>

namespace quarrot {
struct MeshTraits : public OpenMesh::DefaultTraits
{
  typedef glm::dvec3 Point;
  typedef glm::dvec3 Normal;
  typedef double     TexCoord1D;
  typedef glm::dvec2 TexCoord2D;
  typedef glm::dvec3 TexCoord3D;
  typedef int        TextureIndex;
  typedef glm::dvec3 Color;

  VertexAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal);
  HalfedgeAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::PrevHalfedge);
  EdgeAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal);
  FaceAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal);
};

using Mesh = OpenMesh::PolyMesh_ArrayKernelT<MeshTraits>;

void simplify(Mesh& mesh);

}  // namespace quarrot
