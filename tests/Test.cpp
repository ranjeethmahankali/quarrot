#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>

#include <filesystem>

#include <OpenMesh/Core/IO/MeshIO.hh>

#include <Quarrot.h>

using namespace quarrot;

TEST_CASE("Simplify", "[simplify]")
{
  std::filesystem::path fpath = "/home/rnjth94/buffer/parametrization/bimba.obj";
  Mesh                  mesh;
  REQUIRE(OpenMesh::IO::read_mesh(mesh, fpath.string()));
  size_t before = mesh.n_faces();
  simplify(mesh);
  REQUIRE(mesh.n_faces() < before);
}
