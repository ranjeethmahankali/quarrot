#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>

#include <filesystem>

#include <OpenMesh/Core/IO/MeshIO.hh>

#include <quarrot.h>

using namespace quarrot;

TEST_CASE("Debugging", "[debugging]")
{
  std::filesystem::path fpath = "/home/rnjth94/buffer/parametrization/bimba.obj";
  Mesh                  mesh;
  REQUIRE(OpenMesh::IO::read_mesh(mesh, fpath.string()));
  size_t before = mesh.n_faces();
  std::cout << "Before: V " << mesh.n_vertices() << "; E " << mesh.n_edges() << "; F "
            << mesh.n_faces() << std::endl;
  pair_triangles(mesh, 0.1);
  std::cout << "After pairing: V " << mesh.n_vertices() << "; E " << mesh.n_edges()
            << "; F " << mesh.n_faces() << std::endl;
  REQUIRE(OpenMesh::IO::write_mesh(
    mesh, "/home/rnjth94/buffer/parametrization/bimba_paired.obj"));
  REQUIRE(mesh.n_faces() < before);
  subdivide(mesh);
  std::cout << "After subdivision: V " << mesh.n_vertices() << "; E " << mesh.n_edges()
            << "; F " << mesh.n_faces() << std::endl;
  REQUIRE(OpenMesh::IO::write_mesh(
    mesh, "/home/rnjth94/buffer/parametrization/bimba_subdivided.obj"));
  simplify(mesh);
}
