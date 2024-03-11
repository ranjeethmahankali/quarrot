#include <debug.h>
#include <quarrot.h>

#include <filesystem>
#include <fstream>
#include <iomanip>

namespace quarrot {

namespace fs = std::filesystem;

namespace debug {

static fs::path from_label(const std::string& label)
{
  return fs::path("temp/" + label + ".txt");
}

void write_faces(std::span<const FaceH> faces, const std::string& label)
{
  std::ofstream f(from_label(label), std::ios::ate);
  for (FaceH fh : faces) {
    f << fh.idx() << std::endl;
  }
}

void append(double val, const std::string& label)
{
  std::ofstream f(from_label(label), std::ios::app);
  f << std::setprecision(14) << val << std::endl;
}

void clear(const std::string& label)
{
  std::ofstream f(from_label(label), std::ios::ate);
  // Truncalte the file and close the file. Nothing else to do.
}

size_t count_singularities(const Mesh& mesh)
{
  return std::count_if(mesh.vertices_begin(), mesh.vertices_end(), [&mesh](VertH vh) {
    return mesh.valence(vh) != 4;
  });
}

void copy_test(Mesh mesh)
{
  mesh.garbage_collection();
  Mesh copy;
  copy.reserve(mesh.n_vertices(), mesh.n_edges(), mesh.n_faces());
  for (VertH vh : mesh.vertices()) {
    copy.add_vertex(mesh.point(vh));
  }
  std::vector<VertH> fvs;
  for (FaceH fh : mesh.faces()) {
    fvs.clear();
    std::copy(mesh.cfv_begin(fh), mesh.cfv_end(fh), std::back_inserter(fvs));
    copy.add_face(fvs.data(), fvs.size());
  }
}

}  // namespace debug
}  // namespace quarrot
