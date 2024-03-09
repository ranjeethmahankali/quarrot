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

}  // namespace debug
}  // namespace quarrot
