#include <debug.h>
#include <quarrot.h>

#include <filesystem>
#include <fstream>
#include <iomanip>

namespace quarrot {
namespace debug {

void write_faces(std::span<const FaceH> faces, const fs::path& path)
{
  std::ofstream f(path, std::ios::ate);
  for (FaceH fh : faces) {
    f << fh.idx() << std::endl;
  }
}

void append(double val, const fs::path& path)
{
  std::ofstream f(path, std::ios::app);
  f << std::setprecision(14) << val << std::endl;
}

}  // namespace debug
}  // namespace quarrot
