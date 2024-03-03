#pragma once
#include <quarrot.h>

#include <filesystem>
#include <span>

namespace quarrot {
namespace debug {

namespace fs = std::filesystem;

void write_faces(std::span<const FaceH> faces, const fs::path& path);
void append(double val, const fs::path& path);

}  // namespace debug
}  // namespace quarrot
