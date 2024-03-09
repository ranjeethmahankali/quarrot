#pragma once
#include <quarrot.h>

#include <filesystem>
#include <span>

namespace quarrot {

namespace debug {

void   write_faces(std::span<const FaceH> faces, const std::string& label);
void   append(double val, const std::string& label);
void   clear(const std::string& label);
size_t count_singularities(const Mesh& mesh);

}  // namespace debug
}  // namespace quarrot
