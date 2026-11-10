#pragma once

#include <cstdio>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "pcp/core/types.h"

namespace pcp::text {

// Number formatting/parsing helpers. They assume the process stays in the default "C"
// locale (the library never calls setlocale), so output is identical on all platforms.
inline void append_double(std::string& out, double v) {
  char buf[32];
  int n = std::snprintf(buf, sizeof(buf), "%.17g", v);
  out.append(buf, static_cast<size_t>(n));
}
inline void append_int(std::string& out, long long v) {
  char buf[24];
  int n = std::snprintf(buf, sizeof(buf), "%lld", v);
  out.append(buf, static_cast<size_t>(n));
}

// Reads a whole file into memory (binary mode). Throws pcp::Error if it cannot be opened.
std::string read_file(const std::filesystem::path& path);
// Writes the string as-is in binary mode (always '\n' line endings).
void write_file(const std::filesystem::path& path, std::string_view data);

// Whitespace tokenizer that skips '#' comments to end of line.
class Tokenizer {
 public:
  explicit Tokenizer(std::string_view s) : s_(s) {}
  bool next(std::string_view& tok);
  // Throwing helpers; `what` names the expected item for the error message.
  std::string_view expect(const char* what);
  double expect_double(const char* what);
  long long expect_int(const char* what);
  size_t line() const { return line_; }

 private:
  std::string_view s_;
  size_t pos_ = 0;
  size_t line_ = 1;
};

}  // namespace pcp::text
