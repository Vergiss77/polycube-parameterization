#include "pcp/io/text.h"

#include <cerrno>
#include <cstdlib>
#include <fstream>
#include <sstream>

namespace pcp::text {

std::string read_file(const std::filesystem::path& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw Error("cannot open file: " + path.string());
  std::ostringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

void write_file(const std::filesystem::path& path, std::string_view data) {
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out) throw Error("cannot write file: " + path.string());
  out.write(data.data(), static_cast<std::streamsize>(data.size()));
  if (!out) throw Error("failed while writing file: " + path.string());
}

bool Tokenizer::next(std::string_view& tok) {
  const size_t n = s_.size();
  while (pos_ < n) {
    char c = s_[pos_];
    if (c == '#') {
      while (pos_ < n && s_[pos_] != '\n') ++pos_;
    } else if (c == '\n') {
      ++line_;
      ++pos_;
    } else if (c == ' ' || c == '\t' || c == '\r') {
      ++pos_;
    } else {
      break;
    }
  }
  if (pos_ >= n) return false;
  size_t start = pos_;
  while (pos_ < n) {
    char c = s_[pos_];
    if (c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '#') break;
    ++pos_;
  }
  tok = s_.substr(start, pos_ - start);
  return true;
}

std::string_view Tokenizer::expect(const char* what) {
  std::string_view tok;
  if (!next(tok)) throw Error(std::string("unexpected end of file, expected ") + what);
  return tok;
}

double Tokenizer::expect_double(const char* what) {
  std::string tok(expect(what));
  char* end = nullptr;
  errno = 0;
  double v = std::strtod(tok.c_str(), &end);
  if (end != tok.c_str() + tok.size() || errno == ERANGE) {
    throw Error("line " + std::to_string(line_) + ": invalid number '" + tok + "' for " + what);
  }
  return v;
}

long long Tokenizer::expect_int(const char* what) {
  std::string tok(expect(what));
  char* end = nullptr;
  errno = 0;
  long long v = std::strtoll(tok.c_str(), &end, 10);
  if (end != tok.c_str() + tok.size() || errno == ERANGE) {
    throw Error("line " + std::to_string(line_) + ": invalid integer '" + tok + "' for " + what);
  }
  return v;
}

}  // namespace pcp::text
