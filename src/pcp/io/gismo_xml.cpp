#include "pcp/io/gismo_xml.h"

#include <string>

#include "pcp/io/text.h"

namespace pcp {

GismoMultiPatch gismo_topology(const BlockLayout& L, std::vector<TensorBSpline3> volumes) {
  GismoMultiPatch mp;
  mp.volumes = std::move(volumes);
  for (size_t b = 0; b < L.blocks.size(); ++b)
    for (int f = 0; f < 6; ++f) {
      const BlockFace& bf = L.blocks[b].faces[static_cast<size_t>(f)];
      if (!bf.interior()) {
        mp.boundary.emplace_back(static_cast<int>(b), f + 1);
      } else if (bf.neighbor > static_cast<int>(b)) {
        std::array<int, 10> it = {static_cast<int>(b), f + 1, bf.neighbor, (f ^ 1) + 1, 0, 1, 2, 1, 1, 1};
        it[static_cast<size_t>(7 + f / 2)] = 0;
        mp.interfaces.push_back(it);
      }
    }
  return mp;
}

void write_gismo_xml(const std::filesystem::path& path, const GismoMultiPatch& mp) {
  std::string s = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<xml>\n";
  for (size_t g = 0; g < mp.volumes.size(); ++g) {
    const TensorBSpline3& V = mp.volumes[g];
    s += " <Geometry type=\"TensorBSpline3\" id=\"";
    text::append_int(s, static_cast<long long>(g));
    s += "\">\n  <Basis type=\"TensorBSplineBasis3\">\n";
    for (int d = 0; d < 3; ++d) {
      const KnotVector& kv = V.basis[static_cast<size_t>(d)];
      s += "   <Basis type=\"BSplineBasis\" index=\"";
      text::append_int(s, d);
      s += "\">\n    <KnotVector degree=\"";
      text::append_int(s, kv.degree);
      s += "\">";
      for (size_t i = 0; i < kv.knots.size(); ++i) {
        if (i) s += ' ';
        text::append_double(s, kv.knots[i]);
      }
      s += "</KnotVector>\n   </Basis>\n";
    }
    s += "  </Basis>\n  <coefs geoDim=\"3\">\n";
    for (const Eigen::Vector3d& c : V.coefs) {
      for (int k = 0; k < 3; ++k) {
        if (k) s += ' ';
        text::append_double(s, c(k));
      }
      s += '\n';
    }
    s += "  </coefs>\n </Geometry>\n";
  }
  s += " <MultiPatch parDim=\"3\" id=\"";
  text::append_int(s, static_cast<long long>(mp.volumes.size()));
  s += "\">\n  <patches type=\"id_range\">0 ";
  text::append_int(s, static_cast<long long>(mp.volumes.size()) - 1);
  s += "</patches>\n  <interfaces>\n";
  for (const auto& it : mp.interfaces) {
    for (size_t i = 0; i < it.size(); ++i) {
      if (i) s += ' ';
      text::append_int(s, it[i]);
    }
    s += '\n';
  }
  s += "  </interfaces>\n  <boundary>\n";
  for (const auto& [p, side] : mp.boundary) {
    text::append_int(s, p);
    s += ' ';
    text::append_int(s, side);
    s += '\n';
  }
  s += "  </boundary>\n </MultiPatch>\n</xml>\n";
  text::write_file(path, s);
}

namespace {

// Content between <tag ...> and </tag> starting at pos; advances pos past the closing tag.
std::string_view element(std::string_view s, std::string_view tag, size_t& pos, std::string_view* open = nullptr) {
  const std::string o = "<" + std::string(tag), c = "</" + std::string(tag) + ">";
  size_t a = pos;
  while (true) {
    a = s.find(o, a);
    if (a == std::string_view::npos) throw Error("G+Smo XML: missing <" + std::string(tag) + ">");
    const char next = s[a + o.size()];
    if (next == ' ' || next == '>') break;
    a += o.size();
  }
  const size_t gt = s.find('>', a);
  const size_t e = s.find(c, gt);
  if (gt == std::string_view::npos || e == std::string_view::npos)
    throw Error("G+Smo XML: unterminated <" + std::string(tag) + ">");
  if (open) *open = s.substr(a, gt - a);
  pos = e + c.size();
  return s.substr(gt + 1, e - gt - 1);
}

int attribute_int(std::string_view open, std::string_view name) {
  const std::string key = std::string(name) + "=\"";
  const size_t a = open.find(key);
  if (a == std::string_view::npos) throw Error("G+Smo XML: missing attribute " + std::string(name));
  text::Tokenizer t(open.substr(a + key.size(), open.find('"', a + key.size()) - a - key.size()));
  return static_cast<int>(t.expect_int(name.data()));
}

}  // namespace

GismoMultiPatch read_gismo_xml(const std::filesystem::path& path) {
  const std::string data = text::read_file(path);
  const std::string_view s = data;
  GismoMultiPatch mp;
  const size_t mp_pos = s.find("<MultiPatch");
  if (mp_pos == std::string_view::npos) throw Error("G+Smo XML: missing <MultiPatch>");
  size_t pos = 0;
  while (true) {
    const size_t g = s.find("<Geometry", pos);
    if (g == std::string_view::npos || g > mp_pos) break;
    pos = g;
    std::string_view open;
    const std::string_view geo = element(s, "Geometry", pos, &open);
    if (attribute_int(open, "id") != static_cast<int>(mp.volumes.size()))
      throw Error("G+Smo XML: geometry ids must be 0..N-1 in order");
    TensorBSpline3 V;
    size_t gp = 0;
    for (int d = 0; d < 3; ++d) {
      std::string_view kopen;
      const std::string_view knots = element(geo, "KnotVector", gp, &kopen);
      KnotVector& kv = V.basis[static_cast<size_t>(d)];
      kv.degree = attribute_int(kopen, "degree");
      text::Tokenizer t(knots);
      std::string_view tok;
      while (t.next(tok)) kv.knots.push_back(std::stod(std::string(tok)));
    }
    text::Tokenizer t(element(geo, "coefs", gp));
    const size_t n = static_cast<size_t>(V.n(0)) * V.n(1) * V.n(2);
    for (size_t i = 0; i < n; ++i) {
      Eigen::Vector3d c;
      for (int k = 0; k < 3; ++k) c(k) = t.expect_double("coefficient");
      V.coefs.push_back(c);
    }
    std::string_view extra;
    if (t.next(extra)) throw Error("G+Smo XML: too many coefficients");
    mp.volumes.push_back(std::move(V));
  }
  size_t p = mp_pos;
  {
    text::Tokenizer t(element(s, "interfaces", p));
    std::string_view tok;
    while (t.next(tok)) {
      std::array<int, 10> it{};
      it[0] = std::stoi(std::string(tok));
      for (size_t i = 1; i < 10; ++i) it[i] = static_cast<int>(t.expect_int("interface entry"));
      mp.interfaces.push_back(it);
    }
  }
  {
    text::Tokenizer t(element(s, "boundary", p));
    std::string_view tok;
    while (t.next(tok)) {
      const int patch = std::stoi(std::string(tok));
      mp.boundary.emplace_back(patch, static_cast<int>(t.expect_int("boundary side")));
    }
  }
  return mp;
}

}  // namespace pcp
