#include "pcp/polycube/feature_edges.h"

#include <algorithm>
#include <map>

#include "pcp/io/text.h"

namespace pcp {

std::vector<FeatureEdgeEnergy::Edge> read_feature_edges(const std::filesystem::path& path,
                                                        const BoundaryTopology& topo) {
  std::map<std::pair<int, int>, int> lookup;
  for (size_t i = 0; i < topo.pairs.size(); ++i) {
    const FacePair& p = topo.pairs[i];
    lookup[{std::min(p.v0, p.v1), std::max(p.v0, p.v1)}] = static_cast<int>(i);
  }
  const std::string data = text::read_file(path);
  text::Tokenizer tk(data);
  std::vector<FeatureEdgeEnergy::Edge> edges;
  std::string_view tok;
  try {
    while (tk.next(tok)) {
      long long a = std::stoll(std::string(tok));
      long long b = tk.expect_int("second edge vertex");
      std::string_view mode = tk.expect("edge mode");
      auto it = lookup.find({static_cast<int>(std::min(a, b)), static_cast<int>(std::max(a, b))});
      if (it == lookup.end())
        throw Error("line " + std::to_string(tk.line()) + ": (" + std::to_string(a) + ", " +
                    std::to_string(b) + ") is not an edge of the surface");
      FeatureEdgeEnergy::Edge e;
      e.pair = it->second;
      if (mode == "corner")
        e.target = 0.0;
      else if (mode == "flat")
        e.target = 1.0;
      else
        throw Error("line " + std::to_string(tk.line()) + ": unknown edge mode '" +
                    std::string(mode) + "' (expected corner or flat)");
      edges.push_back(e);
    }
  } catch (const std::invalid_argument&) {
    throw Error(path.string() + ": invalid vertex index near line " + std::to_string(tk.line()));
  } catch (const Error& e) {
    throw Error(path.string() + ": " + e.what());
  }
  return edges;
}

}  // namespace pcp
