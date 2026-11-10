// pcp_cli: runs the polycube pipeline (validate -> tetrahedralize -> l1 polycube -> structure
// -> height simplification) on one closed triangle mesh and writes all results into an output directory.
#include <cstdio>
#include <optional>
#include <string>

#include <CLI/CLI.hpp>

#include "pcp/core/version.h"
#include "pcp/io/text.h"
#include "pcp/pipeline/pipeline.h"

int main(int argc, char** argv) {
  CLI::App app{"Polycube construction pipeline (Huang et al. 2014 l1 polycube)"};
  app.set_version_flag("--version", pcp::version());

  std::string input, output, until = "simplify", resume, params_file, feature_edges,
                             linear_solver;
  bool quiet = false, verbose = false, no_rotation = false, no_simplify = false;
  std::optional<double> beta, alpha0, alpha_factor, eps0, eps_factor, threshold, omega,
      tet_size_factor, tet_max_volume, tet_ratio, max_stretch;
  std::optional<int> max_stages;

  app.add_option("input", input, "Input closed triangle mesh (.off)")->required();
  app.add_option("-o,--output", output, "Output directory")->required();
  app.add_option("--until", until, "Last stage to run: validate | tet | l1 | structure | simplify")
      ->capture_default_str();
  app.add_option("--resume-from", resume, "Start from a saved tet mesh (.tmsh) instead of tetrahedralizing");
  app.add_option("--params", params_file, "JSON parameter file (command-line options override it)");
  app.add_option("--beta", beta, "Complexity weight beta");
  app.add_option("--alpha0", alpha0, "Initial l1 weight alpha");
  app.add_option("--alpha-factor", alpha_factor, "alpha multiplier per stage");
  app.add_option("--eps0", eps0, "Initial l1 smoothing epsilon");
  app.add_option("--eps-factor", eps_factor, "epsilon multiplier per stage");
  app.add_option("--max-stages", max_stages, "Maximum number of schedule stages");
  app.add_option("--threshold", threshold, "Polycube error threshold");
  app.add_flag("--no-rotation", no_rotation, "Disable the global orientation search");
  app.add_option("--feature-edges", feature_edges, "Feature edge file (\"a b corner|flat\" per line)");
  app.add_option("--omega", omega, "Feature edge weight");
  app.add_option("--linear-solver", linear_solver, "auto | direct | cg");
  app.add_option("--tet-size-factor", tet_size_factor, "Tet edge length / average surface edge length");
  app.add_option("--tet-max-volume", tet_max_volume, "Maximum tet volume (overrides the size factor)");
  app.add_option("--tet-ratio", tet_ratio, "TetGen radius-edge ratio bound");
  app.add_flag("--no-simplify", no_simplify, "Skip the height simplification (Chen 2019)");
  app.add_option("--max-stretch", max_stretch, "Height simplification: maximum chain stretch ratio (>= 1)");
  auto* q = app.add_flag("-q,--quiet", quiet, "Only print errors");
  app.add_flag("-v,--verbose", verbose, "Print detailed progress")->excludes(q);

  CLI11_PARSE(app, argc, argv);

  try {
    pcp::PipelineParams params;
    if (!params_file.empty())
      params.update_from_json(nlohmann::json::parse(pcp::text::read_file(params_file)));
    if (beta) params.l1.beta = *beta;
    if (alpha0) params.l1.alpha0 = *alpha0;
    if (alpha_factor) params.l1.alpha_factor = *alpha_factor;
    if (eps0) params.l1.eps0 = *eps0;
    if (eps_factor) params.l1.eps_factor = *eps_factor;
    if (max_stages) params.l1.max_stages = *max_stages;
    if (threshold) params.l1.error_threshold = *threshold;
    if (no_rotation) params.l1.global_rotation = false;
    if (omega) params.l1.omega = *omega;
    if (!linear_solver.empty()) params.l1.linear_solver = linear_solver;
    params.l1.linear_options();  // validate
    if (!feature_edges.empty()) params.feature_edges = feature_edges;
    if (tet_size_factor) params.tet.size_factor = *tet_size_factor;
    if (tet_max_volume) params.tet.max_volume = *tet_max_volume;
    if (tet_ratio) params.tet.radius_edge_ratio = *tet_ratio;
    if (no_simplify) params.simplify.enabled = false;
    if (max_stretch) params.simplify.update_from_json({{"max_stretch", *max_stretch}});  // validates

    pcp::PipelineOptions opt;
    opt.input = input;
    opt.output_dir = output;
    opt.until = pcp::parse_stage(until);
    if (!resume.empty()) opt.resume_from = resume;
    opt.verbosity = quiet ? pcp::Verbosity::Quiet
                          : verbose ? pcp::Verbosity::Verbose : pcp::Verbosity::Normal;

    pcp::run_pipeline(opt, params, [](const std::string& s) {
      std::fputs(s.c_str(), stdout);
      std::fputc('\n', stdout);
      std::fflush(stdout);
    });
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
  return 0;
}
