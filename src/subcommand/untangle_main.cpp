#include "algorithms/untangle.hpp"
#include "args.hxx"
#include "odgi.hpp"
#include "subcommand.hpp"
#include "utils.hpp"
#include <omp.h>

namespace odgi {

using namespace odgi::subcommand;

int main_untangle(int argc, char **argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc - 1; ++i) {
        argv[i] = argv[i + 1];
    }
    const std::string prog_name = "odgi untangle";
    argv[0] = (char *)prog_name.c_str();
    --argc;

    args::ArgumentParser parser("Project paths into reference-relative BEDPE (optionally PAF), to decompose paralogy relationships.");
    args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
    args::ValueFlag<std::string> og_in_file(
        mandatory_opts, "FILE",
        "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually "
        "ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format "
        "requires additional time!",
        {'i', "idx"});
    args::Group untangling_opts(parser, "[ Untangling Options ]");
    args::ValueFlag<std::string> _query_path(untangling_opts, "NAME", "Use this query path.",
                       {'q', "query-path"});
    args::ValueFlag<std::string> _target_path(untangling_opts, "NAME", "Use this target (reference) path.",
                       {'r', "target-path"});
    args::ValueFlag<std::string> _query_paths(untangling_opts, "FILE", "Use query paths listed (one per line) in FILE.",
                       {'Q', "query-paths"});
    args::ValueFlag<std::string> _target_paths(untangling_opts, "FILE", "Use target (reference) paths listed (one per line) in FILE.",
                       {'R', "target-paths"});
    args::ValueFlag<uint64_t> merge_dist(untangling_opts, "N", "Merge segments shorter than this length into previous segments.",
                                         {'m', "merge-dist"});
    args::ValueFlag<double> _max_self_coverage(untangling_opts, "N", "Skip mappings with greater than this level of self-coverage.",
                                               {'s', "max-self-coverage"});
    args::ValueFlag<uint64_t> _best_n_mappings(untangling_opts, "N", "Report up to the Nth best target (reference) mapping for each query segment (default: 1).",
                                               {'n', "n-best"});
    args::ValueFlag<double> _jaccard_threshold(untangling_opts, "F", "Report target mappings >= the given jaccard threshold, with 0 <= F <= 1.0 (default: 0.0).",
                                               {'j', "min-jaccard"});
    args::ValueFlag<uint64_t> _cut_every(untangling_opts, "N", "Start a segment boundary every Nbp of the sorted graph (default: 0/off).",
                                         {'e', "cut-every"});
    args::Flag paf_output(untangling_opts, "paf_output", "Emit the output in PAF format.",
                          {'p', "paf-output"});
    args::Flag gene_order_output(untangling_opts, "gene_order_output", "Write each query as a series of target gene segments.",
                                 {'G', "gene-order"});
    args::Flag gggenes_output(untangling_opts, "gggenes_output", "Emit the output in gggenes-compatible tabular format.",
                              {'g', "gggenes-output"});
    args::Flag gggenes_schematic(untangling_opts, "gggenes_schematic", "Emit the output in gggenes-compatible *schematic* tabular format, where each gene is rendered as 100bp.",
                                 {'X', "gggenes-schematic"});
    args::ValueFlag<std::string> input_cut_points(untangling_opts, "FILE", "A text file of node identifiers (one identifier per row) where to start the segment boundaries."
                                                                           "When specified, no further starting points will be added.", {'c', "cut-points-input"});
    args::ValueFlag<std::string> output_cut_points(untangling_opts, "FILE", "Emit node identifiers where segment boundaries started (one identifier per row).",
                                                  {'d', "cut-points-output"});
    args::Group debugging_opts(parser, "[ Debugging Options ]");
    args::Flag make_self_dotplot(debugging_opts, "DOTPLOT", "Render a table showing the positional dotplot of the query against itself.",
                                 {'S', "self-dotplot"});
	args::Group step_index_opts(parser, "[ Step Index Options ]");
	args::ValueFlag<std::string> _step_index(step_index_opts, "FILE", "Load the step index from this *FILE*. The file name usually ends with *.stpidx*. (default: build the step index from scratch with sampling rate 8).",
											 {'a', "step-index"});
	args::ValueFlag<uint64_t> _step_sampling(step_index_opts, "N", "Sampling rate to use for step index construction. (default: 8).",
											 {'D', "step-sampling"});
    args::Group threading(parser, "[ Threading ]");
    args::ValueFlag<uint64_t> nthreads(
        threading, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
    args::Group processing_info_opts(parser, "[ Processing Information ]");
    args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.",
                        {'P', "progress"});
    args::Group program_info_opts(parser, "[ Program Information ]");
    args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi untangle.",
                        {'h', "help"});

    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    if (argc == 1) {
        std::cout << parser;
        return 1;
    }

    if (!og_in_file) {
        std::cerr << "[odgi::untangle] error: please specify an input file from where to load the "
                     "graph via -i=[FILE], --idx=[FILE]."
                  << std::endl;
        return 1;
    }

    const uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;

    graph_t graph;
    assert(argc > 0);
    {
        const std::string infile = args::get(og_in_file);
        if (!infile.empty()) {
            if (infile == "-") {
                graph.deserialize(std::cin);
            } else {
                utils::handle_gfa_odgi_input(infile, "untangle", args::get(progress), num_threads,
                                             graph);
            }
        }
        graph.set_static();
    }

    // path loading
    auto load_paths = [&](const std::string& path_names_file) {
        std::ifstream path_names_in(path_names_file);
        uint64_t num_of_paths_in_file = 0;
        std::vector<bool> path_already_seen;
        path_already_seen.resize(graph.get_path_count(), false);
        std::string line;
        std::vector<path_handle_t> paths;
        while (std::getline(path_names_in, line)) {
            if (!line.empty()) {
                if (graph.has_path(line)) {
                    const path_handle_t path = graph.get_path_handle(line);
                    const uint64_t path_rank = as_integer(path) - 1;
                    if (!path_already_seen[path_rank]) {
                        path_already_seen[path_rank] = true;
                        paths.push_back(path);
                    } else {
                        std::cerr << "[odgi::untangle] error: in the path list there are duplicated path names."
                                  << std::endl;
                        exit(1);
                    }
                }
                ++num_of_paths_in_file;
            }
        }
        path_names_in.close();
        std::cerr << "[odgi::untangle] found " << paths.size() << "/" << num_of_paths_in_file
                  << " paths to consider." << std::endl;
        if (paths.empty()) {
            std::cerr << "[odgi::untangle] error: no path to consider." << std::endl;
            exit(1);
        }
        return paths;
    };

    std::vector<path_handle_t> target_paths;
    std::vector<path_handle_t> query_paths;
    if (_target_path) {
        auto& path_name = args::get(_target_path);
        if (graph.has_path(path_name)) {
            target_paths.push_back(graph.get_path_handle(path_name));
        } else {
            std::cerr << "[odgi::untangle] error: no path "
                      << path_name << " found in graph." << std::endl;
            exit(1);
        }
    } else if (_target_paths) {
        target_paths = load_paths(args::get(_target_paths));
    } else {
        target_paths.reserve(graph.get_path_count());
        graph.for_each_path_handle([&](const path_handle_t path) {
            target_paths.push_back(path);
        });
    }
    if (_query_path) {
        auto& path_name = args::get(_query_path);
        if (graph.has_path(path_name)) {
            query_paths.push_back(graph.get_path_handle(path_name));
        } else {
            std::cerr << "[odgi::untangle] error: no path "
                      << path_name << " found in graph." << std::endl;
            exit(1);
        }
    } else if (_query_paths) {
        query_paths = load_paths(args::get(_query_paths));
    } else {
        query_paths.reserve(graph.get_path_count());
        graph.for_each_path_handle([&](const path_handle_t path) {
            query_paths.push_back(path);
        });
    }

	std::vector<path_handle_t> paths;
	paths.insert(paths.end(), query_paths.begin(), query_paths.end());
	paths.insert(paths.end(), target_paths.begin(), target_paths.end());
	std::sort(paths.begin(), paths.end());
	paths.erase(std::unique(paths.begin(), paths.end()), paths.end());

    algorithms::untangle_output_t output_type = algorithms::untangle_output_t::BEDPE;
    if (args::get(paf_output)) {
        output_type = algorithms::untangle_output_t::PAF;
    } else if (args::get(gene_order_output)) {
        output_type = algorithms::untangle_output_t::ORDER;
    } else if (args::get(gggenes_output)) {
        output_type = algorithms::untangle_output_t::GGGENES;
    } else if (args::get(gggenes_schematic)) {
        output_type = algorithms::untangle_output_t::SCHEMATIC;
    }


    if (make_self_dotplot) {
        for (auto& query : query_paths) {
            algorithms::self_dotplot(graph, query);
        }
    } else {
		if (!_step_index) {
            auto step_sampling = _step_sampling ? args::get(_step_sampling) : 8;
			if (progress) {
				std::cerr
                    << "[odgi::untangle] warning: no step index specified. Building one with a sample rate of " << step_sampling << ". This may take additional time. "
                    "A step index can be provided via -a, --step-index. A step index can be built using odgi stepindex."
                    << std::endl;
			}
			algorithms::step_index_t step_index(graph, paths, num_threads, progress, step_sampling);
			algorithms::untangle(graph,
								 query_paths,
								 target_paths,
								 args::get(merge_dist),
								 (_max_self_coverage ? args::get(_max_self_coverage) : 0),
								 (_best_n_mappings ? args::get(_best_n_mappings) : 1),
								 (_jaccard_threshold ? args::get(_jaccard_threshold) : 0.0),
								 (_cut_every ? args::get(_cut_every) : 0),
                                 output_type,
								 args::get(input_cut_points),
								 args::get(output_cut_points),
								 num_threads,
								 progress,
								 step_index,
								 paths);
		} else {
			algorithms::step_index_t step_index;
			step_index.load(args::get(_step_index));
			algorithms::untangle(graph,
								 query_paths,
								 target_paths,
								 args::get(merge_dist),
								 (_max_self_coverage ? args::get(_max_self_coverage) : 0),
								 (_best_n_mappings ? args::get(_best_n_mappings) : 1),
								 (_jaccard_threshold ? args::get(_jaccard_threshold) : 0.0),
								 (_cut_every ? args::get(_cut_every) : 0),
                                 output_type,
								 args::get(input_cut_points),
								 args::get(output_cut_points),
								 num_threads,
								 progress,
								 step_index,
                                 paths);
		}
    }

    return 0;
}

static Subcommand odgi_untangle("untangle", "Project paths into reference-relative, to decompose paralogy relationships.", PIPELINE, 3, main_untangle);

} // namespace odgi
