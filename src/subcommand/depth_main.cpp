#include "subcommand.hpp"
#include "odgi.hpp"
#include "position.hpp"
#include "args.hxx"
#include "split.hpp"
#include "algorithms/bfs.hpp"
#include <omp.h>

namespace odgi {

    using namespace odgi::subcommand;

    int main_depth(int argc, char **argv) {

        // trick argumentparser to do the right thing with the subcommand
        for (uint64_t i = 1; i < argc - 1; ++i) {
            argv[i] = argv[i + 1];
        }
        std::string prog_name = "odgi depth";
        argv[0] = (char *) prog_name.c_str();
        --argc;

        args::ArgumentParser parser("find the depth of graph as defined by query criteria");
        args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
        args::ValueFlag<std::string> og_file(parser, "FILE", "compute path depths in this graph", {'i', "input"});
        args::ValueFlag<std::string> path_name(parser, "PATH_NAME", "compute the depth of the given path in the graph",
                                               {'r', "path"});
        args::ValueFlag<std::string> path_file(parser, "FILE", "compute depth for the paths listed in FILE",
                                               {'R', "paths"});
        args::ValueFlag<std::string> graph_pos(parser, "[node_id][,offset[,(+|-)]*]*",
                                               "compute the depth at the given node, e.g. 7 or 3,4 or 42,10,+ or 302,0,-",
                                               {'g', "graph-pos"});
        args::ValueFlag<std::string> graph_pos_file(parser, "FILE", "a file with one graph position per line",
                                                    {'G', "graph-pos-file"});
        args::ValueFlag<std::string> path_pos(parser, "[path_name][,offset[,(+|-)]*]*",
                                              "return depth at the given path position e.g. chrQ or chr3,42 or chr8,1337,+ or chrZ,3929,-",
                                              {'p', "path-pos"});
        args::ValueFlag<std::string> path_pos_file(parser, "FILE", "a file with one path position per line",
                                                   {'F', "path-pos-file"});
        args::ValueFlag<std::string> bed_input(parser, "FILE", "a BED file of ranges in paths in the graph",
                                               {'b', "bed-input"});

        // TODO: (or subset covered by the given paths)
        args::Flag graph_depth(parser, "graph-depth",
                               "compute the depth on each node in the graph (or subset covered by the given paths)",
                               {'d', "graph-depth"});
        args::ValueFlag<uint64_t> nthreads(parser, "N", "number of threads to use", {'t', "threads"});

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

        if (!og_file) {
            std::cerr << "[odgi::depth] error: please specify a target graph via -i=[FILE], --idx=[FILE]." << std::endl;
            return 1;
        }

        odgi::graph_t graph;
        assert(argc > 0);
        std::string infile = args::get(og_file);
        if (!infile.empty()) {
            if (infile == "-") {
                graph.deserialize(std::cin);
            } else {
                ifstream f(infile.c_str());
                graph.deserialize(f);
                f.close();
            }
        }

        uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;
        omp_set_num_threads(num_threads);

        // these options are exclusive (probably we should say with a warning)
        std::vector<odgi::pos_t> graph_positions;
        std::vector<odgi::path_pos_t> path_positions;
        std::vector<odgi::path_range_t> path_ranges;

        auto add_graph_pos = [&graph_positions](const odgi::graph_t &graph,
                                                const std::string &buffer) {
            auto vals = split(buffer, ',');
            /*
            if (vals.size() != 3) {
                std::cerr << "[odgi::depth] error: graph position record is incomplete" << std::endl;
                std::cerr << "[odgi::depth] error: got '" << buffer << "'" << std::endl;
                exit(1); // bail
            }
            */
            uint64_t id = std::stoi(vals[0]);
            if (!graph.has_node(id)) {
                std::cerr << "[odgi::depth] error: no node " << id << " in graph" << std::endl;
                exit(1);
            }
            uint64_t offset = 0;
            if (vals.size() >= 2) {
                offset = std::stoi(vals[1]);
                handle_t h = graph.get_handle(id);
                if (graph.get_length(h) < offset) {
                    std::cerr << "[odgi::depth] error: offset of " << offset << " lies beyond the end of node " << id
                              << std::endl;
                    exit(1);
                }
            }
            bool is_rev = false;
            if (vals.size() == 3) {
                is_rev = vals[2] == "-";
            }
            graph_positions.push_back(make_pos_t(id, is_rev, offset));
        };

        auto add_path_pos = [&path_positions](const odgi::graph_t &graph,
                                              const std::string &buffer) {
            auto vals = split(buffer, ',');
            /*
            if (vals.size() != 3) {
                std::cerr << "[odgi::depth] error: path position record is incomplete" << std::endl;
                std::cerr << "[odgi::depth] error: got '" << buffer << "'" << std::endl;
                exit(1); // bail
            }
            */
            auto &path_name = vals[0];
            if (!graph.has_path(path_name)) {
                std::cerr << "[odgi::depth] error: ref path " << path_name << " not found in graph" << std::endl;
                exit(1);
            } else {
                path_positions.push_back({
                                                 graph.get_path_handle(path_name),
                                                 (vals.size() > 1 ? (uint64_t) std::stoi(vals[1]) : 0),
                                                 (vals.size() == 3 ? vals[2] == "-" : false)
                                         });
            }
        };

        auto add_bed_range = [&path_ranges](const odgi::graph_t &graph,
                                            const std::string &buffer) {
            auto vals = split(buffer, '\t');
            /*
            if (vals.size() != 3) {
                std::cerr << "[odgi::depth] error: path position record is incomplete" << std::endl;
                std::cerr << "[odgi::depth] error: got '" << buffer << "'" << std::endl;
                exit(1); // bail
            }
            */
            auto &path_name = vals[0];
            if (!graph.has_path(path_name)) {
                std::cerr << "[odgi::depth] error: ref path " << path_name << " not found in graph" << std::endl;
                exit(1);
            } else {
                uint64_t start = vals.size() == 3 ? (uint64_t) std::stoi(vals[1]) : 0;
                uint64_t end = 0;
                if (vals.size() == 3) {
                    end = (uint64_t) std::stoi(vals[2]);
                } else {
                    graph.for_each_step_in_path(graph.get_path_handle(path_name), [&](const step_handle_t &s) {
                        end += graph.get_length(graph.get_handle_of_step(s));
                    });
                }

                path_ranges.push_back(
                        {
                                {
                                        graph.get_path_handle(path_name),
                                        start,
                                        false
                                },
                                {
                                        graph.get_path_handle(path_name),
                                        end,
                                        false
                                },
                                (vals.size() > 3 && vals[3] == "-"),
                                buffer
                        });
            }
        };

        if (graph_depth) {
            graph.for_each_handle([&](const handle_t& h) {
                add_graph_pos(graph, std::to_string(graph.get_id(h)));
            });
        } else if (graph_pos) {
            // if we're given a graph_pos, we'll convert it into a path pos
            add_graph_pos(graph, args::get(graph_pos));
        } else if (graph_pos_file) {
            std::ifstream gpos(args::get(graph_pos_file).c_str());
            std::string buffer;
            while (std::getline(gpos, buffer)) {
                add_graph_pos(graph, buffer);
            }
        } else if (path_pos) {
            // if given a path pos, we convert it into a path pos in our reference set
            add_path_pos(graph, args::get(path_pos));
        } else if (path_pos_file) {
            // if we're given a file of path depths, we'll convert them all
            std::ifstream refs(args::get(path_pos_file).c_str());
            std::string buffer;
            while (std::getline(refs, buffer)) {
                add_path_pos(graph, buffer);
            }
        } else if (bed_input) {
            std::ifstream bed_in(args::get(bed_input).c_str());
            std::string buffer;
            while (std::getline(bed_in, buffer)) {
                add_bed_range(graph, buffer);
            }
        } else if (path_name) {
            add_bed_range(graph, args::get(path_name));
        } else if (path_file) {
            // for thing in things
            std::ifstream refs(args::get(path_file).c_str());
            std::string path_name;
            while (std::getline(refs, path_name)) {
                add_bed_range(graph, path_name);
            }
        }/* else {
            // using all the paths in the graph
            graph.for_each_path_handle([&](const path_handle_t& path) { paths_to_consider.push_back(path); });
        }*/

        auto get_graph_pos = [](const odgi::graph_t &graph,
                                const path_pos_t &pos) {
            auto path_end = graph.path_end(pos.path);
            uint64_t walked = 0;
            for (step_handle_t s = graph.path_begin(pos.path);
                 s != path_end; s = graph.get_next_step(s)) {
                handle_t h = graph.get_handle_of_step(s);
                uint64_t node_length = graph.get_length(h);
                if (walked + node_length > pos.offset) {
                    return make_pos_t(graph.get_id(h), graph.get_is_reverse(h), pos.offset - walked);
                }
                walked += node_length;
            }

#pragma omp critical (cout)
            std::cerr << "[odgi::depth] warning: position " << graph.get_path_name(pos.path) << ":" << pos.offset
                      << " outside of path" << std::endl;
            return make_pos_t(0, false, 0);
        };

        auto get_offset_in_path = [](const odgi::graph_t &graph,
                                     const path_handle_t &path, const step_handle_t &target) {
            auto path_end = graph.path_end(path);
            uint64_t walked = 0;
            step_handle_t s = graph.path_begin(path);
            for (; s != target; s = graph.get_next_step(s)) {
                handle_t h = graph.get_handle_of_step(s);
                walked += graph.get_length(h);
            }
            assert(s != path_end);
            return walked;
        };

        auto get_graph_node_coverage = [](const odgi::graph_t &graph, const nid_t node_id) {
            uint64_t node_coverage = 0;
            std::set<uint64_t> unique_paths;

            handle_t h = graph.get_handle(node_id);

            graph.for_each_step_on_handle(h, [&](const step_handle_t &occ) {
                ++node_coverage;
                unique_paths.insert(as_integer(graph.get_path(occ)));
            });

            return make_pair(node_coverage, unique_paths.size());
        };

        if (!graph_positions.empty()) {
            std::cout << "#graph_position\tcoverage\tcoverage_uniq" << std::endl;
#pragma omp parallel for schedule(dynamic, 1)
            for (auto &pos : graph_positions) {
                nid_t node_id = id(pos);

                auto coverage = get_graph_node_coverage(graph, node_id);

#pragma omp critical (cout)
                std::cout << node_id << "," << offset(pos) << "," << (is_rev(pos) ? "-" : "+") << "\t"
                          << coverage.first << "\t" << coverage.second << std::endl;
            }
        }

        if (!path_positions.empty()) {
            std::cout << "#path_position\tcoverage\tcoverage_uniq" << std::endl;
#pragma omp parallel for schedule(dynamic, 1)
            for (auto &path_pos : path_positions) {
                pos_t pos = get_graph_pos(graph, path_pos);

                nid_t node_id = id(pos);

                auto coverage = get_graph_node_coverage(graph, node_id);

#pragma omp critical (cout)
                std::cout << (graph.get_path_name(path_pos.path)) << "," << path_pos.offset << ","
                          << (path_pos.is_rev ? "-" : "+") << "\t"
                          << coverage.first << "\t" << coverage.second << std::endl;
            }
        }

        if (!path_ranges.empty()) {
            std::cout << "#path\tstart\tend\tmean_coverage" << std::endl;
#pragma omp parallel for schedule(dynamic, 1)
            for (auto &path_range : path_ranges) {
                uint64_t start = path_range.begin.offset;
                uint64_t end = path_range.end.offset;
                path_handle_t path_handle = path_range.begin.path;

                uint64_t coverage = 0;

                uint64_t walked = 0;
                auto path_end = graph.path_end(path_handle);
                for (step_handle_t cur_step = graph.path_begin(path_handle);
                     cur_step != path_end && walked <= end; cur_step = graph.get_next_step(cur_step)) {
                    handle_t cur_handle = graph.get_handle_of_step(cur_step);
                    walked += graph.get_length(cur_handle);
                    if (walked > start) {
                        bool seen = false;
                        graph.for_each_step_on_handle(cur_handle, [&](const step_handle_t &step) {
                            if (!seen) {
                                if (graph.get_path_handle_of_step(step) == path_handle) {
                                    seen = true;
                                    return;
                                }
                            }
                            ++coverage;
                        });
                    }
                }

#pragma omp critical (cout)
                std::cout << (graph.get_path_name(path_handle)) << "\t" << start << "\t" << end << "\t"
                          << (double) coverage / (double) (end - start) << std::endl;
            }
        }

        return 0;
    }

    static Subcommand odgi_depth("depth", "find the depth of graph as defined by query criteria",
                                 PIPELINE, 3, main_depth);

}
