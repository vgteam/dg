#include "subcommand.hpp"
#include "odgi.hpp"
#include "position.hpp"
#include "args.hxx"
#include "split.hpp"
#include "subgraph/region.hpp"
#include "algorithms/bfs.hpp"
#include "algorithms/path_jaccard.hpp"
#include <omp.h>
#include "utils.hpp"
#include "picosha2.h"

namespace odgi {

using namespace odgi::subcommand;

int main_position(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi position";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("Find, translate, and liftover graph and path positions between graphs. Results are printed to stdout.");
    args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
    args::ValueFlag<std::string> og_target_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "target"});
    args::Group position_opts(parser, "[ Position Options ]");
    args::ValueFlag<std::string> og_source_file(position_opts, "FILE", "Translate positions from this *FILE graph into the target graph using common"
                                                                " *-l, --lift-paths* shared between both graphs (default: use the same"
                                                                " source/target graph). It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'x', "source"});
    args::ValueFlag<std::string> ref_path_name(position_opts, "PATH_NAME", "Translate the given positions into positions relative to this reference path.", {'r', "ref-path"});
    args::ValueFlag<std::string> ref_path_file(position_opts, "FILE", "Use the ref-paths in *FILE* for positional translation.", {'R', "ref-paths"});
    args::ValueFlag<std::string> lift_path_name(position_opts, "PATH_NAME", "Lift positions from *-x, --source* to *-i, --target* via coordinates in"
                                                                            " this path common to both graphs (default: all common paths between"
                                                                            " *-x, --source* and *-i, --target*).", {'l', "lift-path"});
    args::ValueFlag<std::string> lift_path_file(position_opts, "FILE", "Same as in *-l, --lift-paths*, but for all paths in *FILE*.", {'L', "lift-paths"});
    args::ValueFlag<std::string> graph_pos(position_opts, "[node_id][,offset[,(+|-)]*]*", "A graph position, e.g. 42,10,+ or 302,0,-.", {'g', "graph-pos"});
    args::ValueFlag<std::string> graph_pos_file(position_opts, "FILE", "Same as in *-g, --graph-pos*, but for all graph positions in *FILE*.", {'G', "graph-pos-file"});
    args::ValueFlag<std::string> path_pos(position_opts, "[path_name][,offset[,(+|-)]*]*", "A path position, e.g. chr8,1337,+ or chrZ,3929,-.", {'p', "path-pos"});
    args::ValueFlag<std::string> path_pos_file(position_opts, "FILE", "A *FILE* with one path position per line.", {'F', "path-pos-file"});
    args::ValueFlag<std::string> bed_input(position_opts, "FILE", "A BED file of ranges in paths in the graph to lift into the target"
                                                                  " graph *-v, --give-graph-pos* emit graph positions.", {'b', "bed-input"});
	args::ValueFlag<std::string> gff_input(position_opts, "FILE", "A GFF/GTF file with annotation of ranges in paths in the graph to lift into the target"
																  " (sub)graph emitting graph identifiers with annotation. The output is a CSV reading for the visualization within Bandage."
																  " The first column is the node identifier, the second column the annotation. "
																  "If several annotations exist for the same node, they are combined via ';'.", {'E', "gff-input"});
    args::Flag give_graph_pos(position_opts, "give-graph-pos", "Emit graph positions (node, offset, strand) rather than path positions.", {'v', "give-graph-pos"});
    args::Flag all_immediate(position_opts, "all-immediate", "Emit all positions immediately at the given graph/path position.", {'I', "all-immediate"});
    args::ValueFlag<uint64_t> _search_radius(position_opts, "DISTANCE", "Limit coordinate conversion breadth-first search up to DISTANCE bp from each given position (default: 10000).", {'d',"search-radius"});
	args::ValueFlag<uint64_t> _walking_dist(position_opts, "N", "Maximum walking distance in nucleotides for one orientation when finding the best target (reference) range for each query path (default: 10000). Note: If we walked 9999 base pairs and **w, --jaccard-context** is **10000**, we will also include the next node, even if we overflow the actual limit.",
											{'w', "jaccard-context"});
    args::Flag all_positions_of_ref_path(position_opts, "all-positions", "Emit all positions for all nodes in the specified ref-paths.", {"all-positions"});
    args::Group threading_opts(parser, "[ Threading ]");
    args::ValueFlag<uint64_t> threads(threading_opts, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
	args::Group processing_info_opts(parser, "[ Processing Information ]");
	args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
    args::Group program_info_opts(parser, "[ Program Information ");
    args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi position.", {'h', "help"});


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
    if (argc==1) {
        std::cout << parser;
        return 1;
    }

    if (!og_target_file) {
        std::cerr << "[odgi::position] error: please specify a target graph via -i=[FILE], --idx=[FILE]." << std::endl;
        return 1;
    }

	const uint64_t num_threads = args::get(threads) ? args::get(threads) : 1;
	omp_set_num_threads(num_threads);

    odgi::graph_t target_graph;
    assert(argc > 0);
    {
        const std::string infile = args::get(og_target_file);
        if (!infile.empty()) {
            if (infile == "-") {
                target_graph.deserialize(std::cin);
            } else {
                utils::handle_gfa_odgi_input(infile, "position", args::get(progress), num_threads, target_graph);
            }
        }
    }

    bool lifting = false;
    odgi::graph_t source_graph; // will be empty if we don't have different source and target
    if (og_source_file) {
        lifting = true;
        std::string infile = args::get(og_source_file);
        if (infile.size()) {
            if (infile == "-") {
                source_graph.deserialize(std::cin);
            } else {
				utils::handle_gfa_odgi_input(infile, "position", args::get(progress), num_threads, source_graph);
            }
        }
    }

    // todo: load many positions from a file
    // todo: convert a BED file
    // to simplify parallelism, collect our positions when doing so

    // collect our reference paths
    std::vector<path_handle_t> ref_paths;
    if (ref_path_name) {
        std::string path_name = args::get(ref_path_name);
        if (!target_graph.has_path(path_name)) {
            std::cerr << "[odgi::position] error: ref path " << path_name << " not found in graph" << std::endl;
            return 1;
        } else {
            ref_paths.push_back(target_graph.get_path_handle(path_name));
        }
    } else if (ref_path_file) {
        // for thing in things
        std::ifstream refs(args::get(ref_path_file).c_str());
        std::string path_name;
        while (std::getline(refs, path_name)) {
            if (!target_graph.has_path(path_name)) {
                std::cerr << "[odgi::position] error: ref path " << path_name << " not found in graph" << std::endl;
                return 1;
            } else {
                ref_paths.push_back(target_graph.get_path_handle(path_name));
            }
        }
    } else {
        // using all the paths in the graph
        target_graph.for_each_path_handle([&](const path_handle_t& path) { ref_paths.push_back(path); });
    }

    if (ref_paths.size() > 0 && args::get(all_positions_of_ref_path)){
        std::cout << "path\tnode_id\tposition" << std::endl;
        for (auto &path_handle : ref_paths) {
            uint64_t walked = 0;
            const auto path_end = target_graph.path_end(path_handle);
            auto path_name = target_graph.get_path_name(path_handle);
            for (step_handle_t cur_step = target_graph.path_begin(path_handle);
                 cur_step != path_end;
                 cur_step = target_graph.get_next_step(cur_step)) {
                const handle_t cur_handle = target_graph.get_handle_of_step(cur_step);
                std::cout << path_name << "\t" 
                          << target_graph.get_id(cur_handle) << "\t" 
                          << walked << std::endl;
                walked += target_graph.get_length(cur_handle);
            }
        }
    }

	std::unordered_map<std::string, std::tuple<std::string, uint64_t, uint64_t>> path_start_end_pos_map;
	std::string gff_in_file;
	if (gff_input) {
		gff_in_file = args::get(gff_input);
		if (!std::filesystem::exists(gff_in_file)) {
			std::cerr << "[odgi::position] error: the given file \"" << gff_in_file << "\" does not exist. "
																					   "Please specify an existing GFF/GTF file -E=[FILE], --gff-input=[FILE]." << std::endl;
			exit(1);
		}
		target_graph.for_each_path_handle([&](const path_handle_t& path) {
			uint64_t start_pos = 0;
			uint64_t end_pos = 0;
			std::string path_name = target_graph.get_path_name(path);
			std::string short_path_name;
			uint64_t len = 0;
			auto vals = split(path_name, ':');
			// TODO what if we have a subgraph of a subgraph?
			// do we have a subgraph in front of us?
			if (vals.size() > 1) {
				auto pos = split(vals[1], '-');
				start_pos = std::stoi(pos[0]);
				short_path_name = vals[0];
				end_pos = std::stoi(pos[1]);
			} else {
				target_graph.for_each_step_in_path(path, [&](const step_handle_t& step) {
					len += target_graph.get_length(target_graph.get_handle_of_step(step));
				});
				end_pos = len -1; // 0-based
				short_path_name = path_name;
			}
			path_start_end_pos_map[short_path_name] = std::make_tuple(path_name, start_pos, end_pos);
		});
	}

    // for translating source to target
    std::vector<path_handle_t> lift_paths_source;
    std::vector<path_handle_t> lift_paths_target;
    if ((lift_path_name || lift_path_file) && !lifting) {
        std::cerr << "[odgi::position] error: lifting requires a separate source and target graph, specify --source" << std::endl;
        return 1;
    } else if (lifting) {
        if (lift_path_name) {
            std::string path_name = args::get(lift_path_name);
            if (!target_graph.has_path(path_name)
                || !source_graph.has_path(path_name)) {
                std::cerr << "[odgi::position] error: lift path " << path_name << " not found in both source and target graph" << std::endl;
                return 1;
            } else {
                lift_paths_source.push_back(source_graph.get_path_handle(path_name));
                lift_paths_target.push_back(target_graph.get_path_handle(path_name));
            }
        } else if (lift_path_file) {
            // for thing in things
            std::ifstream refs(args::get(lift_path_file).c_str());
            std::string path_name;
            while (std::getline(refs, path_name)) {
                if (!target_graph.has_path(path_name)
                    || !source_graph.has_path(path_name)) {
                    std::cerr << "[odgi::position] error: lift path " << path_name << " not found in both source and target graph" << std::endl;
                    return 1;
                } else {
                    lift_paths_source.push_back(source_graph.get_path_handle(path_name));
                    lift_paths_target.push_back(target_graph.get_path_handle(path_name));
                }
            }
        } else {
            // using the common set of paths between the two graphs
            // make a set intersection
            std::vector<std::string> lift_names_source, lift_names_target, lift_names_common;
            source_graph.for_each_path_handle([&](const path_handle_t& path) {
                                                  lift_names_source.push_back(source_graph.get_path_name(path)); });
            target_graph.for_each_path_handle([&](const path_handle_t& path) {
                                                  lift_names_target.push_back(target_graph.get_path_name(path)); });
            std::sort(lift_names_source.begin(), lift_names_source.end());
            std::sort(lift_names_target.begin(), lift_names_target.end());
            std::set_intersection(lift_names_source.begin(), lift_names_source.end(),
                                  lift_names_target.begin(), lift_names_target.end(),
                                  std::back_inserter(lift_names_common));
            for (auto& path_name : lift_names_common) {
                lift_paths_source.push_back(source_graph.get_path_handle(path_name));
                lift_paths_target.push_back(target_graph.get_path_handle(path_name));
            }
            //target_graph.
            //lift_paths
        }
        if (lift_paths_source.size() != lift_paths_target.size()) {
            std::cerr << "[odgi::position] error: differing number of lift paths in target and source, suggests error" << std::endl;
            return 1;
        }
        if (lift_paths_source.empty() || lift_paths_target.empty()) {
            std::cerr << "[odgi::position] error: no lift paths common to both target and source, cannot proceed" << std::endl;
            std::cerr << "[odgi::position] select a set of common paths as lift paths or ensure that there are paths in common" << std::endl;
            return 1;
        }
    }

    // these options are exclusive (probably we should say with a warning)
    std::vector<odgi::pos_t> graph_positions;
    std::vector<odgi::path_pos_t> path_positions;
    std::vector<odgi::path_range_t> path_ranges;

    // TODO the parsers here should find the last 2 , delimiters and split on them
    
    auto add_graph_pos =
        [&graph_positions](const odgi::graph_t& graph,
                           const std::string& buffer) {
            auto vals = split(buffer, ',');
            /*
            if (vals.size() != 3) {
                std::cerr << "[odgi::position] error: graph position record is incomplete" << std::endl;
                std::cerr << "[odgi::position] error: got '" << buffer << "'" << std::endl;
                exit(1); // bail
            }
            */
            uint64_t id = std::stoi(vals[0]);
            if (!graph.has_node(id)) {
                std::cerr << "[odgi::position] error: no node " << id << " in graph" << std::endl;
                exit(1);
            }
            uint64_t offset = 0;
            if (vals.size() >= 2) {
                offset = std::stoi(vals[1]);
                handle_t h = graph.get_handle(id);
                // offsets are 0-based!
                if ((graph.get_length(h) - 1) < offset) {
                    std::cerr << "[odgi::position] error: offset of " << offset << " lies beyond the end of node " << id << std::endl;
                    exit(1);
                }
            }
            bool is_rev = false;
            if (vals.size() == 3) {
                is_rev = vals[2] == "-";
            }
            graph_positions.push_back(make_pos_t(id, is_rev, offset));
        };

    auto add_path_pos =
        [&path_positions](const odgi::graph_t& graph,
                          const std::string& buffer) {
            auto vals = split(buffer, ',');
            /*
            if (vals.size() != 3) {
                std::cerr << "[odgi::position] error: path position record is incomplete" << std::endl;
                std::cerr << "[odgi::position] error: got '" << buffer << "'" << std::endl;
                exit(1); // bail
            }
            */
            auto& path_name = vals[0];
            if (!graph.has_path(path_name)) {
                std::cerr << "[odgi::position] error: ref path " << path_name << " not found in graph" << std::endl;
                exit(1);
            } else {
                path_positions.push_back({
                        graph.get_path_handle(path_name),
                        (vals.size() > 1 ? (uint64_t)std::stoi(vals[1]) : 0),
                        (vals.size() == 3 ? vals[2] == "-" : false)
                    });
            }
        };

	auto add_gff_range =
			[&path_ranges](const odgi::graph_t& graph,
						   const std::string& buffer,
						   std::unordered_map<std::string, std::tuple<std::string, uint64_t, uint64_t>>& path_start_end_pos_map) {
				if (!buffer.empty() && buffer[0] != '#') {
					auto vals = split(buffer, '\t');
					/*
					if (vals.size() != 3) {
						std::cerr << "[odgi::position] error: path position record is incomplete" << std::endl;
						std::cerr << "[odgi::position] error: got '" << buffer << "'" << std::endl;
						exit(1); // bail
					}
					*/
					auto& path_name = vals[0];
					if (path_start_end_pos_map.count(path_name) == 0) {
						std::cerr << "[odgi::position] error: GFF/GTF path " << path_name << " not found in path_start_end_pos_map!" << std::endl;
						exit(1);
					} else {
						uint64_t start = vals.size() > 2 ? (uint64_t) std::stoi(vals[3]) : 0;
						uint64_t end = 0;
						if (vals.size() > 3) {
							end = (uint64_t) std::stoi(vals[4]);
						} else {
							graph.for_each_step_in_path(graph.get_path_handle(path_name), [&](const step_handle_t &s) {
								end += graph.get_length(graph.get_handle_of_step(s));
							});
						}

						if (start > end) {
							std::cerr << "[odgi::position::add_gff_range] error: wrong input coordinates in row: " << buffer << std::endl;
							exit(1);
						}


						// in the GFF format, the start and end are 1-based.
						// is the GFF entry actually within range?
						// we have to adjust for the given subgraph range which:
						// pos_in_gff - 1 - start_pos_(sub)graph + 1
						uint64_t graph_start_pos = get_long_path_start(path_start_end_pos_map[path_name]);
						uint64_t graph_end_pos = get_long_path_end(path_start_end_pos_map[path_name]);
						// the GFF/GTF can hit certain regions of the graph
						if (start >= graph_end_pos || end <= graph_start_pos) {
							// we do nothing here but return, because the GFF annotation is completely out of scope
							return;
						} else if (start <= graph_start_pos && end <= graph_end_pos) {
							start = 0; // graph_start_pos - graph_start_pos;
							end = end - graph_start_pos - 1;
						} else if (start >= graph_start_pos && end >= graph_end_pos) {
							end = graph_end_pos - graph_start_pos - 1;
							start = start - graph_start_pos - 1;
						} else if (start >= graph_start_pos && end <= graph_end_pos) {
							start = start - graph_start_pos - 1;
							end = end - graph_start_pos - 1;
						} else if (start <= graph_start_pos && end >= graph_end_pos) {
							start = 0;
							end = graph_end_pos - graph_start_pos;
						} else {
							return;
						}

						if (start > end) {
							std::cerr << "[odgi::position::add_gff_range] error: wrong input coordinates in row: " << buffer <<
							"for detected start: " << start << " and end: " << end << std::endl;
							exit(1);
						}

						// we can store things in the path ranges, but then we need to traverse these differently
						path_ranges.push_back(
								{
										{
												graph.get_path_handle(get_long_path_name(path_start_end_pos_map[path_name])),
												start,
												false
										},
										{
												graph.get_path_handle(get_long_path_name(path_start_end_pos_map[path_name])),
												end,
												false
										},
										(vals.size() > 6 && vals[6] == "-"),
										vals[8],
										vals[8]
								});
					}
				}
			};

	if (!gff_input) {
		if (graph_pos) {
			// if we're given a graph_pos, we'll convert it into a path pos
			if (lifting) {
				add_graph_pos(source_graph, args::get(graph_pos));
			} else {
				add_graph_pos(target_graph, args::get(graph_pos));
			}
		} else if (graph_pos_file) {
			std::ifstream gpos(args::get(graph_pos_file).c_str());
			std::string buffer;
			if (lifting) {
				while (std::getline(gpos, buffer)) {
					add_graph_pos(source_graph, buffer);
				}
			} else {
				while (std::getline(gpos, buffer)) {
					add_graph_pos(target_graph, buffer);
				}
			}
		} else if (path_pos) {
			// if given a path pos, we convert it into a path pos in our reference set
			if (lifting) {
				add_path_pos(source_graph, args::get(path_pos));
			} else {
				add_path_pos(target_graph, args::get(path_pos));
			}
		} else if (path_pos_file) {
			// if we're given a file of path positions, we'll convert them all
			std::ifstream refs(args::get(path_pos_file).c_str());
			std::string buffer;
			while (std::getline(refs, buffer)) {
				if (lifting) {
					add_path_pos(source_graph, buffer);
				} else {
					add_path_pos(target_graph, buffer);
				}
			}
		} else if (bed_input) {
			std::ifstream bed_in(args::get(bed_input).c_str());
			std::string buffer;
			while (std::getline(bed_in, buffer)) {
				if (lifting) {
				    add_bed_range(path_ranges, source_graph, buffer);
				} else {
				    add_bed_range(path_ranges, target_graph, buffer);
				}
			}
		}
	} else {
		std::ifstream gff_in(gff_in_file.c_str());
		std::string buffer;
		while (std::getline(gff_in, buffer)) {
			// we have to add to our own vector, because we are doing something different
			add_gff_range(target_graph, buffer, path_start_end_pos_map);
		}
	}

    uint64_t search_radius = _search_radius ? args::get(_search_radius) : 10000;
    uint64_t walking_dist = _walking_dist ? args::get(_walking_dist) : 10000;

    // make an hash set of our ref path ids for quicker lookup
    hash_set<uint64_t> ref_path_set;
    for (auto& path : ref_paths) {
        ref_path_set.insert(as_integer(path));
    }
    hash_set<uint64_t> lift_path_set_source;
    hash_set<uint64_t> lift_path_set_target;
    for (auto& path : lift_paths_source) {
        lift_path_set_source.insert(as_integer(path));
    }
    for (auto& path : lift_paths_target) {
        lift_path_set_target.insert(as_integer(path));
    }

    auto get_graph_pos =
        [](const odgi::graph_t& graph,
           const path_pos_t& pos,
           step_handle_t& step) {
            auto path_end = graph.path_end(pos.path);
            uint64_t walked = 0;
            for (step_handle_t s = graph.path_begin(pos.path);
                 s != path_end; s = graph.get_next_step(s)) {
                handle_t h = graph.get_handle_of_step(s);
                uint64_t node_length = graph.get_length(h);
                if (walked + node_length - 1 >= pos.offset) {
                	step = s;
                    return make_pos_t(graph.get_id(h), graph.get_is_reverse(h), pos.offset - walked);
                }
                walked += node_length;
            }
#pragma omp critical (cout)
            std::cerr << "[odgi::position] warning: position " << graph.get_path_name(pos.path) << ":" << pos.offset << " outside of path. Walked " << walked << std::endl;
            return make_pos_t(0, false, 0);
        };

	auto get_graph_node_ids_annotation =
			[](const odgi::graph_t& graph,
			   const path_range_t& path_range) {
				auto path_end = graph.path_end(path_range.begin.path);
				std::unordered_map<uint64_t , std::set<std::string>> node_annotation_map;
				uint64_t walked = 0;
				uint64_t path_pos_start = path_range.begin.offset;
				uint64_t path_pos_end = path_range.end.offset;
				for (step_handle_t s = graph.path_begin(path_range.begin.path);
					 s != path_end; s = graph.get_next_step(s)) {
					handle_t h = graph.get_handle_of_step(s);
					uint64_t nid = graph.get_id(h);
					uint64_t node_length = graph.get_length(h);
					uint64_t local_min_pos = walked;
					uint64_t local_max_pos = local_min_pos + node_length - 1;
					if ((path_pos_start >= local_min_pos && path_pos_start <= local_max_pos)
						|| (path_pos_end <= local_max_pos && path_pos_end >= local_min_pos)
						|| (path_pos_start <= local_min_pos && path_pos_end >= local_max_pos)) {
						// TODO add to our hashmap of node_id -> hash_set of annotation
						if (node_annotation_map.count(nid) == 0) {
							std::set<std::string> anno_set;
							anno_set.insert(path_range.data);
							node_annotation_map[nid] = anno_set;
							// we just blindly at the value again
						} else {
							std::set<std::string> anno_set = node_annotation_map[nid];
							anno_set.insert(path_range.data);
						}
					} else if (local_min_pos > path_pos_end) {
						// we can return here, nothing more to do
						return node_annotation_map;
					}
					walked += node_length;
				}
				return node_annotation_map;
			};

    auto get_offset_in_path =
        [](const odgi::graph_t& graph,
           const path_handle_t& path, const step_handle_t& target) {
            auto path_end = graph.path_end(path);
            uint64_t walked = 0;
            step_handle_t s = graph.path_begin(path);
            for ( ;  s != target; s = graph.get_next_step(s)) {
                handle_t h = graph.get_handle_of_step(s);
                walked += graph.get_length(h);
            }
            assert(s != path_end);
            return walked;
        };

	auto set_adj_last_node =
			[](const odgi::graph_t& graph,
			   const step_handle_t& ref_hit, const handle_t& h_bfs,
			   const bool& used_bidirectional, const uint64_t& d_bfs, const pos_t& pos,
			   bool& rev_vs_ref, uint64_t& adj_last_node) {
				// this check is confusing, but it's due to us walking the reverse graph from our start position in the BFS
				rev_vs_ref = graph.get_is_reverse(graph.get_handle_of_step(ref_hit)) == graph.get_is_reverse(h_bfs);
				if ((d_bfs == 0) || (d_bfs == graph.get_length(h_bfs) && used_bidirectional)) { // if we're on the start node
					if (rev_vs_ref) {
						// and if the path orientation is the same as our traversal orientation
						// then we need to add the remaining distance from our original offset to the end of node
						// to the final path position offset
						adj_last_node = graph.get_length(h_bfs) - offset(pos);
					} else {
						// otherwise if the original path is in the same orientation
						// then we add the original forward offset to the ref path offset
						adj_last_node = offset(pos);
					}
				} else { // if we're not on the first node
					if (rev_vs_ref) {
						// and we come onto the result in the same orientation
						// it means the ref pos is at the node end
						adj_last_node = 0; // so we have no adjustment
					} else {
						// otherwise, it means the original path is in the same orientation
						// then we need to adjust by the length of this stepb
						// because we enter at node end, but we'll get the graph position for the step
						// at the node beginning
						adj_last_node = graph.get_length(h_bfs);
					}
				}
			};

    // TODO this part needs to include adjustments for in-node offsets vs. where we find the ref path
    // TODO should we always look "backwards" when seeking the ref pos?

    struct lift_result_t {
        int64_t path_offset = 0;
        step_handle_t ref_hit;
        uint64_t walked_to_hit_ref = 0;
        bool is_rev_vs_ref = false;
        bool used_bidirectional = false;
    };

    // get the reference path positions right where we are
    auto get_immediate =
        [&search_radius,&get_offset_in_path](const odgi::graph_t& graph,
                                             const hash_set<uint64_t>& path_set,
                                             const pos_t& pos, std::vector<lift_result_t>& lifts) {
            // unpacking our args
            handle_t h = graph.get_handle(id(pos), is_rev(pos));
            graph.for_each_step_on_handle(
                h, [&](const step_handle_t& s) {
                       auto p = graph.get_path_handle_of_step(s);
                       if (path_set.count(as_integer(p))) {
                           // are we reverse against the reference
                           bool rev_vs_ref = graph.get_is_reverse(graph.get_handle_of_step(s)) != graph.get_is_reverse(h);
                           int adj_node = 0;
                           if (rev_vs_ref) {
                               // and if the path orientation is the same as our traversal orientation
                               // then we need to add the remaining distance from our original offset to the end of node
                               // to the final path position offset
                               adj_node = graph.get_length(h) - offset(pos);
                           } else {
                               // otherwise if the original path is in the same orientation
                               // then we add the original forward offset to the ref path offset
                               adj_node = offset(pos);
                           }
                           lifts.emplace_back();
                           auto& lift = lifts.back();
                           lift.ref_hit = s;
                           lift.path_offset = get_offset_in_path(graph, p, lift.ref_hit) + adj_node;
                           lift.walked_to_hit_ref = 0;
                           lift.is_rev_vs_ref = rev_vs_ref;
                       }
                   });
            return lifts.size() > 0;
        };

    auto get_position =
        [&search_radius,&get_offset_in_path,&walking_dist,&set_adj_last_node](const odgi::graph_t& graph,
                                             const hash_set<uint64_t>& path_set,
                                             const pos_t& pos, lift_result_t& lift,
                                             const step_handle_t target_step_handle,
                                             const bool path_jaccard,
                                             const path_handle_t* target_path=NULL) {
            // unpacking our args
            int64_t& path_offset = lift.path_offset;
            step_handle_t& ref_hit = lift.ref_hit;
            uint64_t& walked_to_hit_ref = lift.walked_to_hit_ref;
            bool& rev_vs_ref = lift.is_rev_vs_ref;
            bool& used_bidirectional = lift.used_bidirectional;
            handle_t start_handle = graph.get_handle(id(pos), is_rev(pos));
            bool found_hit = false;
            uint64_t adj_last_node = 0;
            hash_set<uint64_t> seen;
            uint64_t d_bfs;
            handle_t h_bfs;
            for (auto try_bidirectional : { false, true }) {
                if (try_bidirectional) {
					used_bidirectional = true;
					seen.erase(as_integer(graph.flip(start_handle)));
                }
                odgi::algorithms::bfs(
                    graph,
                    [&](const handle_t& h, const uint64_t& r, const uint64_t& l, const uint64_t& d) {
                        seen.insert(as_integer(h));
                        bool got_hit = false;
                        step_handle_t hit;
                        graph.for_each_step_on_handle(
                            h, [&](const step_handle_t& s) {
                                   auto p = graph.get_path_handle_of_step(s);
                                   if (!got_hit && path_set.count(as_integer(p)) && (target_path==NULL||p==*target_path)) {
                                       //std::cerr << "thought I got a hit" << std::endl;
                                       got_hit = true;
                                       hit = s;
                                       walked_to_hit_ref += l; // how far we came to get to this node
									   d_bfs = d; // we need this for the path jaccard calculations
									   h_bfs = h;
									   set_adj_last_node(graph, s, h, used_bidirectional, d, pos, rev_vs_ref, adj_last_node);
                                   }
                               });
                        if (got_hit) {
                            ref_hit = hit;
                            found_hit = true;
                        }
                    },
                    [&seen](const handle_t& h) { return seen.count(as_integer(h)); },
                    [](const handle_t& l, const handle_t& h) { return false; },
                    [&found_hit](void) { return found_hit; },
                    { graph.flip(start_handle) },
                    { },
                    used_bidirectional,
                    0,
                    search_radius);
                if (found_hit) break; // if we got a hit, don't go bidirectional
            }
            if (found_hit) {
            	if (path_jaccard) {
					std::vector<step_handle_t> query_step_handles;
					path_handle_t ref_hit_path = graph.get_path_handle_of_step(ref_hit);
					graph.for_each_step_on_handle(
							h_bfs,
							[&](const step_handle_t& s) {
								/// we can do these expensive iterations here, because we only have to do it once for each walk
								if (graph.get_path_handle_of_step(s) == ref_hit_path) {
									// collect only the steps for the given target
									query_step_handles.push_back(s);
								}
							});
					// iterate over the node to get the list of canditate query step handles
					std::vector<algorithms::step_jaccard_t> target_jaccard_indices = algorithms::jaccard_indices_from_step_handles(graph,
																																   walking_dist,
																																   target_step_handle,
																																   query_step_handles);
					ref_hit = target_jaccard_indices[0].step;
					set_adj_last_node(graph, ref_hit, h_bfs, used_bidirectional, d_bfs, pos, rev_vs_ref, adj_last_node);
				}

				path_handle_t p = graph.get_path_handle_of_step(ref_hit);
				// TODO ORIENTATION
				path_offset = get_offset_in_path(graph, p, ref_hit) + adj_last_node;
                return true;
            } else {
                path_offset = -1;
                return false;
            }
        };

    if (graph_positions.size()) {
        if (lifting) {
            std::cout << "#source.graph.pos\ttarget.graph.pos\t";
        } else {
            std::cout << "#target.graph.pos\t";
        }
        if (give_graph_pos) {
            std::cout << "target.graph.pos" << std::endl;
        } else if (all_immediate) {
            std::cout << "target.path.pos\tdist.to.ref\tstrand.vs.ref" << std::endl;
        } else {
        	if (ref_path_name || ref_path_file) {
				std::cout << "target.path.pos\tdist.to.ref\tstrand.vs.ref" << std::endl;
        	} else {
				std::cout << "target.path.pos\tdist.to.path\tstrand.vs.ref" << std::endl;
        	}
        }
    }
    // for each position that we want to look up
#pragma omp parallel for schedule(dynamic,1)
    for (auto& _pos : graph_positions) {
        // go to the graph
        // do a little BFS, bounded by our limit
        // now, if we found our hit, print
        // optionally, we will translate from a source graph into a target graph
        lift_result_t source_result;
        //bool ok = false;
        pos_t pos;
        // empty step handle
		step_handle_t step_handle_graph_pos;
        if (lifting) {
            if (get_position(source_graph, lift_path_set_source, _pos, source_result, step_handle_graph_pos, false)) {
                pos = get_graph_pos(target_graph,
                                    { target_graph.get_path_handle(
                                            source_graph.get_path_name(
                                                source_graph.get_path_handle_of_step(
                                                    source_result.ref_hit))),
                                      (uint64_t)source_result.path_offset,
                                      source_result.is_rev_vs_ref },
                                      step_handle_graph_pos);
            } else {
                pos = make_pos_t(0,false,0); // couldn't lift
            }
        } else {
            pos = _pos;
        }
        //= (lifting ? ok = get_position(source_graph, lay_path_set
        lift_result_t result;
        std::vector<lift_result_t> result_v;
        if (id(pos) && give_graph_pos) {
            // force graph position in target
#pragma omp critical (cout)
            {
                if (lifting) {
                    std::cout << id(_pos) << "," << offset(_pos) << "," << (is_rev(_pos) ? "-" : "+") << "\t";
                }
                std::cout << id(pos) << "," << offset(pos) << "," << (is_rev(pos) ? "-" : "+") << "\t"
                          << "\t" << id(pos) << "," << offset(pos) << "," << (is_rev(pos) ? "-" : "+") << std::endl;
            }
        } else if (args::get(all_immediate) && get_immediate(target_graph, ref_path_set, pos, result_v)) {
            bool ref_is_rev = false;
            for (auto& result : result_v) {
                path_handle_t p = target_graph.get_path_handle_of_step(result.ref_hit);
#pragma omp critical (cout)
                {
                    if (lifting) {
                        std::cout << id(_pos) << "," << offset(_pos) << "," << (is_rev(_pos) ? "-" : "+") << "\t";
                    }
                    std::cout << id(pos) << "," << offset(pos) << "," << (is_rev(pos) ? "-" : "+") << "\t"
                              << target_graph.get_path_name(p) << "," << result.path_offset << "," << (ref_is_rev ? "-" : "+") << "\t"
                              << result.walked_to_hit_ref << "\t" << (result.is_rev_vs_ref ? "-" : "+") << std::endl;
                }
            }
        } else if (get_position(target_graph, ref_path_set, pos, result, step_handle_graph_pos, false)) {
            bool ref_is_rev = false;
            path_handle_t p = target_graph.get_path_handle_of_step(result.ref_hit);
#pragma omp critical (cout)
            {
                if (lifting) {
                    std::cout << id(_pos) << "," << offset(_pos) << "," << (is_rev(_pos) ? "-" : "+") << "\t";
                }
                std::cout << id(pos) << "," << offset(pos) << "," << (is_rev(pos) ? "-" : "+") << "\t"
                          << target_graph.get_path_name(p) << "," << result.path_offset << "," << (ref_is_rev ? "-" : "+") << "\t"
                          << result.walked_to_hit_ref << "\t" << (result.is_rev_vs_ref ? "-" : "+") << std::endl;
            }
        }
    }

#pragma omp parallel for schedule(dynamic,1)
    for (auto& path_pos : path_positions) {
        // TODO we need a better input format
        pos_t pos;
		step_handle_t step_handle_graph_pos;
        // handle the lift into the target graph
        if (lifting) {
            lift_result_t source_result;
            pos_t _pos = get_graph_pos(source_graph, path_pos, step_handle_graph_pos);
            if (id(_pos) && get_position(source_graph, lift_path_set_source, _pos, source_result, step_handle_graph_pos, true)) {
                pos = get_graph_pos(target_graph,
                                    { target_graph.get_path_handle(
                                            source_graph.get_path_name(
                                                source_graph.get_path_handle_of_step(
                                                    source_result.ref_hit))),
                                      (uint64_t)source_result.path_offset,
                                      source_result.is_rev_vs_ref },
                                      step_handle_graph_pos);
            } else {
                pos = make_pos_t(0,false,0); // couldn't lift
            }

        } else {
            //path_pos = _path_pos;
            pos = get_graph_pos(target_graph, path_pos, step_handle_graph_pos);
        }
        lift_result_t result;
        //std::cerr << "Got graph pos " << id(pos) << std::endl;
        if (id(pos)) {
            if (give_graph_pos) {
#pragma omp critical (cout)
                std::cout << "#source.path.pos\ttarget.graph.pos" << std::endl
                          << (lifting ? source_graph.get_path_name(path_pos.path) : target_graph.get_path_name(path_pos.path))
                          << "," << path_pos.offset << "," << (path_pos.is_rev ? "-" : "+")
                          << "\t" << id(pos) << "," << offset(pos) << "," << (is_rev(pos) ? "-" : "+") << std::endl;
            } else if (get_position(target_graph, ref_path_set, pos, result, step_handle_graph_pos, true)) {
                bool ref_is_rev = false;
                path_handle_t p = target_graph.get_path_handle_of_step(result.ref_hit);
#pragma omp critical (cout)
                std::cout << "#source.path.pos\ttarget.path.pos\tdist.to.ref\tstrand.vs.ref" << std::endl
                          << (lifting ? source_graph.get_path_name(path_pos.path) : target_graph.get_path_name(path_pos.path)) << ","
                          << path_pos.offset << "," << (path_pos.is_rev ? "-" : "+") << "\t"
                          << target_graph.get_path_name(p) << "," << result.path_offset << "," << (ref_is_rev ? "-" : "+") << "\t"
                          << result.walked_to_hit_ref << "\t" << (result.is_rev_vs_ref ? "-" : "+") << std::endl;
            }
        }
    }

	std::vector<std::unordered_map<uint64_t , std::set<std::string>>> node_annotation_maps;

#pragma omp parallel for schedule(dynamic,1)
    for (auto& path_range : path_ranges) {
		pos_t pos_begin, pos_end;
        // handle the lift into the target graph
		step_handle_t step_handle_graph_pos_begin;
		step_handle_t step_handle_graph_pos_end;
        if (lifting) {
            lift_result_t source_begin_result, source_end_result;
            pos_t _pos_begin = get_graph_pos(source_graph, path_range.begin, step_handle_graph_pos_begin);
            pos_t _pos_end = get_graph_pos(source_graph, path_range.end, step_handle_graph_pos_end);
            if (id(_pos_begin) && get_position(source_graph, lift_path_set_source, _pos_begin, source_begin_result, step_handle_graph_pos_begin, true)
                && id(_pos_end) && get_position(source_graph, lift_path_set_source, _pos_end, source_end_result, step_handle_graph_pos_end, true)) {
                pos_begin = get_graph_pos(target_graph,
                                          { target_graph.get_path_handle(
                                                  source_graph.get_path_name(
                                                      source_graph.get_path_handle_of_step(
                                                          source_begin_result.ref_hit))),
                                            (uint64_t)source_begin_result.path_offset,
                                            source_begin_result.is_rev_vs_ref },
										  step_handle_graph_pos_begin);
                pos_end = get_graph_pos(target_graph,
                                        { target_graph.get_path_handle(
                                                source_graph.get_path_name(
                                                    source_graph.get_path_handle_of_step(
                                                        source_end_result.ref_hit))),
                                          (uint64_t)source_end_result.path_offset,
                                          source_end_result.is_rev_vs_ref },
										step_handle_graph_pos_end);
            } else {
                pos_begin = make_pos_t(0,false,0); // couldn't lift
                pos_end = make_pos_t(0,false,0); // couldn't lift
            }
        } else {
            //path_pos = _path_pos;
			if (gff_input) {
				std::unordered_map<uint64_t , std::set<std::string>> node_annotation_map = get_graph_node_ids_annotation(target_graph, path_range);
#pragma omp critical (node_annotation_maps)
				node_annotation_maps.push_back(node_annotation_map);
			} else {
				pos_begin = get_graph_pos(target_graph, path_range.begin, step_handle_graph_pos_begin);
				pos_end = get_graph_pos(target_graph, path_range.end, step_handle_graph_pos_end);
			}
        }
        if (id(pos_begin) && id(pos_end) && !gff_input) {
            lift_result_t lift_begin;
            lift_result_t lift_end;
            // TODO add a GAF-style path to the record to say where the BED range walks in the graph
            // TODO optionally list out the nodes in this particular range (e.g. those within it in our sort order)
            if (give_graph_pos) {
#pragma omp critical (cout)
                std::cout << path_range.data << "\t"
                          << id(pos_begin) << "," << offset(pos_begin) << "," << (is_rev(pos_begin)?"-":"+") << "\t"
                          << id(pos_end) << "," << offset(pos_end) << "," << (is_rev(pos_end)?"-":"+") << std::endl;
            } else {

                for (const path_handle_t& P: ref_paths) {
                    if (get_position(target_graph, ref_path_set, pos_begin, lift_begin, step_handle_graph_pos_begin, true,&P)
                        && get_position(target_graph, ref_path_set, pos_end, lift_end, step_handle_graph_pos_end, true,&P)) {

                        bool ref_is_rev = false;
                        path_handle_t p_begin = target_graph.get_path_handle_of_step(lift_begin.ref_hit);
                        path_handle_t p_end = target_graph.get_path_handle_of_step(lift_end.ref_hit);
                        // XXX TODO assert these to be equal......
#pragma omp critical (cout)
                        std::cout << path_range.data << "\t"
                              << target_graph.get_path_name(p_begin) << ","
                              << lift_begin.path_offset << ","
                              << (lift_begin.is_rev_vs_ref ? "-" : "+") << "\t"
                              << target_graph.get_path_name(p_end) << ","
                              << lift_end.path_offset << ","
                              << (lift_end.is_rev_vs_ref ? "-" : "+") << "\t"
                              << (lift_begin.is_rev_vs_ref ^ path_range.is_rev ? "-" : "+") << std::endl;
                              //<< walked_to_hit_ref << "\t" << (is_rev_vs_ref ? "-" : "+") << std::endl;
                    }
                }
            }
        }
    }
	if (gff_input) {
		//  clean up duplicates
		std::map<uint64_t , std::set<std::string>> final_node_annotation_map;
		std::cout << "NODE_ID,ANNOTATION,COLOR" << std::endl;
		for (auto& node_annotation_map : node_annotation_maps) {
			for (auto& elem : node_annotation_map) {
				if (final_node_annotation_map.count(elem.first) == 0) {
					final_node_annotation_map[elem.first] = elem.second;
				} else {
					for (auto& anno : elem.second) {
						final_node_annotation_map[elem.first].insert(anno);
					}
				}
			}
		}
		std::string prev_anno = "";
		std::set<std::string> prev_set;
		uint64_t prev_node_id = -1;
		std::map<uint64_t , std::set<std::string>>::iterator it_map;
		for (it_map = final_node_annotation_map.begin(); it_map != final_node_annotation_map.end(); ++it_map) {
			std::cout << std::dec << it_map->first << ",";
			std::string anno = "";
			for(auto it = it_map->second.begin() ; it != it_map->second.end() ; ++it) {
				if(it != it_map->second.begin()) {
					anno += ";";
				}
				anno = anno + *it;
			}

			// do we match with the previous set?
			// do we match with the next set, if available
			if (prev_set != it_map->second
				|| (std::prev(final_node_annotation_map.end())->first == it_map->first)
				|| (std::next(it_map)->second != it_map->second && std::next(it_map) != final_node_annotation_map.end())) {

				cout << anno;
			}

			prev_set = it_map->second;

			// use a sha256 to get a few bytes that we'll use for a color
			picosha2::byte_t hashed[picosha2::k_digest_size];
			picosha2::hash256(anno.begin(), anno.end(), hashed, hashed + picosha2::k_digest_size);

			uint8_t path_r = hashed[24];
			uint8_t path_g = hashed[8];
			uint8_t path_b = hashed[16];

			cout << ",#" << std::setfill('0') << std::setw(6) << std::hex << (path_r << 16 | path_g << 8 | path_b );
			std::cout << std::endl;
		}
	}

    // todo - lift the position into another graph
    // requires an input of target paths in the final graph
    // and optionally the set of paths in common (we can compute this by default) to drive the lift

    return 0;
}

static Subcommand odgi_position("position", "Find, translate, and liftover graph and path positions between graphs.",
                                PIPELINE, 3, main_position);


}
