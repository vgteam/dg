#include "untangle.hpp"

namespace odgi {
namespace algorithms {

using namespace handlegraph;

std::vector<step_handle_t> untangle_cuts(
    const PathHandleGraph& graph,
    const path_handle_t& path,
    const std::string& path_name,
    const path_step_index_t::step_it& _start,
    const path_step_index_t::step_it& _end,
    const step_index_t& step_index,
    const path_step_index_t& self_index,
    const std::function<bool(const handle_t&)>& is_cut) {

    /*
    std::cerr << "untangle_cuts(" << path_name << ", "
              << step_index.get_position(_start) << ", "
              << step_index.get_position(_end) << ")" << std::endl;
    */

    ska::flat_hash_set<const step_handle_t*> seen_fwd_step;
    ska::flat_hash_set<const step_handle_t*> seen_rev_step;
    auto is_seen_fwd_step = [&seen_fwd_step,&self_index](const path_step_index_t::step_it& step) {
        return seen_fwd_step.count(&*step);
    };
    auto is_seen_rev_step = [&seen_rev_step,&self_index](const path_step_index_t::step_it& step) {
        return seen_rev_step.count(&*step);
    };
    auto mark_seen_fwd_step = [&seen_fwd_step,&self_index](const path_step_index_t::step_it& step) {
        bool state = seen_fwd_step.count(&*step);
        seen_fwd_step.insert(&*step);
        return state;
    };
    auto mark_seen_rev_step = [&seen_rev_step,&self_index](const path_step_index_t::step_it& step) {
        bool state = seen_rev_step.count(&*step);
        seen_rev_step.insert(&*step);
        return state;
    };

    std::vector<std::pair<uint64_t, step_handle_t>> cut_points;
    auto clean_cut_points = [&cut_points,&step_index,&graph](void) {
        std::sort(cut_points.begin(),
                  cut_points.end(),
                  [&](const std::pair<uint64_t,step_handle_t>& a,
                      const std::pair<uint64_t,step_handle_t>& b) {
                      return a.first < b.first;
                  });
        std::vector<step_handle_t> c;
        for (auto& x : cut_points) c.push_back(x.second);
        // then take unique positions
        c.erase(std::unique(c.begin(),
                            c.end()),
                c.end());
        return c;
    };

    std::deque<std::pair<path_step_index_t::step_it, path_step_index_t::step_it>> todo;
    todo.push_back(std::make_pair(_start, _end));
    while (!todo.empty()) {
        auto start = todo.front().first;
        auto end = todo.front().second;
        uint64_t start_pos = step_index.get_position(*start, graph);
        uint64_t curr_pos = start_pos;
        uint64_t end_pos = step_index.get_position(*end, graph);
        //std::cerr << "todo: " << start_pos << " " << end_pos << std::endl;
        cut_points.push_back(std::pair(curr_pos, *start));
        todo.pop_front();
        // we go forward until we see a loop, where the other step has position < end_pos and > start_pos
        for (auto step = start; step != end; ++step) {
            //  we take the first and shortest loop we find
            handle_t handle = graph.get_handle_of_step(*step);
            if (mark_seen_fwd_step(step)) {
                curr_pos += graph.get_length(handle);
                continue;
            }
            if (is_cut(handle)) {
                /*
                if (curr_pos != step_index.get_position(step, graph)) {
                    std::cerr << "position error fwd " << curr_pos << " " << step_index.get_position(step, graph) << std::endl;
                    exit(1);
                }
                */
                cut_points.push_back(std::pair(curr_pos, *step));
                //cut_points.push_back(std::make_pair(step_index.get_position(step, graph), step));
            }
            bool found_loop = false;
            path_step_index_t::step_it other;
            uint64_t other_pos = 0;
            auto x = self_index.get_next_step_on_node(graph.get_id(handle), *step);
            if (x.first) {
                other_pos = x.second.second;
                if (other_pos > start_pos
                    && other_pos < end_pos
                    && other_pos > curr_pos
                    && !is_seen_fwd_step(x.second.first)) {
                    other = x.second.first;
                    found_loop = true;
                }
            }
            if (found_loop) {
                //  recurse this function into it, taking start as our current handle other side of the loop as our end
                //  to cut_points we add the start position, the result from recursion, and our end position
                //std::cerr << "Found loop! " << step_index.get_position(step) << " " << step_index.get_position(other) << std::endl;
                todo.push_back(std::make_pair(step, other));
                //  we then step forward to the loop end and continue iterating
                curr_pos = other_pos + graph.get_length(graph.get_handle_of_step(*other));
                step = other;
            } else {
                curr_pos += graph.get_length(handle);
            }
        }
        // TODO this block is the same as the previous one, but in reverse
        // the differences in how positions are managed are subtle, making it hard to fold the
        // forward and reverse version together
        // now we reverse it
        step_handle_t path_begin = graph.path_begin(path);
        if (*end == path_begin || !graph.has_previous_step(*end)) {
            return clean_cut_points();
        }
        //step_handle_t _step = graph.get_previous_step(end);
        curr_pos = step_index.get_position(*end, graph) + graph.get_length(graph.get_handle_of_step(*end));
        //cut_points.push_back(std::make_pair(curr_pos, end));
        //std::cerr << "reversing" << std::endl;
        for (auto step = end; step != start; --step) {
            handle_t handle = graph.get_handle_of_step(*step);
            curr_pos -= graph.get_length(handle);
            if (mark_seen_rev_step(step)) {
                continue;
            }
            //  we take the first and shortest loop we find
            if (is_cut(handle)) {
                /*
                if (curr_pos != step_index.get_position(step, graph)) {
                    std::cerr << "position error rev " << curr_pos << " " << step_index.get_position(step, graph) << std::endl;
                    exit(1);
                }
                */
                cut_points.push_back(std::pair(curr_pos, *step));
            }
            //std::cerr << "rev on step " << graph.get_id(handle) << " " << curr_pos << std::endl;
            bool found_loop = false;
            path_step_index_t::step_it other;
            auto x = self_index.get_prev_step_on_node(graph.get_id(handle), *step);
            uint64_t other_pos = 0;
            if (x.first) {
                other_pos = x.second.second;
                if (other_pos > start_pos
                    && other_pos < end_pos
                    && other_pos < curr_pos
                    && !is_seen_rev_step(x.second.first)) {
                    other = x.second.first;
                    found_loop = true;
                }
            }
            if (found_loop) {
                //  recurse this function into it, taking start as our current handle other side of the loop as our end
                //  to cut_points we add the start position, the result from recursion, and our end position
                todo.push_back(std::make_pair(other, step));
                //  we then step forward to the loop end and continue iterating
                curr_pos = other_pos;
                step = other;
            }
        }
    }
    // and sort
    return clean_cut_points();
}

void write_cuts(
    const PathHandleGraph& graph,
    const path_handle_t& path,
    const std::vector<step_handle_t>& cuts,
    const step_index_t& step_index) {
    auto path_name = graph.get_path_name(path);
    std::cout << "name\tcut" << std::endl;
    for (auto& step : cuts) {
        std::cout << path_name << "\t" << step_index.get_position(step, graph) << std::endl;
    }
}

std::vector<step_handle_t> merge_cuts(
    const std::vector<step_handle_t>& cuts,
    const uint64_t& dist,
    const step_index_t& step_index,
	const PathHandleGraph& graph) {
    std::vector<step_handle_t> merged;
    uint64_t last = 0;
    //std::cerr << "dist is " << dist << std::endl;
    for (auto& step : cuts) {
        auto pos = step_index.get_position(step, graph);
        if (pos == 0 || pos > (last + dist)) {
            merged.push_back(step);
            last = pos;
        }
    }
    // add the end step
    if (cuts.size()) {
        merged.push_back(
            graph.path_end(graph.get_path_handle_of_step(cuts.front())));
    }
    return merged;
}

void self_dotplot(
    const PathHandleGraph& graph,
    const path_handle_t& path) {
    auto step_pos = make_step_index(graph, { path }, 1);
    auto path_name = graph.get_path_name(path);
    std::cout << "name\tfrom\tto" << std::endl;
    graph.for_each_step_in_path(
        path,
        [&](const step_handle_t& step) {
            // TODO this is given by the walk, no need for a hash table lookup
            const auto& curr_pos = step_pos.find(step)->second;
            handle_t handle = graph.get_handle_of_step(step);
            graph.for_each_step_on_handle(
                handle,
                [&](const step_handle_t& s) {
                    if (graph.get_path_handle_of_step(s) == path) {
                        const auto& other_pos = step_pos.find(s)->second;
                        std::cout << path_name << "\t"
                                  << curr_pos << "\t"
                                  << other_pos << std::endl;
                        /*
                        if (other_pos > curr_pos
                            && other_pos > start_pos
                            && other_pos < end_pos
                            && other_pos - start_pos > min_length) {
                        }
                        */
                    }
                });
        });
}

ska::flat_hash_map<step_handle_t, uint64_t> make_step_index(
    const PathHandleGraph& graph,
    const std::vector<path_handle_t>& paths,
    const size_t& num_threads) {
    ska::flat_hash_map<step_handle_t, uint64_t> step_pos;
#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
    for (auto& path : paths) {
        uint64_t pos = 0;
        graph.for_each_step_in_path(
            path,
            [&](const step_handle_t& step) {
#pragma omp critical (step_pos)
                step_pos[step] = pos;
                handle_t handle = graph.get_handle_of_step(step);
                pos += graph.get_length(handle);
            });
#pragma omp critical (step_pos)
        step_pos[graph.path_end(path)] = pos; // record the end position
    }
    return step_pos;
}

void show_steps(
    const PathHandleGraph& graph,
    const ska::flat_hash_map<step_handle_t, uint64_t>& steps) {
    for (auto& pos : steps) {
        auto h = graph.get_handle_of_step(pos.first);
        auto p = graph.get_path_handle_of_step(pos.first);
        auto name = graph.get_path_name(p);
        auto id = graph.get_id(h);
        auto rev = graph.get_is_reverse(h);
        std::cerr << name << " " << id << (rev?"-":"+") << " " << pos.second << std::endl;
    }
}

// compute the reference segmentations
// and map them onto the graph using a static multiset index structure based on two arrays
// we'll build up a big vector of node -> path segment pairings
// then we'll sort them and build an index
segment_map_t::segment_map_t(
    const PathHandleGraph& graph,
    const std::vector<path_handle_t>& paths,
    const step_index_t& step_index,
    const std::function<bool(const handle_t&)>& is_cut,
    const uint64_t& merge_dist,
    const size_t& num_threads,
    const bool& show_progress) {

    std::unique_ptr<algorithms::progress_meter::ProgressMeter> progress;

    if (show_progress) {
        std::cerr << "[odgi::algorithms::untangle] building target segment index" << std::endl;

        progress = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                paths.size(), "[odgi::algorithms::untangle] untangle and merge cuts");
    }

    std::vector<std::pair<uint64_t, int64_t>> node_to_segment;
    std::vector<std::vector<step_handle_t>> all_cuts(paths.size());
#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
    for (uint64_t i = 0; i < paths.size(); ++i) { //auto& path : paths) {
        auto& path = paths[i];
        auto path_name = graph.get_path_name(path);
        //std::cerr << "path " << path_name << std::endl;
        auto self_index = path_step_index_t(graph, path, 1);
        all_cuts[i] =
            merge_cuts(
                untangle_cuts(graph,
                              path,
                              path_name,
                              self_index.path_begin(),
                              self_index.path_back(),
                              step_index,
                              self_index,
                              is_cut
                    ),
                merge_dist,
                step_index,
				graph);

        if (show_progress) {
            progress->increment(1);
        }
    }
    if (show_progress) {
        progress->finish();
    }
    // the index construction must be serial

    // Put fake stuff in the 1-st position to avoid having segments with id 0
    // becahse we can't discriminate +0 and -0 for the strandness
    if (show_progress) {
        progress = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                paths.size(), "[odgi::algorithms::untangle] prepare segment cuts");
    }
    segment_cut.push_back(graph.path_begin(paths[0]));
    segment_length.push_back(0);

    for (uint64_t i = 0; i < paths.size(); ++i) {
        auto& path = paths[i];
        auto& cuts = all_cuts[i];
        //std::cerr << "reference segmentation" << std::endl;
        //write_cuts(graph, path, cuts, step_pos);
        // walk the path to get the segmentation
        uint64_t curr_segment_idx = 0;
        uint64_t segment_idx = segment_cut.size();
        uint64_t* curr_length = nullptr;
        for (step_handle_t step = graph.path_begin(path);
             step != graph.path_end(path);
             step = graph.get_next_step(step)) {
            // if we are at a segment cut
            if (step == cuts[curr_segment_idx]) {
                segment_idx = segment_cut.size();
                segment_cut.push_back(step);
                segment_length.push_back(0);
                curr_length = &segment_length.back();
                ++curr_segment_idx;
            }
            handle_t h = graph.get_handle_of_step(step);
            bool is_rev = graph.get_is_reverse(h);
            node_to_segment.push_back(
                std::make_pair(graph.get_id(h),
                               (is_rev ? -segment_idx : segment_idx)));
            uint64_t node_length = graph.get_length(h);
            *curr_length += node_length;
        }

        if (show_progress) {
            progress->increment(1);
        }
    }

    if (show_progress) {
        progress->finish();
    }
    //std::cerr << "segment_cut.size() " << segment_cut.size() << std::endl;
    //std::cerr << "segment_length.size() " << segment_length.size() << std::endl;

    ips4o::parallel::sort(node_to_segment.begin(),
                          node_to_segment.end(),
                          std::less<>(),
                          num_threads);

    if (show_progress) {
        progress = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                node_to_segment.size(), "[odgi::algorithms::untangle] make the mapping");
    }

    // make the mapping
    uint64_t prev_node = 0;
    for (auto& node_segment : node_to_segment) {
        auto& node_id = node_segment.first;
        if (node_id > prev_node) {
            while (prev_node < node_id) {
                node_idx.push_back(segments.size());
                ++prev_node;
            }
        }
        segments.push_back(node_segment.second);
        prev_node = node_id;

        if (show_progress) {
            progress->increment(1);
        }
    }
    if (show_progress) {
        progress->finish();
    }

    auto max_id = graph.get_node_count();
    while (prev_node < max_id) {
        node_idx.push_back(segments.size());
        ++prev_node;
    }
    node_idx.push_back(segments.size()); // to avoid special casing the last node
}

void segment_map_t::for_segment_on_node(
    uint64_t node_id,
    const std::function<void(const uint64_t& segment_id, const bool& is_rev)>& func) const {
    uint64_t from = node_idx[node_id-1];
    uint64_t to = node_idx[node_id];
    for (uint64_t i = from; i < to; ++i) {
        auto& j = segments[i];
        func(std::abs(j), j < 0);
    }
}

uint64_t segment_map_t::get_segment_length(const uint64_t& segment_id) const {
    return segment_length.at(segment_id);
}

struct isec_t {
    uint64_t len = 0;
    uint64_t inv = 0;
    void incr(const uint64_t& l, const bool& is_inv) {
        len += l;
        inv += (is_inv ? l : 0);
    }
};

std::vector<segment_mapping_t>
segment_map_t::get_matches(
        const PathHandleGraph& graph,
        const step_handle_t& start,
        const step_handle_t& end,
        const uint64_t& query_length) const {
    // collect the target segments that overlap our segment
    // computing the intersection size (in bp) as we go
    // our final metric is jaccard of intersection over total length for each overlapped target
    //path_handle_t query_path = graph.get_path_handle_of_step(start);
    ska::flat_hash_map<uint64_t, isec_t> target_isec;
    ska::flat_hash_map<uint64_t, uint64_t> query_seen;
    for (step_handle_t step = start;
         step != end;
         step = graph.get_next_step(step)) {
        handle_t h = graph.get_handle_of_step(step);
        uint64_t node_id = graph.get_id(h);
        uint64_t node_length = graph.get_length(h);
        bool is_rev = graph.get_is_reverse(h);
        uint64_t query_idx = query_seen[node_id]++;
        ska::flat_hash_map<uint64_t, uint64_t> target_seen;
        for_segment_on_node(
            node_id,
            [&](const uint64_t& segment_id, const bool& segment_rev) {
                // n.b. we do not skip self matches
                uint64_t target_idx = target_seen[segment_id]++;
                if (query_idx == target_idx) {
                    target_isec[segment_id].incr(node_length, is_rev != segment_rev);
                }
            });
    }
    // compute the jaccards
    path_handle_t curr_path = graph.get_path_handle_of_step(start);
    std::vector<segment_mapping_t> jaccards;
    for (auto& p : target_isec) {
        auto& segment_id = p.first;
        path_handle_t segment_path
            = graph.get_path_handle_of_step(get_segment_cut(segment_id));
        auto& isec = p.second.len;
        //auto& inv = p.second.inv;
        bool is_inv = (double)p.second.inv/(double)isec > 0.5;
        // intersection / union
        jaccards.push_back(
            {
                segment_id,
                segment_path == curr_path,
                is_inv,
                (double)isec
                / (double)(get_segment_length(segment_id)
                           + query_length - isec)
            });
    }
    // sort the target segments by their jaccard with the query
    // keep self-matches ahead of equivalent non-self matches
    // and order equals using segment_id
    std::sort(jaccards.begin(), jaccards.end(),
              [](const segment_mapping_t& a,
                 const segment_mapping_t& b) {
                  return std::tie(a.jaccard, a.self_map, a.is_inv, a.segment_id) >
                      std::tie(b.jaccard, b.self_map, b.is_inv, b.segment_id);
              });
    return jaccards;
}

const step_handle_t& segment_map_t::get_segment_cut(
    const uint64_t& idx) const {
    return segment_cut[idx];
}

double self_mean_coverage(
    const PathHandleGraph& graph,
    const path_step_index_t& self_index,
    const path_handle_t& path,
    const step_handle_t& begin,
    const step_handle_t& end) {
    uint64_t sum = 0;
    uint64_t bp = 0;
    for (step_handle_t step = begin;
         step != end;
         step = graph.get_next_step(step)) {
        handle_t handle = graph.get_handle_of_step(step);
        uint64_t len = graph.get_length(handle);
        bp += len;
        sum += len * self_index.n_steps_on_node(graph.get_id(handle));
    }
    return (double)sum / (double)bp;
}

uint64_t query_hits_target_front(
		const PathHandleGraph& graph,
		const path_handle_t& query,
		const atomicbitvector::atomic_bv_t& targets_node_idx) {
	// currently, we don't care about self mapping
	bool tip_reached_target = false;
	step_handle_t cur_step = graph.path_begin(query);
	handle_t cur_h = graph.get_handle_of_step(cur_step);
	uint64_t cur_id = graph.get_id(cur_h);
	while (!tip_reached_target) {
		if (targets_node_idx.test(cur_id)) {
			tip_reached_target = true;
			return cur_id;
		}
		if (graph.has_next_step(cur_step)) {
			cur_step = graph.get_next_step(cur_step);
			cur_h = graph.get_handle_of_step(cur_step);
			cur_id = graph.get_id(cur_h);
		} else {
			tip_reached_target = true;
		}
	}
	return 0;
}

uint64_t query_hits_target_back(
		const PathHandleGraph& graph,
		const path_handle_t& query,
		const atomicbitvector::atomic_bv_t& targets_node_idx) {
	// currently, we don't care about self mapping
	bool tip_reached_target = false;
	step_handle_t cur_step = graph.path_back(query);
	handle_t cur_h = graph.get_handle_of_step(cur_step);
	uint64_t cur_id = graph.get_id(cur_h);
	while (!tip_reached_target) {
		if (targets_node_idx.test(cur_id)) {
			tip_reached_target = true;
			return cur_id;
		}
		if (graph.has_previous_step(cur_step)) {
			cur_step = graph.get_previous_step(cur_step);
			cur_h = graph.get_handle_of_step(cur_step);
			cur_id = graph.get_id(cur_h);
		} else {
			tip_reached_target = true;
		}
	}
	return 0;
}


void map_segments(
    const PathHandleGraph& graph,
    const path_handle_t& path,
    const std::vector<step_handle_t>& cuts,
    const segment_map_t& target_segments,
    const step_index_t& step_index,
    const path_step_index_t& self_index,
    const double& max_self_coverage,
    const uint64_t& n_best,
    const double& min_jaccard,
    const untangle_output_t& output_type,
    ska::flat_hash_map<path_handle_t, uint64_t>& path_to_len) {
    // query name is the first field in our outputs
    std::string query_name = graph.get_path_name(path);
    // helper for building up gene order lists and gggenes plot data
    struct path_range_t {
        path_handle_t target_path;
        uint64_t query_begin;
        uint64_t query_end;
        uint64_t target_begin;
        uint64_t target_end;
        bool is_inv;
    };
    std::vector<path_range_t> gene_order;
    for (uint64_t i = 0; i < cuts.size()-1; ++i) {
        auto& begin = cuts[i];
        auto& end = cuts[i+1];
        auto begin_pos = step_index.get_position(begin, graph);
        auto end_pos = step_index.get_position(end, graph);
        uint64_t length = end_pos - begin_pos;
        // get the self coverage TODO
        double self_coverage = self_mean_coverage(graph, self_index, path, begin, end);
        if (max_self_coverage && self_coverage > max_self_coverage) continue;
        std::vector<segment_mapping_t> target_mapping =
            target_segments.get_matches(graph, cuts[i], cuts[i+1], length);
        uint64_t nth_best = 0;
        //std::cerr << "search-> " << query_name << "\t" << begin_pos << "\t" << end_pos << "\t" << target_mapping.size() << std::endl;
        // todo for each target mapping, up to the Nth, emitting if they're >= min_jaccard
        if (!target_mapping.empty()) {
            for (auto& mapping : target_mapping) {
                ++nth_best;
                if (nth_best > n_best) break;
                auto& jaccard = mapping.jaccard;
                if (jaccard >= min_jaccard) {
                    double dist = -log(2.0 * jaccard / (1. + jaccard));
                    if (dist > 1.0) dist = 1.0;
                    auto& idx = mapping.segment_id; // segment index
                    auto& target_begin = target_segments.get_segment_cut(idx);
                    auto target_begin_pos = step_index.get_position(target_begin, graph);
                    auto target_end_pos = target_begin_pos + target_segments.get_segment_length(idx);
                    path_handle_t target_path = graph.get_path_handle_of_step(target_begin);
                    std::string target_name = graph.get_path_name(target_path);
                    if (output_type == untangle_output_t::PAF){
                        // PAF format
#pragma omp critical (cout)
                        std::cout << query_name << "\t"
                        << path_to_len[path] << "\t"
                        << begin_pos << "\t"
                        << end_pos << "\t"          // Query end (0-based; BED-like; open)
                        << (mapping.is_inv ? "-" : "+") << "\t"
                        << target_name << "\t"
                        << path_to_len[target_path] << "\t"
                        << target_begin_pos << "\t"
                        << target_end_pos << "\t"    // Target end (0-based; BED-like; open)
                        << 0 << "\t"
                        << std::max(target_end_pos - target_begin_pos, end_pos - begin_pos) << "\t"
                        << 255 << "\t"
                        << "id:f:" << ((double) 1.0 - dist) * (double) 100 << "\t"
                        << "jc:f:" << jaccard << "\t"
                        << "sc:f:" << self_coverage << "\t"
                        << "nb:i:" << nth_best << "\t"
                        << std::endl;
                    } else if (output_type == untangle_output_t::ORDER
                               || output_type == untangle_output_t::GGGENES
                               || output_type == untangle_output_t::SCHEMATIC) {
                        if (gene_order.size() && gene_order.back().target_path == target_path
                            && gene_order.back().query_end == begin_pos
                            && gene_order.back().target_end == target_begin_pos
                            && gene_order.back().is_inv == mapping.is_inv) {
                            // extend the last range
                            gene_order.back().query_end = end_pos;
                            gene_order.back().target_end = target_end_pos;
                        } else {
                            gene_order.push_back({ target_path, begin_pos,
                                    end_pos, target_begin_pos, target_end_pos,
                                    mapping.is_inv });
                        }
                    } else if (output_type == untangle_output_t::BEDPE) {
#pragma omp critical (cout)
                        // BEDPE format
                        std::cout << query_name << "\t"
                        << begin_pos << "\t"
                        << end_pos << "\t"              // chrom1 end (1-based)
                        << target_name << "\t"
                        << target_begin_pos << "\t"
                        << target_end_pos << "\t"       // chrom2 end (1-based)
                        << jaccard << "\t"
                        << (mapping.is_inv ? "-" : "+") << "\t"
                        << self_coverage << "\t"
                        << nth_best << std::endl;
                    }

                }
            }
        }
    }
    if (output_type == untangle_output_t::ORDER) {
        std::stringstream ss;
        ss << query_name << "\t";
        for (auto& range : gene_order) {
            ss << graph.get_path_name(range.target_path) << ":"
               << range.target_begin << "-" << range.target_end << ",";
        }
        std::string s = ss.str();
        if (s.size() && s.at(s.size()-1) == ',') { s.pop_back(); }
#pragma omp critical (cout)
        std::cout << s << std::endl;
    }
    if (output_type == untangle_output_t::GGGENES
        || output_type == untangle_output_t::SCHEMATIC) {
        std::stringstream ss;
        if (output_type == untangle_output_t::SCHEMATIC) {
            uint64_t idx = 0;
            for (auto& range : gene_order) {
                range.query_begin = idx;
                idx += 100;
                range.query_end = idx;
                idx += 50;
            }
        }
        for (auto& range : gene_order) {
            ss << query_name << "\t"
               << graph.get_path_name(range.target_path) << "\t"
               << range.query_begin << "\t"
               << range.query_end << "\t"
               << (range.is_inv ? "0" : "1") << std::endl;
        }
#pragma omp critical (cout)
        std::cout << ss.str();
    }
}

// BEDPE (pair-BED) projection of the graph
// that describe nonlinear query : target relationships
void untangle(
    const PathHandleGraph& graph,
    const std::vector<path_handle_t>& queries,
    const std::vector<path_handle_t>& targets,
    const uint64_t& merge_dist,
    const double& max_self_coverage,
    const uint64_t& n_best,
    const double& min_jaccard,
    const uint64_t& cut_every,
    const untangle_output_t& output_type,
    const std::string& cut_points_input,
    const std::string& cut_points_output,
    const size_t& num_threads,
    const bool& show_progress,
	const step_index_t& step_index,
	const std::vector<path_handle_t>& paths) {

    std::unique_ptr<algorithms::progress_meter::ProgressMeter> progress;

    if (show_progress) {
        std::cerr << "[odgi::algorithms::untangle] untangling " << queries.size() << " queries with " << targets.size() << " targets" << std::endl;
    }

    int threads_per = std::max(1, (int)std::floor((double)num_threads/(double)paths.size()));

    // collect all possible cuts
    // we'll use this to drive the subsequent segmentation
    atomicbitvector::atomic_bv_t cut_nodes(graph.get_node_count()+1);

    if (cut_points_input.empty()) {
        if (show_progress) {
            std::cerr << "[odgi::algorithms::untangle] establishing initial cuts for " << paths.size() << " paths" << std::endl;
        }

        if (show_progress) {
            progress = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                    targets.size(), "[odgi::algorithms::untangle] set target nodes");
        }

        // which nodes are traversed by our target paths?
        atomicbitvector::atomic_bv_t target_nodes(graph.get_node_count() + 1);
#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
        for (auto& target : targets) {
            graph.for_each_step_in_path(
                target, [&](const step_handle_t& step) {
                    target_nodes.set(graph.get_id(graph.get_handle_of_step(step)), true);
                });

            if (show_progress) {
                progress->increment(1);
            }
        }
        if (show_progress) {
            progress->finish();
        }

        if (show_progress) {
            progress = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                    targets.size(), "[odgi::algorithms::untangle] untangle and merge cuts");
        }

#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
        for (auto& path : paths) {
            // test path_step_index_t
            auto self_index = path_step_index_t(graph, path, threads_per);
            auto path_name = graph.get_path_name(path);
            std::vector<step_handle_t> cuts
            = merge_cuts(
                    untangle_cuts(graph,
                                  path,
                                  path_name,
                                  self_index.path_begin(),
                                  self_index.path_back(),
                                  step_index,
                                  self_index,
                                  [](const handle_t& h) { return false; }),
                    merge_dist,
                    step_index,
                    graph);
            for (auto& step : cuts) {
                cut_nodes.set(graph.get_id(graph.get_handle_of_step(step)));
            }
            // also add the nodes here where the query path touches the target for the first time
            // we start from the front until we found a target node
            const uint64_t node_id_front = query_hits_target_front(graph, path, target_nodes);
            if (graph.has_node(node_id_front)) {
                cut_nodes.set(node_id_front, true);
            }
            // we start from the back until we found a target node
            const uint64_t node_id_back = query_hits_target_back(graph, path, target_nodes);
            if (graph.has_node(node_id_back)) {
                cut_nodes.set(node_id_back, true);
            }

            if (show_progress) {
                progress->increment(1);
            }
        }

        if (show_progress) {
            progress->finish();
        }

        if (cut_every > 0) {
            if (show_progress) {
                progress = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                        graph.get_node_count(), "[odgi::algorithms::untangle] mark nodes every " + to_string(cut_every) + " bp");
            }

            /*
    #pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
            for (auto& target : targets) {
                uint64_t pos = 0;
                uint64_t last = 0;
                graph.for_each_step_in_path(
                    target, [&](const step_handle_t& step) {
                    });
            }
            */
            // walk along the node space in sorted order
            // marking nodes every cut_every bp
            uint64_t pos = 0;
            uint64_t last = 0;
            uint64_t segment = 0;
            std::vector<uint64_t> node_to_segment(graph.get_node_count()+1);
            graph.for_each_handle(
                [&](const handle_t& h) {
                    auto l = graph.get_length(h);
                    pos += l;
                    if (pos - last > cut_every) {
                        last = pos;
                        // this is possible, but may introduce too many cut points
                        //cut_nodes.set(graph.get_id(h), true);
                        ++segment;
                    }
                    node_to_segment[graph.get_id(h)] = segment;

                    if (show_progress) {
                        progress->increment(1);
                    }
                });

            if (show_progress) {
                progress->finish();
            }

            // todo: split up the graph space into regions of cut_every bp
            // write a map from node ids to segments
            // walk along every path
            // mark cut points the first nodes in each segment that we get to

            if (show_progress) {
                progress = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                        paths.size(), "[odgi::algorithms::untangle] add new cut points");
            }

    #pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
            for (auto& path : paths) {
                uint64_t segment = 0;
                uint64_t last = 0;
                graph.for_each_step_in_path(
                    path, [&](const step_handle_t& step) {
                        auto h = graph.get_handle_of_step(step);
                        auto id = graph.get_id(h);
                        auto segment = node_to_segment[id];
                        if (segment != last) {
                            cut_nodes.set(id, true);
                        }
                        last = segment;
                    });

                if (show_progress) {
                    progress->increment(1);
                }
            }

            if (show_progress) {
                progress->finish();
            }

        }
    } else {
        uint64_t num_cut_points_read = 0;

        if (show_progress) {
            std::cerr << "[odgi::algorithms::untangle] loading input cuts" << std::endl;
        }
        std::ifstream bed_in(cut_points_input.c_str());
        std::string buffer;
        while (std::getline(bed_in, buffer)) {
            if (!buffer.empty()) {
                const uint64_t handle_id = std::stoull(buffer);

                if (!graph.has_node(handle_id)) {
                    std::cerr << "[odgi::algorithms::untangle] error: node identifier " << handle_id << " not found in graph" << std::endl;
                    exit(1);
                }

                cut_nodes.set(
                        handle_id,
                        true);
                ++num_cut_points_read;
            }
        }

        if (num_cut_points_read == 0) {
            std::cerr << "[odgi::algorithms::untangle] error: no cut points loaded" << std::endl;
            exit(1);
        }

        if (show_progress) {
            std::cerr << "[odgi::algorithms::untangle] loaded " << num_cut_points_read << " cuts points" << std::endl;
        }
    }

    //auto step_pos = make_step_index(graph, queries);
    // node to reference segmentation mapping

    segment_map_t target_segments(graph,
                                  targets,
                                  step_index,
                                  [&cut_nodes,&graph](const handle_t& h) {
                                      return cut_nodes.test(graph.get_id(h));
                                  },
                                  merge_dist,
                                  num_threads,
                                  show_progress);

    //show_steps(graph, step_pos);
    //std::cout << "path\tfrom\tto" << std::endl;
    //auto step_pos = make_step_index(graph, queries);

    ska::flat_hash_map<path_handle_t, uint64_t> path_to_len;

    if (output_type == untangle_output_t::PAF) {
        // PAF format
        auto get_path_length = [](const PathHandleGraph &graph, const path_handle_t &path_handle) {
            uint64_t path_len = 0;
            graph.for_each_step_in_path(path_handle, [&](const step_handle_t &s) {
                path_len += graph.get_length(graph.get_handle_of_step(s));
            });
            return path_len;
        };

#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
        for (uint64_t i = 0; i < paths.size(); ++i) {
            auto& path = paths[i];
            const uint64_t path_len = get_path_length(graph, paths[i]);

            // You can't write on such a data structure in parallel
#pragma omp critical (path_to_len)
            path_to_len[path] = path_len;
        }
    } else if (output_type == untangle_output_t::BEDPE) {
        std::cout << "#query.name\tquery.start\tquery.end\tref.name\tref.start\tref.end\tscore\tinv\tself.cov\tnth.best" << std::endl;
    } else if (output_type == untangle_output_t::GGGENES
               || output_type == untangle_output_t::SCHEMATIC) {
        // gggenes format
        std::cout << "molecule\tgene\tstart\tend\tstrand" << std::endl;
    } else if (output_type == untangle_output_t::ORDER) {
        // nothing to do
    }

    if (show_progress) {
        progress = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                queries.size(), "[odgi::algorithms::untangle] untangling " + to_string(queries.size()) + " queries");
    }

#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
    for (auto& query : queries) {
        auto self_index = path_step_index_t(graph, query, threads_per);
        auto query_name = graph.get_path_name(query);
        std::vector<step_handle_t> cuts
            = merge_cuts(
                untangle_cuts(graph,
                              query,
                              query_name,
                              self_index.path_begin(),
                              self_index.path_back(),
                              step_index,
                              self_index,
                              [&cut_nodes,&graph](const handle_t& h) {
                                  return cut_nodes.test(graph.get_id(h));
                              }),
                merge_dist,
                step_index,
				graph);
        map_segments(graph, query, cuts, target_segments,
                     step_index, self_index,
                     max_self_coverage, n_best, min_jaccard,
                     output_type, path_to_len);

        //write_cuts(graph, query, cuts, step_pos);

        if (show_progress) {
            progress->increment(1);
        }
    }

    if (show_progress) {
        progress->finish();
    }

    //self_dotplot(graph, query, step_pos);

    // If requested, write cut points to a file
    if (!cut_points_output.empty()) {
        std::ofstream f(cut_points_output.c_str());
        for (auto x : cut_nodes) {
            f << x << std::endl;
        }
        f.close();
    }
}

}
}
