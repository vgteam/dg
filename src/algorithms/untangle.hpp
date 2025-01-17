#pragma once

#include <handlegraph/path_handle_graph.hpp>
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/util.hpp>
#include <fstream>
#include <vector>
#include <set>
#include <deque>
#include <atomic_bitvector.hpp>
#include "hash_map.hpp"
#include "ips4o.hpp"
#include "stepindex.hpp"

namespace odgi {
namespace algorithms {

using namespace handlegraph;

enum class untangle_output_t {
    BEDPE,
    PAF,
    GGGENES,
    ORDER,
    SCHEMATIC
};

struct segment_mapping_t {
    uint64_t segment_id = 0;
    bool self_map = false;
    bool is_inv = false;
    double jaccard = 0;
};

class segment_map_t {
public:
    // each segment is identified by its starting step
    std::vector<step_handle_t> segment_cut;
    // and a length
    std::vector<uint64_t> segment_length;
    // maps from node id-1 to idx in segments
    std::vector<uint64_t> node_idx;
    // stores segment assignments sorted by node
    // segment ids stored here map into segment_cuts and segment_lengths
    // +/- indicates path orientation
    std::vector<int64_t> segments;
    segment_map_t(const PathHandleGraph& graph,
                  const std::vector<path_handle_t>& paths,
                  const step_index_t& step_index,
                  const std::function<bool(const handle_t&)>& is_cut,
                  const uint64_t& merge_dist,
                  const size_t& num_threads,
                  const bool& show_progress);
    void for_segment_on_node(
        uint64_t node_id,
        const std::function<void(const uint64_t& segment_id, const bool& is_rev)>& func) const;
    uint64_t get_segment_length(const uint64_t& segment_id) const;
    std::vector<segment_mapping_t> get_matches(
        const PathHandleGraph& graph,
        const step_handle_t& start,
        const step_handle_t& end,
        const uint64_t& query_length) const;
    const step_handle_t& get_segment_cut(
        const uint64_t& idx) const;
};

std::vector<step_handle_t> untangle_cuts(
    const PathHandleGraph& graph,
    const path_handle_t& path,
    const std::string& path_name,
    const path_step_index_t::step_it& start,
    const path_step_index_t::step_it& end,
    const step_index_t& step_index,
    const path_step_index_t& self_index,
    const std::function<bool(const handle_t&)>& is_cut);

std::vector<step_handle_t> merge_cuts(
    const std::vector<step_handle_t>& cuts,
    const uint64_t& dist,
    const step_index_t& step_index,
	const PathHandleGraph& graph);

void write_cuts(
    const PathHandleGraph& graph,
    const path_handle_t& path,
    const std::vector<step_handle_t>& cuts,
    const step_index_t& step_index);

void self_dotplot(
    const PathHandleGraph& graph,
    const path_handle_t& path);

ska::flat_hash_map<step_handle_t, uint64_t> make_step_index(
    const PathHandleGraph& graph,
    const std::vector<path_handle_t>& paths,
    const size_t& num_threads);

double self_mean_coverage(
    const PathHandleGraph& graph,
    const path_step_index_t& self_index,
    const path_handle_t& path,
    const step_handle_t& begin,
    const step_handle_t& end);

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
    const ska::flat_hash_map<path_handle_t, uint64_t>& path_to_len);

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
	const std::vector<path_handle_t>& paths);

}
}
