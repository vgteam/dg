#include "node.hpp"

namespace odgi {

uint64_t node_t::sequence_size() const {
    return sequence.size();
}

void node_t::set_sequence(const std::string& seq) {
    sequence = seq;
}

void node_t::set_id(const uint64_t& new_id) {
    id = new_id;
}

const uint64_t& node_t::get_id() const {
    return id;
}

const std::string& node_t::get_sequence() const {
    return sequence;
}

void node_t::set_volatile(void) {
    dynamic = true;
}

void node_t::set_static(void) {
    dynamic = false;
}

// encode an internal representation of an external id (adding if none exists)
uint64_t node_t::encode(const uint64_t& other_id) {
    uint64_t delta = to_delta(other_id);
    uint64_t i = 0;
    uint64_t s = decoding.size();
    while (i < s && decoding.at(i) != delta) {
        ++i;
    }
    if (i == s) {
        decoding.push_back(delta);
    }
    return i;
}

// decode an internal representation of an external id
uint64_t node_t::decode(const uint64_t& idx) const {
    return from_delta(decoding.at(idx));
}

void node_t::for_each_edge(const std::function<bool(uint64_t other_id,
                                                    bool other_rev,
                                                    bool to_curr,
                                                    bool on_rev)>& func) const {
    for (uint64_t i = 0; i < edges.size(); ) {
        uint64_t other_id = edges.at(i++);
        uint8_t packed_edge = edges.at(i++);
        if (!func(other_id,
                  edge_helper::unpack_other_rev(packed_edge),
                  edge_helper::unpack_to_curr(packed_edge),
                  edge_helper::unpack_on_rev(packed_edge))) {
            break;
        }
    }
}

void node_t::add_edge(const uint64_t& other_id,
                      const bool& other_rev,
                      const bool& to_curr,
                      const bool& on_rev) {
    //std::cerr << "add edge " << "relative_id " << relative_id << " edge_type " << edge_type << std::endl;
    auto edge_type = edge_helper::pack(other_rev, to_curr, on_rev);
    edges.push_back(other_id);
    edges.push_back(edge_type);
}

bool node_t::remove_edge(const uint64_t& target_id,
                         const bool& target_rev,
                         const bool& ends_here,
                         const bool& is_rev) {
    for (uint64_t i = 0; i < edges.size(); i+=EDGE_RECORD_LENGTH) {
        uint64_t other_id = edges.at(i);
        if (other_id == target_id) {
            uint8_t packed_edge = edges.at(i+1);
            bool on_rev = edge_helper::unpack_on_rev(packed_edge);
            bool other_rev = edge_helper::unpack_other_rev(packed_edge);
            bool to_curr = edge_helper::unpack_to_curr(packed_edge);
            if (is_rev != on_rev) {
                other_rev ^= 1;
                to_curr ^= 1;
            }
            if (other_rev == target_rev
                && to_curr == ends_here) {
                edges.remove(i);
                edges.remove(i);
                return true;
            }
        }
    }
    return false;
}

void node_t::add_path_step(const uint64_t& path_id, const bool& is_rev,
                           const bool& is_start, const bool& is_end,
                           const uint64_t& prev_id, const uint64_t& prev_rank,
                           const uint64_t& next_id, const uint64_t& next_rank) {
    //std::cerr << "packing " << path_id << " " << is_rev << " " << is_start << " " << is_end << std::endl;
    paths.push_back(path_id);
    paths.push_back(step_type_helper::pack(is_rev, is_start, is_end));
    // we store the smallest possible delta for path starts and ends
    paths.push_back(encode(!is_start ? prev_id : id));
    paths.push_back(prev_rank);
    paths.push_back(encode(!is_end ? next_id : id));
    paths.push_back(next_rank);
}

void node_t::add_path_step(const node_t::step_t& step) {
    add_path_step(step.path_id,
                  step.is_rev,
                  step.is_start,
                  step.is_end,
                  step.prev_id,
                  step.prev_rank,
                  step.next_id,
                  step.next_rank);
}

const std::vector<node_t::step_t> node_t::get_path_steps() const {
    uint64_t n_paths = path_count();
    if (n_paths == 0) return {};
    std::vector<node_t::step_t> steps(n_paths);
    for (uint64_t i = 0; i < n_paths; ++i) {
        steps[i] = get_path_step(i);
    }
    return steps;
}

const node_t::step_t node_t::get_path_step(const uint64_t& rank) const {
    if (rank >= path_count()) assert(false);
    uint64_t i = PATH_RECORD_LENGTH*rank;
    uint64_t t = paths.at(i+1);
    //std::cerr << "got paths at i+2 " << paths.at(i+2) << std::endl;
    return {
        paths.at(i),
        step_type_helper::unpack_is_rev(t),
        step_type_helper::unpack_is_start(t),
        step_type_helper::unpack_is_end(t),
        decode(paths.at(i+2)),
        paths.at(i+3),
        decode(paths.at(i+4)),
        paths.at(i+5),
    };
}

void node_t::for_each_path_step(
    const std::function<bool(uint64_t rank,
                             uint64_t path_id,
                             bool is_rev)>& func) const {
    uint64_t n_paths = path_count();
    for (uint64_t i = 0; i < n_paths; ++i) {
        if (!step_is_del(i) && !func(i, paths.at(PATH_RECORD_LENGTH*i), step_is_rev(i))) {
            break;
        }
    }
}

void node_t::for_each_path_step(const std::function<bool(step_t step)>& func) const {
    uint64_t n_paths = path_count();
    for (uint64_t i = 0; i < n_paths; ++i) {
        if (!step_is_del(i) && !func(get_path_step(i))) {
            break;
        }
    }
}

void node_t::set_path_step(const uint64_t& rank, const step_t& step) {
    if (rank >= path_count()) assert(false);
    uint64_t i = PATH_RECORD_LENGTH*rank;
    paths[i] = step.path_id;
    paths[i+1] = step_type_helper::pack(step.is_rev,
                                        step.is_start,
                                        step.is_end);
    paths[i+2] = encode(!step.is_start ? step.prev_id : id);
    paths[i+3] = step.prev_rank;
    paths[i+4] = encode(!step.is_end ? step.next_id : id);
    paths[i+5] = step.next_rank;
}

void node_t::clear_path_step(const uint64_t& rank) {
    if (rank >= path_count()) assert(false);
    uint64_t i = PATH_RECORD_LENGTH*rank;
    paths[i] = 0;
    paths[i+1] = step_type_helper::pack(false, false, false);
    paths[i+2] = encode(id);
    paths[i+3] = 0;
    paths[i+4] = encode(id);
    paths[i+5] = 0;
}

void node_t::set_step_path_id(const uint64_t& rank, const uint64_t& path_id) {
    paths[PATH_RECORD_LENGTH*rank] = path_id;
}

void node_t::set_step_prev_id(const uint64_t& rank, const uint64_t& prev_id) {
    paths[PATH_RECORD_LENGTH*rank+2] = encode(prev_id);
}

void node_t::set_step_prev_rank(const uint64_t& rank, const uint64_t& prev_rank) {
    paths[PATH_RECORD_LENGTH*rank+3] = prev_rank;
}

void node_t::set_step_next_id(const uint64_t& rank, const uint64_t& next_id) {
    paths[PATH_RECORD_LENGTH*rank+4] = encode(next_id);
}

void node_t::set_step_next_rank(const uint64_t& rank, const uint64_t& next_rank) {
    paths[PATH_RECORD_LENGTH*rank+5] = next_rank;
}

void node_t::set_step_is_rev(const uint64_t& rank, const bool& is_rev) {
    uint64_t idx = PATH_RECORD_LENGTH*rank+1;
    paths[idx] = paths[idx] & ~(1UL) | is_rev;
}

void node_t::set_step_is_start(const uint64_t& rank, const bool& is_start) {
    uint64_t idx = PATH_RECORD_LENGTH*rank+1;
    paths[idx] = paths[idx] & ~(1UL << 1) | is_start;
}

void node_t::set_step_is_end(const uint64_t& rank, const bool& is_end) {
    uint64_t idx = PATH_RECORD_LENGTH*rank+1;
    paths[idx] = paths[idx] & ~(1UL << 2) | is_end;
}

void node_t::set_step_is_del(const uint64_t& rank, const bool& is_del) {
    uint64_t idx = PATH_RECORD_LENGTH*rank+1;
    paths[idx] = paths[idx] & ~(1UL << 3) | is_del;
}

uint64_t node_t::step_path_id(const uint64_t& rank) const {
    return paths.at(PATH_RECORD_LENGTH*rank);
}

uint64_t node_t::step_prev_id(const uint64_t& rank) const {
    return decode(paths.at(PATH_RECORD_LENGTH*rank+2));
}

uint64_t node_t::step_prev_rank(const uint64_t& rank) const {
    return paths.at(PATH_RECORD_LENGTH*rank+3);
}

uint64_t node_t::step_next_id(const uint64_t& rank) const {
    return decode(paths.at(PATH_RECORD_LENGTH*rank+4));
}

uint64_t node_t::step_next_rank(const uint64_t& rank) const {
    return paths.at(PATH_RECORD_LENGTH*rank+5);
}

bool node_t::step_is_rev(const uint64_t& rank) const {
    return step_type_helper::unpack_is_rev(paths.at(PATH_RECORD_LENGTH*rank+1));
}

bool node_t::step_is_start(const uint64_t& rank) const {
    return step_type_helper::unpack_is_start(paths.at(PATH_RECORD_LENGTH*rank+1));
}

bool node_t::step_is_end(const uint64_t& rank) const {
    return step_type_helper::unpack_is_end(paths.at(PATH_RECORD_LENGTH*rank+1));
}

bool node_t::step_is_del(const uint64_t& rank) const {
    return step_type_helper::unpack_is_del(paths.at(PATH_RECORD_LENGTH*rank+1));
}

std::pair<std::map<uint64_t, std::pair<uint64_t, bool>>, // path fronts
          std::map<uint64_t, std::pair<uint64_t, bool>>> // path backs
node_t::flip_paths() {
    const std::vector<node_t::step_t> steps = get_path_steps();
    // remove all path steps
    clear_paths();
    // flip them and replace, recording which path starts and ends should be rewritten
    std::pair<std::map<uint64_t, std::pair<uint64_t, bool>>,
              std::map<uint64_t, std::pair<uint64_t, bool>>> path_start_end_rewrites;
    uint64_t rank = 0;
    for (auto& step : steps) {
        // flip the step
        step_t flipped = { step.path_id, !step.is_rev,
                           step.is_start, step.is_end,
                           step.prev_id, step.prev_rank,
                           step.next_id, step.next_rank };
        if (step.is_start) {
            path_start_end_rewrites.first[step.path_id] = std::make_pair(rank, !step.is_rev);
        }
        if (step.is_end) {
            path_start_end_rewrites.second[step.path_id] = std::make_pair(rank, !step.is_rev);
        }
        add_path_step(flipped);
        ++rank;
    }
    return path_start_end_rewrites;
}

void node_t::remove_path_step(const uint64_t& rank) {
    if (rank >= path_count()) assert(false);
    uint64_t offset = PATH_RECORD_LENGTH*rank;
    for (uint8_t i = 0; i < PATH_RECORD_LENGTH; ++i) {
        paths.remove(offset);
    }
}

node_t::node_t() {
    //edges = dyn::hacked_vector(0,8);
    //decoding = dyn::hacked_vector(0,8);
    //paths = dyn::hacked_vector(0,4);
}

void node_t::clear() {
    sequence.clear(); // not sure this works
    clear_encoding();
    clear_edges();
    clear_paths();
}

void node_t::clear_edges() {
    dyn::hacked_vector null_iv;
    edges = null_iv;
}

void node_t::clear_paths() {
    dyn::hacked_vector null_iv;
    paths = null_iv;
}

void node_t::clear_encoding() {
    dyn::hacked_vector null_iv;
    decoding = null_iv;
}

void node_t::copy(const node_t& other) {
    clear();
    id = other.id;
    dynamic = other.dynamic;
    sequence = other.sequence;
    edges = other.edges;
    decoding = other.decoding;
    paths = other.paths;
}

void node_t::apply_ordering(
    const std::function<uint64_t(uint64_t)>& get_new_id,
    const std::function<bool(uint64_t)>& to_flip) {
    // flip the node sequence if needed
    bool flip = to_flip(id);
    if (flip) {
        reverse_complement_in_place(sequence);
    }
    // rewrite the encoding (affects path storage)
    std::vector<uint64_t> dec_v;
    std::vector<uint64_t> encoding_map; // to compress away deleted nodes
    bool compress_encoding = false;
    uint64_t j = 0;
    for (uint64_t i = 0; i < decoding.size(); ++i) {
        uint64_t old_id = decode(i);
        uint64_t new_id = old_id ? get_new_id(old_id) : 0;
        if (new_id) {
            dec_v.push_back(new_id);
            encoding_map.push_back(j++);
        } else {
            // this means that the node referred to by this entry has been deleted
            // we'll need to rewrite our references to it
            compress_encoding = true;
            encoding_map.push_back(j);
        }
    }
    // update our own id before re-encoding (affects to_delta computation)
    id = get_new_id(id);
    clear_encoding();
    for (auto& other_id : dec_v) {
        decoding.push_back(to_delta(other_id));
    }
    // compress our encoding if requested
    if (compress_encoding) {
        uint64_t n_paths = path_count();
        for (uint64_t i = 0; i < n_paths; ++i) {
            uint64_t q = PATH_RECORD_LENGTH*i;
            paths[q+2] = encoding_map[paths[q+2]];
            paths[q+4] = encoding_map[paths[q+4]];
        }
    }
    // flip path steps if needed
    if (flip) {
        uint64_t n_paths = path_count();
        for (uint64_t i = 0; i < n_paths; ++i) {
            set_step_is_rev(i, !step_is_rev(i));
        }
    }
    // rewrite the edges, reflecting the orientation information we're given
    dyn::hacked_vector new_edges;
    for_each_edge(
        [&](uint64_t other_id,
            bool other_rev,
            bool to_curr,
            bool on_rev) {
            //new_edges.
            auto edge_type = edge_helper::pack(to_flip(other_id)^other_rev,
                                               to_curr,
                                               flip^on_rev);
            new_edges.push_back(get_new_id(other_id));
            new_edges.push_back(edge_type);
            return true;
        });
    edges = new_edges;
}

void node_t::apply_path_ordering(
    const std::function<uint64_t(uint64_t)>& get_new_path_id) {
    uint64_t n_paths = path_count();
    for (uint64_t i = 0; i < n_paths; ++i) {
        if (!step_is_del(i)) {
            set_step_path_id(
                i,
                get_new_path_id(step_path_id(i)));
        }
    }
}

uint64_t node_t::serialize(std::ostream& out) const {
    uint64_t written = 0;
    size_t seq_size = sequence.size();
    out.write((char*)&seq_size, sizeof(size_t));
    written += sizeof(size_t);
    out.write((char*)sequence.c_str(), seq_size*sizeof(char));
    written += seq_size*sizeof(char);
    out.write((char*)&id, sizeof(id));
    written += sizeof(id);
    written += edges.serialize(out);
    written += decoding.serialize(out);
    written += paths.serialize(out);
    return written;
}

void node_t::load(std::istream& in) {
    size_t len = 0;
    in.read((char*)&len, sizeof(size_t));
    sequence.resize(len);
    in.read((char*)sequence.c_str(), len*sizeof(uint8_t));
    in.read((char*)&id, sizeof(id));
    edges.load(in);
    decoding.load(in);
    paths.load(in);
    //display();
}

void node_t::display() const {
    std::cerr << "seq " << sequence << " "
              << "edge_count " << edge_count() << " "
              << "path_count " << path_count();
    std::cerr << " | ";
    if (edge_count()) {
        for (uint64_t i = 0; i < edge_count(); ++i) {
            std::cerr
                << edges.at(i) << ":"
                << edges.at(i+1) << " ";
        }
    }
    std::cerr << " | ";
    if (path_count()) {
        for_each_path_step(
            [&](const step_t& step) {
                std::cerr
                    << step.path_id << ":"
                    << step.is_rev << ":"
                    << step.is_start << ":"
                    << step.is_end << ":"
                    << step.prev_id << ":"
                    << step.prev_rank << ":"
                    << step.next_id << ":"
                    << step.next_rank << " ";
                return true;
            });
    }
    std::cerr << std::endl;
}

}
