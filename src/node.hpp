#pragma once

#include <iostream>
#include <cstdint>
#include <string>
#include <handlegraph/util.hpp>
#include <vector>
#include <map>
#include <utility>
#include <cstring>
#include <cassert>
#include <atomic>
// #include "bmap.hpp"
#include "dynamic.hpp"
#include "varint.hpp"
#include "dna.hpp"

namespace odgi {

using namespace handlegraph;
//using nid_t = handlegraph::nid_t;
const uint8_t EDGE_RECORD_LENGTH = 2;
const uint8_t PATH_RECORD_LENGTH = 6;

/// A node object with the sequence, its edge lists, and paths
class node_t {
    uint64_t id = 0;
    std::atomic_flag lock = ATOMIC_FLAG_INIT;
    bool dynamic = true;
    std::string sequence;
    dyn::hacked_vector edges;
    dyn::hacked_vector decoding;
    dyn::hacked_vector paths;
    // relativistic conversions
    inline uint64_t to_delta(const uint64_t& other_id) const {
        if (other_id > id) {
            return ((other_id - id) << 1) | 1;
        } else if (other_id < id) {
            return ((id - other_id) << 1) | 0;
        } else {
            return 0;
        }
    }
    inline uint64_t from_delta(const uint64_t& delta) const {
        if (delta == 0) {
            return id;
        } else if (delta & 1) {
            return id + (delta >> 1);
        } else {
            return id - (delta >> 1);
        }
    }
    /// edge type conversion
    /// 1 = fwd->fwd, 2 = fwd->rev, 3 = rev->fwd, 4 = rev->rev
    struct edge_helper {
        inline static uint8_t pack(bool other_rev, bool to_curr, bool on_rev) {
            return other_rev | (on_rev << 1) | (to_curr << 2);
        }
        inline static uint8_t unpack_on_rev(uint8_t edge) {
            return edge & (1 << 1);
        }
        inline static uint8_t unpack_other_rev(uint8_t edge) {
            return edge & 1;
        }
        inline static uint8_t unpack_to_curr(uint8_t edge) {
            return edge & (1 << 2);
        }
    };
    struct step_type_helper {
        inline static uint8_t pack(const bool& is_rev, const bool& is_start, const bool& is_end) {
            return is_rev | (is_start << 1) | (is_end << 2);
        }
        inline static bool unpack_is_rev(const uint8_t& type) {
            return type & 1;
        }
        inline static bool unpack_is_start(const uint8_t& type) {
            return type & (1 << 1);
        }
        inline static bool unpack_is_end(const uint8_t& type) {
            return type & (1 << 2);
        }
        inline static bool unpack_is_del(const uint8_t& type) {
            return type & (1 << 3);
        }
    };

public:
    node_t(void); // constructor
    // locking methods
    inline void get_lock(void) {
        if (dynamic) {
            while (lock.test_and_set(std::memory_order_acquire))  // acquire lock
                ; // spin
        }
    }
    inline void clear_lock(void) {
        if (dynamic) {
            lock.clear(std::memory_order_release);
        }
    }
    inline const uint64_t edge_count(void) const { return edges.size()/EDGE_RECORD_LENGTH; }
    inline const uint64_t path_count(void) const { return paths.size()/PATH_RECORD_LENGTH; }
    struct step_t {
        uint64_t path_id;
        bool is_rev;
        bool is_start;
        bool is_end;
        uint64_t prev_id;
        uint64_t prev_rank;
        uint64_t next_id;
        uint64_t next_rank;
    };

    uint64_t encode(const uint64_t& other_id);
    uint64_t decode(const uint64_t& idx) const;

    uint64_t sequence_size(void) const;
    const std::string& get_sequence(void) const;
    void set_sequence(const std::string& seq);
    const uint64_t& get_id(void) const;
    void set_id(const uint64_t& new_id);
    void set_volatile(void);
    void set_static(void);
    void for_each_edge(const std::function<bool(uint64_t other_id,
                                                bool on_rev,
                                                bool other_rev,
                                                bool to_curr)>& func) const;
    void add_edge(const uint64_t& other_id,
                  const bool& on_rev,
                  const bool& to_rev,
                  const bool& to_curr);
    bool remove_edge(const uint64_t& target_id,
                     const bool& target_rev,
                     const bool& ends_here,
                     const bool& is_rev);
    void add_path_step(const uint64_t& path_id, const bool& is_rev,
                       const bool& is_start, const bool& is_end,
                       const uint64_t& prev_id, const uint64_t& prev_rank,
                       const uint64_t& next_id, const uint64_t& next_rank);
    void add_path_step(const node_t::step_t& step);
    const step_t get_path_step(const uint64_t& rank) const;
    const std::vector<step_t> get_path_steps(void) const;
    void set_path_step(const uint64_t& rank, const uint64_t& path_id, const bool& is_rev,
                       const bool& is_start, const bool& is_end,
                       const uint64_t& prev_id, const uint64_t& prev_rank,
                       const uint64_t& next_id, const uint64_t& next_rank);
    void set_path_step(const uint64_t& rank, const step_t& step);
    void clear_path_step(const uint64_t& rank);
    void set_step_path_id(const uint64_t& rank, const uint64_t& path_id);
    void set_step_prev_id(const uint64_t& rank, const uint64_t& prev_id);
    void set_step_prev_rank(const uint64_t& rank, const uint64_t& prev_rank);
    void set_step_next_id(const uint64_t& rank, const uint64_t& next_id);
    void set_step_next_rank(const uint64_t& rank, const uint64_t& next_rank);
    void set_step_is_rev(const uint64_t& rank, const bool& is_rev);
    void set_step_is_start(const uint64_t& rank, const bool& is_start);
    void set_step_is_end(const uint64_t& rank, const bool& is_end);
    void set_step_is_del(const uint64_t& rank, const bool& is_del);
    uint64_t step_path_id(const uint64_t& rank) const;
    uint64_t step_prev_id(const uint64_t& rank) const;
    uint64_t step_prev_rank(const uint64_t& rank) const;
    uint64_t step_next_id(const uint64_t& rank) const;
    uint64_t step_next_rank(const uint64_t& rank) const;
    bool step_is_rev(const uint64_t& rank) const;
    bool step_is_start(const uint64_t& rank) const;
    bool step_is_end(const uint64_t& rank) const;
    bool step_is_del(const uint64_t& rank) const;
    void for_each_path_step(const std::function<bool(uint64_t rank,
                                                     uint64_t path_id,
                                                     bool is_rev)>& func) const;
    void for_each_path_step(const std::function<bool(step_t step)>& func) const;
    std::pair<std::map<uint64_t, std::pair<uint64_t, bool>>, // path fronts and backs
          std::map<uint64_t, std::pair<uint64_t, bool>>> flip_paths(void);
    void remove_path_step(const uint64_t& rank);
    void update_path_last_bytes(void);
    void clear(void);
    void clear_edges(void);
    void clear_paths(void);
    void clear_encoding(void);
    uint64_t serialize(std::ostream& out) const;
    void load(std::istream& in);
    void display(void) const;
    void copy(const node_t& other);
    void apply_ordering(
        const std::function<uint64_t(uint64_t)>& get_new_id,
        const std::function<bool(uint64_t)>& to_flip);
    void apply_path_ordering(
        const std::function<uint64_t(uint64_t)>& get_new_path_id);

};

}
