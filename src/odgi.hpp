//
//  odgi
//
//  graph.hpp
//
//  main dynamic compact graph definition
//

#pragma once

#include <cstdio>
#include <cstdint>
#include <vector>
#include <utility>
#include <functional>
#include <thread>
#include <handlegraph/types.hpp>
#include <handlegraph/iteratee.hpp>
#include <handlegraph/util.hpp>
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/path_handle_graph.hpp>
#include <handlegraph/mutable_handle_graph.hpp>
#include <handlegraph/mutable_path_handle_graph.hpp>
#include <handlegraph/mutable_path_mutable_handle_graph.hpp>
#include <handlegraph/deletable_handle_graph.hpp>
#include <handlegraph/mutable_path_deletable_handle_graph.hpp>
#include <handlegraph/serializable_handle_graph.hpp>
#include "dynamic.hpp"
#include "dynamic_types.hpp"
#include "lockfree_hashtable.hpp"
#include "dna.hpp"
#include "hash_map.hpp"
#include "node.hpp"

#include <omp.h>
#include "atomic_bitvector.hpp"
#include <mutex>

namespace odgi {

using namespace handlegraph;

// Resolve ambiguous nid_t typedef by putting it in our namespace.
using nid_t = handlegraph::nid_t;

class graph_t : public MutablePathDeletableHandleGraph, public SerializableHandleGraph, public RankedHandleGraph {

public:

    graph_t(void) {
        // set up initial delimiters
        path_metadata_h = std::make_unique<lockfree::LockFreeHashTable<uint64_t,
                                                                       path_metadata_t*>>();
        path_name_h = std::make_unique<lockfree::LockFreeHashTable<std::string,
                                                                   path_metadata_t*>>();
        _edge_count = 0;
        _path_count = 0;
        _path_handle_next = 0;
    }

    ~graph_t(void) { clear(); }

    /// Method to check if a node exists by ID
    bool has_node(nid_t node_id) const;

    /// If the handle has been deleted internally (these are removed in optimize())
    bool is_deleted(const handle_t& handle) const;

    /// Look up the handle for the node with the given ID in the given orientation
    handle_t get_handle(const nid_t& node_id, bool is_reverse = false) const;

    /// Get the ID from a handle
    nid_t get_id(const handle_t& handle) const;

    /// Get the orientation of a handle
    bool get_is_reverse(const handle_t& handle) const;

    /// Invert the orientation of a handle (potentially without getting its ID)
    handle_t flip(const handle_t& handle) const;

    /// Get the length of a node
    size_t get_length(const handle_t& handle) const;

    /// Get the sequence of a node, presented in the handle's local forward orientation.
    std::string get_sequence(const handle_t& handle) const;

protected:
    /// Loop over all the handles to next/previous (right/left) nodes. Passes
    /// them to a callback which returns false to stop iterating and true to
    /// continue. Returns true if we finished and false if we stopped early.
    bool follow_edges_impl(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const;

    /// Loop over all the nodes in the graph in their local forward
    /// orientations, in their internal stored order. Stop if the iteratee
    /// returns false. Can be told to run in parallel, in which case stopping
    /// after a false return value is on a best-effort basis and iteration
    /// order is not defined.
    bool for_each_handle_impl(const std::function<bool(const handle_t&)>& iteratee, bool parallel = false) const;

public:

    /// Return the number of nodes in the graph
    /// TODO: can't be node_count because XG has a field named node_count.
    size_t get_node_count(void) const;

    /// Return the smallest ID in the graph, or some smaller number if the
    /// smallest ID is unavailable. Return value is unspecified if the graph is empty.
    nid_t min_node_id(void) const;

    /// Return the largest ID in the graph, or some larger number if the
    /// largest ID is unavailable. Return value is unspecified if the graph is empty.
    nid_t max_node_id(void) const;

    /// Set a minimum id to increment the id space by, used as a hint during construction.
    /// May have no effect on a backing implementation.
    void set_id_increment(const nid_t& min_id);

    /// Increment node ids, using the builtin id increment, assumes we're increasing by a positive value
    void increment_node_ids(nid_t increment);

    ////////////////////////////////////////////////////////////////////////////
    // Additional optional interface with a default implementation
    ////////////////////////////////////////////////////////////////////////////

    /// Get the number of edges on the right (go_left = false) or left (go_left
    /// = true) side of the given handle. The default implementation is O(n) in
    /// the number of edges returned, but graph implementations that track this
    /// information more efficiently can override this method.
    size_t get_degree(const handle_t& handle, bool go_left) const;

    /// Get the locally forward version of a handle
    handle_t forward(const handle_t& handle) const;

    /// A pair of handles can be used as an edge. When so used, the handles have a
    /// canonical order and orientation.
    edge_t edge_handle(const handle_t& left, const handle_t& right) const;

    /// Such a pair can be viewed from either inward end handle and produce the
    /// outward handle you would arrive at.
    handle_t traverse_edge_handle(const edge_t& edge, const handle_t& left) const;

    ////////////////////////////////////////////////////////////////////////////
    // Rank handle interface
    ////////////////////////////////////////////////////////////////////////////
    /// Return the rank of a node (ranks start at 1 and are dense).
    size_t id_to_rank(const nid_t& node_id) const;

    /// Return the node with a given rank.
    nid_t rank_to_id(const size_t& rank) const;

    // If you define node ID ranks you get a default implementation of handle ranks.

    /// Return the rank of a handle (ranks start at 1 and are dense, and each
    /// orientation has its own rank). Handle ranks may not have anything to do
    /// with node ranks.
    size_t handle_to_rank(const handle_t& handle) const;

    /// Return the handle with a given rank.
    handle_t rank_to_handle(const size_t& rank) const;


    ////////////////////////////////////////////////////////////////////////////
    // Path handle interface
    ////////////////////////////////////////////////////////////////////////////

    /// Determine if a path name exists and is legal to get a path handle for.
    bool has_path(const std::string& path_name) const;

    /// Look up the path handle for the given path name.
    /// The path with that name must exist.
    path_handle_t get_path_handle(const std::string& path_name) const;

    /// Look up the name of a path from a handle to it
    std::string get_path_name(const path_handle_t& path_handle) const;

    /// Returns the number of node steps in the path
    size_t get_step_count(const path_handle_t& path_handle) const;

    /// Returns the number of paths stored in the graph
    size_t get_path_count(void) const;

    /// Returns a vector of all steps of a node on paths. Optionally restricts to
    /// steps that match the handle in orientation.
    std::vector<step_handle_t> steps_of_handle(const handle_t& handle,
                                               bool match_orientation = false) const;

protected:

    /// Execute a function on each path in the graph
    bool for_each_path_handle_impl(const std::function<bool(const path_handle_t&)>& iteratee) const;

    /// Enumerate the path steps on a given handle (strand agnostic)
    bool for_each_step_on_handle_impl(const handle_t& handle, const std::function<bool(const step_handle_t&)>& iteratee) const;

public:

    /// Returns the number of node steps on the handle
    size_t get_step_count(const handle_t& handle) const;

    /// Get a node handle (node ID and orientation) from a handle to an step on a path
    handle_t get_handle_of_step(const step_handle_t& step_handle) const;

    /// Get a path handle (path ID) from a handle to an step on a path
    path_handle_t get_path(const step_handle_t& step_handle) const;

    /// Get a handle to the first step in a path.
    /// The path MUST be nonempty.
    step_handle_t path_begin(const path_handle_t& path_handle) const;

    /// Get a handle to a fictitious handle one past the end of the path
    step_handle_t path_end(const path_handle_t& path_handle) const;

    /// Get a handle to the last step, which is arbitrary in the case of a circular path
    step_handle_t path_back(const path_handle_t& path_handle) const;

    /// Get a handle to a fictitious handle one past the start of the path
    step_handle_t path_front_end(const path_handle_t& path_handle) const;

    /// Returns true if the step handle is a front end magic handle
    bool is_path_front_end(const step_handle_t& step_handle) const;

    /// Returns true if the step handle is an end magic handle
    bool is_path_end(const step_handle_t& step_handle) const;

    /// Returns true if the step is not the last step on the path, else false
    bool has_next_step(const step_handle_t& step_handle) const;

    /// Returns true if the step is not the first step on the path, else false
    bool has_previous_step(const step_handle_t& step_handle) const;

    /// Returns a handle to the next step on the path
    step_handle_t get_next_step(const step_handle_t& step_handle) const;

    /// Returns a handle to the previous step on the path
    step_handle_t get_previous_step(const step_handle_t& step_handle) const;

    /// Returns a handle to the path that an step is on
    path_handle_t get_path_handle_of_step(const step_handle_t& step_handle) const;

    /// Returns the 0-based ordinal rank of a step on a path
    size_t get_ordinal_rank_of_step(const step_handle_t& step_handle) const;

    /// Returns true if the given path is empty, and false otherwise
    bool is_empty(const path_handle_t& path_handle) const;

    /// Loop over all the steps along a path, from first through last
    void for_each_step_in_path(const path_handle_t& path, const std::function<void(const step_handle_t&)>& iteratee) const;

    /// Returns true if the path is circular
    bool get_is_circular(const path_handle_t& path_handle) const;

    /// Set if the path is circular or not
    void set_circularity(const path_handle_t& path_handle, bool circular);

    /// Create a new node with the given sequence and return the handle.
    handle_t create_handle(const std::string& sequence);

    /// Create a new node with the given id and sequence, then return the handle.
    handle_t create_handle(const std::string& sequence, const nid_t& id);

    /// Remove the node belonging to the given handle and all of its edges.
    /// Does not update any stored paths.
    /// Invalidates the destroyed handle.
    /// May be called during serial for_each_handle iteration **ONLY** on the node being iterated.
    /// May **NOT** be called during parallel for_each_handle iteration.
    /// May **NOT** be called on the node from which edges are being followed during follow_edges.
    void destroy_handle(const handle_t& handle);

    /// Create an edge connecting the given handles in the given order and orientations.
    /// Ignores existing edges.
    void create_edge(const handle_t& left, const handle_t& right);

    /// Check if an edge exists
    bool has_edge(const handle_t& left, const handle_t& right) const;

    /// Convenient wrapper for create_edge.
    inline void create_edge(const edge_t& edge) {
        create_edge(edge.first, edge.second);
    }

    /// Remove the edge connecting the given handles in the given order and orientations.
    /// Ignores nonexistent edges.
    /// Does not update any stored paths.
    void destroy_edge(const handle_t& left, const handle_t& right);

    /// Convenient wrapper for destroy_edge.
    inline void destroy_edge(const edge_t& edge) {
        destroy_edge(edge.first, edge.second);
    }

    /// Remove all nodes and edges. Does not update any stored paths.
    void clear(void);

    /// Remove all stored paths
    void clear_paths(void);

    /// Swap the nodes corresponding to the given handles, in the ordering used
    /// by for_each_handle when looping over the graph. Other handles to the
    /// nodes being swapped must not be invalidated. If a swap is made while
    /// for_each_handle is running, it affects the order of the handles
    /// traversed during the current traversal (so swapping an already seen
    /// handle to a later handle's position will make the seen handle be visited
    /// again and the later handle not be visited at all).
    void swap_handles(const handle_t& a, const handle_t& b);

    /// Reorder the graph's internal structure to match that given.
    /// Optionally compact the id space of the graph to match the ordering, from 1->|ordering|.
    bool apply_ordering(const std::vector<handle_t>& order, bool compact_ids = false);

    /// Organize the graph for better performance and memory use
    void optimize(bool allow_id_reassignment = true);

    /// is the graph optimized?
    /// definition:
    /// smallest node identifier is 1 and largest node identifier is equal to get_node_count()
    bool is_optimized(void);

    /// Reassign the node ids
    void reassign_node_ids(const std::function<nid_t(const nid_t&)>& get_new_id);

    /// Reorder the graph's paths as given.
    void apply_path_ordering(const std::vector<path_handle_t>& order);

    /// Alter the node that the given handle corresponds to so the orientation
    /// indicated by the handle becomes the node's local forward orientation.
    /// Rewrites all edges pointing to the node and the node's sequence to
    /// reflect this. Invalidates all handles to the node (including the one
    /// passed). Returns a new, valid handle to the node in its new forward
    /// orientation. Note that it is possible for the node's ID to change.
    /// Does not update any stored paths. May change the ordering of the underlying
    /// graph.
    handle_t apply_orientation(const handle_t& handle);

    /// Split a handle's underlying node at the given offsets in the handle's
    /// orientation. Returns all of the handles to the parts. Other handles to
    /// the node being split may be invalidated. The split pieces stay in the
    /// same local forward orientation as the original node, but the returned
    /// handles come in the order and orientation appropriate for the handle
    /// passed in.
    /// Updates stored paths.
    std::vector<handle_t> divide_handle(const handle_t& handle, const std::vector<size_t>& offsets);

    /// Specialization of divide_handle for a single division point
    inline std::pair<handle_t, handle_t> divide_handle(const handle_t& handle, size_t offset) {
        auto parts = divide_handle(handle, std::vector<size_t>{offset});
        return std::make_pair(parts.front(), parts.back());
    }

    handle_t combine_handles(const std::vector<handle_t>& handles);

/**
 * This is the interface for a handle graph with embedded paths where the paths can be modified.
 * Note that if the *graph* can also be modified, the implementation will also
 * need to inherit from MutableHandleGraph, via the combination
 * MutablePathMutableHandleGraph interface.
 * TODO: This is a very limited interface at the moment. It will probably need to be extended.
 */

    /**
     * Destroy the given path. Invalidates handles to the path and its node steps.
     */
    void destroy_path(const path_handle_t& path);

    /**
     * Create a path with the given name. The caller must ensure that no path
     * with the given name exists already, or the behavior is undefined.
     * Returns a handle to the created empty path. Handles to other paths must
     * remain valid.
     */
    path_handle_t create_path_handle(const std::string& name,
                                     bool is_circular = false);

    /**
     * Append a visit to a node to the given path. Returns a handle to the new
     * final step on the path which is appended. Handles to prior
     * steps on the path, and to other paths, must remain valid.
     */
    step_handle_t prepend_step(const path_handle_t& path, const handle_t& to_append);

    /**
     * Append a visit to a node to the given path. Returns a handle to the new
     * final step on the path which is appended. Handles to prior
     * steps on the path, and to other paths, must remain valid.
     */
    step_handle_t append_step(const path_handle_t& path, const handle_t& to_append);

    /**
     * Insert a visit to a node to the given path between the given steps.
     * Returns a handle to the new step on the path which is appended.
     * Handles to prior steps on the path, and to other paths, must remain valid.
     */
    step_handle_t insert_step(const step_handle_t& before, const step_handle_t& after, const handle_t& to_insert);

    /// Set the step to the given handle, possibly re-linking and cleaning up if needed
    step_handle_t set_step(const step_handle_t& step_handle, const handle_t& handle);

    /// Replace the path range with the new segment
    std::pair<step_handle_t, step_handle_t> rewrite_segment(const step_handle_t& segment_begin,
                                                            const step_handle_t& segment_end,
                                                            const std::vector<handle_t>& new_segment);

    /// A helper function to visualize the state of the graph
    void display(void) const;

    /// Convert to GFA (for debugging)
    void to_gfa(std::ostream& out, const bool& emit_node_annotation = false) const;

    /// Magic number header for serialization
    uint32_t get_magic_number(void) const;

    /// Serialize
    void serialize_members(std::ostream& out) const;

    /// Load
    void deserialize_members(std::istream& in);

    void set_number_of_threads(uint64_t num_threads);

    uint64_t get_number_of_threads();

    /// copy the other graph into this one
    void copy(const graph_t& other);

/// These are the backing data structures that we use to fulfill the above functions

    /// Records the handle to node_id mapping
    /// Use the special value "0" to indicate deleted nodes so that
    /// handle references in the id_map and elsewhere are not immediately destroyed

    // lock for the node vector and the following variables
    // TODO use it in create_handle and friends
    std::atomic_flag node_lock = ATOMIC_FLAG_INIT;
    std::vector<node_t*> node_v; // not threadsafe
    node_t& get_node_ref(const handle_t& handle) const;
    const node_t& get_node_cref(const handle_t& handle) const;
    /// Mark deleted nodes here for translating graph ids into internal ranks
    //dyn::hacked_vector deleted_nodes;
    hash_set<uint64_t> deleted_nodes;
    /// efficient id to handle/sequence conversion
    std::atomic<nid_t> _max_node_id = 0;
    std::atomic<nid_t> _min_node_id = 0;
    std::atomic<nid_t> _id_increment = 0;
    uint64_t _num_threads = 1;

    inline void canonicalize_edge(handle_t& left, handle_t& right) const {
        if (number_bool_packing::unpack_bit(left) && number_bool_packing::unpack_bit(right)
            || ((number_bool_packing::unpack_bit(left) || number_bool_packing::unpack_bit(right)) && as_integer(left) > as_integer(right))) {
            std::swap(left, right);
            left = number_bool_packing::toggle_bit(left);
            right = number_bool_packing::toggle_bit(right);
        }
    }

    struct path_metadata_t {
        std::atomic<path_handle_t> handle = as_path_handle(0);
        std::atomic<uint64_t> length;
        std::atomic<step_handle_t> first;
        std::atomic<step_handle_t> last;
        std::string name;
        std::atomic<bool> is_circular;
        std::atomic_flag lock = ATOMIC_FLAG_INIT;
        inline void get_lock(void) {
            while (lock.test_and_set(std::memory_order_acquire))  // acquire lock
                ; // spin
        }
        inline void clear_lock(void) {
            lock.clear(std::memory_order_release);
        }
        void copy(const path_metadata_t& other) {
            handle.store(other.handle);
            length.store(other.length);
            first.store(other.first);
            last.store(other.last);
            name = other.name;
            is_circular.store(other.is_circular);
        }
    };

    /// maps between path identifier and the start, end, and length of the path
    std::unique_ptr<lockfree::LockFreeHashTable<uint64_t, path_metadata_t*>> path_metadata_h;
    std::unique_ptr<lockfree::LockFreeHashTable<std::string, path_metadata_t*>> path_name_h;
    path_metadata_t& get_path_metadata(const path_handle_t& path) const;
    const path_metadata_t& path_metadata(const path_handle_t& path) const;

    /// A helper to record the number of live nodes
    //std::atomic<uint64_t> _node_count; // = 0;

    /// A helper to record the number of live edges
    std::atomic<uint64_t> _edge_count;// = 0;

    /// A helper to record the number of live paths
    std::atomic<uint64_t> _path_count; // = 0;

    /// A helper to record the next path handle id
    std::atomic<uint64_t> _path_handle_next; // = 0;

    /// A helper to record the number of live paths
    //std::atomic<uint64_t> _step_count; // = 0; // TODO

    /// Helper to simplify removal of path handle records
    void destroy_path_handle_records(uint64_t i);

    /// Helper to create the internal records for the step
    step_handle_t create_step(const path_handle_t& path, const handle_t& handle);

    /// Helper to destroy the internal records for the step
    void destroy_step(const step_handle_t& step_handle);

    /// Helper to stitch up partially built paths
    void link_steps(const step_handle_t& from, const step_handle_t& to);

    /// Decrement the step rank references for this step
    void decrement_rank(const step_handle_t& step_handle);

    /// Compact away the deleted nodes info
    //void rebuild_id_handle_mapping(void);

    /// Set the handle sequence
    void set_handle_sequence(const handle_t& handle, const std::string& seq);

    /// get the backing node rank for a given node id
    uint64_t get_node_rank(const nid_t& node_id) const;

    /// set the graph into static mode, which avoids spinlocks
    void set_static(void);

    /// set the graph into volatile mode, which uses (spin)locks for dynamic hogwilding
    void set_volatile(void);

};

//const static uint64_t path_begin_marker = std::numeric_limits<uint64_t>::max();
//const static uint64_t path_end_marker = 2;

} // end dankness
