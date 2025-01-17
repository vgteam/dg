//
//  graph.cpp
//

#include "odgi.hpp"

namespace odgi {

node_t& graph_t::get_node_ref(const handle_t& handle) const {
    return *node_v[number_bool_packing::unpack_number(handle)];
}

const node_t& graph_t::get_node_cref(const handle_t& handle) const {
    return *node_v[number_bool_packing::unpack_number(handle)];
}

/// Method to check if a node exists by ID
bool graph_t::has_node(nid_t node_id) const {
    uint64_t rank = get_node_rank(node_id);
    return (rank >= node_v.size() ? false : node_v[rank] != nullptr);
}

/// If the handle has been deleted internally (these are removed in optimize())
bool graph_t::is_deleted(const handle_t& handle) const {
    return node_v[number_bool_packing::unpack_number(handle)] == nullptr;
}

/// Look up the handle for the node with the given ID in the given orientation
handle_t graph_t::get_handle(const nid_t& node_id, bool is_reverse) const {
    return number_bool_packing::pack(get_node_rank(node_id), is_reverse);
}

/// Get the ID from a handle
nid_t graph_t::get_id(const handle_t& handle) const {
    return number_bool_packing::unpack_number(handle) + 1 + _id_increment;
}

/// get the backing node for a given node id
uint64_t graph_t::get_node_rank(const nid_t& node_id) const {
    return node_id - _id_increment - 1;
}

/// set the id increment, used when the graph starts at a high id to reduce loading costs
void graph_t::set_id_increment(const nid_t& min_id) {
    _id_increment = min_id;
}

/// Increment node ids, using the builtin id increment, assumes we're increasing by a positive value
void graph_t::increment_node_ids(nid_t increment) {
    _id_increment += increment;
}

/// Get the orientation of a handle
bool graph_t::get_is_reverse(const handle_t& handle) const {
    return number_bool_packing::unpack_bit(handle);
}

/// Invert the orientation of a handle (potentially without getting its ID)
handle_t graph_t::flip(const handle_t& handle) const {
    return number_bool_packing::toggle_bit(handle);
}

/// Get the length of a node
size_t graph_t::get_length(const handle_t& handle) const {
    auto& node = get_node_ref(handle);
    node.get_lock();
    auto l = node.sequence_size();
    node.clear_lock();
    return l;
}

/// Get the sequence of a node, presented in the handle's local forward orientation.
std::string graph_t::get_sequence(const handle_t& handle) const {
    auto& node = get_node_ref(handle);
    node.get_lock();
    auto seq = node.get_sequence();
    node.clear_lock();
    return (get_is_reverse(handle) ? reverse_complement(seq) : seq);
}

/// Loop over all the handles to next/previous (right/left) nodes. Passes
/// them to a callback which returns false to stop iterating and true to
/// continue. Returns true if we finished and false if we stopped early.
bool graph_t::follow_edges_impl(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const {
    const node_t& node = get_node_cref(handle);
    nid_t node_id = get_id(handle);
    bool is_rev = get_is_reverse(handle);
    node.for_each_edge(
        [&](nid_t other_id,
            bool other_rev,
            bool to_curr,
            bool on_rev) {
            if (other_id == node_id && on_rev == other_rev) {
                // non-inverting self loop
                // we can go either direction
                to_curr = go_left;
                other_rev = is_rev;
            } else if (is_rev != on_rev) {
                other_rev ^= 1;
                to_curr ^= 1;
            }
            if (!go_left && !to_curr
                || go_left && to_curr) {
                if (!iteratee(get_handle(other_id, other_rev))) {
                    return false;
                }
            }
            return true;
        });
    return true;
}

/// Loop over all the nodes in the graph in their local forward
/// orientations, in their internal stored order. Stop if the iteratee
/// returns false. Can be told to run in parallel, in which case stopping
/// after a false return value is on a best-effort basis and iteration
/// order is not defined.
bool graph_t::for_each_handle_impl(const std::function<bool(const handle_t&)>& iteratee, bool parallel) const {
    if (parallel) {
        volatile bool flag=true;
#pragma omp parallel for
        for (uint64_t i = 0; i < node_v.size(); ++i) {
            handle_t h = number_bool_packing::pack(i,false);
            if (is_deleted(h)) continue;
            if (!flag) continue;
            bool result = iteratee(h);
#pragma omp atomic
            flag &= result;
        }
        return flag;
    } else {
        for (uint64_t i = 0; i < node_v.size(); ++i) {
            handle_t h = number_bool_packing::pack(i,false);
            if (is_deleted(h)) continue;
            if (!iteratee(h)) return false;
        }
        return true;
    }
}

/// Return the number of nodes in the graph
/// TODO: can't be node_count because XG has a field named node_count.
size_t graph_t::get_node_count() const {
    return node_v.size()-deleted_nodes.size();
}

/// Return the smallest ID in the graph, or some smaller number if the
/// smallest ID is unavailable. Return value is unspecified if the graph is empty.
nid_t graph_t::min_node_id() const {
    return _min_node_id + _id_increment;
}

/// Return the largest ID in the graph, or some larger number if the
/// largest ID is unavailable. Return value is unspecified if the graph is empty.
nid_t graph_t::max_node_id() const {
    return _max_node_id + _id_increment;
}

////////////////////////////////////////////////////////////////////////////
// Additional optional interface with a default implementation
////////////////////////////////////////////////////////////////////////////

/// Get the number of edges on the right (go_left = false) or left (go_left
/// = true) side of the given handle. The default implementation is O(n) in
/// the number of edges returned, but graph implementations that track this
/// information more efficiently can override this method.
size_t graph_t::get_degree(const handle_t& handle, bool go_left) const {
    size_t degree = 0;
    follow_edges(handle, go_left, [&degree](const handle_t& h) { ++degree; });
    return degree;
}

/// Get the locally forward version of a handle
handle_t graph_t::forward(const handle_t& handle) const {
    if (get_is_reverse(handle)) {
        return flip(handle);
    } else {
        return handle;
    }
}

/// A pair of handles can be used as an edge. When so used, the handles have a
/// canonical order and orientation.
edge_t graph_t::edge_handle(const handle_t& left, const handle_t& right) const {
    return std::make_pair(left, right);
}

/**
 * This is the interface for a handle graph that stores embedded paths.
 */

////////////////////////////////////////////////////////////////////////////
// Path handle interface that needs to be implemented
////////////////////////////////////////////////////////////////////////////

graph_t::path_metadata_t& graph_t::get_path_metadata(const path_handle_t& path) const {
    graph_t::path_metadata_t* p;
    if (path_metadata_h->Find(as_integer(path), p)) {
        return *p;
    } else {
        assert(false);
        return *p; // won't reach unless assert is disabled
    }
}

const graph_t::path_metadata_t& graph_t::path_metadata(const path_handle_t& path) const {
    graph_t::path_metadata_t* p;
    if (path_metadata_h->Find(as_integer(path), p)) {
        return *p;
    } else {
        assert(false);
        return *p; // won't reach unless assert is disabled
    }
}

/// Determine if a path name exists and is legal to get a path handle for.
bool graph_t::has_path(const std::string& path_name) const {
    graph_t::path_metadata_t* p;
    return path_name_h->Find(path_name, p);
}

/// Look up the path handle for the given path name.
/// The path with that name must exist.
path_handle_t graph_t::get_path_handle(const std::string& path_name) const {
    graph_t::path_metadata_t* p;
    if (path_name_h->Find(path_name, p)) {
        return p->handle;
    } else {
        assert(false);
        return as_path_handle(0); // won't reach unless assert is disabled
    }
}

/// Look up the name of a path from a handle to it
std::string graph_t::get_path_name(const path_handle_t& path_handle) const {
    return path_metadata(path_handle).name;
}

/// Returns the number of node steps in the path
size_t graph_t::get_step_count(const path_handle_t& path_handle) const {
    return path_metadata(path_handle).length;
}

/// Returns the number of paths stored in the graph
size_t graph_t::get_path_count() const {
    return _path_count;
}

/// Execute a function on each path in the graph
bool graph_t::for_each_path_handle_impl(const std::function<bool(const path_handle_t&)>& iteratee) const {
    bool flag = true;
    for (uint64_t i = 1; i <= _path_handle_next; ++i) {
        path_metadata_t* p;
        if (path_metadata_h->Find(i, p)) {
            flag &= iteratee(as_path_handle(i));
        }
    }
    return flag;
}

bool graph_t::for_each_step_on_handle_impl(const handle_t& handle, const std::function<bool(const step_handle_t&)>& iteratee) const {
    uint64_t handle_n = number_bool_packing::unpack_number(handle);
    const node_t& node = get_node_cref(handle);
    bool flag = true;
    node.for_each_path_step(
        [&](uint64_t rank,
            uint64_t path_id,
            bool is_rev) {
            step_handle_t step_handle;
            as_integers(step_handle)[0] = as_integer(number_bool_packing::pack(handle_n, is_rev));
            as_integers(step_handle)[1] = rank;
            flag &= iteratee(step_handle);
            return flag;
        });
    return flag;
}

/// Returns a vector of all steps of a node on paths. Optionally restricts to
/// steps that match the handle in orientation.
std::vector<step_handle_t> graph_t::steps_of_handle(const handle_t& handle,
                                                    bool match_orientation) const {
    std::vector<step_handle_t> res;
    for_each_step_on_handle(handle, [&](const step_handle_t& step) {
            handle_t h = get_handle_of_step(step);
            if (!match_orientation || get_is_reverse(h) == get_is_reverse(handle)) {
                res.push_back(step);
            }
        });
    return res;
}

size_t graph_t::get_step_count(const handle_t& handle) const {
    auto& node = get_node_ref(handle);
    node.get_lock();
    auto count = node.path_count();
    node.clear_lock();
    return count;
}

/// Get a node handle (node ID and orientation) from a handle to an step on a path
handle_t graph_t::get_handle_of_step(const step_handle_t& step_handle) const {
    return as_handle(as_integers(step_handle)[0]);
}

/// Get a path handle (path ID) from a handle to an step on a path
path_handle_t graph_t::get_path(const step_handle_t& step_handle) const {
    auto& node = get_node_ref(get_handle_of_step(step_handle));
    node.get_lock();
    auto p = node.get_path_step(as_integers(step_handle)[1]).path_id;
    node.clear_lock();
    return as_path_handle(p);
}

/// Get a handle to the first step in a path.
/// The path MUST be nonempty.
step_handle_t graph_t::path_begin(const path_handle_t& path_handle) const {
    return get_path_metadata(path_handle).first.load();
}

/// Get a handle to the last step in a path
/// The path MUST be nonempty.
step_handle_t graph_t::path_back(const path_handle_t& path_handle) const {
    return get_path_metadata(path_handle).last.load();
}

/// Get the reverse end iterator
step_handle_t graph_t::path_front_end(const path_handle_t& path_handle) const {
    step_handle_t step;
    as_integers(step)[0] = as_integer(path_handle);
    as_integers(step)[1] = std::numeric_limits<uint64_t>::max()-1;
    return step;
}

/// Get the forward end iterator
step_handle_t graph_t::path_end(const path_handle_t& path_handle) const {
    step_handle_t step;
    as_integers(step)[0] = as_integer(path_handle);
    as_integers(step)[1] = std::numeric_limits<uint64_t>::max();
    return step;
}

/// is the path circular
bool graph_t::get_is_circular(const path_handle_t& path_handle) const {
    return path_metadata(path_handle).is_circular;
}

/// set the circular flag for the path
void graph_t::set_circularity(const path_handle_t& path_handle, bool circular) {
    get_path_metadata(path_handle).is_circular = circular;
}

/// Returns true if the step is not the last step on the path, else false
bool graph_t::has_next_step(const step_handle_t& step_handle) const {
    if (is_path_end(step_handle)) return false;
    auto& node = get_node_ref(get_handle_of_step(step_handle));
    node.get_lock();
    auto b = !node.step_is_end(as_integers(step_handle)[1]);
    node.clear_lock();
    return b;
}

/// Returns true if the step is not the first step on the path, else false
bool graph_t::has_previous_step(const step_handle_t& step_handle) const {
    if (is_path_end(step_handle)) return true;
    auto& node = get_node_ref(get_handle_of_step(step_handle));
    node.get_lock();
    auto b = !node.step_is_start(as_integers(step_handle)[1]);
    node.clear_lock();
    return b;
}

bool graph_t::is_path_front_end(const step_handle_t& step_handle) const {
    return as_integers(step_handle)[1] == std::numeric_limits<uint64_t>::max()-1;
}

bool graph_t::is_path_end(const step_handle_t& step_handle) const {
    return as_integers(step_handle)[1] == std::numeric_limits<uint64_t>::max();
}

/// Returns a handle to the next step on the path
/// Returns the forward end iterator if none exists
step_handle_t graph_t::get_next_step(const step_handle_t& step_handle) const {
    handle_t curr_handle;
    // check if we have a magic path iterator step handle
    if (is_path_front_end(step_handle)) {
        curr_handle = get_handle_of_step(path_begin(as_path_handle(as_integers(step_handle)[0])));
    } else if (is_path_end(step_handle)) {
        return step_handle;
    } else {
        curr_handle = get_handle_of_step(step_handle);
    }
    nid_t curr_id = get_id(curr_handle);
    node_t& node = get_node_ref(curr_handle);
    node.get_lock();
    auto step_rank = as_integers(step_handle)[1];
    if (node.step_is_end(step_rank)) {
        node.clear_lock();
        return path_end(get_path_handle_of_step(step_handle));
    }
    nid_t next_id = node.step_next_id(step_rank);
    auto next_rank = node.step_next_rank(step_rank);
    node.clear_lock();
    handle_t next_handle = get_handle(next_id);
    node_t& next = get_node_ref(next_handle);
    next.get_lock();
    bool next_rev = next.step_is_rev(next_rank);
    next.clear_lock();
    step_handle_t next_step;
    as_integers(next_step)[0] = as_integer(get_handle(next_id, next_rev));
    as_integers(next_step)[1] = next_rank;
    return next_step;
}

/// Returns a handle to the previous step on the path
step_handle_t graph_t::get_previous_step(const step_handle_t& step_handle) const {
    handle_t curr_handle;
    // check if we have a magic path iterator step handle
    if (is_path_front_end(step_handle)) {
        return step_handle;
    } else if (is_path_end(step_handle)) {
        return path_back(as_path_handle(as_integers(step_handle)[0]));
		// curr_handle = get_handle_of_step(path_back(as_path_handle(as_integers(step_handle)[0])));
    } else {
        curr_handle = get_handle_of_step(step_handle);
    }
    nid_t curr_id = get_id(curr_handle);
    node_t& node = get_node_ref(curr_handle);
    node.get_lock();
    auto step_rank = as_integers(step_handle)[1];
    if (node.step_is_start(step_rank)) {
        node.clear_lock();
        return path_front_end(get_path_handle_of_step(step_handle));
    }
    nid_t prev_id = node.step_prev_id(step_rank);
    auto prev_rank = node.step_prev_rank(step_rank);
    node.clear_lock();
    handle_t prev_handle = get_handle(prev_id);
    node_t& prev = get_node_ref(prev_handle);
    prev.get_lock();
    bool prev_rev = prev.step_is_rev(prev_rank);
    prev.clear_lock();
    step_handle_t prev_step;
    as_integers(prev_step)[0] = as_integer(get_handle(prev_id, prev_rev));
    as_integers(prev_step)[1] = prev_rank;
    return prev_step;
}

path_handle_t graph_t::get_path_handle_of_step(const step_handle_t& step_handle) const {
    node_t& node = get_node_ref(get_handle_of_step(step_handle));
    node.get_lock();
    auto path = as_path_handle(node.step_path_id(as_integers(step_handle)[1]));
    node.clear_lock();
    return path;
}

////////////////////////////////////////////////////////////////////////////
// Additional optional interface with a default implementation
////////////////////////////////////////////////////////////////////////////

/// Returns the 0-based ordinal rank of a step on a path
size_t graph_t::get_ordinal_rank_of_step(const step_handle_t& step_handle) const {
    return 0; // no implementation of this in odgi (yet)
}

/// Returns true if the given path is empty, and false otherwise
bool graph_t::is_empty(const path_handle_t& path_handle) const {
    return get_step_count(path_handle) == 0;
}

/**
 * This is the interface for a handle graph that supports modification.
 */
/*
 * Note: All operations may invalidate path handles and step handles.
 */

/// Loop over all the steps along a path, from first through last
void graph_t::for_each_step_in_path(const path_handle_t& path, const std::function<void(const step_handle_t&)>& iteratee) const {
    auto& p = path_metadata(path);
    if (is_empty(path)) return;
    step_handle_t step = path_begin(path);
    step_handle_t end_step = path_back(path);
    bool keep_going = true;
    do {
        iteratee(step);
        // in circular paths, we'll always have a next step, so we always check if we're at our path's last step
        if (step != end_step && has_next_step(step)) {
            step = get_next_step(step);
        } else {
            keep_going = false;
        }
    } while (keep_going);
}

/// Create a new node with the given sequence and return the handle.
handle_t graph_t::create_handle(const std::string& sequence) {
    // get first deleted node to recycle
    if (deleted_nodes.size()) {
        uint64_t id = *deleted_nodes.begin();
        return create_handle(sequence, id);
    } else {
        return create_handle(sequence, node_v.size()+1);
    }
}

/// Create a new node with the given id and sequence, then return the handle.
handle_t graph_t::create_handle(const std::string& sequence, const nid_t& id) {
    assert(sequence.size());
    assert(id > 0);
    assert(!has_node(id));

    if (id > node_v.size()) {
        uint64_t old_size = node_v.size();
        // realloc
        node_v.resize((uint64_t)id, nullptr);
        // mark empty nodes
        // we'll use the last one for our new node
        for (uint64_t i = old_size+1; i <= id; ++i) {
            deleted_nodes.insert(i);
        }
    }
    // update min/max node ids
    _max_node_id = std::max(id, _max_node_id.load());
    if (_min_node_id) {
        _min_node_id = (uint64_t)min(id, _min_node_id.load());
    } else {
        _min_node_id = id;
    }
    // add to node vector
    uint64_t handle_rank = (uint64_t)id-1;
    // set its values
    auto& n = node_v[handle_rank];
    if (n == nullptr) {
        assert(deleted_nodes.count(id));
        deleted_nodes.erase(id);
    }
    n = new node_t();
    auto& node = *n;
    node.set_id(id);
    node.set_sequence(sequence);
    return number_bool_packing::pack(handle_rank, 0);
}

/// Remove the node belonging to the given handle and all of its edges.
/// Does not update any stored paths.
/// Invalidates the destroyed handle.
/// May be called during serial for_each_handle iteration **ONLY** on the node being iterated.
/// May **NOT** be called during parallel for_each_handle iteration.
/// May **NOT** be called on the node from which edges are being followed during follow_edges.
void graph_t::destroy_handle(const handle_t& handle) {
    handle_t fwd_handle = get_is_reverse(handle) ? flip(handle) : handle;
    uint64_t id = get_id(handle);
    if (!has_node(id)) return; // deleted already
    // remove steps in edge lists
    // enumerate the edges
    std::vector<edge_t> edges_to_destroy;
    follow_edges(fwd_handle, false,
                 [&edges_to_destroy,&fwd_handle,this](const handle_t& h) {
                     edges_to_destroy.push_back(make_pair(fwd_handle, h));
                 });
    follow_edges(fwd_handle, true,
                 [&edges_to_destroy,&fwd_handle,this](const handle_t& h) {
                     edges_to_destroy.push_back(make_pair(h, fwd_handle));
                 });
    // and then remove them
    for (auto& edge : edges_to_destroy) {
        destroy_edge(edge);
    }
    // clear the node storage
    auto& node = node_v[number_bool_packing::unpack_number(handle)];
    delete node;
    // remove from the graph
    node = nullptr;
    // add the index to our list of open node slots
    deleted_nodes.insert(id);
}

/*
void graph_t::rebuild_id_handle_mapping() {
    // for each live node, record the id in a new vector
    uint64_t j = 0;
    for (uint64_t i = 0; i < node_v.size(); ++i) {
        auto id = node_v[i].id();
        if (id == 0) continue;
        graph_id_map[id] = j++;
    }
    for (uint64_t i = 0; i < node_v.size(); ++i) {
        auto id = node_v[i].id();
        if (id == 0) {
            deleted_node_bv.remove(i);
        }
    }
    node_v.erase(std::remove_if(node_v.begin(), node_v.end(),
                                [&](const node_t& node) { return node.id() == 0; }),
                 node_v.end());
    _deleted_node_count = 0;
}
*/

/// Create an edge connecting the given handles in the given order and orientations.
/// Ignores existing edges.
void graph_t::create_edge(const handle_t& left_h, const handle_t& right_h) {
    uint64_t left_rank = number_bool_packing::unpack_number(left_h);
    uint64_t right_rank = number_bool_packing::unpack_number(right_h);
    bool create_edge = false;
    auto& left_node = get_node_ref(left_h);
    auto& right_node = get_node_ref(right_h);
    // ordered locking
    if (left_rank < right_rank) {
        left_node.get_lock();
        right_node.get_lock();
    } else if (left_rank != right_rank) {
        right_node.get_lock();
        left_node.get_lock();
    } else {
        left_node.get_lock();
    }
    if (!has_edge(left_h, right_h)) { // do nothing if edge exists
        create_edge = true;
    }
    // insert the edge for each side
    if (create_edge) {
        ++_edge_count;
        left_node.add_edge(get_id(right_h),
                           get_is_reverse(right_h),
                           false,
                           get_is_reverse(left_h));
        left_node.clear_lock();
        // only insert the second side if it's on a different node
        if (left_rank != right_rank) {
            right_node.add_edge(get_id(left_h),
                                get_is_reverse(left_h),
                                true,
                                get_is_reverse(right_h));
        }
    }
    // ordered unlocking
    if (left_rank < right_rank) {
        right_node.clear_lock();
        left_node.clear_lock();
    } else if (left_rank != right_rank) {
        left_node.clear_lock();
        right_node.clear_lock();
    } else {
        left_node.clear_lock();
    }
}

/*
uint64_t graph_t::edge_delta_to_id(uint64_t base, uint64_t delta) const {
    //assert(delta != 0);
    return node_v[base-1].delta_to_id(base, delta);
}
*/

/*
uint64_t graph_t::edge_to_delta(const handle_t& left, const handle_t& right) const {
    int64_t delta = get_id(right) - get_id(left);
    return (delta == 0 ? 1 : (delta > 0 ? 2*abs(delta) : 2*abs(delta)+1));
}
*/

bool graph_t::has_edge(const handle_t& left, const handle_t& right) const {
    bool exists = false;
    follow_edges(left, false, [&right, &exists](const handle_t& next) {
            if (next == right) exists = true;
            return !exists;
        });
    return exists;
}

/// Remove the edge connecting the given handles in the given order and orientations.
/// Ignores nonexistent edges.
/// Does not update any stored paths.
void graph_t::destroy_edge(const handle_t& left_h, const handle_t& right_h) {
    auto& left_node = get_node_ref(left_h);
    auto& right_node =  get_node_ref(right_h);
    bool left_rev = get_is_reverse(left_h);
    bool right_rev = get_is_reverse(right_h);
    _edge_count -=
        left_node.remove_edge(get_id(right_h), right_rev, false, left_rev)
        && right_node.remove_edge(get_id(left_h), left_rev, true, right_rev);
}

/// Remove all nodes and edges. Does not update any stored paths.
void graph_t::clear() {
    suc_bv null_bv;
    _max_node_id = 0;
    _min_node_id = 0;
    _edge_count = 0;
    deleted_nodes.clear();
    for (auto& n : node_v) {
        delete n;
    }
    node_v.clear();
    for_each_path_handle(
        [&](const path_handle_t& p) {
            // remove from both hash tables
            auto s = get_path_name(p);
            delete &get_path_metadata(p);
            path_metadata_h->Delete(as_integer(p));
            path_name_h->Delete(s);
        });
    _path_count = 0;
    _path_handle_next = 0;
}

void graph_t::clear_paths() {
    for_each_handle(
        [&](const handle_t& handle) {
            node_t& node = get_node_ref(handle);
            node.clear_paths();
        });
    for_each_path_handle(
        [&](const path_handle_t& p) {
            // remove from both hash tables
            auto s = get_path_name(p);
            delete &get_path_metadata(p);
            path_metadata_h->Delete(as_integer(p));
            path_name_h->Delete(s);
        });
    _path_count = 0;
    _path_handle_next = 0;
}

/// Swap the nodes corresponding to the given handles, in the ordering used
/// by for_each_handle when looping over the graph. Other handles to the
/// nodes being swapped must not be invalidated. If a swap is made while
/// for_each_handle is running, it affects the order of the handles
/// traversed during the current traversal (so swapping an already seen
/// handle to a later handle's position will make the seen handle be visited
/// again and the later handle not be visited at all).
void graph_t::swap_handles(const handle_t& a, const handle_t& b) {
    //assert(false);
}

void graph_t::optimize(bool allow_id_reassignment) {
    apply_ordering({}, allow_id_reassignment);
}

bool graph_t::is_optimized(void) {
	if (min_node_id() == 1 && max_node_id() == get_node_count()) {
		return true;
	} else {
		return false;
	}
}

void graph_t::reassign_node_ids(const std::function<nid_t(const nid_t&)>& get_new_id) {
    assert(false); // to implement
    /*
    graph_t reassigned;
    // nodes
    for_each_handle(
        [&](const handle_t& handle) {
            reassigned.create_handle(get_sequence(handle), get_new_id(get_id(handle)));
        });
    // edges
    for_each_edge(
        [&](const edge_t& edge) {
            reassigned.create_edge(reassigned.get_handle(get_new_id(get_id(edge.first)), get_is_reverse(edge.first)),
                                   reassigned.get_handle(get_new_id(get_id(edge.second)), get_is_reverse(edge.second)));
        });
    // paths
    for_each_path_handle(
        [&](const path_handle_t& old_path) {
            path_handle_t new_path = reassigned.create_path_handle(get_path_name(old_path));
            for_each_step_in_path(
                old_path,
                [&](const step_handle_t& step) {
                    handle_t old_handle = get_handle_of_step(step);
                    handle_t new_handle = reassigned.get_handle(get_new_id(get_id(old_handle)),
                                                                get_is_reverse(old_handle));
                    reassigned.append_step(new_path, new_handle);
                });
        });
    *this = reassigned;
    */
    // XXXXXX TODO
}

/// Reorder the graph's internal structure to match that given.
/// Optionally compact the id space of the graph to match the ordering, from 1->|ordering|.
bool graph_t::apply_ordering(const std::vector<handle_t>& order_in, bool compact_ids) {
    // get mapping from old to new id
    // if we're given an empty order, just compact the ids based on our ordering
    const std::vector<handle_t>* order;
    std::vector<handle_t> base_order;
    if (order_in.empty()) {
        base_order.reserve(get_node_count());
        for_each_handle([&](const handle_t& handle) {
                base_order.push_back(handle);
            });
        order = &base_order;
    } else {
        order = &order_in;
    }

    // establish id mapping
    std::vector<std::pair<nid_t, bool>> ids;
    // fill even for deleted nodes, which we map to 0
    ids.resize(node_v.size(), std::make_pair(0, false));

    if (compact_ids) {
        for (uint64_t i = 0; i < order->size(); ++i) {
            ids[number_bool_packing::unpack_number(order->at(i))] =
                std::make_pair(i+1,
                               get_is_reverse(order->at(i)));
        }
    } else {
        for (auto handle : *order) {
            ids[number_bool_packing::unpack_number(handle)] =
                std::make_pair(get_id(handle),
                               get_is_reverse(handle));
        }
    }

    // helpers to map from current to new id and orientation
    auto get_new_id =
        [&](uint64_t id) {
            return ids[id - 1].first;
        };
    auto to_flip =
        [&](uint64_t id) {
            return ids[id - 1].second;
        };

    // nodes, edges, and path steps
#pragma omp parallel for schedule(static, 1) num_threads(_num_threads)
    for (uint64_t i = 0; i < node_v.size(); ++i) {
        handle_t h = number_bool_packing::pack(i,false);
        if (!is_deleted(h)) {
            auto& node = get_node_ref(h);
            node.apply_ordering(get_new_id, to_flip);
        }
    }

    // path metadata
#pragma omp parallel for schedule(static, 1) num_threads(_num_threads)
    for (uint64_t i = 1; i <= _path_handle_next; ++i) {
        path_metadata_t* p;
        if (path_metadata_h->Find(i, p)) {
            const auto& path = as_path_handle(i);
            auto& old_meta = path_metadata(as_path_handle(i));
            path_metadata_h->Delete(as_integer(path));
            path_name_h->Delete(p->name);
            p = new path_metadata_t();
            p->handle.store(old_meta.handle); // same by def
            p->length.store(old_meta.length); // same
            // reassign the handle ids
            step_handle_t f = old_meta.first.load();
            handle_t& f_h = as_handle((uint64_t&)as_integers(f)[0]);
            uint64_t f_id = get_id(f_h);
            f_h = number_bool_packing::pack(get_new_id(f_id)-1, // note -1
                                            get_is_reverse(f_h)^to_flip(f_id));
            p->first.store(f);
            step_handle_t l = old_meta.last.load();
            handle_t& l_h = as_handle((uint64_t&)as_integers(l)[0]);
            uint64_t l_id = get_id(l_h);
            l_h = number_bool_packing::pack(get_new_id(l_id)-1, // note -1
                                            get_is_reverse(l_h)^to_flip(l_id));
            p->last.store(l);
            p->name = old_meta.name;
            p->is_circular.store(old_meta.is_circular);
            path_metadata_h->Insert(as_integer(path), p);
            path_name_h->Insert(p->name, p);
            delete &old_meta;
        }
    }

    // now we actually apply the ordering to our node_v, while removing deleted slots
    std::vector<node_t*> new_node_v; //(order->size());
    _min_node_id = 1;
    if (compact_ids) {
        uint64_t j = 0;
        for (uint64_t i = 0; i < node_v.size(); ++i) {
            if (node_v[i] != nullptr) {
                auto h = (*order)[j++];
                new_node_v.push_back(&get_node_ref(h));
            }
        }
        _max_node_id = new_node_v.size();
    } else {
        uint64_t j = 0;
        for (auto n_v : node_v) {
            if (n_v != nullptr) {
                auto h = (*order)[j++];
                new_node_v.push_back(&get_node_ref(h));
            } else {
                new_node_v.push_back(nullptr);
            }
        }
        _max_node_id = new_node_v.size();
    }
    node_v = new_node_v;
    deleted_nodes.clear();

    return true;
}

void graph_t::apply_path_ordering(const std::vector<path_handle_t>& order) {
    std::vector<path_handle_t> curr_to_new(order.size());
    {
        uint64_t i = 0;
        for (auto& p : order) {
            curr_to_new[as_integer(p)-1] = as_path_handle(++i);
        }
    }
    auto get_new_path_handle =
        [&](const path_handle_t& p) {
            return curr_to_new[as_integer(p)-1];
        };
    // now we save our metadata
    std::vector<path_metadata_t*> metadata;
    metadata.reserve(_path_count.load());
    // then we'll apply this to our metadata map in parallel
    for_each_path_handle(
        [&](const path_handle_t& path) {
            auto& p_m = get_path_metadata(path);
            assert(p_m.handle == path);
            // update our internal handle
            p_m.handle.store(get_new_path_handle(path));
            metadata.push_back(&p_m);
            path_metadata_h->Delete(as_integer(path));
        });
    for (auto* m : metadata) {
        path_metadata_h->Insert(as_integer(m->handle), m);
    }
    // and to the nodes in parallel
    auto get_new_path_id =
        [&](const uint64_t& id) {
            return as_integer(curr_to_new[id-1]);
        };
#pragma omp parallel for schedule(static, 1) num_threads(_num_threads)
    for (uint64_t i = 0; i < node_v.size(); ++i) {
        handle_t h = number_bool_packing::pack(i,false);
        if (!is_deleted(h)) {
            auto& node = get_node_ref(h);
            node.apply_path_ordering(get_new_path_id);
        }
    }
}

/// Alter the node that the given handle corresponds to so the orientation
/// indicated by the handle becomes the node's local forward orientation.
/// Rewrites all edges pointing to the node and the node's sequence to
/// reflect this. Invalidates all handles to the node (including the one
/// passed). Returns a new, valid handle to the node in its new forward
/// orientation. Note that it is possible for the node's ID to change.
/// Updates all stored paths. May change the ordering of the underlying
/// graph.
handle_t graph_t::apply_orientation(const handle_t& handle) {
    // do nothing if we're already in the right orientation
    if (!get_is_reverse(handle)) return handle;
    handle_t fwd_handle = flip(handle);
    handle_t rev_handle = handle;
    // store edges
    vector<handle_t> edges_fwd_fwd;
    vector<handle_t> edges_fwd_rev;
    vector<handle_t> edges_rev_fwd;
    vector<handle_t> edges_rev_rev;
    follow_edges(fwd_handle, false, [&](const handle_t& h) {
            edges_fwd_fwd.push_back(h);
        });
    follow_edges(fwd_handle, true, [&](const handle_t& h) {
            edges_fwd_rev.push_back(h);
        });
    follow_edges(rev_handle, false, [&](const handle_t& h) {
            edges_rev_fwd.push_back(h);
        });
    follow_edges(rev_handle, true, [&](const handle_t& h) {
            edges_rev_rev.push_back(h);
        });
    for (auto& h : edges_fwd_fwd) destroy_edge(fwd_handle, h);
    for (auto& h : edges_fwd_rev) destroy_edge(h, fwd_handle);
    for (auto& h : edges_rev_fwd) destroy_edge(rev_handle, h);
    for (auto& h : edges_rev_rev) destroy_edge(h, rev_handle);
    // save the sequence's reverse complement, which we will use to add the new handle
    //const std::string seq = get_sequence(handle);
    // we have the technology. we can rebuild it.
    // replace the handle sequence
    //set_handle_sequence(handle, seq);
    auto& node = get_node_ref(handle);

    // flip the node sequence
    node.set_sequence(get_sequence(handle));

    // we need to flip any stored path fronts and backs as well
    std::pair<std::map<uint64_t, std::pair<uint64_t, bool>>, // path fronts
              std::map<uint64_t, std::pair<uint64_t, bool>>> // path backs
        path_rewrites = node.flip_paths();
    for (auto& r : path_rewrites.first) {
        auto& p = get_path_metadata(as_path_handle(r.first));
        step_handle_t step;
        as_integers(step)[0]
            = as_integer(number_bool_packing::pack(number_bool_packing::unpack_number(handle),
                                                   r.second.second));
        as_integers(step)[1] = r.second.first;
        p.first.store(step);
    }
    for (auto& r : path_rewrites.second) {
        auto& p = get_path_metadata(as_path_handle(r.first));
        step_handle_t step;
        as_integers(step)[0]
            = as_integer(number_bool_packing::pack(number_bool_packing::unpack_number(handle),
                                                   r.second.second));
        as_integers(step)[1] = r.second.first;
        p.last.store(step);
    }
    // reconnect it to the graph
    for (auto& h : edges_fwd_fwd) {
        create_edge(handle, h);
    }
    for (auto& h : edges_fwd_rev) {
        create_edge(h, handle);
    }
    for (auto& h : edges_rev_fwd) {
        create_edge(flip(handle), h);
    }
    for (auto& h : edges_rev_rev) {
        create_edge(h, flip(handle));
    }
    return flip(handle);
}

void graph_t::set_handle_sequence(const handle_t& handle, const std::string& seq) {
    assert(seq.size());
    auto& node = get_node_ref(handle);
    node.get_lock();
    node.set_sequence(seq);
    node.clear_lock();
}

/// Split a handle's underlying node at the given offsets in the handle's
/// orientation. Returns all of the handles to the parts. Other handles to
/// the node being split may be invalidated. The split pieces stay in the
/// same local forward orientation as the original node, but the returned
/// handles come in the order and orientation appropriate for the handle
/// passed in.
/// Updates stored paths.
std::vector<handle_t> graph_t::divide_handle(const handle_t& handle, const std::vector<size_t>& offsets) {
    // convert the offsets to the forward strand, if needed
    std::vector<uint64_t> fwd_offsets = { 0 };
    uint64_t length = get_length(handle);
    if (get_is_reverse(handle)) {
        for (auto& o : offsets) fwd_offsets.push_back(length-o);
    } else {
        for (auto& o : offsets) fwd_offsets.push_back(o);
    }
    std::sort(fwd_offsets.begin(), fwd_offsets.end());
    handle_t fwd_handle = get_is_reverse(handle) ? flip(handle) : handle;
    // break it into the given pieces by building up the new node sequences
    std::string seq = get_sequence(fwd_handle);
    fwd_offsets.push_back(seq.size());
    std::vector<std::string> seqs;
    for (uint64_t i = 0; i < fwd_offsets.size()-1; ++i) {
        seqs.push_back(seq.substr(fwd_offsets[i], fwd_offsets[i+1]-fwd_offsets[i]));
    }
    // make the handles
    std::vector<handle_t> handles;
    for (auto& s : seqs) {
        handles.push_back(create_handle(s));
    }

    // update the backwards references back onto this node
    handle_t last_handle = handles.back();
    follow_edges(fwd_handle, false, [&](const handle_t& h) {
            if (h == flip(fwd_handle)) {
                // destroy_edge(last_handle, h);
                create_edge(last_handle, flip(handles.back()));
            }
    });
    // update the forward references back onto this node
    follow_edges(handle, left, [&](const handle_t& h) {
        if (h == flip(handle)) {
            create_edge(h, handles.front());
        }
    });

    // and record their reverse, for use in path fixup
    std::vector<handle_t> rev_handles;
    for (auto& h : handles) {
        rev_handles.push_back(flip(h));
    }
    std::reverse(rev_handles.begin(), rev_handles.end());


    // connect the pieces head to tail
    for (uint64_t i = 0; i < handles.size()-1; ++i) {
        create_edge(handles[i], handles[i+1]);
    }
    // collect the handle's path context
    vector<step_handle_t> steps;
    for_each_step_on_handle(handle, [&](const step_handle_t& step) {
            steps.push_back(step);
        });
    // reverse the order to allow for safe rewriting
    std::reverse(steps.begin(), steps.end());
    // replace path steps with the new handles
    for (auto& step : steps) {
        handle_t h = get_handle_of_step(step);
        if (get_is_reverse(h)) {
            rewrite_segment(step, step, rev_handles);
        } else {
            rewrite_segment(step, step, handles);
        }
    }
    get_node_ref(handle).clear_paths();
    // collect the context of the handle
    vector<handle_t> edges_fwd_fwd;
    vector<handle_t> edges_fwd_rev;
    vector<handle_t> edges_rev_fwd;
    vector<handle_t> edges_rev_rev;
    follow_edges(fwd_handle, false, [&](const handle_t& h) {
            edges_fwd_fwd.push_back(h);
        });
    follow_edges(fwd_handle, true, [&](const handle_t& h) {
            edges_fwd_rev.push_back(h);
        });
    follow_edges(flip(fwd_handle), false, [&](const handle_t& h) {
            edges_rev_fwd.push_back(h);
        });
    follow_edges(flip(fwd_handle), true, [&](const handle_t& h) {
            edges_rev_rev.push_back(h);
        });
    // connect the ends to the previous context
    for (auto& h : edges_fwd_fwd) create_edge(handles.back(), h);
    for (auto& h : edges_fwd_rev) create_edge(h, handles.front());
    for (auto& h : edges_rev_fwd) create_edge(rev_handles.back(), h);
    for (auto& h : edges_rev_rev) create_edge(h, rev_handles.front());
    // destroy the handles
    bool h_is_rev =  get_is_reverse(handle);
    destroy_handle(fwd_handle);
    destroy_handle(handle);
    return h_is_rev ? rev_handles : handles;
}

handle_t graph_t::combine_handles(const std::vector<handle_t>& handles) {
    std::string seq;
    for (auto& handle : handles) {
        seq.append(get_sequence(handle));
    }
    handle_t combined = create_handle(seq);
    // relink the inbound and outbound nodes
    // get the edge context
    vector<handle_t> edges_fwd_fwd;
    vector<handle_t> edges_fwd_rev;
    vector<handle_t> edges_rev_fwd;
    vector<handle_t> edges_rev_rev;
    follow_edges(handles.back(), false, [&](const handle_t& h) {
            edges_fwd_fwd.push_back(h);
        });
    follow_edges(handles.front(), true, [&](const handle_t& h) {
            edges_fwd_rev.push_back(h);
        });
    // destroy the old handles
    for (auto& handle : handles) {
        destroy_handle(handle);
    }
    // connect the ends to the previous context
    // check that we're not trying to make edges that connect back with the nodes in the component
    // there are three cases
    // self looping, front and rear inverting
    for (auto& h : edges_fwd_fwd) {
        if (h == handles.front()) {
            create_edge(combined, combined);
        } else if (h == flip(handles.back())) {
            create_edge(combined, flip(combined));
        } else {
            create_edge(combined, h);
        }
    }
    for (auto& h : edges_fwd_rev) {
        if (h == handles.back()) {
            create_edge(combined, combined);
        } else if (h == flip(handles.front())) {
            create_edge(flip(combined), combined);
        } else {
            create_edge(h, combined);
        }
    }
    return combined;
}

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
void graph_t::destroy_path(const path_handle_t& path) {
    // select everything with that handle in the path_handle_wt
    std::vector<step_handle_t> path_v;
    for_each_step_in_path(path, [this,&path_v](const step_handle_t& step) {
            path_v.push_back(step);
        });
    // this zeros out the steps
    // final removal of these steps depends on a call to graph_t::optimize
    for (auto& step : path_v) {
        destroy_step(step);
    }
    auto& p = get_path_metadata(path);
    // our length should be 0
    assert(p.length == 0);
    path_metadata_h->Delete(as_integer(p.handle));
    path_name_h->Delete(p.name);
    delete &p;
    --_path_count;
}

/**
 * Create a path with the given name. The caller must ensure that no path
 * with the given name exists already, or the behavior is undefined.
 * Returns a handle to the created empty path. Handles to other paths must
 * remain valid.
 */
path_handle_t graph_t::create_path_handle(const std::string& name, bool is_circular) {
    path_handle_t path = as_path_handle(++_path_handle_next);
    path_metadata_t* _p = new path_metadata_t();
    auto& p = *_p;
    step_handle_t step;
    p.handle = path;
    as_integers(step)[0] = 0;
    as_integers(step)[1] = 0;
    p.first.store(step);
    p.last.store(step);
    p.length.store(0);
    p.name = name;
    p.is_circular = is_circular;
    ++_path_count; // atomic
    path_metadata_h->Insert(as_integer(path), _p);
    path_name_h->Insert(name, _p);
    auto& q = path_metadata(path);
    return path;
}

step_handle_t graph_t::create_step(const path_handle_t& path, const handle_t& handle) {
    // where are we going to insert?
    auto& node = get_node_ref(handle);
    node.get_lock();
    uint64_t rank_on_handle = node.path_count();
    // build our step
    step_handle_t step;
    as_integers(step)[0] = as_integer(handle);
    as_integers(step)[1] = rank_on_handle;
    node.add_path_step(as_integer(path), get_is_reverse(handle),
                       true, true,
                       0, 0, 0, 0);
    node.clear_lock();
    return step;
}

void graph_t::link_steps(const step_handle_t& from, const step_handle_t& to) {
    path_handle_t path = get_path(from);
    assert(path == get_path(to));
    const handle_t& from_handle = get_handle_of_step(from);
    const handle_t& to_handle = get_handle_of_step(to);
    const uint64_t& from_rank = as_integers(from)[1];
    const uint64_t& to_rank = as_integers(to)[1];
    node_t& from_node = get_node_ref(from_handle);
    from_node.get_lock();
    from_node.set_step_next_id(from_rank, get_id(to_handle));
    from_node.set_step_next_rank(from_rank, to_rank);
    from_node.set_step_is_end(from_rank, false);
    from_node.clear_lock();
    node_t& to_node = get_node_ref(to_handle);
    to_node.get_lock();
    to_node.set_step_prev_id(to_rank, get_id(from_handle));
    to_node.set_step_prev_rank(to_rank, from_rank);
    to_node.set_step_is_start(to_rank, false);
    to_node.clear_lock();
}

void graph_t::destroy_step(const step_handle_t& step_handle) {
    // erase reference to this step
    bool has_prev = has_previous_step(step_handle);
    bool has_next = has_next_step(step_handle);
    path_handle_t path = get_path_handle_of_step(step_handle);
    auto& path_meta = get_path_metadata(path);
    if (has_prev) {
        // nothing
    } else if (has_next) {
        auto step = get_next_step(step_handle);
        path_meta.first = step;
    }
    if (has_next) {
        // nothing
    } else if (has_prev) {
        auto step = get_previous_step(step_handle);
        path_meta.last = step;
    }
    // reduce the step count in the path
    --path_meta.length;
    node_t& curr_node = get_node_ref(get_handle_of_step(step_handle));
    curr_node.get_lock();
    curr_node.clear_path_step(as_integers(step_handle)[1]);
    curr_node.clear_lock();
}

step_handle_t graph_t::prepend_step(const path_handle_t& path, const handle_t& to_append) {
    // get the last step
    auto& p = get_path_metadata(path);
    // create the new step
    step_handle_t new_step = create_step(path, to_append);
    if (!p.length) {
        p.last.store(new_step);
    } else {
        step_handle_t first_step = path_begin(path);
        // link it to the last step
        link_steps(new_step, first_step);
    }
    // point to the new last step
    p.first.store(new_step);
    // update our step count
    ++p.length;
    return new_step;
}

step_handle_t graph_t::append_step(const path_handle_t& path, const handle_t& to_append) {
    // get the last step
    auto& p = get_path_metadata(path);
    // create the new step
    step_handle_t new_step = create_step(path, to_append);
    if (!p.length) {
        p.first.store(new_step);
    } else {
        step_handle_t last_step = path_back(path);
        // link it to the last step
        link_steps(last_step, new_step);
    }
    // point to the new last step
    p.last.store(new_step);
    // update our step count
    ++p.length;
    return new_step;
}

/// helper to handle the case where we remove an step from a given path
/// on a node that has other steps from the same path, thus invalidating the
/// ranks used to refer to it
void graph_t::decrement_rank(const step_handle_t& step_handle) {
    // what is the actual rank of this step?
    //std::cerr << "in decrement rank " << get_handle_of_step(step_handle) << ":" << as_integers(step_handle)[1] << std::endl;
    if (has_previous_step(step_handle)) {
        auto step = get_previous_step(step_handle);
        // decrement the rank information
        node_t& step_node = get_node_ref(get_handle_of_step(step));
        step_node.get_lock();
        uint64_t step_rank = as_integers(step)[1];
        step_node.set_step_next_rank(step_rank,
                                     step_node.step_next_rank(step_rank)-1);
        step_node.clear_lock();
    } else {
        // update path metadata
        auto& p = get_path_metadata(get_path(step_handle));
        step_handle_t s = p.first.load();
        --as_integers(s)[1];
        p.first.store(s);
    }
    if (has_next_step(step_handle)) {
        auto step = get_next_step(step_handle);
        node_t& step_node = get_node_ref(get_handle_of_step(step));
        step_node.get_lock();
        uint64_t step_rank = as_integers(step)[1];
        step_node.set_step_prev_rank(step_rank,
                                     step_node.step_prev_rank(step_rank)-1);
        step_node.clear_lock();
    } else {
        // update path metadata
        auto& p = get_path_metadata(get_path(step_handle));
        step_handle_t s = p.last.load();
        --as_integers(s)[1];
        p.last.store(s);
    }
}

// Insert a visit to a node to the given path between the given steps.
step_handle_t graph_t::insert_step(const step_handle_t& before, const step_handle_t& after, const handle_t& to_insert) {
    auto p = rewrite_segment(before, after, { to_insert });
    return get_next_step(p.first);
}

/// reassign the given step to the new handle
step_handle_t graph_t::set_step(const step_handle_t& step_handle, const handle_t& assign_to) {
    return rewrite_segment(step_handle, step_handle, { assign_to }).first;
}

/// Replace the path range with the new segment (range is _inclusive_)
std::pair<step_handle_t, step_handle_t> graph_t::rewrite_segment(const step_handle_t& segment_begin,
                                                                 const step_handle_t& segment_end,
                                                                 const std::vector<handle_t>& new_segment) {
    // collect the steps to replace
    std::vector<step_handle_t> steps;
    //std::string old_seq, new_seq;
    for (step_handle_t step = segment_begin; ; step = get_next_step(step)) {
        steps.push_back(step);
        //old_seq.append(get_sequence(get_handle_of_step(step)));
        if (step == segment_end) break;
    }
    //for (auto& handle : new_segment) new_seq.append(get_sequence(handle));
    // verify that we're making a valid rewrite
    //assert(old_seq == new_seq); // not required, we could be modifying the path
    // find the before and after steps, which we'll link into
    bool is_begin = !has_previous_step(segment_begin);
    bool is_end = !has_next_step(segment_end);
    step_handle_t before = get_previous_step(segment_begin);
    step_handle_t after = get_next_step(segment_end);
    // get the path metadata
    path_handle_t path = get_path(segment_begin);
    auto& path_meta = get_path_metadata(path);
    // step destruction simply zeros out our step data
    // a final removal of the deleted steps requires a call to graph_t::optimize

    for (auto& step : steps) {
        destroy_step(step);
    }
    if (path_meta.length == 0 && new_segment.size() == 0) {
        //std::cerr << "destroyed path" << std::endl;
    }
    // create the new steps
    std::vector<step_handle_t> new_steps;
    for (auto& handle : new_segment) {
        new_steps.push_back(create_step(path, handle));
    }
    // delete the previous steps
    path_meta.length += new_steps.size();
    if (new_steps.size()) {
        // link new steps together
        for (uint64_t i = 0; i < new_steps.size()-1; ++i) {
            link_steps(new_steps[i], new_steps[i+1]);
        }
        if (!is_begin) {
            link_steps(before, new_steps.front());
        } else {
            path_meta.first = new_steps.front();
        }
        if (!is_end) {
            link_steps(new_steps.back(), after);
        } else {
            path_meta.last = new_steps.back();
        }
        return make_pair(new_steps.front(), new_steps.back());
    } else {
        return make_pair(path_front_end(path), path_end(path));
    }
}

////////////////////////////////////////////////////////////////////////////
// Rank handle interface
////////////////////////////////////////////////////////////////////////////
//Return the rank of a node (ranks start at 1 and are dense).
size_t graph_t::id_to_rank(const nid_t &node_id) const {
    return node_id + 1;
}

//Return the node with a given rank.
nid_t graph_t::rank_to_id(const size_t &rank) const {
    return rank - 1;
}

// Return the rank of a handle (ranks start at 1 and are dense, and each
// orientation has its own rank). Handle ranks may not have anything to do
// with node ranks.
size_t graph_t::handle_to_rank(const handle_t &handle) const {
    return number_bool_packing::unpack_number(handle);
}

// Return the handle with a given rank.
handle_t graph_t::rank_to_handle(const size_t &rank) const {
    return number_bool_packing::pack(rank, false);
}

void graph_t::display() const {
    std::cerr << "------ graph state ------" << std::endl;

    std::cerr << "_max_node_id = " << _max_node_id << std::endl;
    std::cerr << "_min_node_id = " << _min_node_id << std::endl;

    //std::cerr << "graph_id_map" << "\t";
    //for (auto& k : graph_id_map) std::cerr << k.first << "->" << k.second << " "; std::cerr << std::endl;
    std::cerr << "node_v" << "\t";
    for (uint64_t i = 0; i < node_v.size(); ++i) {
        if (node_v.at(i) == nullptr) {
            std::cerr << "null | ";
            continue;
        }
        auto& node = *node_v.at(i);
        nid_t node_id = i+1;
        std::cerr << node_id << ":" << node.get_sequence() << " ";
        //const auto& node_edges = node.edges;
        uint64_t j = 0;
        node.for_each_edge(
            [&](nid_t other_id,
                bool other_rev,
                bool to_curr,
                bool on_rev) {
                std::cerr << other_id << ":"
                          << other_rev << ":"
                          << to_curr << ":"
                          << on_rev << ",";
                return true;
            });
        std::cerr << " _ ";
        const std::vector<node_t::step_t> steps = node.get_path_steps();
        for (auto& step : steps) {
            std::cerr << step.path_id << ":"
                      << step.is_rev << ":"
                      << step.is_start << ":"
                      << step.is_end << ":"
                      << step.prev_id << ":"
                      << step.prev_rank << ":"
                      << step.next_id << ":"
                      << step.next_rank << " ";
        }
        std::cerr << " | ";
    }
    std::cerr << std::endl;
    std::cerr << "deleted_nodes" << "\t";
    for (auto& idx : deleted_nodes) std::cerr << idx << " "; std::cerr << std::endl;
    /// Ordered across the nodes in graph_id_iv, stores the path ids (1-based) at each
    /// segment in seq_wt, delimited by 0, one for each path step (node traversal).
    std::cerr << "path_metadata" << "\t";
    for_each_path_handle(
        [&](const path_handle_t& path) {
            auto& p = path_metadata(path);
            std::cerr << as_integer(path) << ":" << p.name << ":"
                      << get_id(as_handle(as_integers(p.first)[0])) << "/" << as_integers(p.first)[1] << "->"
                      << get_id(as_handle(as_integers(p.last)[0])) << "/" << as_integers(p.last)[1] << " ";
        });
    std::cerr << std::endl;

    // not dumped...
    /// Stores path names in their internal order, delimited by '$'
    //dyn::wt_fmi path_name_fmi;
    /// Marks the beginning of each path name
    //dyn::suc_bv path_name_bv;

}

void graph_t::to_gfa(std::ostream& out, const bool& emit_node_annotation) const {
    out << "H\tVN:Z:1.0" << std::endl;
    // for each node
    for_each_handle([&out,&emit_node_annotation, this](const handle_t& h) {
            out << "S\t" << get_id(h) << "\t" << get_sequence(h);
            if (emit_node_annotation) {
                out << "\t"
                << "DP:i:" << get_step_count(h) << "\t"
                << "RC:i:" << get_step_count(h) * get_length(h);
            }
            out << std::endl;

            {
                // use this direct iteration to avoid double counting edges
                // we only consider write the edges relative to their start
                const node_t& node = get_node_cref(h);
                bool is_rev = get_is_reverse(h);
                nid_t node_id = get_id(h);
                node.for_each_edge(
                    [&](nid_t other_id,
                        bool other_rev,
                        bool to_curr,
                        bool on_rev) {
                        if (!to_curr) {
                            out << "L\t" << node_id << "\t"
                                << (on_rev?"-":"+") << "\t"
                                << other_id << "\t"
                                << (other_rev?"-":"+") << "\t"
                                << "0M" << std::endl;
                        }
                        return true;
                    });
            }
        });
    for_each_path_handle([&out,this](const path_handle_t& p) {
            out << "P\t" << get_path_name(p) << "\t";
            auto& path_meta = path_metadata(p);
            uint64_t i = 0;
            for_each_step_in_path(p, [this,&i,&out](const step_handle_t& step) {
                    handle_t h = get_handle_of_step(step);
                    out << get_id(h) << (get_is_reverse(h)?"-":"+");
                    if (has_next_step(step)) out << ",";
                    ++i;
                });
            out << "\t" << "*"; // always put at least a "*" in the overlaps field
            if (get_is_circular(p)) {
                out << "\t" << "TP:Z:circular";
            }
            assert(i == path_meta.length);
            //out << "\t" << "steps:i:" << path_meta.length;
            out << std::endl;
        });
}

uint32_t graph_t::get_magic_number() const {
    return 1988148666ul; // TODO update me
}

void graph_t::serialize_members(std::ostream& out) const {
    //rebuild_id_handle_mapping();
    uint64_t written = 0;
    out.write((char*)&_max_node_id,sizeof(_max_node_id));
    written += sizeof(_max_node_id);
    out.write((char*)&_min_node_id,sizeof(_min_node_id));
    written += sizeof(_min_node_id);
    uint64_t node_count = node_v.size();
    out.write((char*)&node_count,sizeof(node_count));
    written += sizeof(node_count);
    out.write((char*)&_edge_count,sizeof(_edge_count));
    written += sizeof(_edge_count);
    out.write((char*)&_path_count,sizeof(_path_count));
    written += sizeof(_path_count);
    out.write((char*)&_path_handle_next,sizeof(_path_handle_next));
    written += sizeof(_path_handle_next);
    out.write((char*)&_id_increment,sizeof(_id_increment));
    written += sizeof(_id_increment);
    //assert(node_count == node_v.size());
    // hack
    // todo big mess, middle of removal of deleted node bv
    node_t empty_node;
    for (auto& node : node_v) {
        // check if node is null
        if (node == nullptr) {
            written += empty_node.serialize(out);
        } else {
            written += node->serialize(out);
        }
    }
    // there are _path_count of these to write
    uint64_t j = 0;
    for_each_path_handle(
        [&](const path_handle_t& path) {
            auto& m = path_metadata(path);
            out.write((char*)&m.length,sizeof(m.length));
            written += sizeof(m.length);
            out.write((char*)&m.first,sizeof(m.first));
            written += sizeof(m.first);
            out.write((char*)&m.last,sizeof(m.last));
            written += sizeof(m.last);
            size_t k = m.name.size();
            out.write((char*)&k,sizeof(k));
            written += sizeof(k);
            out.write((char*)m.name.c_str(),m.name.size());
            written += k;
            ++j;
        });
    assert(j == _path_count);
}

void graph_t::deserialize_members(std::istream& in) {
    in.read((char*)&_max_node_id,sizeof(_max_node_id));
    in.read((char*)&_min_node_id,sizeof(_min_node_id));
    uint64_t node_count = node_v.size();
    in.read((char*)&node_count,sizeof(node_count));
    in.read((char*)&_edge_count,sizeof(_edge_count));
    in.read((char*)&_path_count,sizeof(_path_count));
    in.read((char*)&_path_handle_next,sizeof(_path_handle_next));
    in.read((char*)&_id_increment,sizeof(_id_increment));
    node_v.resize(node_count,nullptr);
    for (size_t i = 0; i < node_count; ++i) {
        node_v[i] = new node_t;
        auto& node = node_v[i];
        node->load(in);
        if (node->get_id() == 0) {
            // detect which nodes are deleted
            // these must be the only ones with id == 0
            // they have been stored as empty node records
            delete node;
            node = nullptr;
            deleted_nodes.insert(i+1);
        }
    }
    for (size_t j = 0; j < _path_count; ++j) {
        path_metadata_t* _p = new path_metadata_t();
        auto& m = *_p;
        m.handle = as_path_handle(j+1);
        in.read((char*)&m.length,sizeof(m.length));
        in.read((char*)&m.first,sizeof(m.first));
        in.read((char*)&m.last,sizeof(m.last));
        uint64_t s;
        in.read((char*)&s,sizeof(s));
        char n[s+1]; n[s] = '\0';
        in.read(n,s);
        m.name = string(n);
        path_metadata_h->Insert(as_integer(m.handle), _p);
        path_name_h->Insert(m.name, _p);
    }
}


void graph_t::set_number_of_threads(uint64_t num_threads) {
    _num_threads = num_threads;
}

uint64_t graph_t::get_number_of_threads() {
    return _num_threads;
}

void graph_t::copy(const graph_t& other) {
    clear();
    _max_node_id.store(other._max_node_id);
    _min_node_id.store(other._min_node_id);
    _edge_count.store(other._edge_count);
    _path_count.store(other._path_count);
    _path_handle_next.store(other._path_handle_next);
    _id_increment.store(other._id_increment);
    node_v.resize(other.node_v.size());
    for (size_t i = 0; i < other.node_v.size(); ++i) {
        node_v[i] = new node_t;
        auto* node = node_v[i];
        node->copy(other.get_node_cref(as_handle(i)));
    }
    deleted_nodes = other.deleted_nodes;
    // copy the path metadata
    // the paths themselves should have been copied
    // XXX SLOW
    other.for_each_path_handle(
        [&](const path_handle_t& p) {
            auto new_path = create_path_handle(other.get_path_name(p),
                                               other.get_is_circular(p));
            auto& new_path_meta = get_path_metadata(new_path);
            new_path_meta.copy(other.path_metadata(p));
        });
}

void graph_t::set_static(void) {
    for (auto& node : node_v) {
        node->set_static();
    }
}

void graph_t::set_volatile(void) {
    for (auto& node : node_v) {
        node->set_volatile();
    }
}

}
