#include "stepindex.hpp"
#include "progress.hpp"

namespace odgi {
namespace algorithms {

step_index_t::step_index_t() {
	step_mphf = new boophf_step_t();
}



step_index_t::step_index_t(const PathHandleGraph& graph,
                           const std::vector<path_handle_t>& paths,
                           const uint64_t& nthreads,
                           const bool progress,
						   const uint64_t& sample_rate) {
	this->sample_rate = sample_rate;

    // iterate through the paths, recording steps in the structure we'll use to build the mphf
    std::vector<step_handle_t> steps;
	std::unique_ptr<algorithms::progress_meter::ProgressMeter> collecting_steps_progress_meter;
	if (progress) {
		collecting_steps_progress_meter = std::make_unique<algorithms::progress_meter::ProgressMeter>(
				paths.size(), "[odgi::algorithms::stepindex] Collecting Steps Progress:");
	}
	path_len.resize(paths.size());
#pragma omp parallel for schedule(dynamic,1)
    for (auto& path : paths) {
        std::vector<step_handle_t> my_steps;
		uint64_t path_length = 0;
        graph.for_each_step_in_path(
            path, [&](const step_handle_t& step) {
				path_length += graph.get_length(graph.get_handle_of_step(step));
				// sampling
				if (sample_rate == 0 || 0 == utils::modulo(graph.get_id(graph.get_handle_of_step(step)), sample_rate)) {
					my_steps.push_back(step);
				}
			});
		// sampling
		if (sample_rate == 0 || 0 == utils::modulo(graph.get_id(graph.get_handle_of_step(graph.path_end(path))), sample_rate)) {
			my_steps.push_back(graph.path_end(path));
		}

#pragma omp critical (path_len)
		path_len[as_integer(path) - 1] = path_length;
#pragma omp critical (steps_collect)
        steps.insert(steps.end(), my_steps.begin(), my_steps.end());
        if(progress) {
        	collecting_steps_progress_meter->increment(1);
        }
    }
	if (progress) {
		collecting_steps_progress_meter->finish();
	}
    // sort the steps
    ips4o::parallel::sort(steps.begin(), steps.end(), std::less<>(), nthreads);
    // build the hash function (quietly)
    step_mphf = new boophf_step_t(steps.size(), steps, nthreads, 2.0, false, false);
    // use the hash function to record the step positions
    pos.resize(steps.size());
	std::unique_ptr<algorithms::progress_meter::ProgressMeter> building_progress_meter;
	if (progress) {
		building_progress_meter = std::make_unique<algorithms::progress_meter::ProgressMeter>(
				paths.size(), "[odgi::algorithms::stepindex] Building Progress:");
	}
#pragma omp parallel for schedule(dynamic,1)
    for (auto& path : paths) {
        uint64_t offset = 0;
        graph.for_each_step_in_path(
            path, [&](const step_handle_t& step) {
				// sampling
				if (sample_rate == 0 || 0 == utils::modulo(graph.get_id(graph.get_handle_of_step(step)), sample_rate)) {
					pos[step_mphf->lookup(step)] = offset;
				}
				offset += graph.get_length(graph.get_handle_of_step(step));
				});
		// sampling
		if (sample_rate == 0 || 0 == utils::modulo(graph.get_id(graph.get_handle_of_step(graph.path_end(path))), sample_rate)) {
			pos[step_mphf->lookup(graph.path_end(path))] = offset;
		}
        if (progress) {
        	building_progress_meter->increment(1);
        }
    }
	if (progress) {
		building_progress_meter->finish();
	}
}

const uint64_t step_index_t::get_position(const step_handle_t& step, const PathHandleGraph& graph) const {
	// is our step already in a node that we indexed?
	handle_t h = graph.get_handle_of_step(step);
	uint64_t n_id = graph.get_id(h);
	step_handle_t cur_step = step;
	if (this->sample_rate == 0 || 0 == utils::modulo(n_id, this->sample_rate)) {
		return pos[step_mphf->lookup(step)];
	} else {
		// did we hit the first step anyhow?
		if (!graph.has_previous_step(cur_step)) {
			return 0;
		}
		uint64_t walked = 0;
		while (graph.has_previous_step(cur_step)) {
			step_handle_t prev_step = graph.get_previous_step(cur_step);
			handle_t prev_h = graph.get_handle_of_step(prev_step);
			uint64_t prev_n_id = graph.get_id(prev_h);
			walked += graph.get_length(prev_h);
			if (utils::modulo(prev_n_id, this->sample_rate) == 0) {
				return pos[step_mphf->lookup(prev_step)] + walked;
			}
			cur_step = prev_step;
		}
		return walked;
	}
}

const uint64_t step_index_t::get_path_len(const path_handle_t& path) const {
	return path_len[as_integer(path) - 1];
}

void step_index_t::save(const std::string& name) const {
	std::ofstream stpidx_out(name);
	serialize_members(stpidx_out);
	step_mphf->save(stpidx_out);
}

void step_index_t::load(const std::string& name) {
	std::ifstream stpidx_in(name);
	deserialize_members(stpidx_in);
	step_mphf->load(stpidx_in);
}

void step_index_t::serialize_members(std::ostream &out) const {
	serialize_and_measure(out);
}

size_t step_index_t::serialize_and_measure(std::ostream &out, sdsl::structure_tree_node *s, std::string name) const {

	sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
	size_t written = 0;

	// Do the magic number
	std::string sample_rate = std::to_string(this->sample_rate);
	out << "STEP" << sample_rate << "INDEX";
	written += 9;
	written += sample_rate.length();

	// POSITION STUFF
	written += pos.serialize(out, child, "path_position_map");
	// PATH LENGTH STUFF
	written += path_len.serialize(out, child, "path_length_map");

	sdsl::structure_tree::add_size(child, written);
	return written;
}

void step_index_t::deserialize_members(std::istream &in) {
	// simple alias to match an external interface
	load_sdsl(in);
}

void step_index_t::load_sdsl(std::istream &in) {

	if (!in.good()) {
		throw std::runtime_error("[odgi::algorithms::stepindex] error: SDSL step index file does not exist or step index stream cannot be read.");
	}

	// We need to look for the magic value(s)
	char buffer;
	char * step_buffer = new char [4];
	char * index_buffer = new char [4];
	std::string sample_rate = "";
	std::string index = "";

	in.read(step_buffer, 4);
	// https://stackoverflow.com/questions/1195675/convert-a-char-to-stdstring/1195705#1195705
	std::string step(step_buffer, 4);
	if (step == "STEP"){
		// now we need collect all the characters which will form our sample rate
		while(buffer != 'I') {
			in.get(buffer);
			sample_rate += buffer;
		}
		this->sample_rate = std::stoi(sample_rate);
		index += buffer;
		in.read(index_buffer, 4);
		std::string index_buffer_string(index_buffer, 4);
		index += index_buffer_string;
		if (index != "INDEX") {
			throw std::runtime_error("[odgi::algorithms::stepindex] error: SDSL step index file does not have 'INDEX' in its magic value. The file must be malformed.");
		}
	} else  {
		throw std::runtime_error("[odgi::algorithms::stepindex] error: SDSL step index file does not have 'STEP' in its magic value. The file must be malformed.");
	}

	delete[] step_buffer;
	delete[] index_buffer;

	try {
		pos.load(in);
		path_len.load(in);
	} catch (const std::runtime_error &e) {
		// Pass XGFormatErrors through
		throw e;
	} catch (const std::bad_alloc &e) {
		// We get std::bad_alloc generally if we try to read arbitrary data as an xg index.
		std::cerr << "[odgi::algorithms::stepindex] error: SDSL step index input data not in correct format. " << std::endl;
		exit(1);
	} catch (const std::exception &e) {
		// Other things will get re-thrown with a hint.
		std::cerr << "[odgi::algorithms::stepindex] error: SDSL step index file malformed. Is it really on the correct format STEPsample_rateINDEX?" << std::endl;
		throw e;
	}
}

step_index_t::~step_index_t(void) {
    delete step_mphf;
}


// path step index

path_step_index_t::path_step_index_t(const PathHandleGraph& graph,
                                     const path_handle_t& path,
                                     const uint64_t& nthreads) {
    // iterate through the paths, recording steps in the structure we'll use to build the mphf
    {
        std::vector<nid_t> nodes;
        graph.for_each_step_in_path(
            path, [&](const step_handle_t& step) {
                steps.push_back(step);
                nodes.push_back(graph.get_id(graph.get_handle_of_step(step)));
            });
        steps.push_back(graph.path_end(path)); // dangerous...
        // sort the steps, nb. they're unique
        ips4o::parallel::sort(steps.begin(), steps.end(), std::less<>(), nthreads);
        // build the hash function (quietly)
        step_mphf = new boophf_step_t(steps.size(), steps, nthreads, 2.0, false, false);

        ips4o::parallel::sort(nodes.begin(), nodes.end(), std::less<>(), nthreads);
        // then take unique positions
        nodes.erase(std::unique(nodes.begin(),
                                nodes.end()),
                    nodes.end());
        // build the hash function (quietly)
        node_mphf = new boophf_uint64_t(nodes.size(), nodes, nthreads, 2.0, false, false);
        node_count = nodes.size();
        step_count = steps.size();
    }

    {
        // here, we sort steps by handle and then position
        // and build our handle->step list and step->offset maps
        // these are steps sorted by the bbhash of their node id, and then their offset in the path
        std::vector<std::tuple<uint64_t, uint64_t, step_handle_t, step_it>> steps_by_node;
        uint64_t offset = 0;
        for (step_it step = steps.begin(); step != steps.end(); ++step) {
            steps_by_node.emplace_back(node_mphf->lookup(graph.get_id(graph.get_handle_of_step(*step))), offset, *step, step);
            offset += graph.get_length(graph.get_handle_of_step(*step));
        }
        if (offset == 0) {
            std::cerr << "[odgi::algorithms::stepindex] unable to index empty path " << graph.get_path_name(path) << std::endl;
            std::abort();
        }

        node_offset.resize(node_count+1);
        step_offset.resize(step_count+1);
        ips4o::parallel::sort(steps_by_node.begin(), steps_by_node.end(), std::less<>(), nthreads);
        uint64_t last_idx = 0;
        node_offset[0] = 0; // first offset is 0 by definition
        for (auto& node_step : steps_by_node) {
            auto& idx = std::get<0>(node_step);
            auto& offset = std::get<1>(node_step); // just used for sorting
            auto& step = std::get<2>(node_step);
            auto& it_step = std::get<3>(node_step);
            //std::cerr << "idx = " << idx << " " << as_integers(step)[0] << ":" << as_integers(step)[1] << std::endl;
            if (idx != last_idx) {
                node_offset[idx] = node_steps.size();
            }
            step_offset[step_mphf->lookup(step)] = node_steps.size();
            node_steps.push_back(std::make_pair(it_step, offset));
            last_idx = idx;
        }
        if (last_idx != node_count-1) {
            std::cerr << "last_idx vs node count " << last_idx << " " << node_count << std::endl;
            std::cerr << "[odgi::algorithms::stepindex] unexpected mismatch between last_idx and node_count" << std::endl;
            std::abort();
        }
        node_offset[node_count] = node_steps.size();
        step_offset[step_count] = node_steps.size();
    }
}

path_step_index_t::~path_step_index_t(void) {
    delete node_mphf;
    delete step_mphf;
}

uint64_t path_step_index_t::get_node_idx(const nid_t& id) const {
    return node_mphf->lookup(id);
}

uint64_t path_step_index_t::get_step_idx(const step_handle_t& step) const {
    return step_mphf->lookup(step);
}

uint64_t path_step_index_t::n_steps_on_node(const nid_t& id) const {
    auto idx = get_node_idx(id);
    return node_offset[idx+1] - node_offset[idx];
}

std::pair<bool, std::pair<path_step_index_t::step_it, uint64_t>>
path_step_index_t::get_next_step_on_node(const nid_t& id, const step_handle_t& step) const {
    auto node_idx = get_node_idx(id);
    auto curr_steps = node_offset[node_idx];
    auto next_steps = node_offset[node_idx+1];
    auto step_idx = step_offset[get_step_idx(step)];
    bool has_next = step_idx + 1 < next_steps;
    if (has_next) {
        return std::pair(true, node_steps[step_idx+1]);
    } else {
        auto empty_step = steps.end();
        return std::pair(false, std::pair(empty_step, 0));
    }
}

std::pair<bool, std::pair<path_step_index_t::step_it, uint64_t>>
path_step_index_t::get_prev_step_on_node(const nid_t& id, const step_handle_t& step) const {
    auto curr_steps = node_offset[get_node_idx(id)];
    auto step_idx = step_offset[get_step_idx(step)];
    bool has_prev = step_idx > curr_steps;
    if (has_prev) {
        return std::pair(true, node_steps[step_idx-1]);
    } else {
        auto empty_step = steps.end();
        return std::pair(false, std::pair(empty_step, 0));
    }
}

// get the first step in the path
path_step_index_t::step_it path_step_index_t::path_begin(void) const {
    return steps.begin();
}    

// get the last step in the path
path_step_index_t::step_it path_step_index_t::path_back(void) const {
    return std::prev(std::prev(steps.end()));
}


}
}
