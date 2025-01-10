#include "layout.h"
#include <cuda.h>
#include <assert.h>
#include "cuda_runtime_api.h"
#include <openrand/philox.h>


#define CUDACHECK(cmd) do {                         \
  cudaError_t err = cmd;                            \
  if (err != cudaSuccess) {                         \
    printf("Failed: Cuda error %s:%d '%s'\n",       \
        __FILE__,__LINE__,cudaGetErrorString(err)); \
    exit(EXIT_FAILURE);                             \
  }                                                 \
} while(0)

#define NCCLCHECK(cmd) do {                         \
  ncclResult_t res = cmd;                           \
  if (res != ncclSuccess) {                         \
    printf("Failed, NCCL error %s:%d '%s'\n",       \
        __FILE__,__LINE__,ncclGetErrorString(res)); \
    exit(EXIT_FAILURE);                             \
  }                                                 \
} while(0)

namespace cuda {


__device__ double compute_zeta(uint32_t n, double theta) {
    double ans = 0.0;
    for (uint32_t i = 1; i <= n; i++) {
        ans += pow(1.0 / double(i), theta);
    }
    return ans;
}

// this function uses the cuda operation __powf, which is a faster but less precise alternative to the pow operation
__device__ uint32_t cuda_rnd_zipf(openrand::Philox &rnd_state, uint32_t n, double theta, double zeta2, double zetan) {
    double alpha = 1.0 / (1.0 - theta);
    double denominator = 1.0 - zeta2 / zetan;
    if (denominator == 0.0) {
        denominator = 1e-9;
    }
    double eta = (1.0 - __powf(2.0 / double(n), 1.0 - theta)) / (denominator);

    // INFO: curand_uniform generates random values between 0.0 (excluded) and 1.0 (included)
    double u = 1.0 - rnd_state.rand<float>();;
    double uz = u * zetan;

    int64_t val = 0;
    if (uz < 1.0) val = 1;
    else if (uz < 1.0 + __powf(0.5, theta)) val = 2;
    else val = 1 + int64_t(double(n) * __powf(eta * u - eta + 1.0, alpha));

    if (val > n) {
        //printf("WARNING: val: %ld, n: %u\n", val, uint32_t(n));
        val--;
    }
    assert(val >= 0);
    assert(val <= n);
    return uint32_t(val);
}


/**
* @brief: update the coordinates of two visualization nodes in the 2D layout space
* This function is called multiple times in one `gpu_layout_kernel` in order to increase the data reuse. 
* Each time, the warp shuffle intrinsics are used to change the selection of node 2 among the 32 threads in the warp. 
* E.g. Iter : Step Pairs Selected would be: 
*     1: (a0, b0), (a1, b1), (a2, b2), ..., (a31, b31)
*     2: (a0, b9), (a1, b0), (a2, b3), ..., (a31, b4)
*     3: (a0, b1), (a1, b4), (a2, b1), ..., (a31, b10)
*     ...
* `b` is randomly chosen from the 32 threads in the warp. 
* @param n1_pos_in_path: position of node 1 in the current selected path
* @param n1_id: id of node 1
* @param n1_offset: offset of node 1
* @param n2_pos_in_path: position of node 2 in the current selected path
* @param n2_id: id of node 2
* @param n2_offset: offset of node 2
* @param eta: an coefficient used in the update formula
* @param node_data: the data structure that stores the coordinates of all nodes
*/
__device__
void update_pos_gpu(int64_t &n1_pos_in_path, uint32_t &n1_id, int &n1_offset,
                    int64_t &n2_pos_in_path, uint32_t &n2_id, int &n2_offset,
                    double eta, 
                    cuda::node_data_t &node_data) {
    double term_dist = std::abs(static_cast<double>(n1_pos_in_path) - static_cast<double>(n2_pos_in_path));

    if (term_dist < 1e-9) {
        term_dist = 1e-9;
    }

    double w_ij = 1.0 / term_dist;

    double mu = eta * w_ij;
    if (mu > 1.0) {
        mu = 1.0;
    }

    float *x1 = &node_data.nodes[n1_id].coords[n1_offset];
    float *x2 = &node_data.nodes[n2_id].coords[n2_offset];
    float *y1 = &node_data.nodes[n1_id].coords[n1_offset + 1];
    float *y2 = &node_data.nodes[n2_id].coords[n2_offset + 1];
    double x1_val = double(*x1);
    double x2_val = double(*x2);
    double y1_val = double(*y1);
    double y2_val = double(*y2);

    double dx = x1_val - x2_val;
    double dy = y1_val - y2_val;

    if (dx == 0.0) {
        dx = 1e-9;
    }

    double mag = sqrt(dx * dx + dy * dy);
    double delta = mu * (mag - term_dist) / 2.0;
    //double delta_abs = std::abs(delta);

    // TODO implement delta max stop functionality
    double r = delta / mag;
    double r_x = r * dx;
    double r_y = r * dy;
    // TODO check current value before updating
    atomicExch(x1, float(x1_val - r_x));
    atomicExch(x2, float(x2_val + r_x));
    atomicExch(y1, float(y1_val - r_y));
    atomicExch(y2, float(y2_val + r_y)); 
}

__global__ 
void gpu_layout_kernel(int iter, cuda::layout_config_t config, double eta, double *zetas, 
                                   cuda::node_data_t node_data, cuda::path_data_t path_data, int sm_count) {
    uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;

    //curandStateCoalesced_t *thread_rnd_state = &rnd_state[smid];
    openrand::Philox rng(tid, iter); // unique stream for each thread at each iteration
    
    // need upto 4 coin flips per thread, this has 32 bits, 32 flips
    const int coin_flips = rng.rand<int>();

    const bool flip1 = coin_flips & 1;
    const bool flip2 = coin_flips & 2;
    const bool flip3 = coin_flips & 4;
    const bool flip4 = coin_flips & 8;

    __shared__ bool cooling[BLOCK_SIZE / WARP_SIZE]; 
    if (threadIdx.x % WARP_SIZE == 1) {
        cooling[threadIdx.x / WARP_SIZE] = (iter >= config.first_cooling_iteration) || flip1;
    }

    // select path
    __shared__ uint32_t first_step_idx[BLOCK_SIZE / WARP_SIZE]; // BLOCK_SIZE/WARP_SIZE = 1024/32 = 32
    // each thread picks its own path
    uint32_t step_idx = rng.uniform<int>(0, path_data.total_path_steps);

    uint32_t path_idx = path_data.element_array[step_idx].pidx;
    path_t p = path_data.paths[path_idx];

    if (p.step_count < 2) {
        return;
    }
    assert(p.step_count > 1);

    // INFO: curand_uniform generates random values between 0.0 (excluded) and 1.0 (included)
    uint32_t s1_idx = rng.uniform<int>(0, p.step_count);
    assert(s1_idx < p.step_count);
    uint32_t s2_idx;

    if (cooling[threadIdx.x / WARP_SIZE]) {
        bool backward;
        uint32_t jump_space;
        if (s1_idx > 0 && flip2 || s1_idx == p.step_count-1) {
            // go backward
            backward = true;
            jump_space = min(config.space, s1_idx);
        } else {
            // go forward
            backward = false;
            jump_space = min(config.space, p.step_count - s1_idx - 1);
        }
        uint32_t space = jump_space;
        if (jump_space > config.space_max) {
            space = config.space_max + (jump_space - config.space_max) / config.space_quantization_step + 1;
        }

        uint32_t z_i = cuda_rnd_zipf(rng, jump_space, config.theta, zetas[2], zetas[space]);

        s2_idx = backward ? s1_idx - z_i : s1_idx + z_i;
    } else {
        do {
            s2_idx = rng.uniform<int>(0, p.step_count);
        } while (s1_idx == s2_idx);
    }
    assert(s1_idx < p.step_count);
    assert(s2_idx < p.step_count);
    assert(s1_idx != s2_idx);


    uint32_t n1_id = p.elements[s1_idx].node_id;
    int64_t n1_pos_in_path = p.elements[s1_idx].pos;
    bool n1_is_rev = (n1_pos_in_path < 0)? true: false;
    n1_pos_in_path = std::abs(n1_pos_in_path);

    uint32_t n2_id = p.elements[s2_idx].node_id;
    int64_t n2_pos_in_path = p.elements[s2_idx].pos;
    bool n2_is_rev = (n2_pos_in_path < 0)? true: false;
    n2_pos_in_path = std::abs(n2_pos_in_path);

    uint32_t n1_seq_length = node_data.nodes[n1_id].seq_length;
    bool n1_use_other_end = flip3 ? true : false;
    if (n1_use_other_end) {
        n1_pos_in_path += uint64_t{n1_seq_length};
        n1_use_other_end = !n1_is_rev;
    } else {
        n1_use_other_end = n1_is_rev;
    }

    uint32_t n2_seq_length = node_data.nodes[n2_id].seq_length;
    bool n2_use_other_end = flip4 ? true : false;
    if (n2_use_other_end) {
        n2_pos_in_path += uint64_t{n2_seq_length};
        n2_use_other_end = !n2_is_rev;
    } else {
        n2_use_other_end = n2_is_rev;
    }

    int n1_offset = n1_use_other_end? 2: 0;
    int n2_offset = n2_use_other_end? 2: 0;

    // Update Coordinates based on the data of selected nodes: n_pos_in_path, n_id, n_offset
    update_pos_gpu(n1_pos_in_path, n1_id, n1_offset, 
                   n2_pos_in_path, n2_id, n2_offset, 
                   eta, node_data);
}


void gpu_layout(layout_config_t config, const odgi::graph_t &graph, std::vector<std::atomic<double>> &X, std::vector<std::atomic<double>> &Y) {


    std::cout << "===== Use GPU to compute odgi-layout =====" << std::endl;
    // get cuda device property, and get the SM count
    cudaDeviceProp prop;
    CUDACHECK(cudaGetDeviceProperties(&prop, 0));
    int sm_count = prop.multiProcessorCount;

    // create eta array
    double *etas;
    cudaMallocManaged(&etas, config.iter_max * sizeof(double));

    const int32_t iter_max = config.iter_max;
    const int32_t iter_with_max_learning_rate = config.iter_with_max_learning_rate;
    const double w_max = 1.0;
    const double eps = config.eps;
    const double eta_max = config.eta_max;
    const double eta_min = eps / w_max;
    const double lambda = log(eta_max / eta_min) / ((double) iter_max - 1);
    for (int32_t i = 0; i < config.iter_max; i++) {
        double eta = eta_max * exp(-lambda * (std::abs(i - iter_with_max_learning_rate)));
        etas[i] = isnan(eta)? eta_min : eta;
    }

    // create node data structure
    // consisting of sequence length and coords
    uint32_t node_count = graph.get_node_count();
    assert(graph.min_node_id() == 1);
    assert(graph.max_node_id() == node_count);
    assert(graph.max_node_id() - graph.min_node_id() + 1 == node_count);

    cuda::node_data_t node_data;
    node_data.node_count = node_count;
    cudaMallocManaged(&node_data.nodes, node_count * sizeof(cuda::node_t));
    for (int node_idx = 0; node_idx < node_count; node_idx++) {
        //assert(graph.has_node(node_idx));
        cuda::node_t *n_tmp = &node_data.nodes[node_idx];

        // sequence length
        const handlegraph::handle_t h = graph.get_handle(node_idx + 1, false);
        // NOTE: unable store orientation (reverse), since this information is path dependent
        n_tmp->seq_length = graph.get_length(h);

        // copy random coordinates
        n_tmp->coords[0] = float(X[node_idx * 2].load());
        n_tmp->coords[1] = float(Y[node_idx * 2].load());
        n_tmp->coords[2] = float(X[node_idx * 2 + 1].load());
        n_tmp->coords[3] = float(Y[node_idx * 2 + 1].load());
    }


    // create path data structure
    uint32_t path_count = graph.get_path_count();
    cuda::path_data_t path_data;
    path_data.path_count = path_count;
    path_data.total_path_steps = 0;
    cudaMallocManaged(&path_data.paths, path_count * sizeof(cuda::path_t));

    vector<odgi::path_handle_t> path_handles{};
    path_handles.reserve(path_count);
    graph.for_each_path_handle(
        [&] (const odgi::path_handle_t& p) {
            path_handles.push_back(p);
            path_data.total_path_steps += graph.get_step_count(p);
        });
    cudaMallocManaged(&path_data.element_array, path_data.total_path_steps * sizeof(path_element_t));

    // get length and starting position of all paths
    uint64_t first_step_counter = 0;
    for (int path_idx = 0; path_idx < path_count; path_idx++) {
        odgi::path_handle_t p = path_handles[path_idx];
        int step_count = graph.get_step_count(p);
        path_data.paths[path_idx].step_count = step_count;
        path_data.paths[path_idx].first_step_in_path = first_step_counter;
        first_step_counter += step_count;
    }

#pragma omp parallel for num_threads(config.nthreads)
    for (int path_idx = 0; path_idx < path_count; path_idx++) {
        odgi::path_handle_t p = path_handles[path_idx];
        //std::cout << graph.get_path_name(p) << ": " << graph.get_step_count(p) << std::endl;

        uint32_t step_count = path_data.paths[path_idx].step_count;
        uint64_t first_step_in_path = path_data.paths[path_idx].first_step_in_path;
        if (step_count == 0) {
            path_data.paths[path_idx].elements = NULL;
        } else {
            path_element_t *cur_path = &path_data.element_array[first_step_in_path];
            path_data.paths[path_idx].elements = cur_path;

            odgi::step_handle_t s = graph.path_begin(p);
            int64_t pos = 1;
            // Iterate through path
            for (int step_idx = 0; step_idx < step_count; step_idx++) {
                odgi::handle_t h = graph.get_handle_of_step(s);
                //std::cout << graph.get_id(h) << std::endl;

                cur_path[step_idx].node_id = graph.get_id(h) - 1;
                cur_path[step_idx].pidx = uint32_t(path_idx);
                // store position negative when handle reverse
                if (graph.get_is_reverse(h)) {
                    cur_path[step_idx].pos = -pos;
                } else {
                    cur_path[step_idx].pos = pos;
                }
                pos += graph.get_length(h);

                // get next step
                if (graph.has_next_step(s)) {
                    s = graph.get_next_step(s);
                } else if (!(step_idx == step_count-1)) {
                    // should never be reached
                    std::cout << "Error: Here should be another step" << std::endl;
                }
            }
        }
    }

    // cache zipf zetas
    auto start_zeta = std::chrono::high_resolution_clock::now();
    double *zetas;
    uint64_t zetas_cnt = ((config.space <= config.space_max)? config.space : (config.space_max + (config.space - config.space_max) / config.space_quantization_step + 1)) + 1;
    cudaMallocManaged(&zetas, zetas_cnt * sizeof(double));
    double zeta_tmp = 0.0;
    for (uint64_t i = 1; i < config.space + 1; i++) {
        zeta_tmp += dirtyzipf::fast_precise_pow(1.0 / i, config.theta);
        if (i <= config.space_max) {
            zetas[i] = zeta_tmp;
        }
        if (i >= config.space_max && (i - config.space_max) % config.space_quantization_step == 0) {
            zetas[config.space_max + 1 + (i - config.space_max) / config.space_quantization_step] = zeta_tmp;
        }
    }
    auto end_zeta = std::chrono::high_resolution_clock::now();
    uint32_t duration_zeta_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_zeta - start_zeta).count();

    const uint64_t block_size = BLOCK_SIZE;
    uint64_t block_nbr = (config.min_term_updates + block_size - 1) / block_size; 

    // curandState_t *rnd_state_tmp;
    // curandStateCoalesced_t *rnd_state;
    // CUDACHECK(cudaMallocManaged(&rnd_state_tmp, sm_count * block_size * sizeof(curandState_t)));
    // CUDACHECK(cudaMallocManaged(&rnd_state, sm_count * sizeof(curandStateCoalesced_t)));
    // cuda_device_init<<<sm_count, block_size>>>(rnd_state_tmp, rnd_state);
    // CUDACHECK(cudaGetLastError());
    // CUDACHECK(cudaDeviceSynchronize());
    // cudaFree(rnd_state_tmp);

    // one curandStateCoalesced_t for each sm, not each block. So several blocks 
    // share the same curandStateCoalesced_t.

    for (int iter = 0; iter < config.iter_max; iter++) {
        gpu_layout_kernel<<<block_nbr, block_size>>>(iter, config, etas[iter], zetas, node_data, path_data, sm_count);
        // check error
        CUDACHECK(cudaGetLastError());
        CUDACHECK(cudaDeviceSynchronize());
    }

    // copy coords back to X, Y vectors
    for (int node_idx = 0; node_idx < node_count; node_idx++) {
        cuda::node_t *n = &(node_data.nodes[node_idx]);
        // coords[0], coords[1], coords[2], coords[3] are stored consecutively. 
        float *coords = n->coords;
        // check if coordinates valid (not NaN or infinite)
        for (int i = 0; i < 4; i++) {
            if (!isfinite(coords[i])) {
                std::cout << "WARNING: invalid coordiate" << std::endl;
            }
        }
        X[node_idx * 2].store(double(coords[0]));
        Y[node_idx * 2].store(double(coords[1]));
        X[node_idx * 2 + 1].store(double(coords[2]));
        Y[node_idx * 2 + 1].store(double(coords[3]));
        //std::cout << "coords of " << node_idx << ": [" << X[node_idx*2] << "; " << Y[node_idx*2] << "] ; [" << X[node_idx*2+1] << "; " << Y[node_idx*2+1] <<"]\n";
    }

    // free memory
    cudaFree(etas);
    cudaFree(node_data.nodes);
    cudaFree(path_data.paths);
    cudaFree(path_data.element_array);
    cudaFree(zetas);
    // cudaFree(rnd_state);

    return;
}

}