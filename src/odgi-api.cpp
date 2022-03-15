// odgi-api exports the C-API for FFIs
//
// Copyright (c) 2022 Erik Garrison and Pjotr Prins
//
// ODGI is published under the MIT license

#include "odgi-api.h"
#include "version.hpp"


extern "C" {

  using namespace odgi;

  const graph_t* odgi_load_graph(const char *filen) {
    graph_t *graph = new graph_t();
    std::ifstream in(filen);
    graph->deserialize(in);
    return graph;
  }

  void odgi_free_graph(graph_t* graph) {
    delete graph;
  }

  const size_t odgi_get_node_count(const graph_t* graph) {
    return graph->get_node_count();
  }

  const size_t odgi_max_node_id(const graph_t* graph) {
    return graph->max_node_id();
  }

  const size_t odgi_min_node_id(const graph_t* graph) {
    return graph->min_node_id();
  }

  const size_t odgi_get_path_count(const graph_t* graph) {
    return graph->get_path_count();
  }

  const char *odgi_version() {
    return Version::get_version().c_str();
  }

  void odgi_for_each_path_handle(const graph_t *graph, void (*next) (int i, const path_handle_t &path)) {
    int i = 0;
    graph->for_each_path_handle([&](const path_handle_t& path) {
      next(i,path);
      i++;
    });
  }

  const char *odgi_get_path_name(const graph_t *graph, const path_handle_t& path) {
    return graph->get_path_name(path).c_str();
  }
}
