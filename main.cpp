#include "SlashBurn.h"
#include "benchmark.h"
#include "bitmap.h"
#include "command_line.h"
#include "ips4o.hpp"
#include "spray.hpp"
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <iostream>
#include <omp.h>
#include <random>

int main(int argc, char *argv[]) {
  CLApp cli(argc, argv, "slashburn");
  if (!cli.ParseArgs())
    return -1;
  Builder b(cli);
  Graph g = b.MakeGraph();
  int n_neighbour_rounds = 2;
  float percent = cli.percent();
  Bitmap bmap(g.num_nodes());
  uint32_t num_nodes = g.num_nodes();
  uint num_threads = cli.num_threads();
  SlashBurn sb = SlashBurn(g, n_neighbour_rounds, percent, bmap, num_threads);
  sb.write_permutation(cli.out_filename());
  return 0;
}
