#ifndef GAPBS_SB_COMMAND_LINE_H
#define GAPBS_SB_COMMAND_LINE_H

#include <getopt.h>

#include <algorithm>
#include <cinttypes>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>


/*
GAP Benchmark Suite
Class:  CLBase
Author: Scott Beamer

Handles command line argument parsing
 - Through inheritance, can add more options to object
 - For example, most kernels will use CLApp


Adapted by atrostan on 13/10/22.
*/


class CLBase {
  protected:
  int argc_;
  char **argv_;
  std::string name_;
  std::string get_args_ = "f:g:hk:su:mo:d:p:t:";
  std::vector<std::string> help_strings_;

  int scale_ = -1;
  int degree_ = 16;
  float percent_ = 0.0;// pecrent * num_nodes = k <- slashburn complete when num_nodes in gcc == k
  uint num_threads_ = 1;
  std::string filename_ = "";
  std::string out_filename_ = "";
  std::string db_filename_ = "";
  bool symmetrize_ = false;
  bool uniform_ = false;
  bool in_place_ = false;

  void AddHelpLine(char opt, std::string opt_arg, std::string text,
                   std::string def = "") {
    const int kBufLen = 100;
    char buf[kBufLen];
    if (opt_arg != "")
      opt_arg = "<" + opt_arg + ">";
    if (def != "")
      def = "[" + def + "]";
    snprintf(buf, kBufLen, " -%c %-9s: %-54s%10s", opt, opt_arg.c_str(),
             text.c_str(), def.c_str());
    help_strings_.push_back(buf);
  }

  public:
  CLBase(int argc, char **argv, std::string name = "") : argc_(argc), argv_(argv), name_(name) {
    AddHelpLine('h', "", "print this help message");
    AddHelpLine('f', "file", "load graph from file");
    AddHelpLine('s', "", "symmetrize input edge list", "false");
    AddHelpLine('g', "scale", "generate 2^scale kronecker graph");
    AddHelpLine('u', "scale", "generate 2^scale uniform-random graph");
    AddHelpLine('k', "degree", "average degree for synthetic graph",
                std::to_string(degree_));
    AddHelpLine('m', "", "reduces memory usage during graph building", "false");
    AddHelpLine('o', "output", "location of mapfile to output slashburn order");
    AddHelpLine('d', "db", "location of sqlite db to save preprocessing times");
    AddHelpLine('p', "percent", "pecrent * num_nodes = k <- slashburn complete when num_nodes in gcc == k");
    AddHelpLine('t', "num_threads", "number of threads to use");
  }

  bool ParseArgs() {
    signed char c_opt;
    extern char *optarg;// from and for getopt
    while ((c_opt = getopt(argc_, argv_, get_args_.c_str())) != -1) {
      HandleArg(c_opt, optarg);
    }
    if ((filename_ == "") && (scale_ == -1)) {
      std::cout << "No graph input specified. (Use -h for help)" << std::endl;
      return false;
    }
    if (scale_ != -1)
      symmetrize_ = true;
    return true;
  }

  void virtual HandleArg(signed char opt, char *opt_arg) {
    switch (opt) {
      case 'f':
        filename_ = std::string(opt_arg);
        break;
      case 'g':
        scale_ = atoi(opt_arg);
        break;
      case 'h':
        PrintUsage();
        break;
      case 'k':
        degree_ = atoi(opt_arg);
        break;
      case 's':
        symmetrize_ = true;
        break;
      case 'u':
        uniform_ = true;
        scale_ = atoi(opt_arg);
        break;
      case 'm':
        in_place_ = true;
        break;
      case 'o':
        out_filename_ = std::string(opt_arg);
        break;
      case 'd':
        db_filename_ = std::string(opt_arg);
        break;
      case 'p':
        percent_ = atof(opt_arg);
        break;
      case 't':
        num_threads_ = atoi(opt_arg);
        break;
    }
  }

  void PrintUsage() {
    std::cout << name_ << std::endl;
    // std::sort(help_strings_.begin(), help_strings_.end());
    for (std::string h: help_strings_)
      std::cout << h << std::endl;
    std::exit(0);
  }

  int scale() const { return scale_; }

  int degree() const { return degree_; }

  float percent() const { return percent_; }

  uint num_threads() const { return num_threads_; }

  std::string filename() const { return filename_; }

  std::string out_filename() const { return out_filename_; }

  std::string db_filename() const { return db_filename_; }

  bool symmetrize() const { return symmetrize_; }

  bool uniform() const { return uniform_; }

  bool in_place() const { return in_place_; }
};


class CLApp : public CLBase {
  bool do_analysis_ = false;
  int num_trials_ = 16;
  int64_t start_vertex_ = -1;
  bool do_verify_ = false;
  std::string output_path_;

  public:
  CLApp(int argc, char **argv, std::string name) : CLBase(argc, argv, name) {
    get_args_ += "an:r:v:c:";
    AddHelpLine('a', "", "output analysis of last run", "false");
    AddHelpLine('n', "n", "perform n trials", std::to_string(num_trials_));
    AddHelpLine('r', "node", "start from node r", "rand");
    AddHelpLine('v', "", "verify the output of each run", "false");
    AddHelpLine('c', "c", "write output to this file", "false");
  }

  void HandleArg(signed char opt, char *opt_arg) override {
    switch (opt) {
      case 'a':
        do_analysis_ = true;
        break;
      case 'n':
        num_trials_ = atoi(opt_arg);
        break;
      case 'r':
        start_vertex_ = atol(opt_arg);
        break;
      case 'v':
        do_verify_ = true;
        break;
      case 'c':
        output_path_ = opt_arg;
        break;
      default:
        CLBase::HandleArg(opt, opt_arg);
    }
  }

  bool do_analysis() const { return do_analysis_; }

  int num_trials() const { return num_trials_; }

  int64_t start_vertex() const { return start_vertex_; }

  bool do_verify() const { return do_verify_; }
  std::string output_path() const { return output_path_; }
};


class CLIterApp : public CLApp {
  int num_iters_;

  public:
  CLIterApp(int argc, char **argv, std::string name, int num_iters) : CLApp(argc, argv, name), num_iters_(num_iters) {
    get_args_ += "i:";
    AddHelpLine('i', "i", "perform i iterations", std::to_string(num_iters_));
  }

  void HandleArg(signed char opt, char *opt_arg) override {
    switch (opt) {
      case 'i':
        num_iters_ = atoi(opt_arg);
        break;

      default:
        CLApp::HandleArg(opt, opt_arg);
    }
  }

  int num_iters() const { return num_iters_; }
};


class CLPageRank : public CLApp {
  int max_iters_;
  double tolerance_;
  std::string pr_output_path_;// where to persist the graph's PageRank

  public:
  CLPageRank(int argc, char **argv, std::string name, double tolerance,
             int max_iters) : CLApp(argc, argv, name), max_iters_(max_iters), tolerance_(tolerance) {
    get_args_ += "i:t:c:";
    AddHelpLine('i', "i", "perform at most i iterations",
                std::to_string(max_iters_));
    AddHelpLine('t', "t", "use tolerance t", std::to_string(tolerance_));
    AddHelpLine('c', "c", "save the PageRank values as a binary vector to here", pr_output_path_);
  }

  void HandleArg(signed char opt, char *opt_arg) override {
    switch (opt) {
      case 'i':
        max_iters_ = atoi(opt_arg);
        break;
      case 't':
        tolerance_ = std::stod(opt_arg);
        break;
      case 'c':
        pr_output_path_ = opt_arg;
        break;
      default:
        CLApp::HandleArg(opt, opt_arg);
    }
  }

  int max_iters() const { return max_iters_; }

  double tolerance() const { return tolerance_; }

  std::string pr_output_path() const { return pr_output_path_; }
};


template<typename WeightT_>
class CLDelta : public CLApp {
  WeightT_ delta_ = 1;

  public:
  CLDelta(int argc, char **argv, std::string name) : CLApp(argc, argv, name) {
    get_args_ += "d:";
    AddHelpLine('d', "d", "delta parameter", std::to_string(delta_));
  }

  void HandleArg(signed char opt, char *opt_arg) override {
    switch (opt) {
      case 'd':
        if (std::is_floating_point<WeightT_>::value)
          delta_ = static_cast<WeightT_>(atof(opt_arg));
        else
          delta_ = static_cast<WeightT_>(atol(opt_arg));
        break;
      default:
        CLApp::HandleArg(opt, opt_arg);
    }
  }

  WeightT_ delta() const { return delta_; }
};


class CLConvert : public CLBase {
  std::string out_filename_ = "";
  bool out_weighted_ = false;
  bool out_el_ = false;
  bool out_sg_ = false;

  public:
  CLConvert(int argc, char **argv, std::string name)
      : CLBase(argc, argv, name) {
    get_args_ += "e:b:w";
    AddHelpLine('b', "file", "output serialized graph to file");
    AddHelpLine('e', "file", "output edge list to file");
    AddHelpLine('w', "file", "make output weighted");
  }

  void HandleArg(signed char opt, char *opt_arg) override {
    switch (opt) {
      case 'b':
        out_sg_ = true;
        out_filename_ = std::string(opt_arg);
        break;
      case 'e':
        out_el_ = true;
        out_filename_ = std::string(opt_arg);
        break;
      case 'w':
        out_weighted_ = true;
        break;
      default:
        CLBase::HandleArg(opt, opt_arg);
    }
  }

  std::string out_filename() const { return out_filename_; }

  bool out_weighted() const { return out_weighted_; }

  bool out_el() const { return out_el_; }

  bool out_sg() const { return out_sg_; }
};

#endif//GAPBS_SB_COMMAND_LINE_H
