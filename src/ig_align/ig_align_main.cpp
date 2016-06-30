#include "tclap/CmdLine.h"
#include <algorithm>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <string>

extern "C" {
#include "ig_align.h"
}
// Converts string to char array
char *ToCharArray(std::string str) {
  char *c_array = new char[str.length() + 1];
  std::strcpy(c_array, str.c_str());
  return c_array;
}

// Converts string with spaces as delimiters to array of char arrays
void PathArray(char **c_array, std::string str) {
  char *message = ToCharArray(str);
  int i = 0;
  c_array[i] = strtok(message, " ");
  i++;
  while ((c_array[i] = strtok(NULL, " ")) != NULL) {
    i++;
  }
}

int main(int argc, const char *argv[]) {
  try {
    TCLAP::CmdLine cmd("Aligns reads from fasta files", ' ', "1");

    // tclap command line arguments
    TCLAP::ValueArg<std::string> ref_path_opt(
        "r", "ref_path", "The path to the reference file", true, "", "string");
    cmd.add(ref_path_opt);

    TCLAP::ValueArg<int> n_extra_refs_opt("n", "n_extra_refs",
                                          "Number of extra reference files",
                                          false, 0, "int");
    cmd.add(n_extra_refs_opt);

    TCLAP::ValueArg<std::string> extra_ref_paths_opt(
        "x", "extra_ref_paths", "The paths to the extra reference files", false,
        "", "string (space delimiter)");
    cmd.add(extra_ref_paths_opt);

    TCLAP::ValueArg<std::string> qry_path_opt(
        "q", "qry_path", "The path to the query file", true, "", "string");
    cmd.add(qry_path_opt);

    TCLAP::ValueArg<std::string> output_path_opt("o", "output_path",
                                                 "The path to the output file",
                                                 true, "output.sam", "string");
    cmd.add(output_path_opt);

    TCLAP::ValueArg<int> match_opt("m", "match", "Match: default 2", false, 2,
                                   "int");
    cmd.add(match_opt);

    TCLAP::ValueArg<int> mismatch_opt("i", "mismatch", "Mismatch: default 2",
                                      false, 2, "int");
    cmd.add(mismatch_opt);

    TCLAP::ValueArg<int> gap_o_opt("g", "gap_o", "Gap o: default 3", false, 3,
                                   "int");
    cmd.add(gap_o_opt);

    TCLAP::ValueArg<int> gap_e_opt("e", "gap_e", "Gap e: default 1", false, 1,
                                   "int");
    cmd.add(gap_e_opt);

    TCLAP::ValueArg<unsigned> max_drop_opt("d", "max_drop",
                                           "Max drop: default 1000", false,
                                           1000, "int (unsigned)");
    cmd.add(max_drop_opt);

    TCLAP::ValueArg<int> min_score_opt("s", "min_score", "Min score: default 0",
                                       false, 0, "int");
    cmd.add(min_score_opt);

    TCLAP::ValueArg<unsigned> bandwidth_opt("b", "bandwidth",
                                            "Bandwidth: default 150", false,
                                            150, "int (unsigned)");
    cmd.add(bandwidth_opt);

    TCLAP::ValueArg<uint8_t> n_threads_opt(
        "t", "n_threads", "Number of threads: default 1", false, 1, "int");
    cmd.add(n_threads_opt);

    cmd.parse(argc, argv);

    // ref_path
    std::string str_ref_path = ref_path_opt.getValue();
    char *ref_path = ToCharArray(str_ref_path);
    // n_extra_refs
    int n_extra_refs = n_extra_refs_opt.getValue();
    // extra_ref_paths
    std::string paths = extra_ref_paths_opt.getValue();
    char *extra_paths[n_extra_refs + 1];
    PathArray(extra_paths, extra_ref_paths_opt.getValue());
    const char **extra_ref_paths = (const char **)extra_paths;
    // qry_path
    std::string str_qry_path = qry_path_opt.getValue();
    char *qry_path = ToCharArray(str_qry_path);
    // output_path
    std::string str_out_path = output_path_opt.getValue();
    char *output_path = ToCharArray(str_out_path);
    // match
    int match = match_opt.getValue();
    // mismatch
    int mismatch = mismatch_opt.getValue();
    // gap o
    int gap_o = gap_o_opt.getValue();
    // gap e
    int gap_e = gap_e_opt.getValue();
    // max_drop
    unsigned max_drop = max_drop_opt.getValue();
    // min_score
    int min_score = min_score_opt.getValue();
    // bandwidth
    int bandwidth = bandwidth_opt.getValue();
    // n_threads
    uint8_t n_threads = n_threads_opt.getValue();

    ig_align_reads(ref_path, n_extra_refs, extra_ref_paths, qry_path,
                   output_path, match, mismatch, gap_o, gap_e, max_drop,
                   min_score, bandwidth, n_threads, NULL, NULL);

  } catch (TCLAP::ArgException &e) // catch any exception
  {
    std::cerr << "error: " << e.error() << " for arg " << e.argId()
              << std::endl;
  }
  return 0;
}