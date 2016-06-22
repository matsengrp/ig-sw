/*Test for checking if data from a fastq file is the same after being parsed.
Place file "shorttest.fastq" in current working directory or change the path in
the
code. */

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "kseq.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <zlib.h>

// declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)

extern "C" void ReadFile() {
  gzFile fp;
  kseq_t *seq;
  int l;
  FILE *writefp;
  writefp = fopen("readtest.txt", "w");
  fp = gzopen("shorttest.fastq", "r"); // open the file handler
  seq = kseq_init(fp);                 // initialize seq
  while ((l = kseq_read(seq)) >= 0) {  // read sequence
    fprintf(writefp, "name: @%s\n", seq->name.s);
    if (seq->comment.l)
      fprintf(writefp, "comment: %s\n", seq->comment.s);
    fprintf(writefp, "seq: %s\n", seq->seq.s);
    if (seq->qual.l)
      fprintf(writefp, "qual: %s\n", seq->qual.s);
  }
  kseq_destroy(seq); // destroy seq
  gzclose(fp);       // close the file handler
  fclose(writefp);
}

std::vector<std::string> OpenFile(std::string filename, bool label) {
  std::ifstream test_file;
  std::vector<std::string> original_file;
  test_file.open(filename.c_str());
  std::string read_file;
  if (test_file.is_open()) {
    while (test_file.good()) {
      test_file >> read_file;
      if (label) {
        // strip labels
        if ("name:" != read_file && "comment:" != read_file &&
            "seq:" != read_file && "qual:" != read_file)
          original_file.push_back(read_file);
      } else {
        if ("+" != read_file)
          original_file.push_back(read_file);
      }

      if (test_file.eof())
        break;
      read_file = "";
    }
  }
  return original_file;
}

namespace igsw {
TEST_CASE("Issue 3: Comparing parsed and original files") {
  ReadFile();
  std::vector<std::string> parsed_vector = OpenFile("readtest.txt", true);
  REQUIRE("@GTACTCTGGTTTGTC|PRCONS=20080924-IGHM|CONSCOUNT=17|DUPCOUNT=9" ==
          parsed_vector[0]);
  REQUIRE("NNNNNNNNNNNNNNNNNNNNNNNTCCTACGCTGGTGAAACCCACACAGACCCTCACGCTGACCTGCAC"
          "GTTCTCTGGGCTCTCACTCAGTACTAGTGAAGTGGCTGTGGGCTGGATCCGTCAGCCCCCAGGAAAGG"
          "CCCGGGAGTGGCTTGCACTCCTTTATTGGAATGATGATAAGTACTATAGTCCATCTCTGAAGAGCAGG"
          "CTCACCATTACTAAGGACACCTCCGAAAATCAGGTGGTCCTTACAATGACCAACATGGACCCTGTGGA"
          "CACAGGCACATATTATTGTGCACACGACGTGACAAGAAGTCGGTATGGGATGGACGTCTGGGGCCAAG"
          "GGACCACGGTCACCGTCTCATCAGGGA" == parsed_vector[1]);
  REQUIRE("!!!!!!!!!!!!!!!!!!!!!!!HFHHFFHCD:AFFFGHHHHHHEEEHHHHHFHHHF@"
          "EECFFHFFHHHHFFHHFFGHHCFFGHHBFHDFHHHFHFHHFHHFFH;;5@=C=,=+4?DEDED@@"
          "DEEEE@AC*;B,;CE:'.2'8888CCAAA0.*A**?EECCC:?AEFEECA::?EA:*0ACE?:C?C?"
          "1*CFFEFCEFFDFFE=:?EEA*EFFFFFCC?A?B8'?<CFFFECFFCECAA?A::CAFFCCFEFEEE?"
          "A*EEAEEFFFFFFFE?AFEFEFFECEEFFEEE?ECAA>;"
          "EFEEBECAEFFEEEEEFDDDBFFEDEEEFFFFHHHHHHHHHHHHFFHHHHIHFC<"
          "HEHGDGHGEGFFFDH" == parsed_vector[2]);
  REQUIRE("@GGTTTGCAATGGTTT|PRCONS=20080924-IGHG|CONSCOUNT=3|DUPCOUNT=2" ==
          parsed_vector[3]);
  REQUIRE("NNNNNNNNNNNNNNNNNNNNNNNGGGAGGCGTGGTCCAGCCTGGGACGTCCCTGAGACTCTCCTGTGC"
          "AGCGTCTGGATTCACCTTCAGTAGTTATGGCATGCACTGGGTCCGCCAGGCTCCAGGCAAGGGGCTGG"
          "AGTGGGTGGCACTCGCATGGTATGATGGAAGTACTGCATACTATACAAACTCCGTGAAGGGCCGATTC"
          "ACCATCTCCAGAGACAATTCCAAGAACACGCTGTCTTTGCAAATGAACAGCCTGACAGCCGAGGACAC"
          "GGCTGTGTATTACTGTGCGAGAGGCCACATCCCCTATGCCTACAATTACCTTTTTGACTACTGGGGCC"
          "AGGGAACCCTGGTCACCGTCTCCTCAGCCC" == parsed_vector[4]);
  REQUIRE("!!!!!!!!!!!!!!!!!!!!!!!"
          "oioooolooiihnmmoooooooomoaloomnhmiklnilmmmlnnnomlkkkomimomkmmmmogmkk"
          "kkkkkkkkiihiiihgiiiiiWfiZdiiicegigigNeii]^biiVV_eeiLi\\eai_"
          "eeKNeiiigmmkgikmkgiikggcJc^kkii^"
          "kcgejammmmmmmmkekmmmkki\\mmkmkmmmmmiimmmmkimmmmmkii^"
          "jkmkkmmmmmmmmmmmmmimmmmmmmkmmmmmmmmmkmkmmmmmmmmmmmmmmmmkmmmmmmmmmmkk"
          "mkmmmmmmmmooomoomokfoooooooooooooooqqponpqqpnloqqqqqqpoooonklpqlooon"
          "oo" == parsed_vector[5]);
  REQUIRE("@GCGTAGCACCCGTGG|PRCONS=20080924-IGHG|CONSCOUNT=2|DUPCOUNT=1" ==
          parsed_vector[6]);
  REQUIRE("NNNNNNNNNNNNNNNNNNNNNNNNCCTACGCTGGTGAAACCCACACAGACCCTCACGCTGACCTGCAC"
          "CTTCTCTGGGTTCTCACTCAGCACTGGTGGAGTGGGTGTGGGCTGGATCCGTCAGCCCCCAGGAAAGG"
          "CCCTGGAGTGGCTTGCACTCATTTATTGGAATGATGATAAGCGCTACAGCCCATCTCTGAAGAGCAGG"
          "CTCACCATCACCAAGGACACCTCCAAAAACCAGGTGGTCCTTACAATGACCAACATGGACCCTGTGGA"
          "CACAGCCACATATTACTGTGCACACAGGGGATACGGTCACCTTGCAGCAACTGGTACTGAGTGGTTCG"
          "ACCCCTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAGCCT" == parsed_vector[7]);
  REQUIRE("!!!!!!!!!!!!!!!!!!!!!!!!HHooloomoeckkkV`"
          "hlloooomVjooiooommhiikooomjnoooooooikmolmjl`"
          "jhmlhklkkmmfihkokkmollihhikii^_hiighii\\gic^f_f_Ycfgfcgiiiiiiieageg["
          "NXLCPJZag3IZXZINH]gcigZeZXM\\aaNkkkimki`iekegiigciecikgSN_"
          "amiicUeikmigegaikggmikmg\\kkkmmmmmkmmmkmkmmikcmki`"
          "ikmmkkkkkgkkmmmkgmkkmmmmmimmkkiimmmmkgmkikmmkiiPmkg]["
          "Y\\\\kkjkmkmmmmmmmmmmmlmmoooolooooooooopqpqpqqpqppoqqqpoookjddhhemml"
          "npolmi" == parsed_vector[8]);
}

TEST_CASE("Trivial pass", "[trivial]") { REQUIRE(1 == 1); }
}
