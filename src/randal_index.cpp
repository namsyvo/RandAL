/* 
 *   Copyright (c) 2012 Memphis-CS Genome Assembly Group
 *   Copyright (c) 2010 Yasuo Tabei
 * 
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 * 
 *   1. Redistributions of source code must retain the above Copyright
 *      notice, this list of conditions and the following disclaimer.
 *
 *   2. Redistributions in binary form must reproduce the above Copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *
 *   3. Neither the name of the authors nor the names of its contributors
 *      may be used to endorse or promote products derived from this
 *      software without specific prior written permission.
*/

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <limits.h>

#include <map>
#include <set>
#include <cmath>
#include <vector>
#include <algorithm>

#include "fmIndex.hpp"

using namespace std;

char *fname, *toname, oname[100], oname1[100];
int percent = 25;

void usage();
void parse_parameters (int argc, char **argv);

void preprocess()
{	
	
	string r   = "";
	string ref = "";
	bool flag = true;
	ifstream infile;
	ofstream outfile;
	infile.open(fname);
	outfile.open(oname);

	while (getline(infile, r))
	{
		if (flag)
		{
			flag = false;
		}
		else
		{			
			outfile<<r;
		}
	}
	infile.close();
	outfile.close();

	ifstream infile1;
	infile1.open(oname);
	getline(infile1, ref);

	reverse(ref.begin(), ref.end());
	
	outfile.open(oname1);
	outfile << ref;
	infile1.close();
	outfile.close();
}

int main(int argc, char **argv) {
  parse_parameters(argc, argv);
  preprocess();
  fmIndex f;
  char index1[100], index2[100];
  strcpy(index1, fname);
  strcat(index1, ".bw");
  strcpy(index2, fname);
  strcat(index2, ".fw");


  f.buildFmIndex(oname, percent);
  ofstream os(index1);
  f.save(os);

  f.buildFmIndex(oname1, percent);
  ofstream os1(index2);
  f.save(os1);

  //f.buildFmIndex(fname, percent);
  //ofstream os(oname);
  //f.save(os);

  return 0;
}

void usage(){
  std::cerr << std::endl
       << "Usage: randal-index [OPTIONS]... INFILE INDEXFILE" << std::endl << std::endl
       << "       where [OPTIONs] is a list of zero or more optional arguments" << std::endl
       << "             INFILE       is the name of an input file" << std::endl
       << "             INDEXFILE    is the name of index-file (optional)" << std::endl
       << "       Additional arguments (input and output files may be specified):"  << std::endl
       << "             -percent [val: 0<= val <= 100]: percentage of sampling suffix array positions" << std::endl
       << "             (default: " <<  percent << ")" << std::endl
       << std::endl;
  exit(0);
}

void parse_parameters (int argc, char **argv){
  if (argc == 1) usage();
  int argno;
  for (argno = 1; argno < argc; argno++){
    if (argv[argno][0] == '-'){
      if (!strcmp (argv[argno], "-percent")) {
	if (argno == argc - 1) std::cerr << "Must specify a float value after -percent" << std::endl;
	percent = atof(argv[++argno]);
	if (percent < 0 || percent > 100) {
	  usage();
	}
      }
      else {
	usage();
      }
    } else {
      break;
    }
  }
  if (argno > argc)
    usage();

  if (argc = 3)
  {
     fname = argv[argno];
     strcpy(oname, fname);
     strcat(oname, ".ref");
     strcpy(oname1, fname);
     strcat(oname1, ".rev");
  }
  else
  {
     toname = argv[argno+1];
     strcpy(oname, toname);
     strcat(oname, ".ref");
     strcpy(oname1, toname);
     strcat(oname1, ".rev");
  }
}
