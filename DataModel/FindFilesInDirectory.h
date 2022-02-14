// -*- mode:c++; tab-width:4; -*-
/* vim:set noexpandtab tabstop=4 wrap */
#ifndef FindFilesInDirectory_H
#define FindFilesInDirectory_H

#include "Algorithms.h"
#include "Constants.h"

#include <iostream>
#include <map>

int FindFilesInDirectory(std::string inputdir, std::string pattern, std::vector<std::string> &matches, bool case_sensitive=false, int max_subdir_depth=0, bool use_regex=false, std::vector<std::string>* filenames=nullptr, std::vector<std::vector<std::string>>* output_submatches=nullptr, bool verbose=false);

#endif
