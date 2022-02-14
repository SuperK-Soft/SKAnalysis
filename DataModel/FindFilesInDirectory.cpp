// -*- mode:c++; tab-width:4; -*-
/* vim:set noexpandtab tabstop=4 wrap */
#include "FindFilesInDirectory.h"
#include <exception>
//#include <regex>  // in g++4.8 (in use on sukap/SK container) std::regex is broken
/*
#include <boost/regex.hpp> // since std::regex doesn't work
#include <boost/regex/pattern_except.hpp>
*/

int FindFilesInDirectory(std::string inputdir, std::string pattern, std::vector<std::string> &matches, bool case_sensitive, int max_subdir_depth, bool use_regex, std::vector<std::string>* filenames, std::vector<std::vector<std::string>>* output_submatches, bool verbose){
	// Scan directory for all files matching a glob pattern or regex
	// optionally strip out just the filenames
	// optionally do regex submatch extraction
	// XXX NOTE: regex escapes \n etc will need to be doubled: \\n!
	// https://cs.brown.edu/~jwicks/boost/libs/regex/doc/introduction.htmls
	// =====================================================================
	
	// find doesn't work with relative paths, so convert to absolute
	if(verbose) std::cout<<"converting directory \""<<inputdir<<"\" to absolute path"<<std::endl;
	std::string absdir = GetStdoutFromCommand(std::string("readlink -f ")+inputdir);
	/// strip off the trailing newline
	absdir.erase(absdir.find_last_not_of(" \t\n\015\014\013")+1);
	if(verbose) std::cout<<"absolute directory \""<<absdir<<"\""<<std::endl;
	
	// build command
	std::string lscommand="";
	std::string case_char = case_sensitive ? " -" : " -i";
	std::string depth_string = 
		(max_subdir_depth>0) ? std::string("") : (std::string(" -maxdepth ")+std::to_string(max_subdir_depth));
	if(!use_regex){
		lscommand = "find " + absdir + depth_string + case_char +"name " + pattern;
	} else {
		lscommand = std::string("find ") + absdir + depth_string + " -regextype egrep "
		+ case_char + "regex '.*?/" + pattern + "'";
	}
	if(verbose) std::cout<<"search command: '"<<lscommand<<"'"<<std::endl;
	
	// find matching files
	std::string fileliststring = GetStdoutFromCommand(lscommand);
	if(verbose) std::cout<<"got filestring of length "<<fileliststring.size()<<" chars"<<std::endl;
	
	// convert returned results (full file paths) to vector
	matches.clear();
	std::stringstream ssl;
	ssl << fileliststring;
	std::string nextfilestring;
	// n.b. if we're extracting submatches we need filenames,
	// so if we weren't given a container for them, make one
	bool local_filenames=false;
	if((not filenames) && output_submatches){
		if(verbose) std::cout<<"allocating local container for basenames"<<std::endl;
		local_filenames=true;
		filenames = new std::vector<std::string>{};
	}
	// ok now do the conversion
	if(verbose) std::cout<<"converting to vector of lines"<<std::endl;
	while(getline(ssl,nextfilestring)){
		matches.push_back(nextfilestring);
		// if requested also populate the vector of stripped filenames
		std::size_t last_char_loc = nextfilestring.find_last_of("/\\");
		std::string fname = nextfilestring.substr(last_char_loc+1);
		if(filenames) filenames->push_back(fname);
	}
	if(verbose) std::cout<<"got "<<matches.size()<<" lines"<<std::endl;
	
	// if the user passed a vector to hold submatches,
	// re-run the regex filenames and retrieve submatches
/*  // std::regex is not working in g++4.8 https://gcc.gnu.org/bugzilla/show_bug.cgi?id=53631
	// i guess i'll leave this here as reference if you're g++4.9 or newer...
	if(output_submatches){
		output_submatches->clear();
		if(verbose) std::cout<<"searching for regex submatches with pattern "<<pattern<<std::endl;
		std::regex theexpression;
		try{
			//std::regex theexpression (pattern.c_str()); // TODO consider std::wregex for match to wide chars?
			if(case_sensitive){
				theexpression.assign(pattern.c_str(),std::regex::extended);
			} else {
				theexpression.assign(pattern.c_str(),std::regex::extended|std::regex_constants::icase);
				// TODO should we use std::regex::perl? boost had issues with extended
			}
		} catch (std::regex_error& e){
			if(constants::regex_err_strings.count(e.code())){
				std::cout<<constants::regex_err_strings.at(e.code())<<std::endl;
			} else {
				std::cout<<"unknown std::refex_error?? code: "<<e.code()<<std::endl;
			}
			continue;
		} catch(std::exception& e){
			std::cerr<<"Warning: FindFilesInDirectory regex matching threw: "<<e.what()<<std::endl;
			continue;
		} catch(...){
			std::cerr<<"Warning: FindFilesInDirectory unhandled regex matching exception!"<<std::endl;
			continue;
		}
		if(verbose) std::cout<<"made the expression..."<<std::endl;
		for(auto&& afile : *filenames){
			// pre-allocat the submatch vector, to maintain synchonization in case of failure
			if(verbose) std::cout<<"adding container for submatches for this line"<<std::endl;
			output_submatches->push_back(std::vector<std::string>{});
			// declare something to catch the submatches
			std::smatch submatches; // TODO consider std::wsmatch to allow wide (japanese) chars?
			if(verbose) std::cout<<"doing regex match on '"<<afile<<"'"<<std::endl;
			// use 'std::regex_search' to allow matching of incomplete sections
			// which can be retrieved via smatch.prefix and .suffix
			try{
				std::regex_match (afile, submatches, theexpression);
			} catch(const std::out_of_range& oor){
				std::cout<<"Warning: FindFilesInDirectory regex matching threw out of range error!"<<std::endl;
				continue;
			} catch (std::regex_error& e){
				if(constants::regex_err_strings.count(e.code())){
					std::cout<<constants::regex_err_strings.at(e.code())<<std::endl;
				} else {
					std::cout<<"unknown std::refex_error?? code: "<<e.code()<<std::endl;
				}
				continue;
			} catch(std::exception& e){
				std::cerr<<"Warning: FindFilesInDirectory regex matching threw: "<<e.what()<<std::endl;
				continue;
			} catch(...){
				std::cerr<<"Warning: FindFilesInDirectory unhandled regex matching exception!"<<std::endl;
				continue;
			}
			std::cout<<"checking validity of match"<<std::endl;
			if(not submatches.ready()){
				std::cerr<<"Error: FindFilesInDirectory submatches not ready after match!"<<std::endl;
				if(verbose) std::cout<<"breaking match scan to avoid descnch"<<std::endl;
				continue;
			}
			if(submatches.empty()){
				std::cerr<<"Warning: FindFilesInDirectory extraction found no submatches!"<<std::endl;
				continue;
			}
			// iterate over matches and put them into a temp vector
			if(verbose) std::cout<<"iterating over submatches"<<std::endl;
			// first submatch (index 0) is the whole match, so skip it
			for(int match_i=1; match_i<submatches.size(); ++match_i){
				output_submatches->back().push_back(submatches.str(match_i));
			}
			if(verbose) std::cout<<"extracted "<<output_submatches->back().size()<<" submatches"<<std::endl;
		}
	}
*/
	
/*
	// try that again with boost::regex
	if(output_submatches){
		output_submatches->clear();
		if(verbose) std::cout<<"searching for regex submatches with pattern "<<pattern<<std::endl;
		bool valid_pattern=true;
		boost::regex theexpression;
		try{
			if(case_sensitive){
				theexpression.assign(pattern,boost::regex::perl);
			} else {
				theexpression.assign(pattern,boost::regex::perl|boost::regex::icase);
			}
		} catch( std::exception& e){
			std::cerr<<"Warning: FindFilesInDirectory regex matching threw: "<<e.what()<<std::endl;
			valid_pattern=false;
		} catch(...){
			std::cerr<<"Warning: FindFilesInDirectory unhandled regex matching exception!"<<std::endl;
			valid_pattern=false;
		}
		
		if(valid_pattern){
			if(verbose) std::cout<<"made the expression..."<<std::endl;
			
			// loop over filenames and extract matches
			for(auto&& afile : *filenames){
				// pre-allocate the submatch vector, to maintain synchonization in case of failure
				if(verbose) std::cout<<"adding container for submatches for this line"<<std::endl;
				output_submatches->push_back(std::vector<std::string>{});
				// declare something to catch the submatches
				boost::smatch submatches; // TODO consider std::wsmatch to allow wide (japanese) chars?
				if(verbose) std::cout<<"doing regex match on '"<<afile<<"'"<<std::endl;
				// use 'boost::regex_search' to allow matching of incomplete sections
				// which can be retrieved via smatch.prefix and .suffix
				try{
					boost::regex_match (afile, submatches, theexpression);
				} catch (boost::regex_error& e){
					if(constants::bregex_err_strings.count(e.code())){
						std::cerr<<constants::bregex_err_strings.at(e.code())<<std::endl;
					} else {
						std::cerr<<"unknown boost::refex_error?? code: "<<e.code()<<std::endl;
					}
					continue;
				} catch(std::exception& e){
					std::cerr<<"Warning: FindFilesInDirectory regex matching threw: "<<e.what()<<std::endl;
					continue;
				} catch(...){
					std::cerr<<"Warning: FindFilesInDirectory unhandled regex matching exception!"<<std::endl;
					continue;
				}
				if(verbose) std::cout<<"checking validity of match"<<std::endl;
				if(submatches.empty()){
					std::cerr<<"Warning: FindFilesInDirectory extraction found no submatches!"<<std::endl;
					continue;
				}
				// iterate over matches and put them into a temp vector
				if(verbose) std::cout<<"iterating over submatches"<<std::endl;
				// first submatch (index 0) is the whole match, so skip it
				for(int match_i=1; match_i<submatches.size(); ++match_i){
					output_submatches->back().push_back(submatches[match_i]);
				}
				if(verbose) std::cout<<"extracted "<<output_submatches->back().size()<<" submatches"<<std::endl;
			}
		}
	}
*/
	
	if(verbose) std::cout<<"cleaning up local basename container"<<std::endl;
	if(local_filenames) delete filenames;
	
	if(verbose) std::cout<<"returning results of "<<matches.size()<<" matching files"<<std::endl;
	return matches.size();
}

