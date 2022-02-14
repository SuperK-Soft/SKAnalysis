#ifndef MYTOOLMultiThread_H
#define MYTOOLMultiThread_H

#include <string>
#include <iostream>

#include "Tool.h"

/**
 * \struct MyToolMultiThread_args
 *
 * This is a struct to place data you want your thread to acess or exchange with it. The idea is the datainside is only used by the thread and so will be thread safe
 *
 * $Author: B.Richards $
 * $Date: 2019/05/28 10:44:00 $
 */

struct MyToolMultiThread_args:Thread_args{

  MyToolMultiThread_args();
  ~MyToolMultiThread_args();
  bool busy;
  std::string message;

};

/**
 * \class MyToolMultiThread
 *
 * This is a template for a Tool That employs mulitple threads. Please fill out the descripton and author information.
 *
 * $Author: B.Richards $
 * $Date: 2019/05/28 10:44:00 $
 */

class MyToolMultiThread: public Tool {


 public:

  MyToolMultiThread(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Executre function used to perform Tool perpose. 
  bool Finalise(); ///< Finalise funciton used to clean up resorces.


 private:

  static void Thread(Thread_args* arg); ///< Function to be run by the thread in a loop. Make sure not to block in it
  Utilities* m_util; ///< Pointer to utilities class to help with threading
  std::vector<MyToolMultiThread_args*> args; ///< Vector of thread args (also holds pointers to the threads)

  int m_freethreads; ///< Keeps track of free threads

};


#endif
