#include "MyToolMultiThread.h"

MyToolMultiThread_args::MyToolMultiThread_args():Thread_args(){}

MyToolMultiThread_args::~MyToolMultiThread_args(){}


MyToolMultiThread::MyToolMultiThread():Tool(){}


bool MyToolMultiThread::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;

  int threadcount=0;
  if(!m_variables.Get("Threads",threadcount)) threadcount=4;

  m_util=new Utilities();

  for(int i=0;i<threadcount;i++){
    MyToolMultiThread_args* tmparg=new MyToolMultiThread_args();   
    tmparg->busy=0;
    tmparg->message="";
    args.push_back(tmparg);
    std::stringstream tmp;
    tmp<<"T"<<i; 
    m_util->CreateThread(tmp.str(), &Thread, args.at(i));
  }
  
  m_freethreads=threadcount;
  
    
  
  return true;
}


bool MyToolMultiThread::Execute(){

  for(unsigned int i=0; i<args.size(); i++){
    if(args.at(i)->busy==0){
      std::cout<<"reply="<<args.at(i)->message<<std::endl;
      args.at(i)->message="Hi";
      args.at(i)->busy=1;
      break;
    }

  }

  m_freethreads=0;
  for(unsigned int i=0; i<args.size(); i++){
    if(args.at(i)->busy==0) m_freethreads++;
  }

  std::cout<<"free threads="<<m_freethreads<<":"<<args.size()<<std::endl;
  
  sleep(1);
  
  return true;
}


bool MyToolMultiThread::Finalise(){

  for(unsigned int i=0;i<args.size();i++) m_util->KillThread(args.at(i));
    
  args.clear();
  
  delete m_util;
  m_util=0;

  return true;
}

void MyToolMultiThread::Thread(Thread_args* arg){

  MyToolMultiThread_args* args=reinterpret_cast<MyToolMultiThread_args*>(arg);

  if(!args->busy) usleep(100);
  else{ 

    args->message="Hello";
    sleep(10);

    args->busy=0;
  }

}
