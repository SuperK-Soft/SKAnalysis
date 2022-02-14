#include <Utilities.h>

Utilities::Utilities(){ 
  Threads.clear();
}


Thread_args* Utilities::CreateThread(std::string ThreadName,  void (*func)(Thread_args*), Thread_args* args){
  
  if(Threads.count(ThreadName)==0){
    
    if(args==0) args = new Thread_args();  
    
    args->ThreadName=ThreadName;
    args->func=func;
    args->running=true;
    
    pthread_create(&(args->thread), NULL, Utilities::Thread, args);
    
    Threads[ThreadName]=args;
    
}

  else args=0;

  return args;
  
}


void *Utilities::Thread(void *arg){
  
  Thread_args *args = static_cast<Thread_args *>(arg);

  while (!args->kill){
    
    if(args->running) args->func(args );
    else usleep(100);
  
  }
  
  pthread_exit(NULL);

}


bool Utilities::KillThread(Thread_args* &args){
  
  bool ret=false;
  
  if(args){
    
    args->running=false;
    args->kill=true;
    
    pthread_join(args->thread, NULL);
    //delete args;
    //args=0;    
    
  }
  
  return ret;   
  
}

bool Utilities::KillThread(std::string ThreadName){
  
  return KillThread(Threads[ThreadName]);

}

