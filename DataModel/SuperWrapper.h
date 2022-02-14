#ifndef SUPERMANAGERACCESOR_H
#define SUPERMANAGERACCESOR_H
#include "SuperManager.h"

// Robber, for robbing
template<typename Tag, typename Tag::type M>
struct Rob {
  friend typename Tag::type get(Tag) {
    return M;
  }
};

// tag used to access SuperManager::ListOfTreeManagers
struct A_f {
  typedef std::map<int,TreeManager*> SuperManager::*type;
  friend type get(A_f);
};

template struct Rob<A_f, &SuperManager::ListOfTreeManagers>;

TreeManager* GetTreeManager(int LUN);

#endif
