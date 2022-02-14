#include "SuperWrapper.h"

TreeManager* GetTreeManager(int LUN){
	SuperManager& smgr = *SuperManager::GetManager();
	std::map<int,TreeManager*> mgrs = smgr.*get(A_f());
	if(mgrs.count(LUN)) return mgrs.at(LUN);
	else return nullptr;
};
