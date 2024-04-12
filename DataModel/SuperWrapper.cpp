#include "SuperWrapper.h"

TreeManager* GetTreeManager(int LUN){
	SuperManager& smgr = *SuperManager::GetManager();
	std::map<int,TreeManager*> mgrs = smgr.*get(A_f());
	if(mgrs.count(LUN)) return mgrs.at(LUN);
	else return nullptr;
};

int RemoveDuplicateLuns(){
	SuperManager& smgr = *SuperManager::GetManager();
	std::map<int,TreeManager*>& mgrs = smgr.*get(A_f());
	std::unordered_set<TreeManager*> unique_managers;
	std::vector<int> duplicate_luns;
	for(std::pair<const int, TreeManager*>& mgr : mgrs){
		if(unique_managers.count(mgr.second)){
			// this TreeManager is already associated with another LUN
			// mark the duplicate LUN for removal
			duplicate_luns.push_back(mgr.first);
		} else {
			unique_managers.emplace(mgr.second);
		}
	}
	// remove the duplicates so we don't have double frees when the SuperManager destructor runs
	for(auto&& adup : duplicate_luns){
		mgrs.erase(adup);
	}
	return duplicate_luns.size();
}

bool AddDuplicateLun(int existinglun, int newlun){
	// add newlun as a duplicate handle to the same file as existinglun
	SuperManager& smgr = *SuperManager::GetManager();
	std::map<int,TreeManager*>& mgrs = smgr.*get(A_f());
	if(mgrs.count(existinglun)==0){
		std::cerr<<"Lun "<<existinglun<<" not found, cannot add duplicate"<<std::endl;
		return false;
	}
	mgrs.emplace(newlun, mgrs.at(existinglun));
	return true;
}

bool EraseLun(int lun){
	// erase a lun from the SuperManager
	// bear in mind this will mean skroot_* functions that take a lun
	// will not be able to find it! Make sure you delete the pointed-to TreeManager
	// if it's not a duplicate to prevent leaking memory
	SuperManager& smgr = *SuperManager::GetManager();
	std::map<int,TreeManager*>& mgrs = smgr.*get(A_f());
	if(mgrs.count(lun)==0){
		std::cerr<<"Lun "<<lun<<" not found"<<std::endl;
		return false;
	}
	mgrs.erase(lun);
	return true;
}
