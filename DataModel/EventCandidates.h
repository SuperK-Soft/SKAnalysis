#ifndef EVENTCANDIDATES_HH
#define EVENTCANDIDATES_HH

#include <memory>

#include "Candidate.h"
#include "Cluster.h"

class EventCandidates : public Cluster<Candidate>
{
    public:
        ~EventCandidates();
        
        void Print();
        void Clear() {
            Cluster::Clear();
            for (std::map<std::string, std::vector<float>*>::iterator pair=featureVectorMap.begin(); pair!=featureVectorMap.end(); ++pair){
                pair->second->clear();
            }
        }
        void FillVectorMap();
        void RegisterFeatureNames(std::vector<std::string> keyList)
        { 
            for (std::vector<std::string>::iterator key=keyList.begin(); key!=keyList.end(); ++key)
                RegisterFeatureName(*key);
        };
        void RegisterFeatureName(const std::string& key) 
        { 
            if (!featureVectorMap.count(key))
                featureVectorMap[key] = new std::vector<float>; 
        }

        std::map<std::string, std::vector<float>*> featureVectorMap;
};

#endif
