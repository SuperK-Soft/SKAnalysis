#include "Candidate.h"

Candidate::Candidate(const unsigned int iHit)
: hitID(iHit) {}

Candidate::~Candidate() = default;

//std::map<std::string, float> Candidate::GetFeatureMap()
//{
//    std::map<std::string ,float> map;
//
//    for (auto const& pair: storeMap) {
//        map[pair.first] = std::stof(pair.second);
//    }
//
//    return map;
//}
