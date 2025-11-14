#ifndef SOLAR_RELIC_H
#define SOLAR_RELIC_H
struct SolarRelic {
        // existing info
        int out_entry_num=0;
        uint32_t nevsk=0;
        uint64_t ticks=0;
        std::array<double,3> vtx;
        // new info
        std::vector<int32_t> matched_ev_nums;
        std::vector<uint64_t> matched_in_entry_nums;
        std::vector<uint64_t> matched_out_entry_nums;
        std::vector<double> matched_tdiffs;
        std::vector<double> matched_dists;
        std::vector<int> rejected_by;
        bool rejected=false;
        bool written=false;
};
#endif
