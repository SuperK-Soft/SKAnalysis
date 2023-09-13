#ifndef NEUTRON_INFO_H
#define NEUTRON_INFO_H

struct NeutronInfo {
  double bs_goodness;
  double bs_dirks;
  double bsn50;
  std::array<float, 3> bs_vertex = {};
  std::array<float, 3> bs_dir {};
};

#endif
