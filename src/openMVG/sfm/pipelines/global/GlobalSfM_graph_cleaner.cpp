// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "GlobalSfM_graph_cleaner.hpp"

#include "openMVG/sfm/sfm.hpp"
#include "openMVG/graph/graph.hpp"
#include "openMVG/multiview/rotation_averaging.hpp"
#include "openMVG/stl/stlMap.hpp"

#include "third_party/histogram/histogram.hpp"

namespace openMVG{
namespace globalSfM{
  
using namespace openMVG::rotation_averaging;

////////////////////////////////////////////////////////////////////////////////
//                                    Run                                     //
////////////////////////////////////////////////////////////////////////////////

rotation_averaging::RelativeRotations_map GlobalSfM_Graph_Cleaner::Run() const
{

  FindCycles();
  RotationRejection(5.0f);
  
  return map_relatives;

}

////////////////////////////////////////////////////////////////////////////////
//                                Find Cycles                                 //
////////////////////////////////////////////////////////////////////////////////

void GlobalSfM_Graph_Cleaner::FindCycles() const {
  // vec_cycles
  // mutable rotation_averaging::RelativeRotations_map       map_relatives;

  Pair_Set pair_set;
  for(RelativeRotations_map::const_iterator iter = map_relatives.begin(); iter != map_relatives.end(); ++iter)
      pair_set.insert(iter->first);
  
  std::vector<graphUtils::Triplet> vec_triplets = graphUtils::tripletListing(pair_set);
  
  vec_cycles.clear();
  vec_cycles.reserve(vec_triplets.size());
  
  for ( std::vector<graphUtils::Triplet>::iterator iter=vec_triplets.begin(); iter!=vec_triplets.end(); ++iter ) {
    graphUtils::Triplet triplet_i = *iter;    
    std::vector<Pair> cycle_cst;
    cycle_cst.push_back( std::make_pair(triplet_i.i, triplet_i.j) );
    cycle_cst.push_back( std::make_pair(triplet_i.j, triplet_i.k) );
    cycle_cst.push_back( std::make_pair(triplet_i.k, triplet_i.i) );
    vec_cycles.push_back(Cycle(cycle_cst));
  }
}
  
////////////////////////////////////////////////////////////////////////////////
//                             Rotation Rejection                             //
////////////////////////////////////////////////////////////////////////////////
///  Reject edges of the view graph that do not produce triplets with tiny
///  angular error once rotation composition have been computed.
void GlobalSfM_Graph_Cleaner::RotationRejection(const double max_angular_error) const
{
//  mutable Cycles cycle_vec;
//  mutable rotation_averaging::RelativeRotations_map map_relatives;
  
  const size_t edges_start_count = map_relatives.size();  
  RelativeRotations_map map_relatives_validated;

  //--- -- -  -           Detection of rotation outliers           -  - -- ---//
  
  Cycles vec_cycle_validated;
  
  std::vector<float> vec_errToIdentityPerCycle;
  vec_errToIdentityPerCycle.reserve(vec_cycles.size());
  // Compute for each length 3 cycles: the composition error
  // Error to identity rotation.
  for (size_t i = 0; i < vec_cycles.size(); ++i)
  {
    Cycle & cycle = vec_cycles[i];
    // Compute the consistency error
    const float angularErrorDegree = cycle.errToIdentity(map_relatives);
    vec_errToIdentityPerCycle.push_back(angularErrorDegree);

    if (angularErrorDegree < max_angular_error)
    {
      vec_cycle_validated.push_back(cycle);
      addCycleToMap( cycle, map_relatives_validated );      
    }
  }
  map_relatives.swap(map_relatives_validated);

  //--- -- -  -         Display statistics about cleaning          -  - -- ---//
  
  std::cout << "\nStatistics about rotation triplets:" << std::endl;
  minMaxMeanMedian<float>(vec_errToIdentityPerCycle.begin(),vec_errToIdentityPerCycle.end());

  std::sort(vec_errToIdentityPerCycle.begin(),vec_errToIdentityPerCycle.end());

  Histogram<float> histo(0.0f, *max_element(vec_errToIdentityPerCycle.begin(), vec_errToIdentityPerCycle.end()), 20);
  histo.Add(vec_errToIdentityPerCycle.begin(), vec_errToIdentityPerCycle.end());
  std::cout << histo.ToString() << std::endl;

  {
    std::cout << "\nTriplets filtering based on composition error on unit cycles\n";
    std::cout << "#Triplets before: " << vec_cycles.size() << "\n"
    << "#Triplets after: " << vec_cycle_validated.size() << std::endl;
  }

    vec_cycles.clear();
    vec_cycles = vec_cycle_validated;

  const size_t edges_end_count = map_relatives.size();
  std::cout << "\n #Edges removed by triplet inference: " << edges_start_count - edges_end_count << std::endl;
} // function RotationRejection
  
  
} // namespace globalSfM
} // namespace openMVG
