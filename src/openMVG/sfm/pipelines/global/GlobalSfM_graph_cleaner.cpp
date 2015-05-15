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

RelativeInfo_Map GlobalSfM_Graph_Cleaner::Run()
{

  FindCycles();
  RotationRejection(5.0f);
      
  return relatives_Rt;

}

////////////////////////////////////////////////////////////////////////////////
//                                Find Cycles                                 //
////////////////////////////////////////////////////////////////////////////////

void GlobalSfM_Graph_Cleaner::FindCycles() {
      
  Pair_Set pair_set;
  for(RelativeInfo_Map::const_iterator iter = relatives_Rt.begin(); iter != relatives_Rt.end(); ++iter)
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
void GlobalSfM_Graph_Cleaner::RotationRejection(const double max_angular_error)
{
    
  const size_t edges_start_count = relatives_Rt.size();  
  RelativeInfo_Map map_relatives_validated;

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
    const float angularErrorDegree = cycle.errToIdentity(relatives_Rt);
    vec_errToIdentityPerCycle.push_back(angularErrorDegree);

    if (angularErrorDegree < max_angular_error)
    {
      vec_cycle_validated.push_back(cycle);
      addCycleToMap( cycle, map_relatives_validated );
    }
  }
  relatives_Rt.swap(map_relatives_validated);

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

  const size_t edges_end_count = relatives_Rt.size();
  std::cout << "\n #Edges removed by triplet inference: " << edges_start_count - edges_end_count << std::endl;
} // function RotationRejection
  




////////////////////////////////////////////////////////////////////////////////
//                               MISCELLANEOUS                                //
////////////////////////////////////////////////////////////////////////////////

void GlobalSfM_Graph_Cleaner::disp_graph(const string str) const {
  for (Graph::NodeIt iter(g); iter!=lemon::INVALID; ++iter){
    Vec3 foo = position_GroundTruth[nodeMapIndex.at(iter)];
    std::cout << "\\node[n" << str << "] at (" << foo(0) <<  ","<< foo(2) << ")" << " (" << nodeMapIndex.at(iter) << ") " << "{" << nodeMapIndex.at(iter) << "}; ";
  }
  std::cout << std::endl;
  for (Graph::ArcIt edge(g); edge!=lemon::INVALID; ++edge){
    IndexT s = nodeMapIndex.at(g.source(edge));    IndexT t = nodeMapIndex.at(g.target(edge));
    if (s < t) {
      std::cout << s << " -> " << t << std::endl;
      std::cout << "\\draw[e" << str << " = " << CohenrencyMap.at(std::make_pair(s,t)) << "] (" << s << ") -- (" << t << "); ";
    }
  }
}

  
} // namespace globalSfM
} // namespace openMVG
