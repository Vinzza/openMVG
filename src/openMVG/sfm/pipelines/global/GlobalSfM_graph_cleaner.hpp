// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef GLOBALSFM_GRAPH_CLEANER_HPP
#define GLOBALSFM_GRAPH_CLEANER_HPP


namespace openMVG{
namespace globalSfM{


} // namespace globalSfM
} // namespace openMVG

#include "openMVG/sfm/sfm.hpp"
#include "openMVG/graph/graph.hpp"
#include "openMVG/multiview/rotation_averaging_common.hpp"



namespace openMVG{
namespace globalSfM{

////////////////////////////////////////////////////////////////////////////////
//                                   Cycle                                    //
////////////////////////////////////////////////////////////////////////////////
struct Cycle
{
  Cycle(const std::vector<Pair> &  ccycle)
    :cycle(ccycle)
  { }

  std::vector<Pair> cycle;
  

  bool contain(const Pair & edge) const {
    return (std::find(cycle.begin(), cycle.end(), edge) != cycle.end());
  }
  
  void reverse() {
    std::reverse( cycle.begin()+1, cycle.end() );
    for ( std::vector<Pair>::iterator iter = cycle.begin(); iter != cycle.end(); ++iter ){
      Pair & p = *iter;
      p = std::make_pair( p.second, p.first );
    }
  }
  
  float errToIdentity( const rotation_averaging::RelativeRotations_map map_relatives ){
    Mat3 rot_To_Identity = Mat3::Identity();
    
    for ( std::vector<Pair>::iterator iter=cycle.begin(); iter!=cycle.end(); ++iter ){
      const Pair p = *iter;
      Mat3 RIJ;
      if (map_relatives.find(p) != map_relatives.end())
	RIJ = map_relatives.at(p).Rij;
      else
	RIJ = map_relatives.at( make_pair(p.second,p.first) ).Rij.transpose();
      rot_To_Identity = RIJ * rot_To_Identity;
    }
    return static_cast<float>(R2D(getRotationMagnitude(rot_To_Identity)));
  }
  
  ////////// // // /  /    /       /          /       /    /  / // // //////////
  
  friend bool operator==(const Cycle& c1, const Cycle& c2)  {
    std::vector<Pair> cycle2 = c2.cycle;
    Pair edge0 = c1.cycle[0];
    std::vector<Pair>::iterator p = std::find( cycle2.begin(), cycle2.end(), edge0 );
    if ( cycle2.end() != p ) {
      std::rotate( cycle2.begin(), p, cycle2.end() );
      return cycle2 == c1.cycle;
    }
    p = std::find( cycle2.begin(), cycle2.end(), std::make_pair( edge0.second, edge0.first ) );
    if ( cycle2.end() != p ){
      std::rotate( cycle2.begin(), p, cycle2.end() );
      c2.reverse();
      return cycle2 == c1.cycle;
    }
    return false;
  }

  
  friend bool operator!=(const Cycle& c1, const Cycle& c2)  {
    return !(c1 == c2);
  }

  
  friend std::ostream & operator<<(std::ostream & os, const Cycle & c)
  {
    for (int i = 0; i < c.cycle.size(); ++i){
      os << c.cycle[i] << " ";
    }
    os << std::endl;
    return os;
  }
  
}; // Struct Cycle
  
////////////////////////////////////////////////////////////////////////////////
//                          Global SfM Graph Cleaner                          //
//////////////////////////////////////////////////////////////////////////////// 
class GlobalSfM_Graph_Cleaner
{
  void TripletRotationRejection(
    const double max_angular_error,
    std::vector< graphUtils::Triplet > & vec_triplets,
    rotation_averaging::RelativeRotations & relativeRotations
  ) const; 
};


} // namespace globalSfM
} // namespace openMVG

#endif // GLOBALSFM_GRAPH_CLEANER_HPP
