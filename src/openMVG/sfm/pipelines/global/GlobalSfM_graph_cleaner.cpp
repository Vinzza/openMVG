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
//                         Triplet Rotation Rejection                         //
////////////////////////////////////////////////////////////////////////////////
/// Reject edges of the view graph that do not produce triplets with tiny
///  angular error once rotation composition have been computed.
void GlobalSfM_Graph_Cleaner::TripletRotationRejection(
  const double max_angular_error,
  std::vector< graphUtils::Triplet > & vec_triplets,
  RelativeRotations & relativeRotations) const
{
  const size_t edges_start_count = relativeRotations.size();

  RelativeRotations_map map_relatives = getMap(relativeRotations);
  RelativeRotations_map map_relatives_validated;

  // DETECTION OF ROTATION OUTLIERS
  std::vector< graphUtils::Triplet > vec_triplets_validated;

  std::vector<float> vec_errToIdentityPerTriplet;
  vec_errToIdentityPerTriplet.reserve(vec_triplets.size());
  // Compute for each length 3 cycles: the composition error
  // Error to identity rotation.
  for (size_t i = 0; i < vec_triplets.size(); ++i)
  {
    const graphUtils::Triplet & triplet = vec_triplets[i];
    const IndexT I = triplet.i, J = triplet.j , K = triplet.k;

    //-- Find the three rotations
    const Pair ij(I,J);
    const Pair ji(J,I);

    Mat3 RIJ;
    if (map_relatives.find(ij) != map_relatives.end())
      RIJ = map_relatives.at(ij).Rij;
    else
      RIJ = map_relatives.at(ji).Rij.transpose();

    const Pair jk(J,K);
    const Pair kj(K,J);

    Mat3 RJK;
    if (map_relatives.find(jk) != map_relatives.end())
      RJK = map_relatives.at(jk).Rij;
    else
      RJK = map_relatives.at(kj).Rij.transpose();

    const Pair ki(K,I);
    const Pair ik(I,K);

    Mat3 RKI;
    if (map_relatives.find(ki) != map_relatives.end())
      RKI = map_relatives.at(ki).Rij;
    else
      RKI = map_relatives.at(ik).Rij.transpose();

    const Mat3 Rot_To_Identity = RIJ * RJK * RKI; // motion composition
    const float angularErrorDegree = static_cast<float>(R2D(getRotationMagnitude(Rot_To_Identity)));
    vec_errToIdentityPerTriplet.push_back(angularErrorDegree);

    if (angularErrorDegree < max_angular_error)
    {
      vec_triplets_validated.push_back(triplet);

      if (map_relatives.find(ij) != map_relatives.end())
        map_relatives_validated[ij] = map_relatives.at(ij);
      else
        map_relatives_validated[ji] = map_relatives.at(ji);

      if (map_relatives.find(jk) != map_relatives.end())
        map_relatives_validated[jk] = map_relatives.at(jk);
      else
        map_relatives_validated[kj] = map_relatives.at(kj);

      if (map_relatives.find(ki) != map_relatives.end())
        map_relatives_validated[ki] = map_relatives.at(ki);
      else
        map_relatives_validated[ik] = map_relatives.at(ik);
    }
  }
  map_relatives.swap(map_relatives_validated);

  // update to keep only useful triplets
  relativeRotations.clear();
  relativeRotations.reserve(map_relatives.size());
  std::transform(map_relatives.begin(), map_relatives.end(), std::back_inserter(relativeRotations), std::RetrieveValue());
  std::transform(map_relatives.begin(), map_relatives.end(), std::inserter(used_pairs, used_pairs.begin()), std::RetrieveKey());

  // Display statistics about rotation triplets error:
  std::cout << "\nStatistics about rotation triplets:" << std::endl;
  minMaxMeanMedian<float>(vec_errToIdentityPerTriplet.begin(), vec_errToIdentityPerTriplet.end());

  std::sort(vec_errToIdentityPerTriplet.begin(), vec_errToIdentityPerTriplet.end());

  Histogram<float> histo(0.0f, *max_element(vec_errToIdentityPerTriplet.begin(), vec_errToIdentityPerTriplet.end()), 20);
  histo.Add(vec_errToIdentityPerTriplet.begin(), vec_errToIdentityPerTriplet.end());
  std::cout << histo.ToString() << std::endl;

  {
    std::cout << "\nTriplets filtering based on composition error on unit cycles\n";
    std::cout << "#Triplets before: " << vec_triplets.size() << "\n"
    << "#Triplets after: " << vec_triplets_validated.size() << std::endl;
  }

  vec_triplets.clear();
  vec_triplets = vec_triplets_validated;

  const size_t edges_end_count = relativeRotations.size();
  std::cout << "\n #Edges removed by triplet inference: " << edges_start_count - edges_end_count << std::endl;
}
  
  
} // namespace globalSfM
} // namespace openMVG
