// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// (Z)

#include <cstdlib>
#include <random>

#include "software/SfM/io_regions_type.hpp"

#include "openMVG/system/timer.hpp"
using namespace openMVG;

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include "openMVG/sfm/pipelines/global/sfm_global_engine_relative_motions.hpp"
#include "openMVG/sfm/pipelines/global/GlobalSfM_graph_cleaner.hpp"

int main(int argc, char **argv)
{
  using namespace std;
  std::cout << std::endl
    << "-----------------------------------------------------------\n"
    << "Graph Cleaning test:\n"
    << "-----------------------------------------------------------\n" << std::endl;
 
  openMVG::Timer timer;
  
  //---------------------------------------
  // Global SfM reconstruction process
  //---------------------------------------
  
  std::default_random_engine generator;
  std::uniform_real_distribution<double> unif_distrib(0.0,1.0);
  
  int n_views = 20;
  int total_edge = n_views * (n_views) / 2;
  int n_wrong_edges = int( .3 * total_edge );  // 10%
  double min_angular_error = D2R(30);
  
  RelativeInfo_Map map_relative;
  vector<Mat3> R_global;  R_global.reserve(n_views);
  vector<Vec3> C_global;  C_global.reserve(n_views);
  // -------------------- Test graph computation -------------------- //
  std::cout << "Graph Construction\nNumber Views: " << n_views
	    << "\nWrong correspondences: "<< n_wrong_edges << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;  
  
  // Create global transformation
  for( int i = 0; i < n_views; i++ ){    
    
    Vec3 lookdir, center_pos_i;
    Mat3 Rotation_i;
    const double theta0 = i * 2 * M_PI / n_views;       
    const double theta = 2 * M_PI * unif_distrib(generator);

    center_pos_i << sin(theta0), 0.0, cos(theta0); // Y axis UP
    C_global.push_back(center_pos_i);
    Rotation_i << cos(theta), 0.0, -sin(theta), 0.0, 1.0, 0.0, sin(theta), 0.0, cos(theta);
    R_global.push_back(Rotation_i);  // Y axis UP
  }
  // Build edges
  int edge_left = n_wrong_edges;

  for( int i = 0; i < n_views-1; i++ ){
    for( int j = i+1; j < n_views; j++ ){

      const IndexT id_i = i;      const IndexT id_j = j;
      const Pair Pair_ij = std::make_pair( id_i, id_j );
         
      Mat3 Rij = R_global[i] * R_global[j].transpose();     //   Rij = R(j -> i) = Ri * Rj^T
      Vec3 tij = R_global[i] * (C_global[j] - C_global[i]); //   tij = ti - Rij*tj = Ri * ( (-Ci) - (-Cj) )
      if ( rand() % total_edge < edge_left ){
	Mat3 R_error;
	const double theta = min_angular_error + (2 * M_PI - 2 * min_angular_error) * unif_distrib(generator);
	R_error << cos(theta), 0.0, -sin(theta), 0.0, 1.0, 0.0, sin(theta), 0.0, cos(theta);
	Rij = R_error*Rij;
	edge_left -= 1;
	std::cout << "(" << i << "," << j << "){" << R2D(theta) <<'}'  << std::endl;
      }
      map_relative[Pair_ij] = std::make_pair( Rij, tij );
      total_edge -= 1;
    }
  }
  
  globalSfM::GlobalSfM_Graph_Cleaner graph_cleaner(map_relative);
  graph_cleaner.set_position_groundtruth(C_global);
  std::cout << "-----------------------------------------------------------\nDisplay" << std::endl;
  RelativeInfo_Map old_relatives_Rt = graph_cleaner.Run();
  graph_cleaner.disp_graph("base");

  std::cout << " Total time took (s): " << timer.elapsed() << std::endl;

  return EXIT_SUCCESS;
}