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
  
  openMVG::system::Timer timer;
  
  
  int n_views = 20;
  double wrong_edges_p = .1;
  double total_edge_p = 1;
  double min_angular_error = 30;
  double max_noise = 2;  
  double error_thres = 5.;
  int initial_tree_size = 5.;
  
  CmdLine cmd;
  
  cmd.add( make_option('n', n_views, "nombreDeVues") );
  cmd.add( make_option('a', total_edge_p, "nombreAretes") );
  cmd.add( make_option('f', wrong_edges_p, "faussesCorrespondances") );  
  cmd.add( make_option('b', max_noise, "bruit") );
  cmd.add( make_option('e', min_angular_error, "erreurMin") );
  cmd.add( make_option('s', initial_tree_size, "initialTreeSize") );
  cmd.add( make_option('t', error_thres, "seuilErreur") );
  
  cmd.process(argc, argv);
  
  min_angular_error = D2R(min_angular_error);
  max_noise = D2R(max_noise);
  
  std::cout << std::endl
    << "-----------------------------------------------------------\n"
    << "Graph Cleaning test:\n"
    << "-----------------------------------------------------------\n"
    << "Paramètres : "<< std::endl
    << "Nombre de vues : " << n_views << std::endl
    << "Pourcentage d'arêtes : " << 100* total_edge_p << std::endl
    << "Pourcentage de mauvaises correspondances : " << 100* wrong_edges_p << "%" << std::endl
    << "Bruit : " << R2D(max_noise) << std::endl
    << "Erreur angulaire minimum des mauvaises correspondances : " << R2D(min_angular_error) << std::endl
    << "Seuil d'erreur : " << error_thres << std::endl
    << "Taille d'arbre initial : " << initial_tree_size << std::endl
    << "-----------------------------------------------------------"<<std::endl;
    
  
  
  //---------------------------------------
  // Global SfM reconstruction process
  //---------------------------------------
  
  int total_edge = n_views * (n_views - 1) / 2;
  int n_total_edge = int( total_edge_p * total_edge );
  int n_wrong_edges = int( wrong_edges_p * n_total_edge );  // 10%
  
  std::default_random_engine generator;
  std::uniform_real_distribution<double> unif_distrib(0.0,1.0);
    
  RelativeInfo_Map map_relative;
  vector<Mat3> R_global;  R_global.reserve(n_views);
  vector<Vec3> C_global;  C_global.reserve(n_views);
  std::set<Pair> wrong_edges;
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
  int foo;

  for( IndexT i = 0; i < n_views-1; i++ ){
    for( IndexT j = i+1; j < n_views; j++ ){

      const IndexT id_i = i;      const IndexT id_j = j;
      const Pair Pair_ij = std::make_pair( id_i, id_j );
         
      Mat3 Rij = R_global[i] * R_global[j].transpose();     //   Rij = R(j -> i) = Ri * Rj^T
      Vec3 tij = R_global[i] * (C_global[j] - C_global[i]); //   tij = ti - Rij*tj = Ri * ( (-Ci) - (-Cj) )
      foo = rand() % total_edge;
      if ( foo < n_wrong_edges ){
	Mat3 R_error;
	const double theta = min_angular_error + (2 * M_PI - 2 * min_angular_error) * unif_distrib(generator);
	R_error << cos(theta), 0.0, -sin(theta), 0.0, 1.0, 0.0, sin(theta), 0.0, cos(theta);
	Rij = R_error*Rij;
	wrong_edges.insert(std::make_pair(i,j));
	map_relative[Pair_ij] = std::make_pair( Rij, tij );
	n_wrong_edges -= 1;
	n_total_edge -= 1;
      } else if (foo < n_total_edge ) {	
	Mat3 R_error;
	const double theta = 2 * max_noise * unif_distrib(generator) - max_noise;
	R_error << cos(theta), 0.0, -sin(theta), 0.0, 1.0, 0.0, sin(theta), 0.0, cos(theta);
	Rij = R_error*Rij;
	map_relative[Pair_ij] = std::make_pair( Rij, tij );
	n_total_edge -= 1;
      }
      total_edge -= 1;
    }
  }
    
  sfm::GlobalSfM_Graph_Cleaner graph_cleaner(map_relative, error_thres, initial_tree_size);
  graph_cleaner.set_position_groundtruth(C_global);
  graph_cleaner.set_wrong_edges(wrong_edges);
  std::cout << "-----------------------------------------------------------\n---- Cleaning Graph" << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  RelativeInfo_Map old_relatives_Rt = graph_cleaner.run();
  
  std::cout << " Total time took (s): " << timer.elapsed() << std::endl;

  return EXIT_SUCCESS;
}