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
#include "openMVG/sfm/pipelines/global/GlobalSfM_graph_cleaner2.hpp"

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
  int initial_tree_size = 5;
  int ring_size = 0;
  double error_cst = 5;
  int number_tree = 1;
  
  int nb_cluster = 1;
  int min_ring_size = ring_size;
  double ring_ratio = 1;
  double tree_error_thres = .3;
  double max_error_thres = 20;
  
  int max_descent_size = 10000;
  double edge_inconsistency_thres = 10000;
  double pre_edge_inconsistency_thres = 5;
  
  
  double asuma = 1;
  
  CmdLine cmd;
  
  cmd.add( make_option('N', n_views, "nombreDeVues") );
  cmd.add( make_option('a', total_edge_p, "nombreAretes") );
  cmd.add( make_option('f', wrong_edges_p, "faussesCorrespondances") );  
  cmd.add( make_option('b', max_noise, "bruit") );
  cmd.add( make_option('e', min_angular_error, "erreurMin") );
  
  cmd.add( make_option('R', ring_size, "ringSize") );
  cmd.add( make_option('r', min_ring_size, "minRingSize") );
  cmd.add( make_option('n', nb_cluster, "nbCluster") );
  cmd.add( make_option('S', ring_ratio, "ringRatio") );
  
  cmd.add( make_option('s', initial_tree_size, "initialTreeSize") );
  cmd.add( make_option('t', error_thres, "seuilErreur") );
  cmd.add( make_option('c', error_cst, "errorCst") );
  cmd.add( make_option('z', max_descent_size, "maxDescentSize") );
  cmd.add( make_option('y', edge_inconsistency_thres, "edgeInconsistencyThres") );
  cmd.add( make_option('Y', pre_edge_inconsistency_thres, "preEdgeInconsistencyThres") );
  cmd.add( make_option('x', number_tree, "numberTree") );
  
  cmd.add( make_option('T', tree_error_thres, "treeErrorThres") );
  
  
  cmd.add( make_option('A', asuma, "alphaCoef") );
  
  
  
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
    << "Taille de l'anneau : [" << nb_cluster << "x] " << ring_size << " (" << 100*ring_ratio << "%) " << min_ring_size  << " (" << 100*(1-ring_ratio) << "%) " << std::endl
    << std::endl
    << "Seuil d'erreur : " << error_thres << std::endl
    << "Constante d'erreur : " << error_cst << std::endl
    << "Taille d'arbre initial : " << initial_tree_size << std::endl
    << "Descente maximale : " << max_descent_size << std::endl
    << "Seuil d'arête d'incrémentation : " << edge_inconsistency_thres << std::endl
    << "Nombre d'arbres : " << number_tree << std::endl
    << "Asuma Coef : " << asuma << std::endl
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
    //std::cout << i << " :: " << 57.2957795131*theta << " : rotglob(" << R2D(getRotationAngle(R_global[i])) << ")" << std::endl;

  }
  // Build edges
  int foo;

  for( IndexT i = 0; i < n_views-1; i++ ){
    for( IndexT j = i+1; j < n_views; j++ ){

      const IndexT id_i = i;      const IndexT id_j = j;
      const Pair Pair_ij = std::make_pair( id_i, id_j );
         
      Mat3 Rij = R_global[j] * R_global[i].transpose();     //   Rij = R(i -> j) = Rj * Ri^T
      Vec3 tij = R_global[j] * (C_global[i] - C_global[j]); //   tij = tj - Rij*ti = Rj * ( (-Cj) - (-Ci) )
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
	//std::cout << i << "-" << j << " :: " << 57.2957795131*theta << " : rotloc(" << R2D(getRotationAngle(Rij)) << ")" << std::endl;
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
  
  // Ring creation
  int r_count; // width of the current position
  bool c_foo=true;  // true if the current position is in a cluster
  int r_length = int((ring_ratio * double(n_views)) / double(nb_cluster)); // length of the cluster left
  std::cout << "LENGTH = " << ring_ratio * double(n_views) << " / " << nb_cluster << " = " << (ring_ratio * double(n_views)) / double(nb_cluster) << " = " << r_length << std::endl;
  
  if( ring_size >= 0 ){
    for( IndexT i = 0; i < n_views; i++ ){

      r_length = r_length-1;
      
      if(c_foo)
	r_count = min(ring_size,max(r_length,min_ring_size));
      else
	r_count = min_ring_size;
      
      std::cout << i << " : (" << r_length << ") " << r_count;if(c_foo){std::cout<<" BIG"<<std::endl;}else{std::cout<<" SMALL"<<std::endl;}
      
      if( r_length <= 0 ){
	if(c_foo){r_length = int(((1-ring_ratio) * double(n_views)) / double(nb_cluster)); c_foo=false;}
	else{r_length = int((ring_ratio * double(n_views)) / double(nb_cluster)); c_foo=true;}	
      }
      
      for( IndexT k = 1; k < r_count+1; k++ ){
	const int j = (i + k) % n_views;
	const IndexT id_i = i;      const IndexT id_j = j;
	const Pair Pair_ij = std::make_pair( id_i, id_j );
	  
	Mat3 Rij = R_global[j] * R_global[i].transpose();     //   Rij = R(i -> j) = Rj * Ri^T
	Vec3 tij = R_global[j] * (C_global[i] - C_global[j]); //   tij = tj - Rij*ti = Rj * ( (-Cj) - (-Ci) )
	Mat3 R_error;
	const double theta = 2 * max_noise * unif_distrib(generator) - max_noise;
	R_error << cos(theta), 0.0, -sin(theta), 0.0, 1.0, 0.0, sin(theta), 0.0, cos(theta);
	Rij = R_error*Rij;
	map_relative[Pair_ij] = std::make_pair( Rij, tij );
	
	if( wrong_edges.find(std::make_pair(i,j)) != wrong_edges.end() ){ wrong_edges.erase(std::make_pair(i,j));}
	if( wrong_edges.find(std::make_pair(j,i)) != wrong_edges.end() ){ wrong_edges.erase(std::make_pair(j,i));}
	if( map_relative.find(std::make_pair(j,i)) != map_relative.end() ){ map_relative.erase(std::make_pair(j,i));}
      }
    }
  }

  /* / // DEBUG
  for(RelativeInfo_Map::const_iterator iter = map_relative.begin(); iter != map_relative.end(); ++iter) {
    const IndexT id_i = iter->first.first;      const IndexT id_j = iter->first.second;
    Mat3 Rij = iter->second.first;
    Mat Rijb = R_global[id_i] * R_global[id_j].transpose() * Rij;
    std::cout << id_i << "->" << id_j << " : rotglob(" << R2D(getRotationAngle(R_global[id_j] * R_global[id_i].transpose())) << ") : rotloc(" << R2D(getRotationAngle(Rij)) << ") : erreur =  " << R2D(getRotationMagnitude(Rijb)) << " i.e. " << R2D(getRotationAngle(Rijb)) << std::endl;
  }  
  /**/ // DEBUG

  
  
  sfm::GlobalSfM_Graph_Cleaner graph_cleaner(map_relative, error_thres, initial_tree_size);
    // graph_cleaner.set_error_valid_thres( val );
    // graph_cleaner.set_default_max_error( val );
    // graph_cleaner.set_nb_steps_skip( val );
    // graph_cleaner.set_max_error_thres( val );
    graph_cleaner.set_error_cst( error_cst );
    graph_cleaner.set_max_descent_size( max_descent_size );
    graph_cleaner.set_edge_inconsistency_thres( edge_inconsistency_thres );
    graph_cleaner.set_number_tree( number_tree );
    graph_cleaner.set_tree_error_thres( tree_error_thres );
    graph_cleaner.set_pre_edge_inconsistency_thres( pre_edge_inconsistency_thres );
    
    double foop = error_thres/(2*3.1415);
    
    cout.precision(3);
    std::cout << "\n\n\nTABLEAUX (=> Rapport)";
    for( int nbneg = 0; nbneg<=20; nbneg++ ){
      std::cout << nbneg << "&";
      for( int nbpos = 0; nbpos<=10; nbpos++ ){
//	double foo = initial_tree_size * ( 1 + ((max_error_thres*nbneg + (error_thres/2)*nbpos)/error_cst) ) * ( max(nbneg-nbpos,0) + 1 )/( (nbpos + 1)*(nbpos + .01) );
//	double foo = initial_tree_size * ( 1 + ((max_error_thres*nbneg + (error_thres/2)*nbpos)/error_cst) ) * ( nbneg + 1 )/( (nbpos + 1)*(nbpos + .01) );
	double bar = max_error_thres*nbneg + error_thres*nbpos;
	double bar2 = error_thres*nbneg;
	double foo = graph_cleaner.consistency_error_edge(bar, nbneg, nbpos);
	double foo2 = graph_cleaner.consistency_error_edge(bar2, nbneg, nbpos);
	
	if( foo < .04 ) { std::cout << "\\textcolor{red!25}{$" << foo << "/" << foo2 << "$}";}
	else if( foo < .2 ) { std::cout << "\\textcolor{red!50}{$" << foo << "/" << foo2 << "$}";}
	else if( foo < 1 ) { std::cout << "\\textcolor{red}{$" << foo << "/" << foo2 << "$}";}
	else if( foo < 5 ){ std::cout << "\\textcolor{red!50!black}{$" << foo << "/" << foo2 << "$}";}
	else if( foo < 25 ){ std::cout << "\\textcolor{red!25!black}{$" << foo << "/" << foo2 << "$}";}
	else{ std::cout << "$" << foo << "/" << foo2 << "$";}
	if(nbpos<10){std::cout << "&";} else {std::cout << "\\\\";}

      }
      std::cout << "\n";
    }
    cout.precision(6);
  // err = ( 2 / double(min_count) ) * size * ( 1 + (tree_err/error_cst) ) * ( max(nbneg-nbpos,0) + 1 )/( (nbpos + 1)*(nbpos + .01) );
    
  graph_cleaner.set_position_groundtruth(C_global);
  graph_cleaner.set_wrong_edges(wrong_edges);
  std::cout << "-----------------------------------------------------------\n---- Cleaning Graph" << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  RelativeInfo_Map old_relatives_Rt = graph_cleaner.run();
  
  std::cout << " Total time took (s): " << timer.elapsed() << std::endl;

  return EXIT_SUCCESS;
}