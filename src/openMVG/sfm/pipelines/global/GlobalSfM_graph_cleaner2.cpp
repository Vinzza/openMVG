// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "GlobalSfM_graph_cleaner2.hpp"

#include "openMVG/sfm/sfm.hpp"
#include "openMVG/graph/graph.hpp"
#include "openMVG/multiview/rotation_averaging.hpp"
#include "openMVG/stl/stlMap.hpp"

#include "third_party/histogram/histogram.hpp"

namespace openMVG{
namespace sfm{
  
using namespace openMVG::rotation_averaging;




void Tree::initialise_transformation( const RelativeInfo_Map & relatives_Rt){
  std::list<IndexT> markedIndexT; // (TODO)-> : list<Pair> : the edges we still need to use
  global_transformation.clear();
  // root initialisation
  IndexT root = tree.begin()->first;
  markedIndexT.push_back(root);
  Transformation T_i;
  global_transformation[root] = T_i;
      
  while (!markedIndexT.empty()){
    IndexT s = markedIndexT.front();
    for (Tree::iterator iter = tree.begin(); iter != tree.end(); ++iter) {
      IndexT a = iter->first;	IndexT b = iter->second;
      if (s == a && global_transformation.find(b) == global_transformation.end()){
	
	T_i.load_transformation( relatives_Rt, s, b );
	T_i.compose_right( global_transformation.at(s) );
	global_transformation[b] = T_i;
	markedIndexT.push_back(b);
	
      } else if(s == b && global_transformation.find(a) == global_transformation.end()) {
	
	T_i.load_transformation( relatives_Rt, s, a );
	T_i.compose_right( global_transformation.at(s) );
	global_transformation[a] = T_i;
	markedIndexT.push_back(a);
	
      }
    }
    markedIndexT.pop_front();
  }
}

////////////////////////////////////////////////////////////////////////////////
//                                    Run                                     //
////////////////////////////////////////////////////////////////////////////////

RelativeInfo_Map GlobalSfM_Graph_Cleaner::run()
{

  RelativeInfo_Map new_relatives_Rt;
  return new_relatives_Rt;
  
} // function run

////////////////////////////////////////////////////////////////////////////////
//                              TREE GENERATORS                               //
////////////////////////////////////////////////////////////////////////////////

Tree GlobalSfM_Graph_Cleaner::generate_Random_Tree( const int size ) const{
    if (size > adjacency_map.size()) { throw "Tree to large for the graph."; }
    
    std::vector<IndexT> random_nodes;	random_nodes.reserve(adjacency_map.size());
    int count = 1;

    // Select the first node
    Adjacency_map::const_iterator it = adjacency_map.begin();
    std::advance(it, rand() % adjacency_map.size());
    IndexT r_node = it->first;
    
    Tree tree(r_node);
    
    random_nodes.push_back(r_node);
    while ( count < size ) {

      int r_position = rand() % random_nodes.size();
      r_node = random_nodes[r_position]; // Random node from the tree

      const std::set<IndexT> & adj_set = adjacency_map.at(r_node);
      std::vector<IndexT> adj_vec(adj_set.begin(), adj_set.end()); 

      bool b = true;
      std::random_shuffle(adj_vec.begin(), adj_vec.end());            
      // starting from r_node, we search for a new adjacent node which is not already in tree_nodes
      for ( std::vector<IndexT>::const_iterator iter=adj_vec.begin(); b && iter!=adj_vec.end(); ++iter ) {
	if( not tree.contains(*iter) ){
	  b = false;
	  tree.insert(std::make_pair(r_node, *iter));
	  random_nodes.push_back(*iter);
	  count += 1;	  
	}
      }
      // if (b) {// r_node is only connected to node already in the tree
	// random_nodes.erase(r_position); <-(TODO)
      // }
    }
    return tree;
  } // function generate_Random_Tree
  
////////////////////////////////////////////////////////////////////////////////

double GlobalSfM_Graph_Cleaner::tree_Consistency_Error( Tree & tree, int & nbpos, int & nbneg, int & min_count ) const {
    // Local 'Global' transformation to speedup the computation
    double consistency_error = 0.; nbpos = 0; nbneg = 0;
    Transformation T_i;//    Vec3 t_i;    Mat3 R_i;
    
    tree.initialise_transformation( relatives_Rt );
    
    std::set<IndexT> tree_node = tree.get_indexT_set();
    min_count = tree_node.size();
    
    for(std::set<IndexT>::const_iterator iter = tree_node.begin(); iter != tree_node.end(); ++iter) {
      const std::set<IndexT> adj_set = adjacency_map.at(*iter);
      int count = 0;
      for(std::set<IndexT>::const_iterator foo = adj_set.begin(); foo != adj_set.end(); ++foo) {
	if( tree.contains(*foo) ){
	  count += 1;
	  if( relatives_Rt.find(std::make_pair(*iter,*foo)) != relatives_Rt.end()
			   && not tree.contains(std::make_pair(*iter,*foo))){
	    
	    // ICI ICI ICI ICI ICI
	    
	    T_i = globalTransformation.at(*iter);
	    T_i.compose_left(relatives_Rt.at(std::make_pair(*iter,*foo)), std::make_pair(*iter,*foo));
	    T_i.compose_left_rev( globalTransformation.at(*foo) );
	    const double error = T_i.error();
	    consistency_error += min( error, max_error_thres );
	    if (error < error_thres)
	      nbpos += 1;
	    else
	      nbneg += 1;
	  }
	}
      } // for (foo) in (adj_set)
      min_count = min( min_count, count );
    } // foo (iter) in (tree_node)
   
    if( (nbpos == 0 && nbneg == 0) || min_count <= 1 ){ return default_max_error; }
    return consistency_error;
  } // function tree_Consistency_Error
  
////////////////////////////////////////////////////////////////////////////////

  Tree GlobalSfM_Graph_Cleaner::generate_Consistent_Tree( const int size, double & tree_error ) const {
    const int nb_iter = 10000;
    Tree best_tree;
    double tree_err, err;
    tree_error = 10000.;
    
    int nbpos, nbneg, min_count;    
    for( int foo = 1; foo < nb_iter; foo++ ){
      Tree tree_i = generate_Random_Tree(size);
      tree_err = tree_Consistency_Error( tree_i, nbpos, nbneg, min_count );
      
      std::set<int> consist_set;
      for(Tree::iterator iter = tree_i.begin(); iter != tree_i.end(); ++iter)
	consist_set.insert(get_Pair_Coherence(*iter));
      
      double err = consistency_error_tree( tree_err, nbneg, nbpos, min_count, consist_set.size());
      
      if( err < tree_error ){
	best_tree = tree_i;
	tree_error = err;
	if(dbtalk){std::cout << foo << " - error: " << err << "(" << tree_err << " : " << nbpos << "/" << nbneg << " <" << min_count << "> : " << consist_set.size() << ") : ";
	  for (Tree::const_iterator iter = tree_i.begin(); iter != tree_i.end(); ++iter)
	    std::cout << iter->first << " " << iter->second << " ";
	  std::cout << std::endl;}
      }
    }
    return best_tree;
  } // generate_Consistent_Tree
  
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
  
} // namespace sfm
} // namespace openMVG
