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
namespace sfm{
  
using namespace openMVG::rotation_averaging;




void Tree::initialise_transformation( const RelativeInfo_Map & relatives_Rt){
  std::list<IndexT> markedIndexT; // (TODO)-> : list<Pair> : the edges we still need to use
  global_transformation.clear();
  // root initialisation
  IndexT root = *indexT_set.begin();
  markedIndexT.push_back(root);
  Transformation T_i;
  global_transformation[root] = T_i;
      
  while (!markedIndexT.empty()){
    IndexT s = markedIndexT.front();
    for (std::set<Pair>::iterator iter = pair_set.begin(); iter != pair_set.end(); ++iter) {
      IndexT a = iter->first;	IndexT b = iter->second;
      if (s == a && global_transformation.find(b) == global_transformation.end()){
	T_i.load_transformation( relatives_Rt, a, b );
	T_i.compose_right( global_transformation.at(a) );
	global_transformation[b] = T_i;
	markedIndexT.push_back(b);
	
      } else if(s == b && global_transformation.find(a) == global_transformation.end()) {
	
	T_i.load_transformation( relatives_Rt, b, a );
	T_i.compose_right( global_transformation.at(b) );
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

  double tree_error;
  std::cout << "------------\nConstruction of the initial..." << std::endl;
  Tree coherent_tree = generate_Consistent_Tree(initial_tree_size, tree_error);
  
  std::set<Pair> foo = coherent_tree.get_pair_set();
  for (std::set<Pair>::const_iterator iter = foo.begin(); iter != foo.end(); ++iter)
    tikzfile << "\\draw[itree](" << iter->first << ")--(" << iter->second << ");"<<std::endl;


  std::cout << "------------\nEnlargement of the tree..." << std::endl;
  increase_Tree(coherent_tree);
  
  foo = coherent_tree.get_pair_set();
  for (std::set<Pair>::const_iterator iter = foo.begin(); iter != foo.end(); ++iter)
    tikzfile << "\\draw[tree](" << iter->first << ")--(" << iter->second << ");"<<std::endl;
  
  
  
  tikzfile << "\n\\def\\componentchoice{0}\n";
  disp_Graph("base");
  
  
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
      for(std::set<IndexT>::const_iterator iter_t = adj_set.begin(); iter_t != adj_set.end(); ++iter_t) {
	if( tree.contains(*iter_t) ){
	  count += 1;
	  if( relatives_Rt.find(std::make_pair(*iter,*iter_t)) != relatives_Rt.end()
			   && not tree.contains(std::make_pair(*iter,*iter_t))){
	    Pair foo = std::make_pair(*iter,*iter_t);
	    const double error = tree.transformation_error( Transformation(relatives_Rt.at(foo), foo) );

	    consistency_error += min( error, proba_error_negativethres );
	    if (error < proba_error_positivethres)
	      nbpos += 1;
	    else
	      nbneg += 1;
	  }
	} // if (tree) contains (iter_t)
      } // for (iter_t) in (adj_set)
      min_count = min( min_count, count );
    } // for (iter) in (tree_node)
   
    if( (nbpos == 0 && nbneg == 0) || min_count <= 1 ){ return default_max_error; }
    
    consistency_error = ( 1 + (consistency_error/proba_error_negativethres) )
	    * ( max(nbneg-nbpos,0) + 1 )/( (nbpos + 1)*(nbpos + .01) )
	    * ( 2 / double(min_count) );
    return consistency_error;
  } // function tree_Consistency_Error
  
////////////////////////////////////////////////////////////////////////////////

  Tree GlobalSfM_Graph_Cleaner::generate_Consistent_Tree( const int size, double & best_tree_error ) const {
    Tree best_tree;
    double tree_err;
    best_tree_error = default_max_error;
    
    int nbpos, nbneg, min_count;    
    for( int foo = 1; foo < initial_tree_ransac_nb; foo++ ){
      
      Tree tree_i = generate_Random_Tree(size);
      tree_err = tree_Consistency_Error( tree_i, nbpos, nbneg, min_count );  

      if( tree_err < best_tree_error ){
	best_tree = tree_i;
	best_tree_error = tree_err;
	if(dbtalk>0){std::cout << foo << " - error: " << tree_err << "(" << nbpos << "/" << nbneg << " : " << min_count << ") : ";
	  tree_i.print();
	  std::cout << std::endl;
	}
      }
    }
    return best_tree;
  } // generate_Consistent_Tree

////////////////////////////////////////////////////////////////////////////////
//                              TREE INCREASING                               //
////////////////////////////////////////////////////////////////////////////////

  std::map<IndexT,int> GlobalSfM_Graph_Cleaner::compute_distance(
		const Pair & pair,
		const std::set<IndexT> & tree_node_set) const {
    std::map< IndexT, int > distance;
    std::list<IndexT> markedIndexT;
    
    for(std::set<IndexT>::const_iterator iter = tree_node_set.begin(); iter != tree_node_set.end(); ++iter) {
      distance[*iter] = 0;
      markedIndexT.push_back(*iter);
    }

    while (!markedIndexT.empty()){
      IndexT s = markedIndexT.front();
      const std::set<IndexT> adj_s = adjacency_map.at(s);      
      for(std::set<IndexT>::const_iterator iter = adj_s.begin(); iter != adj_s.end(); ++iter) {
	IndexT t = *iter;
	if( distance.find(t) == distance.end()
	 && (pair.first != s || pair.second != t)
	 && (pair.first != t || pair.second != s) ) {
	  distance[t] = distance.at(s)+1;
	  markedIndexT.push_back(t);
	}
      }      
      markedIndexT.pop_front();
    }
    return distance;
}

////////////////////////////////////////////////////////////////////////////////

double GlobalSfM_Graph_Cleaner::edge_consistency_probability(
	      const IndexT & source_node, const IndexT & target_node,
	      Transformation T,
	      const Tree & tree,
	      const std::map< IndexT, int > & distance,
	      int & nb_total, bool & edge_skip, double & nbpos_s, double & nb_s) const {

  if(dbtalk>0){std::cout << "\n" << target_node << "[" << distance.at(target_node) << "]"; } /////////////////////////////////

  if( nb_total > stop_descent_number ){ edge_skip = true; return proba_neutral; }
  
  T.compose_left( Transformation( relatives_Rt, source_node, target_node ) );
  
  if( tree.contains(target_node) ){
    if(dbtalk>0){std::cout << "("; } /////////////////////////////////
    nb_total += 1;
    const double foo = cycle_probability( tree.transformation_error( T ) );
    if(dbtalk>0){std::cout << "path="; T.print_path(); std::cout << " : proba=" << foo << ")"; std::cout.flush();}
    return foo;
  }
  
  double nb = 0, nbpos = 0;
  
  const std::set<IndexT> adj_t = adjacency_map.at(target_node);
  if(dbtalk>0){std::cout << "<"; } /////////////////////////////////
  for( std::set<IndexT>::const_iterator iter = adj_t.begin(); iter != adj_t.end(); ++iter ) {      
    if( *iter != source_node && distance.at(*iter) < distance.at(target_node) ){
      nb += 1;
      nbpos += edge_consistency_probability( target_node, *iter, T, tree, distance, nb_total, edge_skip, nbpos_s, nb_s );
    }
  }
  
  // double error = pow(1+formular_asuma_coef/(1-proba_error_pos_thres_over2pi),nb-nbpos)
  //             / pow(1+formular_asuma_coef/proba_error_pos_thres_over2pi,nbpos);
  
  const double error = formule_asuma( nbpos, nb-nbpos );
  
  if(dbtalk>0){std::cout << "\n  (" << nbpos << "/" << nb << ") p=" << error << ">"; } /////////////////////////////////
  
  nbpos_s = nbpos;  nb_s = nb;
  
  return error; //1/(1+error);

}

////////////////////////////////////////////////////////////////////////////////

void GlobalSfM_Graph_Cleaner::increase_Tree( Tree & tree ) const {
    int nb_total; double nbpos, nb;
    double proba, max_proba = 0;
    Pair best_pair;
    bool edge_skip;
    IndexT new_node, source_node;

    // Set the initial set of edges which could be add to the graph.
    std::map<Pair, int> edge_indic;
    for(RelativeInfo_Map::const_iterator iter = relatives_Rt.begin(); iter != relatives_Rt.end(); ++iter) {
      Pair pair = iter->first;
      IndexT s = pair.first;
      IndexT t = pair.second;
      if( ( not tree.contains(s) && tree.contains(t) )
       || ( tree.contains(s) && not tree.contains(t) ) ){
	edge_indic[ pair ] = 0;
      }
    }

    while( edge_indic.size() > 0 ) {
      max_proba = 0;

      if(dbtalk>0){std::cout << "\n\nTree_nodes (" << tree.size() << ") = "; tree.print_node(); std::cout << std::endl;}
      
      for(std::map<Pair, int>::iterator iter = edge_indic.begin(); iter != edge_indic.end(); ++iter) {
	Pair pair = iter->first;


	// Edge selection
	if(tree.contains(pair.first) && not tree.contains(pair.second))
	  {new_node = pair.second;	source_node = pair.first;}
	else if (tree.contains(pair.second) && not tree.contains(pair.first))
	  {new_node = pair.first;	source_node = pair.second;}
	else { edge_indic.erase(pair); continue;}

	if( iter->second > 0 ) { iter->second -= 1; continue; }
		
	  
	edge_skip = false;
	nb_total = 0;
	if(dbtalk>0){std::cout << "(" << source_node << "-" << new_node << "), "; } /////////////////////////////////
	
	std::map<IndexT,int> distance = compute_distance( best_pair, tree.get_indexT_set() );
	
	dbtalk--;
	proba = edge_consistency_probability( source_node, new_node, Transformation(), tree, distance, nb_total, edge_skip, nbpos, nb );
	dbtalk++;
	
	if( edge_skip ) { iter->second = 10; } // The computation was to slow
	
	if ( proba > max_proba ){ // We keep the extremum
	  if(dbtalk>0){std::cout << "\n(" << source_node << "-" << new_node << ") = " << nbpos << "/" << nb-nbpos << " => " << proba << std::endl; } /////////////////////////////////
	  max_proba = proba;
	  best_pair = std::make_pair(source_node, new_node);
	
	  /* If the edge is to good, we take it without checking the other edges
	    if ( error < error_valid_thres ){
	      if(dbtalk){std::cout << "Break" << std::endl;} // DBTALK
	      break;
	    }*/
	  } // if

      } // for
      
      if (max_proba < tree_increment_stop_thres){ break; }
	
	tree.insert(best_pair, Transformation( relatives_Rt, best_pair.first, best_pair.second ) );
	
	// We compute the distance to the tree, to ignore the non biconnexe vertexes
	std::map<IndexT,int> distance = compute_distance( best_pair, tree.get_indexT_set() );
	const std::set<IndexT> adj_newnode = adjacency_map.at(best_pair.second);
	
	    if(dbtalk>0){std::cout << "\nEdge: " << best_pair.first << "->" << best_pair.second << " : new edges = "; }
	
	for(std::set<IndexT>::const_iterator iter = adj_newnode.begin(); iter != adj_newnode.end(); ++iter) {
	  if( distance.find(*iter) != distance.end() && distance.at(*iter) != 0 ){ // We dont want to consider the evaluated pair
	    edge_indic[std::make_pair(best_pair.second, *iter)] = 0;
	    if(dbtalk>0){std::cout << " " << best_pair.second << "-" << *iter << "[" << distance.at(*iter) << "]"; }
	  }
	}
	
    } // while
    
  } // function increase_Tree


  
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
  
} // namespace sfm
} // namespace openMVG
