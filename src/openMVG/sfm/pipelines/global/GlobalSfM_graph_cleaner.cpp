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

RelativeInfo_Map GlobalSfM_Graph_Cleaner::run()
{
  
  std::cout << "Construction of the initial..." << std::endl;
  Tree coherent_tree = generate_Consistent_Tree(7);
  
  if(dbtalk){
  for (globalSfM::Tree::const_iterator iter = coherent_tree.begin(); iter != coherent_tree.end(); ++iter)
    std::cout << "\\draw[itree](" << iter->first << ")--(" << iter->second << ");"<<std::endl;}
  std::cout << "Enlargement of the tree..." << std::endl;
  increase_Tree( coherent_tree );  
  update_Coherence( coherent_tree );
  
  std::cout << "-----------------------------------------------------------\nPrint Tree" << std::endl;
  for (globalSfM::Tree::const_iterator iter = coherent_tree.begin(); iter != coherent_tree.end(); ++iter)
    std::cout << "\\draw[tree](" << iter->first << ")--(" << iter->second << "); ";
  std::cout << std::endl;
  
  /*
  std::cout << "-----------------------------------------------------------\nTriplets Inference" << std::endl;
  Cycles cycles = findCycles();
  rotationRejection(5.0f, cycles); */
  return relatives_Rt;
  
} // function run

////////////////////////////////////////////////////////////////////////////////
//                                Find Cycles                                 //
////////////////////////////////////////////////////////////////////////////////

Cycles GlobalSfM_Graph_Cleaner::findCycles() {
  Pair_Set pair_set;
  for(RelativeInfo_Map::const_iterator iter = relatives_Rt.begin(); iter != relatives_Rt.end(); ++iter)
      pair_set.insert(iter->first);
  
  std::vector<graphUtils::Triplet> vec_triplets = graphUtils::tripletListing(pair_set);
  
  Cycles vec_cycles;
  vec_cycles.reserve(vec_triplets.size());
  
  for ( std::vector<graphUtils::Triplet>::iterator iter=vec_triplets.begin(); iter!=vec_triplets.end(); ++iter ) {
    graphUtils::Triplet triplet_i = *iter;    
    std::vector<Pair> cycle_cst;
    cycle_cst.push_back( std::make_pair(triplet_i.i, triplet_i.j) );
    cycle_cst.push_back( std::make_pair(triplet_i.j, triplet_i.k) );
    cycle_cst.push_back( std::make_pair(triplet_i.k, triplet_i.i) );
    vec_cycles.push_back(Cycle(cycle_cst));
  }
  return vec_cycles;
} // function findCycles
  
////////////////////////////////////////////////////////////////////////////////
//                             Rotation Rejection                             //
////////////////////////////////////////////////////////////////////////////////
///  Reject edges of the view graph that do not produce triplets with tiny
///  angular error once rotation composition have been computed.
void GlobalSfM_Graph_Cleaner::rotationRejection(const double max_angular_error,
						Cycles & vec_cycles)
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
      for (int i = 0; i < cycle.cycle.size()-1; ++i) {
	change_Cohenrence(cycle.cycle[i], cycle.cycle[i+1]);
      }
    }
  }
  //relatives_Rt.swap(map_relatives_validated);

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
  const size_t edges_end_count = map_relatives_validated.size();
  std::cout << "\n #Edges removed by triplet inference: " << edges_start_count - edges_end_count << std::endl;
} // function rotationRejection
  


////////////////////////////////////////////////////////////////////////////////
//                               MISCELLANEOUS                                //
////////////////////////////////////////////////////////////////////////////////

  void GlobalSfM_Graph_Cleaner::disp_Graph(const string str) const {
    for (Graph::NodeIt iter(g); iter!=lemon::INVALID; ++iter){
      Vec3 foo = position_GroundTruth[nodeMapIndex.at(iter)];
      std::cout << "\\node[n" << str << "] at (" << foo(0) <<  ","<< foo(2) << ")" << " (" << nodeMapIndex.at(iter) << ") " << "{" << nodeMapIndex.at(iter) << "}; ";
    }
    std::cout << std::endl;
    std::cout << "-----------------------------------------------------------\n"
              << "TikZ Graph Display\n" 
	      << "-----------------------------------------------------------" << std::endl;
    
    for(Adjacency_map::const_iterator iter = adjacency_map.begin();
	iter != adjacency_map.end(); ++iter) {
      IndexT s = iter->first;
      const std::set<IndexT> & indexT_set = iter->second;
    
      for (std::set<IndexT>::const_iterator iterT = indexT_set.begin();
	  iterT != indexT_set.end(); ++iterT) {
	IndexT t = *iterT;
	if (s < t) {
	  std::cout << "\\zdraw["<< str << "]";
	  if (CohenrenceMap.find(std::make_pair(s,t)) != CohenrenceMap.end())
	    std::cout << "[" << CohenrenceMap.at(std::make_pair(s,t)) << "]";
	  else
	    std::cout << "[" << CohenrenceMap.at(std::make_pair(t,s)) << "]";
	  if (wrong_edges.find(std::make_pair(s,t)) != wrong_edges.end() || wrong_edges.find(std::make_pair(t,s)) != wrong_edges.end() )
	    std::cout << "[0]";
	  else
	    std::cout << "[1]";
	  std::cout << " (" << s << ") -- (" << t << "); ";
	}
      }
    }
    /*
    for (Graph::ArcIt edge(g); edge!=lemon::INVALID; ++edge){
      IndexT s = nodeMapIndex.at(g.source(edge));    IndexT t = nodeMapIndex.at(g.target(edge));
      if (s < t) {
	std::cout << "\\zdraw["<< str << "]";
	if (CohenrenceMap.find(std::make_pair(s,t)) != CohenrenceMap.end())
	  std::cout << "[" << CohenrenceMap.at(std::make_pair(s,t)) << "]";
	else
	  std::cout << "[" << CohenrenceMap.at(std::make_pair(t,s)) << "]";
	if (wrong_edges.find(std::make_pair(s,t)) != wrong_edges.end() || wrong_edges.find(std::make_pair(t,s)) != wrong_edges.end() )
	  std::cout << "[0]";
	else
	  std::cout << "[1]";
	std::cout << " (" << s << ") -- (" << t << "); ";
      }
    } */
    std::cout << "\n----------------------------------------------------" << std::endl;
  } // function disp_Graph

    
    
////////////////////////////////////////////////////////////////////////////////
//                                   TREES                                    //
////////////////////////////////////////////////////////////////////////////////
  
  // typedef std::set<Pair> Tree;
     
  void GlobalSfM_Graph_Cleaner::sequential_Tree_Reconstruction(
	const Tree & tree,
	std::map<IndexT,Transformation> & globalTransformation) const{    

    std::list<IndexT> markedIndexT; // (TODO)-> : list<Pair> : the edges we still need to use
    // root initialisation
    IndexT root = tree.begin()->first;
    markedIndexT.push_back(root);
    Transformation T_i;
    Vec3 t_i;    t_i << 0.,0.,0.;
    const Mat3 R_i = Mat3::Identity();
    globalTransformation[root] = std::make_pair(R_i, t_i);
        
    while (!markedIndexT.empty()){
      IndexT s = markedIndexT.front();
      for (Tree::iterator iter = tree.begin(); iter != tree.end(); ++iter) {
	IndexT a = iter->first;	IndexT b = iter->second;
	if (s == a && globalTransformation.find(b) == globalTransformation.end()){
	  T_i = get_transformation(std::make_pair(s,b));
	  globalTransformation[b] = compose_transformation( T_i , globalTransformation.at(s) );
	  markedIndexT.push_back(b);

	} else if(s == b && globalTransformation.find(a) == globalTransformation.end()) {
	  T_i = get_transformation(std::make_pair(s,a));
	  globalTransformation[a] = compose_transformation( T_i , globalTransformation.at(s) );
	  markedIndexT.push_back(a);
	}
      }
      markedIndexT.pop_front();
    }
  } // function sequential_Tree_Reconstruction
  
////////////////////////////////////////////////////////////////////////////////
    
  double GlobalSfM_Graph_Cleaner::tree_Consistency_Error(const Tree & tree, int & nbpos, int & nbneg) const {
    // Local 'Global' transformation to speedup the computation
    double consistency_error = 0.; nbpos = 0; nbneg = 0;
    Transformation p_i;    Vec3 t_i;    Mat3 R_i;
    
    //             Computation of the Local 'Global' Transformation             //
    std::map<IndexT,Transformation> globalTransformation;
    sequential_Tree_Reconstruction(tree, globalTransformation);
    
    //                      Compute the consistency error                       //
    for(RelativeInfo_Map::const_iterator iter = relatives_Rt.begin();
	iter != relatives_Rt.end(); ++iter) {
      const openMVG::relativeInfo & rel = *iter;
      if ( globalTransformation.find(rel.first.first) != globalTransformation.end()
	  && globalTransformation.find(rel.first.second) != globalTransformation.end()
	  && tree.find(rel.first) == tree.end()) {
	R_i = rel.second.first; // R_i == R_(first -> second) == R_s * R_f^T
	p_i = globalTransformation.at(rel.first.first);
	R_i = R_i * p_i.first;  // R_i == R_s * R'_f^T * R'_f
	p_i = globalTransformation.at(rel.first.second);
	R_i = R_i * p_i.first.transpose();
	// R_i ==  R_s * R'_f^T * R'_f * R'_s^T == (R_s * R_f^T) * (R'_s * R'_f^T)^T
	const double error = R2D(getRotationMagnitude(R_i));
	consistency_error += error;
	if (error < error_thres)
	  nbpos += 1;
	else
	  nbneg += 1;
	//std::cout << "(" << rel.first.first << "," << rel.first.second << ") : "<< error << " " << std::endl;
      }
    }
    return consistency_error;
  } // function tree_Consistency_Error
  
////////////////////////////////////////////////////////////////////////////////
  
  Tree GlobalSfM_Graph_Cleaner::generate_Random_Tree( const int size ) const{
    // Assertion
    if (size > adjacency_map.size()) { throw "Tree to large for the graph."; }
    
    Tree tree;
    std::vector<IndexT> random_nodes;	random_nodes.reserve(adjacency_map.size());
    std::set<IndexT> tree_nodes;
    int count = 1;

    // Select the first node
    Adjacency_map::const_iterator it = adjacency_map.begin();
    std::advance(it, rand() % adjacency_map.size());
    IndexT r_node = it->first;
    
    tree_nodes.insert(r_node);
    random_nodes.push_back(r_node);
    while ( count < size ) {
      int r_position = rand() % random_nodes.size();
      r_node = random_nodes[r_position];
      const std::set<IndexT> & adj_set = adjacency_map.at(r_node);
      std::vector<IndexT> adj_vec(adj_set.begin(), adj_set.end()); 
      bool b = true;
      
      std::random_shuffle(adj_vec.begin(), adj_vec.end());
            
      // starting from r_node, we search for a new adjacent node which is not already in tree_nodes
      for ( std::vector<IndexT>::const_iterator iter=adj_vec.begin(); b && iter!=adj_vec.end(); ++iter ) {
	if( tree_nodes.find(*iter) == tree_nodes.end() ){
	  b = false;
	  tree_nodes.insert(*iter);
	  tree.insert(std::make_pair(r_node, *iter));
	  random_nodes.push_back(*iter);
	  count += 1;
	}
      }
      if (b) {// r_node is only connected to node already in the tree
	//random_nodes.erase(r_position); <-(TODO)
      }
    }
    return tree;
  } // function generate_Random_Tree
  
////////////////////////////////////////////////////////////////////////////////
  
  Tree GlobalSfM_Graph_Cleaner::generate_Consistent_Tree( const int size ) const {
    // TODO : make ACRANSAC
    const int nb_iter = 10000;
    Tree best_tree;
    double tree_err, err;
    double best_error = 10000.;
    
    int nbpos, nbneg;    
    for( int foo = 1; foo < nb_iter; foo++ ){
      globalSfM::Tree tree_i = generate_Random_Tree(size);
      tree_err = tree_Consistency_Error( tree_i, nbpos, nbneg );
      err = tree_err * (nbneg + 1)/(nbpos+1); // <-(TODO)
      
      if( err < best_error ){
	best_tree = tree_i;
	best_error = err;
	if(dbtalk){std::cout << foo << "- error: " << err << "(" << tree_err << " : " << nbpos << "/" << nbneg << ")" << std::endl;}
      }
    }
    return best_tree;
  } // generate_Consistent_Tree
  
////////////////////////////////////////////////////////////////////////////////
  
  void GlobalSfM_Graph_Cleaner::update_Coherence( const Tree & tree ) {
    Transformation p_i;    Vec3 t_i;    Mat3 R_i;
    std::map<IndexT,Transformation> globalTransformation;
    sequential_Tree_Reconstruction(tree, globalTransformation);
    
    Pair ref_pair = *tree.begin();
    
    for(RelativeInfo_Map::const_iterator iter = relatives_Rt.begin();
	iter != relatives_Rt.end(); ++iter) {
      const openMVG::relativeInfo & rel = *iter;
      if ( globalTransformation.find(rel.first.first) != globalTransformation.end()
	  && globalTransformation.find(rel.first.second) != globalTransformation.end()) {
	
	R_i = rel.second.first; // R_i == R_(first -> second) == R_s * R_f^T
	p_i = globalTransformation.at(rel.first.first);
	R_i = R_i * p_i.first;  // R_i == R_s * R'_f^T * R'_f
	p_i = globalTransformation.at(rel.first.second);
	R_i = R_i * p_i.first.transpose();
	// R_i ==  R_s * R'_f^T * R'_f * R'_s^T == (R_s * R_f^T) * (R'_s * R'_f^T)^T
      
	const double error = R2D(getRotationMagnitude(R_i));
	if (error < error_thres)	{change_Cohenrence( rel.first, ref_pair );}
      }
    }
    
  } // function update_Coherence
  
////////////////////////////////////////////////////////////////////////////////  
  
  double GlobalSfM_Graph_Cleaner::edge_Consistency_Error(
	      const Pair & pair,
	      const Tree & tree,
	      const std::set<IndexT> & tree_nodes,
	      int & nbpos,
	      int & nbneg) const{
    double consistency_error = 0.;
    nbneg = nbpos = 0;
    std::map< IndexT, int > distance;
    std::list<IndexT> markedIndexT;
    IndexT new_node, tree_node; // new_node is the node in pair which is not in tree_nodes.
    if(tree_nodes.find(pair.first) != tree_nodes.end() && tree_nodes.find(pair.second) == tree_nodes.end()) {
      new_node = pair.second;
      tree_node = pair.first;
    }
    else if (tree_nodes.find(pair.second) != tree_nodes.end() && tree_nodes.find(pair.first) == tree_nodes.end()) {
      new_node = pair.first;
      tree_node = pair.second;
    }
    else
      throw "The increasing edge must link the tree and the other part of the graph.";
    
    if(dbtalk){std::cout << "noeud ajouté (" << new_node << ") et noeud d'origine (" << tree_node << ")" << std::endl;} // DBTALK
    
    // Computation of the distance to the tree (without the node new_node)
    for(std::set<IndexT>::const_iterator iter = tree_nodes.begin(); iter != tree_nodes.end(); ++iter) {
      distance[*iter] = 0;
      markedIndexT.push_back(*iter);
    }
    while (!markedIndexT.empty()){
      IndexT s = markedIndexT.front();
      const std::set<IndexT> adj_s = adjacency_map.at(s);      
      for(std::set<IndexT>::const_iterator iter = adj_s.begin(); iter != adj_s.end(); ++iter) {
	IndexT t = *iter;
	if( distance.find(t) == distance.end() && t != new_node ) {
	  distance[t] = distance.at(s)+1;
	  markedIndexT.push_back(t);
	}
      }      
      markedIndexT.pop_front();
    }
    // Compute the globalTransformation
    std::map<IndexT,Transformation> globalTransformation; // TODO : Precompute this
    sequential_Tree_Reconstruction(tree, globalTransformation);

    Transformation T_pair = get_transformation(std::make_pair(tree_node, new_node));
    Transformation T_pt;
    
    // For every edge starting from new_node, we find a shortest path that go
    // back to the tree and we estimate the lack of consistency for this path.
    const std::set<IndexT> adj_newnode = adjacency_map.at(new_node);
    for(std::set<IndexT>::const_iterator iter = adj_newnode.begin(); iter != adj_newnode.end(); ++iter) {
      IndexT s = new_node;	IndexT t = *iter;
      if(dbtalk){std::cout << s << "[XX]" << " - " << t << "[" <<  distance.at(t) << "] ";} // DBTALK
      if( t != tree_node ){ // We dont want to consider the evaluated pair
	T_pt = compose_transformation(get_transformation(std::make_pair(s, t)), T_pair);

	// build the track back to the tree
	while ( distance.at(t) != 0 ) {
	  const std::set<IndexT> adj_t = adjacency_map.at(s);
	  for(std::set<IndexT>::const_iterator iter_s = adj_t.begin(); iter_s != adj_t.end(); ++iter_s) {
	    if(dbtalk){std::cout << "(" << *iter_s << ")";} // DBTALK
	    if( t != new_node && distance.at(*iter_s) < distance.at(t) ){  // lazy evaluation important : distance[new_node] doesn't exist...
	      s = t;
	      t = *iter_s;
	      T_pt = compose_transformation(get_transformation(std::make_pair(s, t)), T_pt);
	      if(dbtalk){std::cout << "- " << t << "[" <<  distance.at(t) << "] ";} // DBTALK
	      break;
	    }
	  }
	  
	}
	// here, t belongs to the tree.
	T_pt = compose_transformation( inverse_transformation(globalTransformation.at(t)) ,T_pt );
	T_pt = compose_transformation( globalTransformation.at(tree_node) ,T_pt );
	const double error = R2D(getRotationMagnitude(T_pt.first));	
	if(dbtalk){std::cout << "  >>  " << tree_node << "->" << t << " (err = " << error << ")";} // DBTALK	  
	consistency_error += error;
	if (error < error_thres)
	  nbpos += 1;
	else
	  nbneg += 1;
      }
      if(dbtalk){std::cout << std::endl;} // DBTALK
    }
    
    return consistency_error;
  } // function edge_Consistency_Error
  
////////////////////////////////////////////////////////////////////////////////
  
  
  void GlobalSfM_Graph_Cleaner::increase_Tree( Tree & tree ) const {
    std::set<IndexT> tree_nodes;
    int nbpos, nbneg;
    double consistency_error, error;
    double min_error = 0;
    Pair best_pair;
    
    if(dbtalk){std::cout << "Construction de tree_nodes :";} // DBTALK
    
    for(Tree::iterator iter = tree.begin(); iter != tree.end(); ++iter) {
      tree_nodes.insert(iter->first);
      tree_nodes.insert(iter->second);
    }
    if(dbtalk){std::cout << "\nOk\nIteration sur les voisins :" << std::endl;} // DBTALK
    
    while(min_error < 10000) {
      min_error = 10000;
      for(RelativeInfo_Map::const_iterator iter = relatives_Rt.begin(); iter != relatives_Rt.end(); ++iter) {
	Pair pair = iter->first;
	IndexT s = pair.first;
	IndexT t = pair.second;
	if( ( tree_nodes.find(s) == tree_nodes.end()
	      && tree_nodes.find(t) != tree_nodes.end() )
	  ||( tree_nodes.find(s) != tree_nodes.end()
	      && tree_nodes.find(t) == tree_nodes.end() ) ){
	  consistency_error = edge_Consistency_Error( pair, tree, tree_nodes, nbpos, nbneg);
	  error = consistency_error * ( nbneg + 1 ) / ( nbpos + 1 ); // <-(TODO)
		  
	  if ( error < min_error ){
	    min_error = error;
	    best_pair = pair;
	    if(dbtalk){std::cout << s << "," << t << " : " << error << " = " << consistency_error << "(" << nbpos << "/" << nbneg << ")" << std::endl;} // DBTALK
	  } // if
	} // if
      } // for
      if (min_error < 10000){
	if(dbtalk){std::cout << " OK : on prend " << best_pair.first << " et " << best_pair.second << std::endl;} 
	tree_nodes.insert(best_pair.first);
	tree_nodes.insert(best_pair.second);
	tree.insert(best_pair);
      }
    } // while
    
  } // function increase_Tree
  
  
} // namespace globalSfM
} // namespace openMVG
