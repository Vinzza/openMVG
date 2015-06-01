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

////////////////////////////////////////////////////////////////////////////////
//                                    Run                                     //
////////////////////////////////////////////////////////////////////////////////

RelativeInfo_Map GlobalSfM_Graph_Cleaner::run()
{
  
  std::cout << "------------\nConstruction of the initial..." << std::endl;
  dbtalk = true;
  Tree coherent_tree = generate_Consistent_Tree(7);
  
  for (Tree::const_iterator iter = coherent_tree.begin(); iter != coherent_tree.end(); ++iter)
    tikzfile << "\\draw[itree](" << iter->first << ")--(" << iter->second << ");"<<std::endl;
    
  std::cout << "------------\nEnlargement of the tree..." << std::endl;  
  increase_Tree( coherent_tree );
  std::cout << "------------\nFind coherent edges..." << std::endl;  
  update_Coherence( coherent_tree );
  dbtalk = false;
  
  std::cout << "------------\nExtraction of the principal coherent component..." << std::endl;  
  int max_comp = principale_Coherence();
  
  tikzfile << "\n\\def\\componentchoice{" << max_comp << "}\n";
  disp_Graph("base");
  
  for (Tree::const_iterator iter = coherent_tree.begin(); iter != coherent_tree.end(); ++iter)
    tikzfile << "\\draw[tree](" << iter->first << ")--(" << iter->second << "); ";
  tikzfile << "\n";
  
  RelativeInfo_Map new_relatives_Rt;
  
  for (std::map< Pair, int >::const_iterator iter = CohenrenceMap.begin(); iter != CohenrenceMap.end(); ++iter) {
    if (iter->second == max_comp)
      new_relatives_Rt[iter->first] = relatives_Rt.at(iter->first);
  }
  
   
  /*
  std::cout << "-----------------------------------------------------------\nTriplets Inference..." << std::endl;
  Cycles cycles = findCycles();
  rotationRejection(5.0f, cycles);*/
  return new_relatives_Rt;
  
} // function run

////////////////////////////////////////////////////////////////////////////////
//                                Find Cycles                                 //
////////////////////////////////////////////////////////////////////////////////

Cycles GlobalSfM_Graph_Cleaner::findCycles() {
  Pair_Set pair_set;
  for(RelativeInfo_Map::const_iterator iter = relatives_Rt.begin(); iter != relatives_Rt.end(); ++iter)
      pair_set.insert(iter->first);
  
  std::vector<graph::Triplet> vec_triplets = graph::tripletListing(pair_set);
  
  Cycles vec_cycles;
  vec_cycles.reserve(vec_triplets.size());
  
  for ( std::vector<graph::Triplet>::iterator iter=vec_triplets.begin(); iter!=vec_triplets.end(); ++iter ) {
    graph::Triplet triplet_i = *iter;    
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
//    for(Adjacency_map::const_iterator iter = adjacency_map.begin(); iter != adjacency_map.end(); ++iter) {
//      Vec3 foo = position_GroundTruth[iter->first];
//      tikzfile << "\\node[n" << str << "] at (" << foo(0) <<  ","<< foo(2) << ")" << " (" << iter->first << ") " << "{" << iter->first << "};\n";
//    }
    for(Adjacency_map::const_iterator iter = adjacency_map.begin();
	iter != adjacency_map.end(); ++iter) {
      IndexT s = iter->first;
      const std::set<IndexT> & indexT_set = iter->second;
    
      for (std::set<IndexT>::const_iterator iterT = indexT_set.begin();
	  iterT != indexT_set.end(); ++iterT) {
	IndexT t = *iterT;
	if (s < t) {
	  tikzfile << "\\zdraw["<< str << "]";
	  if (CohenrenceMap.find(std::make_pair(s,t)) != CohenrenceMap.end())
	    tikzfile << "[" << CohenrenceMap.at(std::make_pair(s,t)) << "]";
	  else
	    tikzfile << "[" << CohenrenceMap.at(std::make_pair(t,s)) << "]";
	  if (wrong_edges.find(std::make_pair(s,t)) != wrong_edges.end() || wrong_edges.find(std::make_pair(t,s)) != wrong_edges.end() )
	    tikzfile << "[0]";
	  else
	    tikzfile << "[1]";
	  tikzfile << " (" << s << ") -- (" << t << "); ";
	}
      }
    }
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
    globalTransformation[root] = T_i;
    
    //Vec3 t_i;    t_i << 0.,0.,0.;
    //const Mat3 R_i = Mat3::Identity();    
    //globalTransformation[root] = std::make_pair(R_i, t_i);
        
    while (!markedIndexT.empty()){
      IndexT s = markedIndexT.front();
      for (Tree::iterator iter = tree.begin(); iter != tree.end(); ++iter) {
	IndexT a = iter->first;	IndexT b = iter->second;
	if (s == a && globalTransformation.find(b) == globalTransformation.end()){
	  T_i = get_transformation(std::make_pair(s,b));
	  T_i.compose_right( globalTransformation.at(s) );
	  globalTransformation[b] = T_i;      // <-(TODO:L)
	  //globalTransformation[b] = compose_transformation( T_i , globalTransformation.at(s) );
	  markedIndexT.push_back(b);

	} else if(s == b && globalTransformation.find(a) == globalTransformation.end()) {
	  T_i = get_transformation(std::make_pair(s,a));
	  T_i.compose_right( globalTransformation.at(s) );
	  globalTransformation[a] = T_i;      // <-(TODO:L)
	  //globalTransformation[a] = compose_transformation( T_i , globalTransformation.at(s) );
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
    Transformation T_i;//    Vec3 t_i;    Mat3 R_i;
    
    //             Computation of the Local 'Global' Transformation             //
    std::map<IndexT,Transformation> globalTransformation;
    sequential_Tree_Reconstruction(tree, globalTransformation);
    
    //                      Compute the consistency error                       //
    for(RelativeInfo_Map::const_iterator iter = relatives_Rt.begin();
	iter != relatives_Rt.end(); ++iter) {
      const openMVG::relativeInfo & rel = *iter;
      if ( globalTransformation.find(rel.first.first) != globalTransformation.end()
	  && globalTransformation.find(rel.first.second) != globalTransformation.end()
	  && tree.find(rel.first) == tree.end() && tree.find(std::make_pair(rel.first.second, rel.first.first)) == tree.end()) {

	T_i = globalTransformation.at(rel.first.first);
	T_i.compose_left(rel.second, rel.first);
	T_i.compose_left_rev( globalTransformation.at(rel.first.second) );

	const double error = T_i.error();
	consistency_error += error;
	if (error < error_thres)
	  nbpos += 1;
	else
	  nbneg += 1;
      }
    }
    if( nbpos == 0 && nbneg == 0){ return 1000; }
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
      Tree tree_i = generate_Random_Tree(size);
      tree_err = tree_Consistency_Error( tree_i, nbpos, nbneg );
      err = (5+tree_err) * (nbneg + 1)/(nbpos+1); // <-(TODO)
      
      if( err < best_error ){
	best_tree = tree_i;
	best_error = err;
	if(dbtalk){std::cout << foo << " - error: " << err << "(" << tree_err << " : " << nbpos << "/" << nbneg << ") : ";
	  for (Tree::const_iterator iter = tree_i.begin(); iter != tree_i.end(); ++iter)
	    std::cout << iter->first << " " << iter->second << " ";
	  std::cout << std::endl;}
      }
    }
    return best_tree;
  } // generate_Consistent_Tree
  
////////////////////////////////////////////////////////////////////////////////
  
  void GlobalSfM_Graph_Cleaner::update_Coherence( const Tree & tree ) {
    Transformation T_i;    Vec3 t_i;    Mat3 R_i;
    std::map<IndexT,Transformation> globalTransformation;
    sequential_Tree_Reconstruction(tree, globalTransformation);
    
    Pair ref_pair = *tree.begin();
    
    for(RelativeInfo_Map::const_iterator iter = relatives_Rt.begin();
	iter != relatives_Rt.end(); ++iter) {
      const openMVG::relativeInfo & rel = *iter;
      if ( globalTransformation.find(rel.first.first) != globalTransformation.end()
	  && globalTransformation.find(rel.first.second) != globalTransformation.end()) {
	
	T_i = globalTransformation.at(rel.first.first);
	T_i.compose_left(rel.second, rel.first);
	T_i.compose_left_rev( globalTransformation.at(rel.first.second) );
	//R_i = rel.second.first; // R_i == R_(first -> second) == R_s * R_f^T
	//T_i = globalTransformation.at(rel.first.first);
	//R_i = R_i * T_i.first;  // R_i == R_s * R'_f^T * R'_f
	//T_i = globalTransformation.at(rel.first.second);
	//R_i = R_i * T_i.first.transpose();
	// R_i ==  R_s * R'_f^T * R'_f * R'_s^T == (R_s * R_f^T) * (R'_s * R'_f^T)^T
      
	const double error = T_i.error(); // R2D(getRotationMagnitude(R_i));
	if(dbtalk){std::cout << "Erreur chemin = ";    for (std::list<IndexT>::const_iterator it=T_i.path.begin(); it!=T_i.path.end(); ++it){ std::cout << ' ' << *it;} std::cout << " : " << error << std::endl;}
	if (error < error_thres)	{change_Cohenrence( rel.first, ref_pair );}
      }
    }
    
  } // function update_Coherence
  
////////////////////////////////////////////////////////////////////////////////
  
  double GlobalSfM_Graph_Cleaner::edge_Descent_Error(
	      const IndexT & to,
	      Transformation T,
	      const std::map< IndexT, int > & distance,
	      const std::map<IndexT,Transformation> & globalTransformation,
	      int & nbpos_i,
	      int & nbneg_i  ) const {
    double error = 10000.;
    double foo_error;
    std::list<std::pair<bool,Pair>> paths;
    std::pair<bool,Pair> path_foo;
    IndexT s, t;
    Transformation T_foo;
    
    if ( distance.at(to) == 0 ){
      T.compose_left_rev(globalTransformation.at(to));
      error = T.error();
      if (error < error_thres)
	nbpos_i += 1;
      else
	nbneg_i += 1;
      //if(dbtalk){std::cout << "  Path0 = "; T.print_path(); std::cout << " : err = " << error << std::endl;}
      return error;
    }
    
    //if(dbtalk){std::cout << "  " << to << "(+"; std::cout.flush();}

    const std::set<IndexT> adj_t = adjacency_map.at(to);
    for( std::set<IndexT>::const_iterator iter = adj_t.begin(); iter != adj_t.end(); ++iter ) {      
      if( distance.at(to) > distance.at(*iter) ){
	paths.push_front(std::make_pair(true, std::make_pair(to,*iter)));
	//if(dbtalk){std::cout << " " << *iter; std::cout.flush();}
      }
    }
    //if(dbtalk){std::cout << ")"; std::cout.flush();}
    
    while( !paths.empty() ){

      path_foo = paths.front();
      s = path_foo.second.first;
      t = path_foo.second.second;
      paths.pop_front();
      //if(dbtalk){std::cout << "  " << t; std::cout.flush();}
      
      if( path_foo.first ){ // We need to continue the path
	//if(dbtalk){std::cout << "(+"; std::cout.flush();}

	if( 0 == distance.at(t) ){ // The vertex is in the tree
	  //if(dbtalk){std::cout << "/)"; std::cout.flush();}
	  T_foo = T;
	  T_foo.compose_left(get_transformation(path_foo.second));
	  T_foo.compose_left_rev(globalTransformation.at(t));
	  foo_error = T_foo.error();
	  if (foo_error < error_thres)
	    nbpos_i += 1;
	  else
	    nbneg_i += 1;
	  if (foo_error < error) { error = foo_error; }
	  //if(dbtalk){std::cout << "  Path = "; T_foo.print_path(); std::cout << " : err = " << foo_error << std::endl;}

	} else { // The vertex is not in the tree : we must continue
	  
	  T.compose_left(get_transformation(path_foo.second));  
	  paths.push_front(std::make_pair(false, std::make_pair(t,s)));

	  const std::set<IndexT> adj_t = adjacency_map.at(t);
	  for( std::set<IndexT>::const_iterator iter = adj_t.begin(); iter != adj_t.end(); ++iter ) {      
	    if( distance.at(t) > distance.at(*iter) ){
	      paths.push_front(std::make_pair(true, std::make_pair(t,*iter)));
	      //if(dbtalk){std::cout << " " << *iter; std::cout.flush();}
	    }
	  }
	  //if(dbtalk){std::cout << ")"; std::cout.flush();}
	}
	
      } else { // We need the go back the path
	
	T.compose_left(get_transformation(path_foo.second));
	//if(dbtalk){std::cout << "(-)"; std::cout.flush();}
	
      }
      
    }
    
    return error;
  }
  
////////////////////////////////////////////////////////////////////////////////
  
  double GlobalSfM_Graph_Cleaner::edge_Consistency_Error(
	      const Pair & pair,
	      const Tree & tree,
	      const std::set<IndexT> & tree_nodes,
	      int & nbpos,
	      int & nbneg) const{
    double consistency_error = 0.;
    nbneg = nbpos = 0;
    int nbneg_i, nbpos_i;
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
    
    //if(dbtalk){std::cout << "a"; std::cout.flush();} // DBTALK
    
    // Computation of the distance to the tree (without the node new_node) TODO : precompute the distance for efficiency
    for(std::set<IndexT>::const_iterator iter = tree_nodes.begin(); iter != tree_nodes.end(); ++iter) {
      distance[*iter] = 0;
      markedIndexT.push_back(*iter);
      // if(dbtalk){std::cout << *iter << "=0 : ";} // DBTALK
    }
    // if(dbtalk){std::cout << std::endl;} // DBTALK
    while (!markedIndexT.empty()){
      IndexT s = markedIndexT.front();
      const std::set<IndexT> adj_s = adjacency_map.at(s);      
      for(std::set<IndexT>::const_iterator iter = adj_s.begin(); iter != adj_s.end(); ++iter) {
	IndexT t = *iter;
	if( distance.find(t) == distance.end() && t != new_node ) {
	  distance[t] = distance.at(s)+1;
	  markedIndexT.push_back(t);
	  // if(dbtalk){std::cout << t << "=" << distance.at(t)<<" : ";} // DBTALK
	}
      }      
      markedIndexT.pop_front();
    }
    distance[new_node] = adjacency_map.size() + 1;
    //if(dbtalk){std::cout << "b"; std::cout.flush();} // DBTALK
    
    // Compute the globalTransformation
    std::map<IndexT,Transformation> globalTransformation; // TODO : Precompute this
    sequential_Tree_Reconstruction(tree, globalTransformation);

    Transformation T_pair = compose_transformation(
			get_transformation(std::make_pair(tree_node, new_node)),
			globalTransformation.at(tree_node));
    Transformation T_pt;
    
    //if(dbtalk){std::cout << "c"; std::cout.flush();} // DBTALK
    
    // For every edge starting from new_node, we find a shortest path that go
    // back to the tree and we estimate the lack of consistency for this path.    
    const std::set<IndexT> adj_newnode = adjacency_map.at(new_node);
    for(std::set<IndexT>::const_iterator iter = adj_newnode.begin(); iter != adj_newnode.end(); ++iter) {
      if( *iter != tree_node ){ // We dont want to consider the evaluated pair
	// if(dbtalk){std::cout << "error(" << *iter << ")" << std::endl;} // DBTALK
	nbpos_i = nbneg_i = 0;
	//if(dbtalk){std::cout << "."; std::cout.flush();} // DBTALK
	const double error = edge_Descent_Error(
		  *iter,
		  compose_transformation(get_transformation(std::make_pair(new_node,*iter)), T_pair),
		  distance,globalTransformation,nbpos_i,nbneg_i);
	//if(dbtalk){std::cout << ","; std::cout.flush();} // DBTALK
	
	//if(dbtalk){std::cout << " " << new_node << " -> " << *iter << " : " <<  error << "(" << nbpos_i << "/" << nbneg_i << ")" << std::endl;} // DBTALK
	
	consistency_error += error;
	if (error < error_thres)
	  nbpos += 1;
	else
	  nbneg += 1;
      }
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
            
    for(Tree::iterator iter = tree.begin(); iter != tree.end(); ++iter) {
      tree_nodes.insert(iter->first);
      tree_nodes.insert(iter->second);
    }

    while(min_error < 10000) {
      min_error = 10000;
      if(dbtalk){std::cout << "\nTree_nodes = "; // DBTALK
      for(std::set<IndexT>::const_iterator iter = tree_nodes.begin(); iter != tree_nodes.end(); ++iter)
	      std::cout << *iter << " "; // DBTALK
      std::cout << std::endl;}
      for(RelativeInfo_Map::const_iterator iter = relatives_Rt.begin(); iter != relatives_Rt.end(); ++iter) {
	Pair pair = iter->first;
	IndexT s = pair.first;
	IndexT t = pair.second;
	if( ( tree_nodes.find(s) == tree_nodes.end() && tree_nodes.find(t) != tree_nodes.end() )
	  ||( tree_nodes.find(s) != tree_nodes.end() && tree_nodes.find(t) == tree_nodes.end() ) ){
	  
	  consistency_error = edge_Consistency_Error( pair, tree, tree_nodes, nbpos, nbneg);
	  error = (5+consistency_error) * ( nbneg + 1 ) / ( nbpos + 1 ); // <-(TODO)
	  
	  if ( error < min_error ){
	    if(dbtalk){std::cout << s << "," << t << " : " << error << " = " << consistency_error << "(" << nbpos << "/" << nbneg << ")" << std::endl;} // DBTALK
	    min_error = error;
	    best_pair = pair;
	  } // if
	} // if
      } // for
      if (min_error < 10000){
	tree_nodes.insert(best_pair.first);
	tree_nodes.insert(best_pair.second);
	tree.insert(best_pair);
      }
    } // while
    
  } // function increase_Tree
  
  
} // namespace sfm
} // namespace openMVG
