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
 for( int k = 0; k < number_tree; k++ ){
  dbtalk = true;

  std::cout << "------------\nConstruction of the initial..." << std::endl;
  Tree coherent_tree = generate_Consistent_Tree(initial_tree_size);
  
      for (Tree::const_iterator iter = coherent_tree.begin(); iter != coherent_tree.end(); ++iter)
	tikzfile << "\\draw[itree](" << iter->first << ")--(" << iter->second << ");"<<std::endl;
      
  std::cout << "------------\nEnlargement of the tree..." << std::endl;  
  increase_Tree( coherent_tree );
  
  std::cout << "------------\nFind coherent edges..." << std::endl;  
  update_Coherence( coherent_tree );
  dbtalk = false;
     
  for (Tree::const_iterator iter = coherent_tree.begin(); iter != coherent_tree.end(); ++iter)
    tikzfile << "\\draw[tree](" << iter->first << ")--(" << iter->second << "); ";
  tikzfile << "\n";
 }
  /* *//*

  Cycles cycles = findCycles();
  rotationRejection(5.0f, cycles); /* */
  
  std::cout << "------------\nExtraction of the principal coherent component..." << std::endl;  
  int max_comp = principale_Coherence();
  
  tikzfile << "\n\\def\\componentchoice{" << max_comp << "}\n";
  disp_Graph("base");

  
  RelativeInfo_Map new_relatives_Rt;

  dbtalk = true;
  if(dbtalk){ std::cout << "\nEdges removed :"; }
  for (std::map< Pair, int >::const_iterator iter = CohenrenceMap.begin(); iter != CohenrenceMap.end(); ++iter) {
    if (iter->second == max_comp)
      new_relatives_Rt[iter->first] = relatives_Rt.at(iter->first);
    else if( dbtalk )
      std::cout << " (" << iter->first.first << "," << iter->first.second << ")[" << iter->second << "]";
  }
  dbtalk = false;
  
  std::cout << "-------------------------\nINFO:" << std::endl;
  std::cout << relatives_Rt.size() - new_relatives_Rt.size() << " edges removed." << std::endl;
  std::cout << new_relatives_Rt.size() << " edges left" << std::endl;
  std::cout << "-------------------------" << std::endl;


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
	consistency_error += min( error, max_error_thres );
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
    const int nb_iter = 10000;
    Tree best_tree;
    double tree_err, err;
    double best_error = 10000.;
    
    int nbpos, nbneg;    
    for( int foo = 1; foo < nb_iter; foo++ ){
      Tree tree_i = generate_Random_Tree(size);
      tree_err = tree_Consistency_Error( tree_i, nbpos, nbneg );
      
      std::set<int> consist_set;
      for(Tree::iterator iter = tree_i.begin(); iter != tree_i.end(); ++iter)
	consist_set.insert(get_Pair_Coherence(*iter));
      if( consist_set.size() > 1 ){
	err = size * ( 1 + (tree_err/error_cst) ) * ( nbneg + 1 )/( (nbpos + 1)*(nbpos + .01) ); // <-(TODO) & (à mettre dans la fonction de calcul d'erreur)
      } else { continue; }
      
      if( err < best_error ){
	best_tree = tree_i;
	best_error = err;
	if(dbtalk){std::cout << foo << " - error: " << err << "(" << tree_err << " : " << nbpos << "/" << nbneg << " : " << consist_set.size() << ") : ";
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
      
	const double error = T_i.error();
	if(dbtalk){std::cout << rel.first.first << "-" << rel.first.second <<" : Erreur chemin = ";    for (std::list<IndexT>::const_iterator it=T_i.path.begin(); it!=T_i.path.end(); ++it){ std::cout << ' ' << *it;} std::cout << " : " << error << std::endl;}
	if (error < error_thres) 
	  change_Cohenrence( rel.first, ref_pair );
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
    double error = default_max_error;
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
	  
	  // DEBUGZZ
	  // if(dbtalk){std::cout << "{" << foo_error << " |"; T_foo.print_path(); std::cout << "}"; std::cout.flush();} // DBTALKZZ
	  // DEBUGZZ
	  
	  if (foo_error < error) { error = foo_error; }
	  if (nbpos_i + nbneg_i >= descent_max) { return error;}
	  //if(dbtalk){std::cout << "  Path = "; T_foo.print_path(); std::cout << " : err = " << foo_error << " nbpos/neg = " << nbpos_i << "/" << nbneg_i << std::endl;}

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

  std::map< IndexT, int > GlobalSfM_Graph_Cleaner::compute_distance(
		const IndexT & new_node,
		const std::set<IndexT> & tree_node_set) const {
    std::map< IndexT, int > distance;
    std::list<IndexT> markedIndexT;
    
    for(std::set<IndexT>::const_iterator iter = tree_node_set.begin(); iter != tree_node_set.end(); ++iter) {
      if( *iter != new_node ){
	distance[*iter] = 0;
	markedIndexT.push_back(*iter);
      }
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
    distance[new_node] = adjacency_map.size() + 1;
    return distance;
}
  
////////////////////////////////////////////////////////////////////////////////
// tree_node -> new_node : the evaluated edge
// tree_node_set  : the tree
  
  double GlobalSfM_Graph_Cleaner::edge_Consistency_Error(
	      const IndexT & tree_node, const IndexT & new_node,
	      const std::set<IndexT> & tree_node_set,
	      const std::map<IndexT,Transformation> & globalTransformation,
	      int & nbpos, int & nbneg, bool & edge_skip) const{
    double consistency_error = 0.;
    nbneg = nbpos = 0;
    edge_skip = false;
    int nbneg_i, nbpos_i;
        
    std::map< IndexT, int > distance = compute_distance( new_node, tree_node_set );
        
    Transformation T_pair = compose_transformation(
			get_transformation(std::make_pair(tree_node, new_node)),
			globalTransformation.at(tree_node));
    Transformation T_pt;
          
    const std::set<IndexT> adj_newnode = adjacency_map.at(new_node);
    for(std::set<IndexT>::const_iterator iter = adj_newnode.begin(); iter != adj_newnode.end(); ++iter) {
      if( (*iter != tree_node) && (distance.find(*iter) != distance.end()) ){ // We dont want to consider the evaluated pair
	//if(dbtalk){std::cout << "path(" << *iter << ")=";} // DBTALKZZ
	nbpos_i = nbneg_i = 0;
	//if(dbtalk){std::cout << "."; std::cout.flush();} // DBTALK
	double error = edge_Descent_Error(
		  *iter,
		  compose_transformation(get_transformation(std::make_pair(new_node,*iter)), T_pair),
		  distance,globalTransformation,nbpos_i,nbneg_i);
	error = (error + error_cst) * (nbneg_i+1) / (nbpos_i+1);
	// error = ( 1 + (error/error_cst) ) * ( nbneg_i + 1 )/( (nbpos_i + 1)*(nbpos_i + .01) );
	
	if( nbpos_i + nbneg_i >= descent_max ){ edge_skip = true; }
	
	//if(dbtalk){std::cout << " || Total_error=" << error << std::endl;} // DBTALKZZ
	
	//if(dbtalk){std::cout << " " << new_node << " -> " << *iter << " : " <<  error << "(" << nbpos_i << "/" << nbneg_i << ")" << std::endl;} // DBTALK
	
	consistency_error += min( error, max_error_thres );
	if (error < error_thres){
	  //if(dbtalk){std::cout << "+(" << *iter << ")["<< nbpos << "]" << std::endl;} // DBTALKZZ
	  nbpos += 1;
	} //AA
	else
	  nbneg += 1;
      }
    }
    
    return consistency_error;
  } // function edge_Consistency_Error
  
////////////////////////////////////////////////////////////////////////////////
  
  
  void GlobalSfM_Graph_Cleaner::increase_Tree( Tree & tree ) const {
    std::set<IndexT> tree_node_set;
    int nbpos, nbneg;
    double consistency_error, error;
    double min_error = 0;
    Pair best_pair;
    Transformation T_foo;
    bool edge_skip;

    // Build the tree_node_set
    for(Tree::iterator iter = tree.begin(); iter != tree.end(); ++iter) {
      tree_node_set.insert(iter->first);
      tree_node_set.insert(iter->second);
    }

    // Set the initial set of edges which could be add to the graph.
    std::map<Pair, int> edge_indic;
    for(RelativeInfo_Map::const_iterator iter = relatives_Rt.begin(); iter != relatives_Rt.end(); ++iter) {
      Pair pair = iter->first;
      IndexT s = pair.first;
      IndexT t = pair.second;
      if( ( tree_node_set.find(s) == tree_node_set.end() && tree_node_set.find(t) != tree_node_set.end() )
	||( tree_node_set.find(s) != tree_node_set.end() && tree_node_set.find(t) == tree_node_set.end() ) ){
	edge_indic[ pair ] = 0;
      }
    }

    // Compute the globalTransformation
    std::map<IndexT,Transformation> globalTransformation;
    sequential_Tree_Reconstruction(tree, globalTransformation);

    ///////////////////////////////////////////////////////////////////////////

    while( edge_indic.size() > 0 ) {
      min_error = default_max_error;

	  if(dbtalk){std::cout << "\n\nTree_nodes (" << tree_node_set.size() << ") = "; // DBTALK
	  for(std::set<IndexT>::iterator iter = tree_node_set.begin(); iter != tree_node_set.end(); ++iter)
		  std::cout << *iter << " "; // DBTALK
	  std::cout << std::endl;}
      
      for(std::map<Pair, int>::iterator iter = edge_indic.begin(); iter != edge_indic.end(); ++iter) {
	Pair pair = iter->first;
	IndexT new_node, tree_node;

	// Edge selection
	if(tree_node_set.find(pair.first) != tree_node_set.end() && tree_node_set.find(pair.second) == tree_node_set.end()) {
	  new_node = pair.second;	  tree_node = pair.first;	}
	else if (tree_node_set.find(pair.second) != tree_node_set.end() && tree_node_set.find(pair.first) == tree_node_set.end()) {
	  new_node = pair.first;	  tree_node = pair.second;	}
	else {	  edge_indic.erase(pair);	  continue;		}

	if( iter->second > 0 ) {
	  iter->second -= 1;
	} else {
	  
	      if(dbtalk){std::cout << tree_node << "," << new_node << ":"; std::cout.flush();}

	  // We compute the consistency error of the edge
	  consistency_error = edge_Consistency_Error( tree_node, new_node, tree_node_set, globalTransformation, nbpos, nbneg, edge_skip);
	  //error = ( error_cst + consistency_error ) * ( nbneg + 1 ) / ( nbpos + 1 );
	  error = ( 1 + (consistency_error/error_cst) ) * ( nbneg + 1 )/( (nbpos + 1)*(nbpos + .01) );
	  
	  if( edge_skip ) { iter->second = nb_steps_skip; } // The computation was to slow
	  
	  //if(dbtalk){std::cout << "\n" << tree_node << "," << new_node << " : " << error << " = " << consistency_error << "(" << nbpos << "/" << nbneg << ")" << std::endl;} // DBTALKZZ

	  if ( error < min_error ){ // We keep the extremum
	    //if(dbtalk){std::cout << "ZZZ\n";} // DBTALKZZ
	    min_error = error;
	    best_pair = std::make_pair(tree_node, new_node);
	    if(dbtalk){std::cout << "\n" << tree_node << "," << new_node << " : " << error << " = " << consistency_error << "(" << nbpos << "/" << nbneg << ")\n";} // DBTALKzz
	    // If the edge is to good, we take it without checking the other edges
	    if ( error < error_valid_thres ){
	      if(dbtalk){std::cout << "Break" << std::endl;} // DBTALK
	      break;
	    }
	  } // if

	} // if
      } // for
      
      if (min_error < edge_inconsistency_thres){ // We add the pair to the tree
	tree_node_set.insert(best_pair.first);
	tree_node_set.insert(best_pair.second);
	tree.insert(best_pair);
	
	// We compute the global transformation for the new node
	if( globalTransformation.find(best_pair.second) == globalTransformation.end() ){
	  T_foo = globalTransformation.at(best_pair.first);
	  T_foo.compose_left( get_transformation(best_pair) );
	  globalTransformation[ best_pair.second ] = T_foo;
	} else {
	  T_foo = globalTransformation.at(best_pair.second);
	  T_foo.compose_left( get_transformation(std::make_pair(best_pair.second, best_pair.first)) );
	  globalTransformation[ best_pair.first ] = T_foo;
	}
	// We compute the distance to the tree, to ignore the non biconnexe vertexes
	std::map<IndexT,int> distance = compute_distance( best_pair.second, tree_node_set );
	const std::set<IndexT> adj_newnode = adjacency_map.at(best_pair.second);
	
	    if(dbtalk){std::cout << "\nEdge: " << best_pair.first << "->" << best_pair.second << " : new edges = "; }
	
	for(std::set<IndexT>::const_iterator iter = adj_newnode.begin(); iter != adj_newnode.end(); ++iter) {
	  if( distance.find(*iter) != distance.end() && distance.at(*iter) != 0 ){ // We dont want to consider the evaluated pair
	    edge_indic[std::make_pair(best_pair.second, *iter)] = 0;
	    if(dbtalk){std::cout << " " << best_pair.second << "-" << *iter << "[" << distance.at(*iter) << "]"; }
	  }
	}
      } else {   break;   }
    } // while
    
  } // function increase_Tree
  
  
} // namespace sfm
} // namespace openMVG
