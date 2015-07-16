// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef GLOBALSFM_GRAPH_CLEANER_HPP
#define GLOBALSFM_GRAPH_CLEANER_HPP


#include "openMVG/sfm/sfm.hpp"
#include "openMVG/graph/graph.hpp"
#include "openMVG/multiview/rotation_averaging_common.hpp"

// #include <lemon/list_graph.h>

namespace openMVG{
namespace sfm{

  
////////////////////////////////////////////////////////////////////////////////
//                               TRANSFORMATION                               //
//////////////////////////////////////////////////////////////////////////////// 

struct Transformation {
public:

  Transformation():R(Mat3::Identity()),t(Vec3::Zero()) {}
  
  Transformation( std::pair<Mat3,Vec3> nt, Pair np ): R(nt.first), t(nt.second) {
    path.push_back(np.first);    path.push_back(np.second);
  }
  Transformation( Mat3 nr, Vec3 nt, std::list<IndexT> np ): R(nr), t(nt), path(np) {}

  // Constructor for the inverse transformation
  Transformation( std::pair<Mat3,Vec3> nt, Pair np, string foo ):
	  R(nt.first.transpose()),
	  t(- nt.first.transpose()*nt.second) {
      path.push_back(np.second);    path.push_back(np.first);
  }
  
  Transformation( const RelativeInfo_Map & relatives_Rt, IndexT ns, IndexT nt ){
    if (relatives_Rt.find(std::make_pair(ns,nt)) != relatives_Rt.end()){
      const std::pair<Mat3, Vec3> foo = relatives_Rt.at(std::make_pair(ns,nt));
      R = foo.first;      t = foo.second;
    } else {
      const std::pair<Mat3, Vec3> foo = relatives_Rt.at(std::make_pair(nt,ns));
      R = foo.first.transpose();      t = - R*foo.second;
    }
    path.clear();       path.push_back(ns);    path.push_back(nt);
  }
  
  void inverse(){
      R = R.transpose();
      t = - R*t;
      path.reverse();
  }
  
  
  void load_transformation( const RelativeInfo_Map & relatives_Rt, IndexT ns, IndexT nt ){
    if (relatives_Rt.find(std::make_pair(ns,nt)) != relatives_Rt.end()){
      const std::pair<Mat3, Vec3> foo = relatives_Rt.at(std::make_pair(ns,nt));
      R = foo.first;      t = foo.second;
    } else {
      const std::pair<Mat3, Vec3> foo = relatives_Rt.at(std::make_pair(nt,ns));
      R = foo.first.transpose();      t = - R*foo.second;
    }
    path.clear();       path.push_back(ns);    path.push_back(nt);
  }
  

  // (TODO:L)-> Clean the path when compose
  void compose_left( const Transformation & T ){
    R = T.R * R;    t = T.R * t + T.t;
    path.insert(path.end(), T.path.begin(), T.path.end());
  }
  void compose_left( const std::pair<Mat3,Vec3> & nt, const Pair & np ){
    R = nt.first * R;    t = nt.first * t + nt.second;
    path.push_back(np.first);	path.push_back(np.second);
  }
  void compose_right( const Transformation & T ){
    R = R * T.R;    t = R * T.t + t;
    std::list<IndexT> foo = T.path;
    foo.insert( foo.end(), path.begin(), path.end() );
    path = foo;
  }
  void compose_right( const std::pair<Mat3,Vec3> & nt, const Pair & np ){
    R = R * nt.first;    t = R * nt.second + t;
    path.push_front(np.second);	path.push_front(np.first);
  }
  void compose_left_rev( const Transformation & T ){
    R = T.R.transpose() * R;    t = T.R.transpose() * (t - T.t);
    std::list<IndexT> foo = T.path;
    foo.reverse();
    path.insert(path.end(), foo.begin(), foo.end());
  }
  
  void clean_cycle_path(){
    IndexT foo1, bar1, foo2, bar2;
    foo1 = path.front();    path.pop_front();    foo2 = path.front();    path.pop_front();
    bar1 = path.back();    path.pop_back();    bar2 = path.back();    path.pop_back();
      
    if( foo1 != bar1 )
      throw "The path is not a cycle. We can't compute the consistency error."; 
    while( foo2 == bar2 ){
      if( path.empty() ){ return; }
      foo1 = path.front();    path.pop_front();      foo2 = path.front();    path.pop_front();
      if( path.empty() ){ path.push_front(foo2); path.push_front(foo1); return; }
      bar1 = path.back();    path.pop_back();      bar2 = path.back();    path.pop_back();
    }
    path.push_front(foo2);    path.push_front(foo1);    path.push_back(bar2);    path.push_back(bar1);
  }

  void clean_path(){

    std::list<IndexT> cpath;
    IndexT foo1, foo2, foo3, foo4;
    
    if( path.empty() ){ return; }
    foo1 = path.front();    path.pop_front();    foo2 = path.front();    path.pop_front();
    if( path.empty() ){ path.push_back(foo1); path.push_back(foo2); return; }
    foo3 = path.front();    path.pop_front();    foo4 = path.front();    path.pop_front();

    if( foo2 != foo3 ){ throw "The path is not a real path..."; }
    while( !path.empty() ){
      if( foo1 != foo4 ){
	cpath.push_back(foo1);	cpath.push_back(foo2);
	foo1=foo3;		foo2=foo4;
      } else {
	if( cpath.empty() ){
	  if( path.empty() ){ path = cpath; return; }
	  foo1 = path.front();    path.pop_front();      foo2 = path.front();    path.pop_front();	  
	}else{
	  foo2 = cpath.back();    cpath.pop_back();    foo1 = cpath.back();    cpath.pop_back();
	}
      }
      if( path.empty() ){ cpath.push_back(foo1); cpath.push_back(foo2); path = cpath; return;}
      foo3 = path.front();    path.pop_front();      foo4 = path.front();    path.pop_front();
    } // while
    cpath.push_back(foo1);	cpath.push_back(foo2);		cpath.push_back(foo3);	cpath.push_back(foo4);
    
    path = cpath;
  }
  
  double error() {
    clean_cycle_path();
    clean_path();
    const int L = path.size();
    if (L == 0){ return 0.; }
    return (2.44948974278 * R2D(getRotationMagnitude(R)) / sqrt(L));  // 2.44948974278 = sqrt(3*2)
  }
  
  void print_path() const {
    for (std::list<IndexT>::const_iterator it=path.begin(); it!=path.end(); ++it){ std::cout << ' ' << *it;}
  }
  
  Mat3 get_R() const { return R; }
  Vec3 get_t() const { return t; }
  std::list<IndexT> get_path() const { return path; }
  Pair get_path_ends() const { return std::make_pair( path.front(), path.back() ); }
  
  private:
  
    Mat3 R;
    Vec3 t;
    std::list<IndexT> path;  
  
};

////////////////////////////////////////////////////////////////////////////////
//                                    TREE                                    //
//////////////////////////////////////////////////////////////////////////////// 


struct Tree {

  public:
    Tree(){}
    Tree(const IndexT & i){
	indexT_set.insert(i);
	global_transformation[i] = Transformation();
	global_transformation_needed = false;
    }
    
    void print() const {
      for (std::set<Pair>::const_iterator iter = pair_set.begin(); iter != pair_set.end(); ++iter)
	std::cout << iter->first << "-" << iter->second << " ";
    }
    void print_node() const {
      for (std::set<IndexT>::const_iterator iter = indexT_set.begin(); iter != indexT_set.end(); ++iter)
	std::cout << *iter << " ";
    }
    
    int size() const { return indexT_set.size(); };
    
    void insert( const Pair & p ){
	pair_set.insert(p);
	indexT_set.insert(p.first);
	indexT_set.insert(p.second);
	global_transformation_needed = true;
    }
    void insert( const Pair & p, Transformation T ){
	pair_set.insert(p);	indexT_set.insert(p.first);	indexT_set.insert(p.second);
	if( global_transformation.find(p.first) != global_transformation.end() ){
	  T.compose_right( global_transformation.at(p.first) );
	  global_transformation[p.second] = T;
	} else {
	  T.inverse();
  	  T.compose_right( global_transformation.at(p.second) );
	  global_transformation[p.first] = T;
	}
	
    }
    bool contains( const Pair & p ) const {
      return ( (pair_set.find(p)!=pair_set.end()) ||
               (pair_set.find( std::make_pair(p.second,p.first) )!=pair_set.end()) );
    }
    bool contains( const IndexT & i ) const {
      return ( indexT_set.find(i) != indexT_set.end() );
    }
        
    double transformation_error( Transformation T ) const {
      Pair st = T.get_path_ends();
      if( contains(st.first) && contains(st.second) ){
	T.compose_right( global_transformation.at(st.first) );
	T.compose_left_rev( global_transformation.at(st.second) );
	const double err = T.error();
	return err;
      }
      throw "The transformation's error can't be computed.";
    }
    
    void initialise_transformation( const RelativeInfo_Map & relatives_Rt );
    
    std::set<IndexT> get_indexT_set() const { return indexT_set; }
    std::set<Pair> get_pair_set() const { return pair_set; }
    std::map<IndexT,Transformation> get_global_transformation() const { return global_transformation; }
    
  private:
    std::set<Pair> pair_set;
    std::set<IndexT> indexT_set;
    std::map<IndexT,Transformation> global_transformation;
    bool global_transformation_needed;
  
};

////////////////////////////////////////////////////////////////////////////////
//                              CONSISTENCY_MAP                               //
//////////////////////////////////////////////////////////////////////////////// 


struct Consistency_Map {

  public:
    Consistency_Map(){}
    
    Consistency_Map( const RelativeInfo_Map & map_relat ){
      for(RelativeInfo_Map::const_iterator iter = map_relat.begin(); iter != map_relat.end(); ++iter) {
	edge_map[iter->first] = std::make_pair(.5,0);
	global_edge_map[iter->first] = std::make_pair(.5,0);
	nbtree[iter->first] = 0;
      }
    }
    
    void erase( const Pair & p ){
	edge_map.erase(p);
	global_edge_map.erase(p);
	nbtree.erase(p);
    }
    
    void initialize( const Pair & p ){
	edge_map[p] = std::make_pair(.5,0);
	global_edge_map[p] = std::make_pair(.5,0);
	nbtree[p] = 0;
    }
    
    double get_coherence( const Pair & edge ) const {
      if( global_edge_map.find(edge) != global_edge_map.end() ){
	  return global_edge_map.at(edge).first;
      } else {
	  return global_edge_map.at(std::make_pair(edge.second, edge.first)).first;
      }
    }
    
    int get_nbtree( const Pair & edge ) const {
      if( nbtree.find(edge) != nbtree.end() ){
	  return min(nbtree.at(edge),5);
      } else {
	  return min(nbtree.at(std::make_pair(edge.second, edge.first)),5);
      }
    }
    
    void update_coherence( const Pair & edge, double coherence, double & returned_coherence ){
      Pair fooedge;
      
      if( edge_map.find(edge) != edge_map.end() )
	fooedge = edge;
      else
	fooedge = std::make_pair(edge.second, edge.first);

      const std::pair<double,int> fooloc = edge_map.at(fooedge);
      const std::pair<double,int> fooglob = global_edge_map.at(fooedge);
      
      edge_map[fooedge] = std::make_pair( (fooloc.first * fooloc.second + coherence)/(fooloc.second+1) , fooloc.second+1 );
      returned_coherence = (fooloc.first * fooloc.second + coherence)/(fooloc.second+1);
      if( fooglob.second > 0 ){
	const double foo = .5 + .5*fooglob.first;
	returned_coherence = foo * returned_coherence + (1-foo) * fooglob.first;
	if( fooglob.first > .9 ){
	  returned_coherence = .5 * returned_coherence;
	}
      }
    }
    
    void update_coherence( const Pair & edge, double coherence ){
      if( edge_map.find(edge) != edge_map.end() ){
	const std::pair<double,int> foo = edge_map.at(edge);
	edge_map[edge] = std::make_pair( (foo.first * foo.second + coherence)/(foo.second+1) , foo.second+1 );
      } else {
	const std::pair<double,int> foo = edge_map.at(std::make_pair(edge.second, edge.first));      
	edge_map[std::make_pair(edge.second, edge.first)] = std::make_pair( (foo.first * foo.second + coherence)/(foo.second+1) , foo.second+1 );
      }
    }
    
    void update_tree( const Tree & tree ){
      const std::set<Pair> pair_set = tree.get_pair_set();
      for(std::set<Pair>::const_iterator iter = pair_set.begin(); iter != pair_set.end(); ++iter) {
	if( pair_set.find(*iter) != pair_set.end() )
	  nbtree[*iter] += 1;
	else
	  nbtree[std::make_pair(iter->second, iter->first)] += 1;
      }
    }
    
    void push_coherence(){
      for(Edge_Consist_Map::iterator iter = edge_map.begin(); iter != edge_map.end(); ++iter) {
	const Pair edge = iter->first;
	const std::pair<double,int> foo = global_edge_map.at(edge);
	global_edge_map[edge] = std::make_pair( (foo.first * foo.second + iter->second.first)/(foo.second+1), foo.second+1 );
	iter->second = std::make_pair( .5, 0 );
      }
    }
    
    double mix_proba( const Pair & edge, const double & proba ) const {      
      if( global_edge_map.find(edge) != global_edge_map.end() ){
	if( 0 == global_edge_map.at(edge).second )
	  return proba;
	else
	  return .5*proba + .5*global_edge_map.at(edge).first;
      }
      else {
	if( 0 == global_edge_map.at(std::make_pair(edge.second, edge.first)).second )
	  return proba;
	else
	  return .5*proba + .5*global_edge_map.at(std::make_pair(edge.second, edge.first)).first;
      }
    }
    
    void clear_local(){
      for(Edge_Consist_Map::iterator iter = global_edge_map.begin(); iter != global_edge_map.end(); ++iter)
	edge_map[iter->first] = std::make_pair(.5,0);
    }
    
    void clear(){
      for(Edge_Consist_Map::iterator iter = global_edge_map.begin(); iter != global_edge_map.end(); ++iter) {
	edge_map[iter->first] = std::make_pair(.5,0);
	global_edge_map[iter->first] = std::make_pair(.5,0);
	nbtree[iter->first] = 0;
      }
    }
    
  private:
    typedef std::map<Pair, std::pair<double,int>> Edge_Consist_Map;
    Edge_Consist_Map edge_map;
    Edge_Consist_Map global_edge_map;
    std::map<Pair, int> nbtree;
    
};

////////////////////////////////////////////////////////////////////////////////
//                          Global SfM Graph Cleaner                          //
//////////////////////////////////////////////////////////////////////////////// 

class GlobalSfM_Graph_Cleaner
{  
  public:
    
  ////////// // // /  /    /       /          /       /    /  / // // //////////
  //                               CONSTRUCTORS                               //
  ////////// // // /  /    /       /          /       /    /  / // // //////////
  
    
    GlobalSfM_Graph_Cleaner(const RelativeInfo_Map & map_relat) : relatives_Rt(map_relat) {
      // Build SetNodes
      std::set<IndexT> setNodes;
      std::cout << "Graph Initialisation:" << std::endl;
      for(RelativeInfo_Map::const_iterator iter = map_relat.begin(); iter != map_relat.end(); ++iter) {
	setNodes.insert(iter->first.first);	setNodes.insert(iter->first.second);      }
      
      // Add the edge
      for(RelativeInfo_Map::const_iterator iter = map_relat.begin(); iter != map_relat.end(); ++iter) {
	const Pair p = iter->first;
	// Initialisation of PairMapEdge
	adjacency_map[p.first].insert(p.second);
	adjacency_map[p.second].insert(p.first);
	edge_consistency_map.initialize(p);
      }  // Initialisation of conherencMap
      
	  
      initial_tree_size = 5;
      initial_tree_ransac_nb = 1000;
      
      default_max_error = 10000;
      
      proba_error_positivethres = 2.5;
      proba_error_negativethres = 7.5; 
      proba_error_pos_thres_over2pi = proba_error_positivethres / 180;
      formular_asuma_coef = .1;
      
      formule_normalization_coef = log(1+formular_asuma_coef/proba_error_pos_thres_over2pi)
				  /log(1+formular_asuma_coef/(1-proba_error_pos_thres_over2pi));
      
      proba_neutral = 1/(1+formule_normalization_coef);
      stop_descent_number = 1000;
      
      coherence_thres = .5;
      number_tree = 5;
      
      tree_increment_stop_thres = 0;
      dbtalk = 1;

      std::cout << "Proba_neutral = " << proba_neutral << std::endl;
      
      tikzfile.open ("graphe.tex");
      tikzfile << "\\documentclass[tikz]{standalone}\n"
	       << "\\def\\zdraw[#1][#2][#3]{\\ifnum#2>\\componentchoice\\ifnum#3=1\\draw[tp]\\else\\draw[fp]\\fi\\else\\ifnum#3=1\\draw[fn]\\else\\draw[tn]\\fi\\fi}"
	       << "\n\\begin{document}\n\\begin{tikzpicture}[\nscale=10,nbase/.style={circle,draw},\nbase/.style={thick},\ntp/.style={green!50!black,line width=.5},"
	       << "\ntn/.style={red!50!black,line width=.5,opacity=.1},\nfp/.style={line width=1,red},\nfn/.style={line width=.5,green!50!black,dotted},"
	       << "\nftree/.style={green!20!yellow,line width=10, opacity=.3},"
	       << "\ntree/.style={blue,line width=3,opacity=.3,dashed},\nitree/.style={blue!50!red,line width=4,opacity=.3}]\n";      
    }
    // \\def\\zdraw[#1][#2][#3]{\n\\ifnum#2>\\componentchoice\\def\\tbar{}\\else\\def\\tbar{dashed}\\fi\n\\ifnum#3=1\n\\pgfmathsetmacro{\\tikzfoo}{.3 + 0.7*#2/100}\n\\draw[line width = .5, green!50!black, opacity = \\tikzfoo, \\tbar]\n\\else\\pgfmathsetmacro{\\tikzfoo}{.05 + 0.95*#2/100}\n\\draw[line width=1,red, opacity = \\tikzfoo, \\tbar]\n\\fi}
    
    
    ~GlobalSfM_Graph_Cleaner(){
      tikzfile << "\\end{tikzpicture}\n\\end{document}";
      tikzfile.close();
    };
    
    
  ////////// // // /  /    /       /          /       /    /  / // // //////////
  //                             CHANGE PARAMETER                             //
  ////////// // // /  /    /       /          /       /    /  / // // //////////
    
  void set_initial_tree_size( int v ){ initial_tree_size = v; }
  void set_initial_tree_ransac_nb( int v ){ initial_tree_ransac_nb = v; }
  void set_default_max_error( double v ){ default_max_error = v; }
  void set_proba_error_positivethres( double v ){ proba_error_positivethres = v; }
  void set_proba_error_negativethres( double v ){ proba_error_negativethres = v; }
  void set_stop_descent_number( int v ){ stop_descent_number = v; }
  void set_number_tree( int v ){ number_tree = v; }
  void set_formular_asuma_coef( double v ){ formular_asuma_coef = v; }
  void set_tree_increment_stop_thres( double v ){ tree_increment_stop_thres = v; }
  
  void set_dbtalk( int v ){ dbtalk = v; }
  
  ////////// // // /  /    /       /          /       /    /  / // // //////////
  //                                  OTHERS                                  //
  ////////// // // /  /    /       /          /       /    /  / // // //////////
    
    RelativeInfo_Map run();

    
    /* DEBUG DEBUG DEBUG */    
    void set_position_groundtruth( const std::vector<Vec3> & v ){
      position_GroundTruth = v;
      for(Adjacency_map::const_iterator iter = adjacency_map.begin(); iter != adjacency_map.end(); ++iter) {
	Vec3 foo = position_GroundTruth[iter->first];
	tikzfile << "\\node[nbase] at (" << foo(0) <<  ","<< foo(2) << ")" << " (" << iter->first << ") " << "{" << iter->first << "};\n";
      }
    }
    
    void set_wrong_edges( const std::set<Pair> & s ){
      wrong_edges = s;
    }

    void disp_Graph(const string str) const{
      for(Adjacency_map::const_iterator iter = adjacency_map.begin();
	  iter != adjacency_map.end(); ++iter) {
	IndexT s = iter->first;
	const std::set<IndexT> & indexT_set = iter->second;
      
	for (std::set<IndexT>::const_iterator iterT = indexT_set.begin();
	    iterT != indexT_set.end(); ++iterT) {
	  IndexT t = *iterT;
	  if (s < t) {
	    tikzfile << "\\zdraw["<< str << "][" << round( 100 * edge_consistency_map.get_coherence(std::make_pair(s,t)) ) << "]";
	    if (wrong_edges.find(std::make_pair(s,t)) != wrong_edges.end() || wrong_edges.find(std::make_pair(t,s)) != wrong_edges.end() )
	      tikzfile << "[0]";
	    else
	      tikzfile << "[1]";
	    tikzfile << " (" << s << ") -- (" << t << "); ";
	  } // if ( s < t )
	} // for (iterT) in [indexT_set]
      } // for (iter) in [adjacency_map]
    }

    int count_negatives(){  
      int count = 0;
      for(RelativeInfo_Map::const_iterator iter = relatives_Rt.begin(); iter != relatives_Rt.end(); ++iter){
	if( wrong_edges.find(iter->first) != wrong_edges.end() ){count++;}
      }
      return count;
    }
    
    double formule_asuma( double nbpos, double nbneg ) const {
      const double error = pow( 1+formular_asuma_coef/(1-proba_error_pos_thres_over2pi),nbneg)
                         / pow( 1+formular_asuma_coef/  proba_error_pos_thres_over2pi  ,nbpos);
      if( 0 == nbneg && nbpos >= 2 )
	return 1 /( 1 + .2 * formule_normalization_coef * error );
      else
	return 1 /( 1 + formule_normalization_coef * error );
    }
    
    /* DEBUG DEBUG DEBUG */    
    
  private:
    
  ////////// // // /  /    /       /          /       /    /  / // // //////////
  //                                 VARIABLE                                 //
  ////////// // // /  /    /       /          /       /    /  / // // //////////
    
    typedef std::map<IndexT, std::set<IndexT>> Adjacency_map;

    RelativeInfo_Map relatives_Rt;
    Adjacency_map adjacency_map;
    Consistency_Map edge_consistency_map;
    
    /* DEBUG DEBUG DEBUG */
    mutable int dbtalk;
    std::vector<Vec3> position_GroundTruth;
    std::set<Pair> wrong_edges;
    mutable ofstream tikzfile;
    /* DEBUG DEBUG DEBUG */    
        
  ////////// // // /  /    /       /          /       /    /  / // // //////////
  //                                PARAMETERS                                //
  ////////// // // /  /    /       /          /       /    /  / // // //////////
  
  double default_max_error;
  
  int initial_tree_size;  
  int initial_tree_ransac_nb;
    
  double proba_error_positivethres;
  double proba_error_negativethres;  
  int stop_descent_number;
  double tree_increment_stop_thres;
  
  double formule_normalization_coef;
  double formular_asuma_coef;
  
  double coherence_thres;
  int number_tree;
  
  
  double proba_error_pos_thres_over2pi;
  double proba_neutral;
    
  
  ////////// // // /  /    /       /          /       /    /  / // // //////////
  //                                 FUNCTION                                 //
  ////////// // // /  /    /       /          /       /    /  / // // //////////
  
  double cycle_probability( double error ) const {
    if( error < proba_error_positivethres ){ return 1; }
    else if( error < proba_error_negativethres )
      return (proba_error_negativethres-error)/(proba_error_negativethres-proba_error_positivethres);
    else { return 0; }
  }
  
  //
  
  Tree generate_Random_Tree( const int size ) const;
  double tree_Consistency_Error( Tree & tree, double & nbpos, double & nb, int & min_count, int & nbtree ) const;
  Tree generate_Consistent_Tree( const int size, double & tree_error ) const;

  //
  
  std::map<IndexT,int> compute_distance(
		const Pair & pair,
		const std::set<IndexT> & tree_node_set) const;
  
  double edge_consistency_probability(
	      const IndexT & source_node, const IndexT & target_node,
	      Transformation T,
	      const Tree & tree,
	      const std::map< IndexT, int > & distance,
	      int & nb_total, bool & edge_skip, double & nbpos_s, double & nb_s);
	      
  void increase_Tree( Tree & tree );

  Tree final_tree();
  //
  
  void tree_update_coherence( const Tree & tree ){
    edge_consistency_map.update_tree( tree );
    for(RelativeInfo_Map::const_iterator iter = relatives_Rt.begin(); iter != relatives_Rt.end(); ++iter) {
      if ( tree.contains( iter->first.first ) && tree.contains(iter->first.second) ) {
	double foo = tree.transformation_error( Transformation( iter->second, iter->first ) );
	foo = cycle_probability( foo );
	if(dbtalk>1){std::cout << iter->first.first << "-" << iter->first.second << " : " << foo << "(";}
	edge_consistency_map.update_coherence( iter->first, foo, foo );
	if(dbtalk>1){std::cout << foo << ")" << std::endl;}
      }
    }
  }
  
  void clean_graph( const double & threshold ){
    IndexT source, target;
    if(dbtalk){ std::cout << "Edges removed:"; }
    
    std::set<Pair> removedPair;
     
    for(RelativeInfo_Map::const_iterator iter = relatives_Rt.begin(); iter != relatives_Rt.end(); ++iter) {
      source = iter->first.first;     target = iter->first.second;
      if ( edge_consistency_map.get_coherence(iter->first) < threshold ){
	if(dbtalk){ std::cout << " (" << iter->first.first << "," << iter->first.second << ")[" << edge_consistency_map.get_coherence(iter->first) << "]"; }
	removedPair.insert(iter->first);
      }
    }
    for(std::set<Pair>::const_iterator iter = removedPair.begin(); iter != removedPair.end(); ++iter) {
      relatives_Rt.erase(*iter);
      adjacency_map[iter->first].erase(iter->second);
      adjacency_map[iter->second].erase(iter->first);
      edge_consistency_map.erase(*iter);
    }
  }
  
  ////////// // // /  /    /       /          /       /    /  / // // //////////
};


} // namespace sfm
} // namespace openMVG

#endif // GLOBALSFM_GRAPH_CLEANER_HPP
