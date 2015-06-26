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

  Mat3 R;
  Vec3 t;
  std::list<IndexT> path;  
  
  void inverse(){
      R = R.transpose();
      t = - R*t;*
      path.reverse();
  }
  
  
  void load_transformation( const RelativeInfo_Map & relatives_Rt, IndexT s, IndexT t ){
    if (relatives_Rt.find(std::make_pair(s,t)) != relatives_Rt.end()){
      const std::pair<Mat3, Vec3> foo = relatives_Rt.at(std::make_pair(s,t));
      R = foo.first;      t = foo.second;
      path.clear();       path.push_back(s);    path.push_back(t);
    } else {
      const std::pair<Mat3, Vec3> foo = relatives_Rt.at(std::make_pair(t,s));
      R = foo.first.transpose();      t = foo.first.transpose()*foo.second;
      path.clear();       path.push_back(t);    path.push_back(s);
    }
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
        
    if( foo2 != foo3 )
      throw "The path is not a real path...";
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
    }
    cpath.push_back(foo1);	cpath.push_back(foo2);		cpath.push_back(foo3);	cpath.push_back(foo4);
    
    path = cpath;
  }
  
  double error(){
    clean_cycle_path();
    clean_path();
    const int L = path.size();
    // std::cout << "Erreur chemin = ";    for (std::list<IndexT>::const_iterator it=path.begin(); it!=path.end(); ++it){ std::cout << ' ' << *it;}
    // std::cout << std::endl;
    if (L == 0)
      return 0.; 
    return (2.44948974278 * R2D(getRotationMagnitude(R)) / sqrt(L));  // 2.44948974278 = sqrt(3*2)
  }
  
  void print_path(){
    for (std::list<IndexT>::const_iterator it=path.begin(); it!=path.end(); ++it){ std::cout << ' ' << *it;}
  }
};

////////////////////////////////////////////////////////////////////////////////
//                                    TREE                                    //
//////////////////////////////////////////////////////////////////////////////// 


struct Tree {

  public:
    Tree(const IndexT & i){
	indexT_set.insert(i);
	global_transformation[i] = Transformation();
    }
    
    void insert( const Pair & p ){
	pair_set.insert(p);
	indexT_set.insert(p.first);
	indexT_set.insert(p.second);
    }
    void insert( const Pair & p, Transformation T ){
	pair_set.insert(p);	indexT_set.insert(p.first);	indexT_set.insert(p.second);
	if( global_transformation.find(p) != global_transformation.end() ){
	  T.compose_right( global_transformation.at(p.first) );
	  global_transformation[p.second] = T;
	} else {
	  T.inverse();
  	  T.compose_right( global_transformation.at(p.second) );
	  global_transformation[p.first] = T;
	}
	
    }
    bool contains( const Pair & p ){
      return ( (pair_set.find(p)!=pair_set.end()) ||
               (pair_set.find( std::make_pair(p.second,p.first) )!=pair_set.end()) );
    }
    bool contains( const indexT & i ){
      return ( indexT_set.find(i) != pair_set.end() );
    }
    
    void initialise_transformation( const RelativeInfo_Map & relatives_Rt);
    
    std::set<IndexT> get_indexT_set(){ return indexT_set; }
    std::set<Pair> get_pair_set(){ return pair_set; }
    std::map<IndexT,Transformation> get_global_transformation(){ return global_transformation; }

  private:
    std::set<Pair> pair_set;
    std::set<IndexT> indexT_set;
    std::map<IndexT,Transformation> global_transformation;
  
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
      std::cout << "Graph initialisation:" << std::endl;
      for(RelativeInfo_Map::const_iterator iter = map_relat.begin(); iter != map_relat.end(); ++iter) {
	setNodes.insert(iter->first.first);	setNodes.insert(iter->first.second);      }
      
      // Add the edge
      for(RelativeInfo_Map::const_iterator iter = map_relat.begin(); iter != map_relat.end(); ++iter) {
	const Pair p = std::make_pair(iter->first.first, iter->first.second);
	// Initialisation of PairMapEdge
	adjacency_map[p.first].insert(p.second);
	adjacency_map[p.second].insert(p.first);
	edge_consistency_map[p] = .5;
      }  // Initialisation of conherencMap

      tikzfile.open ("graphe.tex");
      tikzfile << "\\documentclass[tikz]{standalone}\n"
	       << "\\def\\zdraw[#1][#2][#3]{\\ifnum#2=\\componentchoice\\ifnum#3=1\\draw[tp]\\else\\draw[fp]\\fi\\else\\ifnum#3=1\\draw[fn]\\else\\draw[tn]\\fi\\fi}"
	       << "\n\\begin{document}\n\\begin{tikzpicture}[\nscale=10,nbase/.style={circle,draw},\nbase/.style={thick},\ntp/.style={green!50!black,line width=.5},"
	       << "\ntn/.style={red!50!black,line width=.5,opacity=.1},\nfp/.style={line width=1,red},\nfn/.style={line width=.5,green!50!black,dotted},"
	       << "\ntree/.style={blue,line width=3,opacity=.3,dashed},\nitree/.style={blue!50!red,line width=4,opacity=.3}]\n";      
    }
    
    
    
    ~GlobalSfM_Graph_Cleaner(){
      tikzfile << "\\end{tikzpicture}\n\\end{document}";
      tikzfile.close();
    };
    
    
  ////////// // // /  /    /       /          /       /    /  / // // //////////
  //                             CHANGE PARAMETER                             //
  ////////// // // /  /    /       /          /       /    /  / // // //////////
    
    
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
	  } // if ( s < t )
	} // for (iterT) in [indexT_set]
      } // for (iter) in [adjacency_map]
    }
    /* DEBUG DEBUG DEBUG */    
    
  private:
    
  ////////// // // /  /    /       /          /       /    /  / // // //////////
  //                                 VARIABLE                                 //
  ////////// // // /  /    /       /          /       /    /  / // // //////////
    
    typedef std::map<IndexT, std::set<IndexT>> Adjacency_map;

    RelativeInfo_Map relatives_Rt;
    Adjacency_map adjacency_map;
    std::map<Pair,double> edge_consistency_map;
    
    
    
    /* DEBUG DEBUG DEBUG */    
    mutable int dbtalk;
    std::vector<Vec3> position_GroundTruth;
    std::set<Pair> wrong_edges;
    mutable ofstream tikzfile;
    /* DEBUG DEBUG DEBUG */    
    
  ////////// // // /  /    /       /          /       /    /  / // // //////////
  //                                PARAMETERS                                //
  ////////// // // /  /    /       /          /       /    /  / // // //////////
  
  ////////// // // /  /    /       /          /       /    /  / // // //////////
  //                                 FUNCTION                                 //
  ////////// // // /  /    /       /          /       /    /  / // // //////////
  
  Tree generate_Random_Tree( const int size ) const;
  double tree_Consistency_Error( const Tree & tree, int & nbpos, int & nbneg, int & min_count ) const;
  Tree generate_Consistent_Tree( const int size, double & tree_error ) const;

    
    
  ////////// // // /  /    /       /          /       /    /  / // // //////////
};


} // namespace sfm
} // namespace openMVG

#endif // GLOBALSFM_GRAPH_CLEANER_HPP
