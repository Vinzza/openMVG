// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef GLOBALSFM_GRAPH_CLEANER_HPP
#define GLOBALSFM_GRAPH_CLEANER_HPP


namespace openMVG{
namespace sfm{


} // namespace sfm
} // namespace openMVG

#include "openMVG/sfm/sfm.hpp"
#include "openMVG/graph/graph.hpp"
#include "openMVG/multiview/rotation_averaging_common.hpp"

// #include <lemon/list_graph.h>

namespace openMVG{
namespace sfm{

////////////////////////////////////////////////////////////////////////////////
//                                   Cycle                                    //
////////////////////////////////////////////////////////////////////////////////
struct Cycle
{
  Cycle(const std::vector<Pair> & ccycle)
    :cycle(ccycle)
  { }

  std::vector<Pair> cycle;
  

  bool contain(const Pair & edge) const {
    return (std::find(cycle.begin(), cycle.end(), edge) != cycle.end());
  }
  
  void reverse() {
    std::reverse( cycle.begin()+1, cycle.end() );
    for ( std::vector<Pair>::iterator iter = cycle.begin(); iter != cycle.end(); ++iter ){
      Pair & p = *iter;
      p = std::make_pair( p.second, p.first );
    }
  }
  
  float errToIdentity( const RelativeInfo_Map & relatives_Rt ) const{
    Mat3 rot_To_Identity = Mat3::Identity();
    
    for ( std::vector<Pair>::const_iterator iter=cycle.begin(); iter!=cycle.end(); ++iter ){
      const Pair p = *iter;
      Mat3 RIJ;
      if (relatives_Rt.find(p) != relatives_Rt.end())
	RIJ = relatives_Rt.at(p).first;
      else
	RIJ = relatives_Rt.at( make_pair(p.second,p.first) ).first.transpose();
      rot_To_Identity = RIJ * rot_To_Identity;
    }
    return static_cast<float>(R2D(getRotationMagnitude(rot_To_Identity)));
  }
  
  ////////// // // /  /    /       /          /       /    /  / // // //////////

  friend bool operator==(const Cycle& c1, Cycle& c2)  {
    std::vector<Pair> cycle2 = c2.cycle;
    Pair edge0 = c1.cycle[0];
    std::vector<Pair>::iterator p = std::find( cycle2.begin(), cycle2.end(), edge0 );
    if ( cycle2.end() != p ) {
      std::rotate( cycle2.begin(), p, cycle2.end() );
      return cycle2 == c1.cycle;
    }
    p = std::find( cycle2.begin(), cycle2.end(), std::make_pair( edge0.second, edge0.first ) );
    if ( cycle2.end() != p ){
      std::rotate( cycle2.begin(), p, cycle2.end() );
      c2.reverse();
      return cycle2 == c1.cycle;
    }
    return false;
  }
  
  friend bool operator!=(const Cycle& c1, Cycle& c2)  {
    return !(c1 == c2);
  }

  
  friend std::ostream & operator<<(std::ostream & os, const Cycle & c)
  {
    os << c.cycle[0].first;
    for (int i = 0; i < c.cycle.size(); ++i){
      os << "--" << c.cycle[i].second << " ";
    }
    os << std::endl;
    return os;
  }
  
}; // Struct Cycle


typedef std::vector< Cycle > Cycles;
  
////////////////////////////////////////////////////////////////////////////////
//                          Global SfM Graph Cleaner                          //
//////////////////////////////////////////////////////////////////////////////// 

typedef lemon::ListGraph Graph;
typedef std::set<Pair> Tree;
typedef std::pair<Mat3, Vec3> Transformationz;

struct Transformation {
public:

  Transformation():R(Mat3::Identity()),t(Vec3::Zero()),length(0) {}
  Transformation(std::pair<Mat3,Vec3> np): R(np.first), t(np.second), length(1) {}
  Transformation(Mat3 nr, Vec3 nt, int nl): R(nr), t(nt), length(nl) {}

  // Constructor for the inverse transformation
  Transformation(std::pair<Mat3,Vec3> np, string foo): R(np.first.transpose()), t(- np.first.transpose()*np.second), length(1) {}

  Mat3 R;
  Vec3 t;
  int length;  
  
  void compose_left( const Transformation & T ){
    R = T.R * R;    t = T.R * t + T.t;    length = length + T.length;
  }
  void compose_left( const std::pair<Mat3, Vec3> & T ){
    R = T.first * R;    t = T.first * t + T.second;    length = length + 1;
  }
  void compose_right( const Transformation & T ){
    R = R * T.R;    t = R * T.t + t;    length = length + T.length;
  }
  void compose_right( const std::pair<Mat3, Vec3> & T ){
    R = R * T.first;    t = R * T.second + t;    length = length + 1;
  }
  void compose_left_rev( const Transformation & T ){
    R = T.R.transpose() * R;    t = T.R.transpose() * (t - T.t);    length = length + T.length;
  }
  
  double error(){
    return R2D(getRotationMagnitude(R));
  }
};

class GlobalSfM_Graph_Cleaner
{  /*
  CONSTRUCTORS:
    GlobalSfM_Graph_Cleaner( RelativeInfo_Map )

  UTILS:
    inverse_transformation( T ) : return the inverse of the transformation
    get_transformation( Pair ) : return the local transformation of the Pair    
    compose_transformation( T1, T2 ) : return the composition : T1 o T2
    
    change_Cohenrence( int, int )
    - - - - - - - - -( Pair, Pair ) : set the two Pair to the same consistent value
    update_Coherence( tree ) : update the CohenrenceMap
    principale_Coherence() : return the int corresponding to the biggest coherent part
                     
    sequential_Tree_Reconstruction( tree, globalTransformation ) : build the globalTransformation sequentially with the tree
    
  PROCESS:
    run() : run the cleaning process
    
    addCycleToMap( Cycle, RelativeInfo_Map ) : add the cycle to the RelativeInfo_Map
    findCycles() : return a vector of Cycles (only triplets for now)
    rotationRejection( max_angular_error, Cycles ) : keep only the edges which are in a consistent cycle.
    
    generate_Random_Tree( size(int) ) : return a random tree with a given size
    tree_Consistency_Error( tree, nbpos, nbneg ) : return the error and set nbneg/nbpos to the number of incoherent/coherent cycles.
    generate_Consistent_Tree( size(int) ) : return a consistent tree with a given size
    edge_Consistency_Error( pair, tree, tree_nodes  ) : return the consistency error made by the pair according to the tree
    increase_Tree( tree ) : add a consistent edge to the tree
    
  DEBUG:
    set_position_groundtruth( vector<Vect3> ) : set the position groundtruth for the disp_Graph function
    set_wrong_edges( set<Pair> ) : set the groundtruth wrond_edges_set for the disp_Graph function
    disp_Graph( string ) : display the graph in a TikZ format               
    
  */
  public:
    
  ////////// // // /  /    /       /          /       /    /  / // // //////////
  //                               CONSTRUCTORS                               //    
  ////////// // // /  /    /       /          /       /    /  / // // //////////
    GlobalSfM_Graph_Cleaner(const RelativeInfo_Map & map_relat, const double & error_thres = 5.)
    :relatives_Rt(map_relat), error_thres(error_thres) {
      // Build SetNodes
      std::set<IndexT> setNodes;
      std::cout << "Graph initialisation:" << std::endl;
      for(RelativeInfo_Map::const_iterator iter = map_relat.begin(); iter != map_relat.end(); ++iter) {
	setNodes.insert(iter->first.first);	setNodes.insert(iter->first.second);      }
      std::cout << "Nodes:";
      // Build association beetween Node and indexT
      for (std::set<IndexT>::const_iterator iter = setNodes.begin(); iter != setNodes.end(); ++iter) {
	indexMapNode[*iter] = g.addNode();      nodeMapIndex[ indexMapNode[*iter] ] = *iter;
	std::cout << " " << g.id(indexMapNode[*iter]);
      }      
      // Add the edge
      int count = 0;
      std::cout << "\nEdges:";
      for(RelativeInfo_Map::const_iterator iter = map_relat.begin(); iter != map_relat.end(); ++iter) {
	const Pair p = std::make_pair(iter->first.first, iter->first.second);
	// Initialisation of PairMapEdge
	pairMapEdge[p] = g.addEdge(indexMapNode[p.first], indexMapNode[p.second]);
	adjacency_map[p.first].insert(p.second);
	adjacency_map[p.second].insert(p.first);
	std::cout << " (" << p.first << "," << p.second << ") ";      
	CohenrenceMap[p] = count; count+=1; }  // Initialisation of conherencMap
      std::cout << "\ninitialisation end\n" << std::endl;
      
      // DEBUG
      dbtalk = true;
      tikzfile.open ("graphe.tex");
      tikzfile << "\\documentclass[tikz]{standalone}\n"
	       << "\\def\\zdraw[#1][#2][#3]{\\ifnum#2=\\componentchoice\\ifnum#3=1\\draw[tp]\\else\\draw[fp]\\fi\\else\\ifnum#3=1\\draw[fn]\\else\\draw[tn]\\fi\\fi}"
	       << "\n\\begin{document}\n\\begin{tikzpicture}[\nscale=10,nbase/.style={circle,draw},\nbase/.style={thick},\ntp/.style={green!50!black,line width=.5},"
	       << "\ntn/.style={red!50!black,line width=.5,opacity=.1},\nfp/.style={line width=1,red},\nfn/.style={line width=.5,green!50!black,dotted},"
	       << "\ntree/.style={blue,line width=3,opacity=.3,dashed},\nitree/.style={blue!70!green,line width=3,opacity=.5,dashed}]\n";      
    }
    
    ~GlobalSfM_Graph_Cleaner(){
      tikzfile << "\\end{tikzpicture}\n\\end{document}";
      tikzfile.close();
    };
  ////////// // // /  /    /       /          /       /    /  / // // //////////
      
    RelativeInfo_Map run();
    
  ////////// // // /  /    /       /          /       /    /  / // // //////////
  //                              MISCELLANEOUS                               //
  ////////// // // /  /    /       /          /       /    /  / // // //////////
    
    void change_Cohenrence( const int a, const int b ){ // TODO : improve speed
      if ( a != b ) {
	for (std::map<Pair,int>::iterator iter = CohenrenceMap.begin(); iter != CohenrenceMap.end();  ++iter) {
	  if ( iter->second == a ) { iter->second = b; }
	}
      }
    }
    int principale_Coherence() const{ // TODO : improve speed
      std::map<int,int> foo;
      for (std::map<Pair,int>::const_iterator iter = CohenrenceMap.begin(); iter != CohenrenceMap.end();  ++iter) {
	if( foo.find(iter->second) == foo.end() )
	  foo[iter->second] = 1;
	else
	  foo[iter->second] = foo[iter->second] + 1;
      }
      int maxint=0, maxkey;
      for (std::map<int,int>::iterator iter = foo.begin(); iter != foo.end();  ++iter) {
	if( iter->second > maxint ){
	  maxint = iter->second;
	  maxkey = iter->first;
	}
      }
      return maxkey;
    }
    void change_Cohenrence( const Pair a, const Pair b ){
      int ca;
      if (CohenrenceMap.find(a) != CohenrenceMap.end())	{ ca = CohenrenceMap.at(a); }
      else { ca = CohenrenceMap.at(std::make_pair(a.second, a.first)); }

      if (CohenrenceMap.find(b) != CohenrenceMap.end())
	change_Cohenrence( ca, CohenrenceMap.at(b));
      else
	change_Cohenrence( ca, CohenrenceMap.at(std::make_pair(b.second, b.first)));
    }
    
  ////////// // // /  /    /       /          /       /    /  / // // //////////
  //                             PRIVATE FONCTION                             //    
  ////////// // // /  /    /       /          /       /    /  / // // //////////
  //private: TODO TODO TODO TODO
    // Relative transformation Map
    RelativeInfo_Map relatives_Rt;
    // Conherence Map code the cohenrent set of edges. 
    std::map< Pair, int > CohenrenceMap;
    
    // Graph
    Graph g;
    // Map to switch from IndexT to Node and Node to IndexT
    std::map<Graph::Node, IndexT> nodeMapIndex;
    std::map<IndexT, Graph::Node> indexMapNode;
    // Map to switch from Pair to Edge
    std::map<Pair, Graph::Edge> pairMapEdge;

    // Adjacency map
    typedef std::map<IndexT, std::set<IndexT>> Adjacency_map;
    Adjacency_map adjacency_map;
    
    //
    double error_thres;
    
  ////////// // // /  /    /       /          /       /    /  / // // //////////
    
    Transformation inverse_transformation( const Transformation & T ) const{
	return Transformation( T.R.transpose(), - T.R.transpose()*T.t, T.length );
    }
    
    Transformation get_transformation( const Pair & pair ) const{
      if (relatives_Rt.find(pair) != relatives_Rt.end())
	return Transformation(relatives_Rt.at(pair));
      else
	return Transformation(relatives_Rt.at( make_pair(pair.second,pair.first) ), "inverse");
    }
    
    Transformation compose_transformation( const Transformation & T1,
					   const Transformation & T2 ) const{
      return Transformation( T1.R * T2.R, T1.R * T2.t + T1.t, T1.length + T2.length );
    }
    
  ////////// // // /  /    /       /          /       /    /  / // // //////////
    
    Cycles findCycles();
    
    void rotationRejection(const double max_angular_error, Cycles & vec_cycles);    
    
    void addCycleToMap( Cycle & c, RelativeInfo_Map & relatives_Rt_new ){
      for ( std::vector<Pair>::iterator iter = c.cycle.begin(); iter != c.cycle.end(); ++iter ){
	const Pair ij = *iter;	
	if (relatives_Rt.find(ij) != relatives_Rt.end())
	  relatives_Rt_new[ij] = relatives_Rt.at(ij);
	else {
	  const Pair ji = std::make_pair( ij.second, ij.first );
	  relatives_Rt_new[ji] = relatives_Rt.at(ji);
	}
      }
    }
    
    ////////////////////////////////////////////////////////////////////////////
    //                                 TREES                                  //
    ////////////////////////////////////////////////////////////////////////////


    Tree generate_Random_Tree( const int size ) const;
    void sequential_Tree_Reconstruction(
		const Tree & tree,
		std::map<IndexT,Transformation> & globalTransformation) const;
    double tree_Consistency_Error( const Tree & tree, int & nbpos, int & nbneg ) const;
    
    Tree generate_Consistent_Tree( const int size ) const;

    void update_Coherence( const Tree & tree );
    
    double edge_Descent_Error( const IndexT & t,
	  Transformation T, const std::map< IndexT, int > & distance,
	  const std::map<IndexT,Transformation> & globalTransformation,
	  int & nbpos_i, int & nbneg_i) const;
    double edge_Consistency_Error( const Pair & pair, const Tree & tree, const std::set<IndexT> & tree_nodes, int & nbpos, int & nbneg  ) const;
    
    void increase_Tree( Tree & tree ) const;
    
  ////////// // // /  /    /       /          /       /    /  / // // //////////
  //                              DEBUG FUNCTION                              //
  ////////// // // /  /    /       /          /       /    /  / // // //////////
    mutable bool dbtalk;
    std::vector<Vec3> position_GroundTruth;
    std::set<Pair> wrong_edges;
    mutable ofstream tikzfile;
    
  // Display for TikZ
    void disp_Graph(const string str) const;
    
  public:
    
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
  ////////// // // /  /    /       /          /       /    /  / // // //////////
};


} // namespace sfm
} // namespace openMVG

#endif // GLOBALSFM_GRAPH_CLEANER_HPP
