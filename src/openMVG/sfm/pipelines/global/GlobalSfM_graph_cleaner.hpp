// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef GLOBALSFM_GRAPH_CLEANER_HPP
#define GLOBALSFM_GRAPH_CLEANER_HPP


namespace openMVG{
namespace globalSfM{


} // namespace globalSfM
} // namespace openMVG

#include "openMVG/sfm/sfm.hpp"
#include "openMVG/graph/graph.hpp"
#include "openMVG/multiview/rotation_averaging_common.hpp"

// #include <lemon/list_graph.h>

namespace openMVG{
namespace globalSfM{

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


class GlobalSfM_Graph_Cleaner
{
  public:
    
  ////////// // // /  /    /       /          /       /    /  / // // //////////
  //                               CONSTRUCTORS                               //    
  ////////// // // /  /    /       /          /       /    /  / // // //////////
    GlobalSfM_Graph_Cleaner(const RelativeInfo_Map & map_relat)
    :relatives_Rt(map_relat) {
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
	pairMapEdge[p] = g.addEdge(indexMapNode[p.first], indexMapNode[p.second]);
        std::cout << " (" << p.first << "," << p.second << ") ";      
	CohenrencyMap[p] = count; count+=1; }
      std::cout << "\ninitialisation end\n" << std::endl;
    }
  ////////// // // /  /    /       /          /       /    /  / // // //////////
    
    RelativeInfo_Map Run();
    
  ////////// // // /  /    /       /          /       /    /  / // // //////////
  //                              MISCELLANEOUS                               //
  ////////// // // /  /    /       /          /       /    /  / // // //////////
    // Display for TikZ
    void disp_graph(const string str) const;
    
  ////////// // // /  /    /       /          /       /    /  / // // //////////
  //                              DEBUG FUNCTION                              //
  ////////// // // /  /    /       /          /       /    /  / // // //////////
    void set_position_groundtruth( const std::vector<Vec3> & v ){
      position_GroundTruth = v;
    }
    
  void change_consistency( const int a, const int b){
    for (std::map<Pair,int>::iterator iter = CohenrencyMap.begin(); iter != CohenrencyMap.end();  ++iter) {
      if ( iter->second == a ) { iter->second = b; }
    }    
  }
    
  ////////// // // /  /    /       /          /       /    /  / // // //////////
  //                             PRIVATE FONCTION                             //    
  ////////// // // /  /    /       /          /       /    /  / // // //////////
  private:
    
    Cycles vec_cycles;
    RelativeInfo_Map relatives_Rt;
    std::map< Pair, int > CohenrencyMap;
    
    Graph g;
    std::map<Graph::Node, IndexT> nodeMapIndex;
    std::map<IndexT, Graph::Node> indexMapNode;
    std::map<Pair, Graph::Edge> pairMapEdge;
    
    // debug
    std::vector<Vec3> position_GroundTruth;
    
  ////////// // // /  /    /       /          /       /    /  / // // //////////
    
    void FindCycles();
    void RotationRejection(const double max_angular_error);    
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
	  
    typedef std::set<Pair> Tree;

    Tree generate_Random_Tree( const Graph & g, const int size ) const{
      Tree tree;
      // TODO
      return tree;
    };

    double tree_consistency_error( const Tree & t, const RelativeInfo_Map & relatives_Rt ) const{
      // TODO
      return 0;
    }
      
  ////////// // // /  /    /       /          /       /    /  / // // //////////
  
};


} // namespace globalSfM
} // namespace openMVG

#endif // GLOBALSFM_GRAPH_CLEANER_HPP
