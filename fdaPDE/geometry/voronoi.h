// This file is part of fdaPDE, a C++ library for physics-informed
// spatial and functional data analysis.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __VORONOI_H__
#define __VORONOI_H__

#include <unordered_map>
#include <map>
#include <unordered_set>
#include <vector>
#include <list>

#include "../utils/symbols.h"
#include "utils.h"

#include <cmath>
#include <unsupported/Eigen/SparseExtra>

using UInt = int;

namespace fdapde {
  namespace core {

    // adaptor adapting a (Delaunay) triangulation to its dual Vornoi diagram
    template <typename Triangulation> class Voronoi;
    
    template <> class Voronoi<Triangulation<3, 3>> {
    public:
      static constexpr UInt local_dim = Triangulation<3,3>::local_dim;
      static constexpr UInt embed_dim = Triangulation<3,3>::embed_dim;
      
      Voronoi() = default;
      Voronoi(const Triangulation<3, 3> & mesh) : mesh_(&mesh) {   // constructs voronoi diagram from Delanoy triangulation	
	// Set the sites
	this->nodes_ = mesh_->nodes();
	
	// Build the actual unrestricted voronoi
	this->build_unrestricted_voronoi();
	
	// clip the surface with respect to the Voronoi, and update the cells
	this->cut_surface_to_cells();
	
	// Finally make a cycle over the Voronoi cells to check which are boundary
	this->nodes_markers_.resize(this->nodes_.rows());
	for(auto iter = this->cells_.begin(); iter != this->cells_.end(); iter++){
	  if((iter->second).on_boundary()){
	    this->nodes_markers_.set(iter->first);
	  }
	}
      }
      
      // Structure that stores a segment in the Voronoi<3,3> notation
      struct segment {
	SVector<3> p1, p2;
	UInt neighbor;
      };

      class ConvexPolygon{
      public:
	ConvexPolygon()=default;
	ConvexPolygon(Triangulation<2,3>::cell_iterator face, const DMatrix<double> & nodes_ , std::unordered_map<int,int> surface_to_mesh, SVector<3> tet_mean){
	  // NOTE: should take an iterator
	  // To be implemented, ask Palummo; Main problem is to retrieve the adiacent polygons information on the surface!! INDEED WE CAN DO THIS USIG THE SURFACE METHOD; NO MORE!!, simply will be referred to the tringles of the sruface (the methods that will use this need to take this into account) (We will also need a transformation to internal map at the end, since the closeness is local)
	  // I need to deepen more this prolem
	  // Idea: before use the indices of surface mesh boundary, only at the last, build your map!!
	  // Build the list of points
	  std::vector<SVector<3>> points;
	  std::vector<UInt> adiacent_faces;
	  
	  UInt p0,p1,p2;
	  
	  auto edge_0 = face->edge(0);
	  auto edge_1 = face->edge(1);
	  auto edge_2 = face->edge(2);
	  
	  p0 = surface_to_mesh.at(edge_0.node_ids()[0]);
	  p1 = surface_to_mesh.at(edge_0.node_ids()[1]);
	  
	  auto neighbors  = edge_0.adjacent_cells(); 
	    for(auto neigh : neighbors){
	      if(neigh != face->id()){
		adiacent_faces.push_back(-1-neigh); // CHECK THAT THERE EXIST ONLY A CLOSE FRIEND
	      }
	    }
	  
	  if(edge_1.node_ids()[0] == edge_0.node_ids()[0]){
	    p2 = surface_to_mesh.at(edge_1.node_ids()[1]);
	    neighbors  = edge_2.adjacent_cells(); // Needs to be edge 2 because the edges are not sorted!
	    for(auto neigh : neighbors){
	      if(neigh != face->id()){
		adiacent_faces.push_back(-1-neigh); // CHECK THAT THERE EXIST ONLY A CLOSE FRIEND
	      }
	    }
	    
	    neighbors  = edge_1.adjacent_cells(); // Now we fiinish the adiacent polygons with the adiacenzy of edge 1
	    for(auto neigh : neighbors){
	      if(neigh != face->id()){
		adiacent_faces.push_back(-1-neigh); // CHECK THAT THERE EXIST ONLY A CLOSE FRIEND
	      }
	    }
	  }
	  
	  if(edge_1.node_ids()[1] == edge_0.node_ids()[0]){
	    p2 = surface_to_mesh.at(edge_1.node_ids()[0]);
	    neighbors  = edge_2.adjacent_cells(); // Needs to be edge 2 because the edges are not sorted!
	    for(auto neigh : neighbors){
	      if(neigh != face->id()){
		adiacent_faces.push_back(-1-neigh); // CHECK THAT THERE EXIST ONLY A CLOSE FRIEND
	      }
	    }
	    
	    neighbors  = edge_1.adjacent_cells(); // Now we fiinish the adiacent polygons with the adiacenzy of edge 1
	    for(auto neigh : neighbors){
	      if(neigh != face->id()){
		adiacent_faces.push_back(-1-neigh); // CHECK THAT THERE EXIST ONLY A CLOSE FRIEND
	      }
	    }
	  }
	  
	  if(edge_1.node_ids()[0] == edge_0.node_ids()[1]){
	    p2 = surface_to_mesh.at(edge_1.node_ids()[1]);
	    neighbors  = edge_1.adjacent_cells(); // Needs to be edge 1 because the edges are sorted!
	    for(auto neigh : neighbors){
	      if(neigh != face->id()){
		adiacent_faces.push_back(-1-neigh); // CHECK THAT THERE EXIST ONLY A CLOSE FRIEND
	      }
	    }
	    
	    neighbors  = edge_2.adjacent_cells(); // Now we fiinish the adiacent polygons with the adiacenzy of edge 2
	    for(auto neigh : neighbors){
	      if(neigh != face->id()){
		adiacent_faces.push_back(-1-neigh); // CHECK THAT THERE EXIST ONLY A CLOSE FRIEND
	      }
	    }
	  }
	  
	  
	  if(edge_1.node_ids()[1] == edge_0.node_ids()[1]){
	    p2 = surface_to_mesh.at(edge_1.node_ids()[0]);
	    neighbors  = edge_1.adjacent_cells(); // Needs to be edge 2 because the edges are not sorted!
	    for(auto neigh : neighbors){
	      if(neigh != face->id()){
		adiacent_faces.push_back(-1-neigh); // CHECK THAT THERE EXIST ONLY A CLOSE FRIEND
	      }
	    }
	    
	    neighbors  = edge_2.adjacent_cells(); // Now we fiinish the adiacent polygons with the adiacenzy of edge 1
	    for(auto neigh : neighbors){
	      if(neigh != face->id()){
		adiacent_faces.push_back(-1-neigh); // CHECK THAT THERE EXIST ONLY A CLOSE FRIEND
	      }
	    }
	  }
	  
	  
	  // manual insertion
	  //for(auto i=0; i<3; i++){ // WRONG
	  //  auto edge = face->edge(i);
	  //  auto neighbors  = edge.adjacent_cells(); 
	  //  for(auto neigh : neighbors){
	  //    if(neigh != face->id()){
	//	adiacent_faces.push_back(-1-neigh); // CHECK THAT THERE EXIST ONLY A CLOSE FRIEND
	 //     }
	  //  }
	    //adiacent_faces.push_back(-1-neigh); // Negative index to comply with notation of ConvexPolygons
	  //}
	    
	  //UInt p0 = surface_to_mesh.at(face->node_ids()[0]);
	  //UInt p1 = surface_to_mesh.at(face->node_ids()[1]);
	  //UInt p2 = surface_to_mesh.at(face->node_ids()[2]);
      
	  //iter++;
	  
	  //UInt p2 = iter->node_ids()(0);
	  //if(p2==p1 || p2==p0){
	  //	UInt p2 = iter->node_ids()(1);
	  //}
	  
	  points.push_back(nodes_.row(p0));
	  points.push_back(nodes_.row(p1));
	  points.push_back(nodes_.row(p2));
	  
	  SVector<3> v1 = nodes_.row(p1)-nodes_.row(p0);
	  SVector<3> v2 = nodes_.row(p2)-nodes_.row(p0);
	  
	  SVector<3> prova = v1.cross(v2);
	  SVector<3> diff = tet_mean - face->supporting_plane().project(tet_mean);
	  
	  double boh = prova.dot(diff);
	  
	  if((v1.cross(v2)).dot(tet_mean - face->supporting_plane().project(tet_mean))<0){ // I use the tetrahedron mean to avoid the possible wrong ordering of the normal of the supporting_plane
	    std::reverse(points.begin(), points.end());
	    std::reverse(adiacent_faces.begin(), adiacent_faces.end());
	  }
	 
	  this->points_ = points;
	  this->adiacent_polygons_ = adiacent_faces;
	  this->adiacent_cell_ = -1;
	  this->closed_=true;
	  SVector<3> mean = points_[0];
	  UInt den =this->points_.size(); 

	  for(std::size_t i = 1; i<this->points_.size(); i++){
	    mean = mean + this->points_[i];
	  }
	  // Build the centroid: Here the polygon is always closed, so it is just the mean of the points
	  this->centroid_ = mean/den;

	  // RECHEK HERE
	}
	void set_points(std::vector<SVector<3>> set_points){
	  this->points_ = set_points;
	  
	  if(set_points.size() > 0){
	    UInt den = set_points.size();
	    SVector<3> mean = set_points[0];

	    for(std::size_t i = 1; i<set_points.size(); i++){
	      mean = mean + set_points[i];
	    }

	    // Build the centroid
	    centroid_ = mean/den; // Note: the centroid may fall out of the volume, and therefore we need to perform further checks afterwards.
	  
	  }
	
	  if(adiacent_polygons_.size() <  this->points_.size()){
	    this->closed_ = false;
	  }else{
	    this->closed_ = true;
	  }
	}
	void set_adiacent_polygons(std::vector<UInt> adiacent_polygons){
	  adiacent_polygons_ = adiacent_polygons;
	  if(adiacent_polygons_.size() <  this->points_.size()){
	    this->closed_ = false;
	  }else{
	    this->closed_ = true;
	  }
	} // positive indexes point internal faces neighbor; Otherwise Surface neighborus
	void set_adiacent_cell(UInt adiacent_cell){adiacent_cell_= adiacent_cell;}
	double measure() const { // since the polygon is convex, I just need to compute the triangle with respect the first point and any successive couple of points. If the polygon is open, the measure is set to -1
	  if(this->closed()){
	    double result=0;
	    SVector<3> v1;
	    SVector<3> v2;
	    for(std::size_t i=1; i<this->points_.size()-1; i++){ //from (1,2) to (n-1, n) // CHECK THE FORMULA FOR TRIANGLE IN 3D space!!!!
	      v1 = this->points_[i]-this->points_[0];
	      v2 = this->points_[i+1]-this->points_[i];
	      result += (v1.cross(v2)).norm();
	    }
	    return result/2;
	  }else{
	    return -1; // NO measure is possible if polygon open
	  }
	}
	bool closed() const {return this->closed_;}
      
	const std::vector<SVector<3>> & points() const {return this->points_;}; // Orderded points of the poligon, with end/init not repated. The polygon is counterclockwise ordered, such that the cell is indicated by the normal vector generated by three ordered points
      
	const std::vector<UInt> & adiacent_polygons() const {return adiacent_polygons_;}; // Size = points_.rows(). Each element i represents the index of the Polygon adiacent to the segment points(i,...) - points(i+1,...). If the segment is adiacent to infinity, contains -1
      
	UInt adiacent_cell() const {
	  return adiacent_cell_;
	}

	HyperPlane<2,3> supporting_plane() const {
	UInt i=0;
	SVector<3> cross = (this->points_[i+1]-this->points_[i]).cross(this->points_[i+2]-this->points_[i+1]);
	while(i+2<this->points_.size() && cross.norm() < 0.0000000001){
	  i++;
	  cross = (this->points_[i+1]-this->points_[i]).cross(this->points_[i+2]-this->points_[i+1]);
	}
	  return HyperPlane<2,3>(this->points_[i], this->points_[i+1], this->points_[i+2]); // Build the hyperplane from any three nonlinear points of ConvexPolygon (at least 3 points always present)
	}
      
	bool check_Direction (SMatrix<3,2> segment) const {
	  double value = (segment.col(0)-this->centroid_).cross(segment.col(1)-segment.col(0)).dot((this->points_[0]-this->centroid_).cross(this->points_[1] - this->points_[0]));
	  return value > 0;
	}
      
	// This routine is of core importance. Implements the sutherland algorithm (adapted to the current polygon structure).
	// Input: - cut plane: the plane upon which we need to cut the Polygon of interesr
	//        - plane_index: index of the cutting plane (cutting polygon) inside the voronoi cell, used for adiacency operations (Yan aglorithm and clipping_by_surface routine)
	bool cut_by_plane(HyperPlane<2,3> cut_plane, UInt cut_plane_index){ //DMatrix<double> cut_points
        
	  // NOTE: THE OUTPUT POLYGON MAY BE EMPTY? WE NEED TO HANDLE WELL THIS CASE!!! ---> to be checked in superior routine
	  // Note: a point is inside if signed_distance >= 0, i.e. points on the border are inside. Moreover, we need to take special care if the intersection of segments with plane is one of the two points:   
	  bool result = false;
        
	  std::vector<SVector<3>> O_points;       // Output vertices
	  std::vector<UInt> O_adiacent_polygons; // Output connectivity vector
        
	  UInt n = this->points_.size();
        
	  // Prepare data structures
	  //O_points.resize(0);
	  //O_adiacent_polygons.resize(0);
        
	  // we need to take into account the maximum size
	  // UInt size=0;
        
	  double offset_0 = 0; // Offset is just the distance with sign with respect to a point
	  double offset_1 = 0; // Offset is just the distance with sign with respect to a point
        
	  double dist_0=0;
	  double dist_1=0;
        
	  double intersection=0;
        
	  double tol = 0.0000000001; // To avoid cancellation errors ---> Armonize with fdaPDE
        
	  for(auto i = 0; i<n; i++){
        
	    // Extract the current couple under analysis
	    SVector<3> P0 = this->points_[i];
	    SVector<3> P1 = P0;
          
	    // Take care of the last point
	    if(i==(n-1)){
	      P1 = this->points_[0];
	    }else{
	      P1 = this->points_[i+1];
	    }  
          
	    // distance of current point
	    offset_0 = (P0 - cut_plane.project(P0)).dot(cut_plane.normal());
	    dist_0 = std::abs(offset_0);
          
	    // distance of next point
	    offset_1 = (P1 - cut_plane.project(P1)).dot(cut_plane.normal());
	    dist_1 = std::abs(offset_1);
          
	    if((offset_0 < 0 && offset_1>0) || (offset_0>0 && offset_1 <0) || (dist_0 < tol && offset_1 < 0) || (dist_1 < tol && offset_0 < 0)){
	      result=true;
	    }
          
	    // Case both inside, current on plane and second inside, current inside and the other on plane
	    if((offset_0 > tol && offset_1 > tol) || (dist_0 < tol && offset_1 > tol) || (offset_0 > tol && dist_1 < tol)){
	      O_points.push_back(P0);
	      O_adiacent_polygons.push_back(adiacent_polygons_[i]);
	      //size++;
	    }
          
	    // Current in, next out
	    if(offset_0 > tol && offset_1 < -tol){
	      // Find intersection
	      intersection = dist_0 / (dist_0 + dist_1);
	      O_points.push_back(P0);
	      O_adiacent_polygons.push_back(adiacent_polygons_[i]); // i
	      //size++;
	      // Add the novel vertex
	      O_points.push_back(P0+intersection*(P1-P0)); // Add the novel intersection point 
	      O_adiacent_polygons.push_back(cut_plane_index); // i+1
	      //size++;
	    }
          
	    // Current out, next in
	    if(offset_0 < -tol && offset_1 > tol){
	      // Find intersection
	      intersection = dist_0 / (dist_0 + dist_1);
	      O_points.push_back(P0+intersection*(P1-P0)); // Add the novel intersection point
	      O_adiacent_polygons.push_back(adiacent_polygons_[i]); // i
	      //size++;
	    }
          
	    // Case current point is on the plane but the other is out or both on the plane
	    if((dist_0 < tol && offset_1 < -tol) || (dist_0<tol && dist_1<tol) ){ // Current on plane, next out of bound
	      // Find intersection
	      O_points.push_back(P0);
	      O_adiacent_polygons.push_back(cut_plane_index); // i
	      //size++;
	    }  
        
	    // In the remaining cases: do nothing!
	  }
        
	  // Update the vertices and the adiacency structure
	  this->points_ = O_points;
	  //O_adiacent_polygons.resize(size);
	  this->adiacent_polygons_ = O_adiacent_polygons;
        
	  if(this->points_.size()==this->adiacent_polygons_.size()){
	    this->closed_=true;
	  }
	  
	  return result; // There has been a cut, but we still need to handle the case of emplty lists!!!
        
	}

	void update_to_surface(std::list<SVector<3>> surface_limit, std::list<UInt> neighbors, const ConvexPolygon * poly_begin, const ConvexPolygon * poly_end, SVector<3> site){ //, const Voronoi* v_){
	  // To find where to place the list of segments, I just need to cycle and cut by the extremes of the
	  // the idea is to first compute where is the intersection between the last segment and the sequence (entering), and a the same time computing where is the last intersection between last segment.
	  // Then, we beed to keep just the inner points and the points of the line!!
	  // Finally er send back the so modified
	
	  HyperPlane<2,3> plane_begin = poly_begin->supporting_plane();
	  HyperPlane<2,3> plane_end = poly_end->supporting_plane();

	  //MOdify sutherland algorithm to find where to cut the elements
	  std::vector<SVector<3>> O_points; // Output vertices
	  std::vector<UInt> O_adiacent_polygons; // Output connectivity vector
	  
	  if(surface_limit.size() == neighbors.size()){ // We need to drop all of the points altoghether
	    auto iter_surf = surface_limit.begin();
	    auto iter_adiacent = neighbors.begin();
	  
	    for(UInt i =  0; (i < surface_limit.size()) && (i < neighbors.size()); i++ ){
	        O_points.push_back(*iter_surf);
	        O_adiacent_polygons.push_back(*iter_adiacent);
	        iter_surf++;
	        iter_adiacent++;
	     }
	     
	     // compute the normal to the new plane and check that it is consistent with notation (pointing to the site)
	     // Build the supporting plane
	     HyperPlane<2,3> current_plane(O_points[0], O_points[1], O_points[2]);
	     
	     SVector<3> mean = O_points[0];
	     // Reset centroid_ 
	     for(UInt l = 1; l < O_points.size(); l++){
	       mean = mean + O_points[l];
	     }
	     this->centroid_ = mean / O_points.size();
	     
	     if((current_plane.normal()).dot(site-current_plane.project(site)) < 0){ // The points we are substituting are not well ordered, reverse them // Here there may be a problem if the surface is locally really complex, because the good normal is still not guaranteed (few cells)
	       std::reverse(O_points.begin(), O_points.end());
	       std::reverse(O_adiacent_polygons.begin(), O_adiacent_polygons.end());
	     }
	     
	     this->points_ = O_points;
	     this->adiacent_polygons_ = O_adiacent_polygons;
	     
	     this->closed_=true;
	     
	     return;
	  }

	  UInt n = this->points_.size();
        
	  // we need to take into account the maximum size
	  UInt size=0;
        
	  double offset_0_begin=0;
	  double offset_1_begin=0;
	
	  double offset_0_end=0;
	  double offset_1_end=0;
        
	  double tol = 0.0000000001; // To avoid cancellation errors +

	  UInt index_begin = 0; // warning: may be used uninitialized
	  UInt index_end = 0;

	  for(auto i = 0; i<n; i++){
        
	    // Extract the current couple under analysis
	    SVector<3> P0 = this->points_[i];
	    SVector<3> P1 = P0;
          
	    // Take care of the last point
	    if(i==(n-1)){
	      P1 = this->points_[0];
	    }else{
	      P1 = this->points_[i+1];
	    }  
          
	    // distances from begin
	    offset_0_begin = (P0 - plane_begin.project(P0)).dot(plane_begin.normal());
	    offset_1_begin = (P1 - plane_begin.project(P1)).dot(plane_begin.normal());
          
	    // distances from end
	    offset_0_end = (P0 - plane_end.project(P0)).dot( plane_end.normal());
	    offset_1_end = (P1 - plane_end.project(P1)).dot(plane_end.normal());
          
	    // Current in or on border, next out from last segment of surface segments
	    if((offset_0_begin > tol) && (offset_1_begin < tol)){
	      index_end = i;
	    }
          
	    // Current out (or on border), next in from last segment of surface segments
	    if((offset_0_end < tol) && (offset_1_end > tol)){
	      if(i<n-1){
		index_begin = i+1; 
	      }else{
		index_begin = 0;
	      }
	    
	    }
        
	  }
	  
	  // Check that the cell limit is well ordered with respect to the ConvexPolygon --> we need to compute a pseudo centroid_ that is guaranteed to be inside the polygon
	  SVector<3> point_inside = this->points_[index_end]; // this point is guaranteed to be inside
	  
	  auto iter_limit = surface_limit.begin();
	  SVector<3> point_0, point_1;
	  point_0 = *iter_limit;
	  iter_limit++;
	  point_1 = *iter_limit;
	  
	  bool well_positioned = ((points_[1] - points_[0]).cross(points_[2]-points_[1])).dot((point_0-point_inside).cross(point_1-point_0)) > 0; //if greater than zero, then the original segment was ok and consequently all the cell limit is ok
	  
	  if(!well_positioned){ // Reverse everything, recompute index begin, index_end
	    std::reverse(surface_limit.begin(), surface_limit.end());
	    std::reverse(neighbors.begin(), neighbors.end());
	    
	    plane_begin = poly_end->supporting_plane(); // Switch the supporting planes
	    plane_end = poly_begin->supporting_plane(); // Switch the supporting planes
	    
	    // Recompute the index_begin and index_end, they CANNOT be inferred by the wrong ones!!!
	    for(auto i = 0; i<n; i++){
        
	    // Extract the current couple under analysis
	    SVector<3> P0 = this->points_[i];
	    SVector<3> P1 = P0;
          
	    // Take care of the last point
	    if(i==(n-1)){
	      P1 = this->points_[0];
	    }else{
	      P1 = this->points_[i+1];
	    }  
          
	    // distances from begin
	    offset_0_begin = (P0 - plane_begin.project(P0)).dot(plane_begin.normal());
	    offset_1_begin = (P1 - plane_begin.project(P1)).dot(plane_begin.normal());
          
	    // distances from end
	    offset_0_end = (P0 - plane_end.project(P0)).dot( plane_end.normal());
	    offset_1_end = (P1 - plane_end.project(P1)).dot(plane_end.normal());
          
	    // Current in or on border, next out from last segment of surface segments
	    if((offset_0_begin > tol) && (offset_1_begin < tol)){
	      index_end = i;
	    }
          
	    // Current out (or on border), next in from last segment of surface segments
	    if((offset_0_end < tol) && (offset_1_end > tol)){
	      if(i<n-1){
		index_begin = i+1; 
	      }else{
		index_begin = 0;
	      }
	    
	    }
        
	  }
	  }
	  
	  // Note: since the cell limit may have been ordered in a bad way, it may happen that we need to reverse the problem. If the ordering was right, the first point that we keep in should be inside the domain, not out. We check this by locate
	  //if(v_->locate_in_mesh(this->points_[index_begin])[0] == -1){ // The point is out, the cell limit was badly sorted
	  //  std::reverse(surface_limit.begin(), surface_limit.end());
	  //  std::reverse(neighbors.begin(), neighbors.end());
	  //  UInt aux_begin = index_begin;
	  //  UInt aux_end = index_end;
	  //  if(aux_begin == 0){
	  //    index_end = this->points_.size()-1;
	  //  }else{
	  //    index_end = aux_begin-1;
	  //  }
	    
	  //  if(aux_end == this->points_.size()-1){
	  //    index_begin = 0;
	  //  }else{
	  //    index_begin = aux_end+1;
	  //  }
	  //}
	
	  if(index_begin <= index_end){
	    //size = (index_end - index_begin + 1) + surface_limit.size();
	    //O_points.resize(size);
	    //O_adiacent_polygons.resize(size);
	    for(UInt i = index_begin; (i< index_end +1) && (i< this->points_.size()) && (i< this->adiacent_polygons_.size()); i++){ // Extra checks required by Valgrind
	      O_points.push_back(this->points_[i]);
	      O_adiacent_polygons.push_back(adiacent_polygons_[i]);
	    }
	  }else{
	    //size = (index_begin - index_end + 1) + surface_limit.size();
	    //O_points.resize(size);
	    //O_adiacent_polygons.resize(size);
	    for(UInt i = index_begin; (i< this->points_.size()) && (i< this->adiacent_polygons_.size()); i++){
	      O_points.push_back(this->points_[i]);
	      O_adiacent_polygons.push_back(adiacent_polygons_[i]);
	    }

	    for(UInt i = 0; (i< index_end+1) && (i< this->points_.size()) && (i< this->adiacent_polygons_.size()); i++){
	      O_points.push_back(this->points_[i]);
	      O_adiacent_polygons.push_back(adiacent_polygons_[i]);
	    }
	  }

	  auto iter_surf = surface_limit.begin();
	  auto iter_adiacent = neighbors.begin();
	  
	  for(UInt i =  0; (i < surface_limit.size()) && (i < neighbors.size()); i++ ){
	      O_points.push_back(*iter_surf);
	      O_adiacent_polygons.push_back(*iter_adiacent);
	      iter_surf++;
	      iter_adiacent++;
	    }
	    
	  // Add the last point, since neighbors has always a lower dimension with respect to surface_limit.size()
	  O_points.push_back(*iter_surf); // Iter surf is already at the last position
	    
	  // Add the neighbor related to the last surface segment, missing from surfaces (it is the one related to the segment entering in the volume, lost when inserting points)
	  if(index_begin > 0){
	    O_adiacent_polygons.push_back(adiacent_polygons_[index_begin-1]);
	  }else{
	    O_adiacent_polygons.push_back(adiacent_polygons_.back());
	  }

	  //if(index_begin < index_end){
	  //  for(UInt i =  index_end-index_begin+1; i < size; i++ ){
	  //    O_points[i] = *iter_surf;
	  //    O_adiacent_polygons[i] = *iter_adiacent;
	  //    iter_surf++;
	  //    iter_adiacent++;
	  //  }
	  //}else{
	  //  for(UInt i = this->points_.size() - index_begin+index_end+1; i < size; i++ ){
	  //    O_points[i] = *iter_surf;
	  //    O_adiacent_polygons[i] = *iter_adiacent;
	  //    iter_surf++;
	  //    iter_adiacent++;
	  //  }
	  //}

	
        
	  // Update the vertices and the adiacency structure
	  this->points_ = O_points;
	  this->adiacent_polygons_ = O_adiacent_polygons;
        
	  if(this->points_.size()==this->adiacent_polygons_.size()){
	    this->closed_=true;
	  }
	  
	  return;
	} 
  

      private:
	std::vector<SVector<3>> points_; // Orderded points of the poligon, with end/init not repated. The polygon is counterclockwise ordered, such that the cell is indicated by the normal vector generated by three ordered points
	std::vector<UInt> adiacent_polygons_; // Size = points.rows(). Each element i represents the index of the Polygon adiacent to the segment points(i,...) - points(i+1,...). If the segment is adiacent to infinity, contains -1
	bool closed_;
	UInt adiacent_cell_; // index of the adiacent cell in the voronoi or -1 if it is not adiacent to anything, after clipping to surface
	SVector<3> centroid_;
      };
      
      // Now I have a set of unsorted segments. I need to transform this into a ConvexPolygon, that is a list of points, ordered counter-clockwise with respect to the centroid
      // To do so, I define the following auxiliary routine
      ConvexPolygon build_convex_polygon(std::vector<segment> & original_segments, UInt adiacent_cell, SVector<3> site){ // Note: in this case site is right, because we are building just the convexPolygons of the inner cells, therefore, site is not on the Polygon!!!
	double tol = 0.000000001;
	
	// This Convex Polygon will contain the result
	ConvexPolygon result;
	
	std::vector<segment> segments;
	for(UInt i = 0; i < original_segments.size(); i++){
	  if((original_segments[i].p1 - original_segments[i].p2).norm()>tol){ // Drop the segments too small (it may happen for "Perfect square triangulations")
	    segments.push_back(original_segments[i]);
	  }
	}
	
	SVector<3> mean;
	mean[0]=0;
	mean[1]=0;
	mean[2]=0;
	UInt den = 2*segments.size();
	
	if(den < 4){ // Degenerate polygon (perfect shape triangulation)
	  std::vector<SVector<3>> dummy_points;
	  //dummy_points.resize(0);
	  result.set_points(dummy_points);
	  
	  return result; // The polygon will be dropped
	}else if(((segments[1].p1-segments[1].p2).cross(segments[0].p1-segments[0].p2)).norm()<tol){
	  // The Polygon is degenerate (at least the first two segments), we drop it
	  std::vector<SVector<3>> dummy_points;
	  //dummy_points.resize(0);
	  result.set_points(dummy_points);
	  
	  return result; // The polygon will be dropped
	}

	for(UInt i = 0; i< segments.size(); i++){
	  mean = mean + segments[i].p1 + segments[i].p2;
	}

	// Build the centroid
	mean = mean/den; 

	// Auxiliary structure used to build the sorted set of points. The key will be the angle with respect the centroid - first point of the polygon
	std::map<double,segment> final_polygon, auxiliary_map;

	SVector<3> v1 = segments[0].p1 - mean;
	SVector<3> v2 = segments[0].p2 - segments[0].p1;

	HyperPlane<2,3> plane(mean, segments[0].p1, segments[0].p2);
	SVector<3> normal = plane.normal();
	if((site - plane.project(site)).dot(normal) < 0){
	  normal = -normal;
	}

	if((v1.cross(v2)).dot(normal) > 0){
	  final_polygon[M_PI] = segments[0]; // Atan2 goes from -pi to pi
	}else{
	  segment seg  = segments[0];
	  seg.p2 = segments[0].p1;
	  seg.p1 = segments[0].p2;
	  final_polygon[M_PI] = seg;
	  v1 =  segments[0].p2 - mean;
	  v2 = -v2;
	}
	
	// normalize v1
	v1 = v1/v1.norm();
	  
	for(UInt i = 1; i<segments.size(); i++){
	  SVector<3> v3 = segments[i].p1 - mean;
	  SVector<3> v4 = segments[i].p2 - segments[i].p1;
	  
	  // normalize (v3)
	  v3= v3/v3.norm();

	  if((v3.cross(v4)).dot(normal) > 0){ // Segment is already well positioned
	    // Compute the angle in rad (counterclockwise) w.r.t. the original axis
	    double angle = 0;
	    //if(normal.dot(v1.cross(v3)) > 0){
	    angle = M_PI + atan2(normal.dot(v1.cross(v3)),v1.dot(v3)); // Chck with palummo!!
	    //}else{
	    //angle = atan2(-normal.dot(v1.cross(v3)),v1.dot(v3)); // Chck with palummo!!
	    //}
	    final_polygon[angle] = segments[i];
	  }else{
	    segment seg  = segments[i];
	    seg.p2 = segments[i].p1;
	    seg.p1 = segments[i].p2;
	    v3 =  segments[i].p2 - mean;
	    v4 = -v4;

	    //double angle = atan2(normal.dot(v3.cross(v1)),v3.dot(v1)); // Chck with palummo!!
	    double angle = 0;
	    //if(normal.dot(v1.cross(v3)) > 0){
	    angle = M_PI + atan2(normal.dot(v1.cross(v3)),v1.dot(v3)); // Chck with palummo!!
	    //}else{
	    //angle = atan2(-normal.dot(v1.cross(v3)),v1.dot(v3)); // Chck with palummo!!
	    //}
	    final_polygon[angle] = seg;
	  }
	}
	
	// Drop the double segments due to numerical errors
	double tolerance_angle = 0.00000001;
	auto iter_0 = final_polygon.begin();
	auto iter_1 = iter_0;
	iter_1++;
	std::vector<double> removal_set; // Set of keys to be removed due to duplication induced by numerical tolerances
	while(iter_1!=final_polygon.end()){
	  if(iter_1->first - iter_0->first < tolerance_angle){
	    removal_set.push_back(iter_1->first);
	  }
	  iter_0++;
	  iter_1++;
	}
	
	for(UInt k = 0; k < removal_set.size(); k++ ){
	  final_polygon.erase(removal_set[k]); // Delete the duplicates#31 0x00007ffff7b5419d in ?? () from /usr/lib/R/lib/libR.so

	}
	
	if(final_polygon.rbegin()->first - final_polygon.begin()->first > 2*M_PI - tol){ // Also remove the last segmetn if it is too close to the first!
	  final_polygon.erase(final_polygon.rbegin()->first);
	}
	
	// Sort the final polygon map to handle well the case of infinite segments without adiacency
	iter_0 = final_polygon.begin();
	iter_1 = iter_0;
	iter_1++; //at least two segments!
	bool found=false;
	while(iter_1!= final_polygon.end() && found == false){
	  if((iter_0->second.p2 - iter_1->second.p1).norm() > 0.000005){ // we found two infinite segments not matching, we need to sort again the map with respect to the iter_1
	    found=true;
	  }else{
	    iter_0++;
	    iter_1++;
	  }
	}

        if(found==true){
          double counter = 0;
          
          for(auto iter = iter_1; iter!= final_polygon.end(); iter++){
            auxiliary_map[counter] = iter->second;
            counter = counter +1;
          }

          for(auto iter = final_polygon.begin(); iter!= iter_1; iter++){
            auxiliary_map[counter] = iter->second;
            counter = counter +1;
          }
          
          final_polygon = auxiliary_map;
        }
        
        
	// Finally build the actual polygon
	std::vector<SVector<3>> points;
	std::vector<UInt> adiacent;

	//points.resize(0);
	//adiacent.resize(0);
	points.push_back(final_polygon.begin()->second.p1);
	for(auto iter = final_polygon.begin(); iter!= final_polygon.end(); iter++){
	  points.push_back(iter->second.p2);
	  adiacent.push_back(iter->second.neighbor);
	}

	// Closed or open
	if((points.front()-points.back()).norm()<0.00005){ // Check tol with Palummo!!!
	  points.pop_back();
	}

        // Fill the result
	result.set_points(points);
	result.set_adiacent_polygons(adiacent);
	result.set_adiacent_cell(adiacent_cell); 
	
	return result;
      };
      
      void build_unrestricted_voronoi(){ // This method will be invoked by the constructor to iniztialize the cells of the voronoi.
	// Idea: preallocate the space;
	//       then make a cycle across the faces of the triangulation. For each we have a segment (either infinite or finite)
	//       add the segment to the polygons associated to the edges
	//       finally sort-build the polygons (convex) by centroid;
	
	// This will contain the segments for each polygon for VoronoiCells
	std::vector<std::vector<std::vector<segment>>> segments; // This structure will store the segments of each polygon, one for every edge insisitng on a node. Final dim: n_nodes x n_edges_on_node x n_tets_insisiting on edge
	std::vector<std::map<UInt,UInt>> edges_to_poly; // Tis map is need to report the global voronoi index of an edge to the local index used to store the associated polygon in the cell, dim n_nodes
	std::vector<std::vector<UInt>> close_cell; // This structure contains the cell-cell adiacency info of a Polygon, dim n_nodes x n_edges_on_node
	std::vector<bool> analyzed_faces; // This structure keeps trace of the faces analyzed in the construction, so to save computational time, dim: n_faces in triangulation
	
	// Sizes and preallocation of space
	DVector<UInt> n_edges_to_cell; // Aux structure

	UInt n = mesh_->n_nodes();
	segments.resize(n);
	edges_to_poly.resize(n);
	close_cell.resize(n);
	n_edges_to_cell.resize(n);

	for(auto i =0; i<n;i++){
	  // at least one inside
	  n_edges_to_cell[i]=0;
	}

	// Fill the maps;
        for(auto edge = mesh_->edges_begin(); edge != mesh_->edges_end(); ++edge){ // NB each point is counted 
	  auto P0 = edge->node_ids()(0); //
	  auto P1 = edge->node_ids()(1); //
	  UInt edge_id = edge->id(); // This exists;

	  edges_to_poly[P0].insert({edge_id,n_edges_to_cell[P0]});
	  n_edges_to_cell[P0]++;

	  edges_to_poly[P1].insert({edge_id,n_edges_to_cell[P1]});
	  n_edges_to_cell[P1]++;
	}

	for(UInt i=0;i<n;i++){
	  segments[i].resize(n_edges_to_cell[i]);
	  close_cell[i].resize(n_edges_to_cell[i]);
	  //for (UInt j = 0; j < n_edges_to_cell[i]; j++){
	  //  segments[i][j].resize(0);
	  //}
	}
	
	analyzed_faces.resize(mesh_->n_faces());
	for(std::size_t i = 0; i < analyzed_faces.size(); i++){
	  analyzed_faces[i]=false;
	}


	// Main cycle: we cycle over the faces, and add the appropriate segments to the nodes.
	// Since we need to employ the information of the tetrahedron used, we need make a cycle based on those.

	// To be modified: this M is too large, we need something really smaller (because in pratice it will be needed to be half the cell: it could be the circumradius!!!!)
	double M = 0;
	//((nodes_.col(0)).maxCoeff() - (nodes_.col(0)).minCoeff())*((nodes_.col(0)).maxCoeff() - (nodes_.col(0)).minCoeff()) + 
	//(((nodes_.col(1)).maxCoeff() - (nodes_.col(1)).minCoeff())*((nodes_.col(1)).maxCoeff() - (nodes_.col(1)).minCoeff()) +
	//((nodes_.col(2)).maxCoeff() - (nodes_.col(2)).minCoeff())*((nodes_.col(2)).maxCoeff() - (nodes_.col(2)).minCoeff()) ; // max distance in mesh

	for(auto tet = mesh_->cells_begin(); tet != mesh_->cells_end(); ++tet ){ // Note: we work on the faces of each tetrahedron
	  
	  M = tet->circumradius(); // NB part of the inner cell may not touch the border. However, in general works well becaus we only use it as plane and then clip to surface puts everything in order
	  
	  SVector<3> tet_mean;
	  tet_mean[0]=0;
	  tet_mean[1]=0;
	  tet_mean[2]=0;
	  for(UInt i=0; i<4;i++){
	    UInt index = tet->node_ids()[i];
	    SVector<3> node = nodes_.row(index);
	    tet_mean = tet_mean + node;
	  }
	  tet_mean = tet_mean/4;
	  
	  for(auto face = tet->faces_begin(); face != tet->faces_end(); ++face ){ // for each face, do something only if the face hasn't been analyzed yet
	    
	    if(analyzed_faces[face->id()] == false){ // Only if the face hasn't been analyzed yet
	      // analyzed_faces[face->id()] == true; // warning: unused-value
	      segment seg;

	      if(face->on_boundary()){
		auto face_plane = face->supporting_plane(); 

                // Check the normal is outgoing or not
                SVector<3> normal = face_plane.normal();
                if((tet_mean-face_plane.project(tet_mean)).dot(normal) > 0){ // Normal is ingoing,I want outgoing 
                  normal = -normal;
                }

		seg.p1 = tet->circumcenter();
		seg.p2 = tet->circumcenter() + normal*M; ////// Correct when normal is outgoing
	      }else{
	        auto neighbors = face->adjacent_cells();
	        UInt neighbor = neighbors[1];
	        if(neighbor == tet->id()){
	          neighbor = neighbors[0];
	        }
		auto tet_1 = mesh_->cells_begin() + neighbor; // get the other tet

		seg.p1 = tet->circumcenter();
		seg.p2 = tet_1->circumcenter();
	      }

	      // Now, save the segments in each position
	      auto e0 = face->edge(0);
	      auto e1 = face->edge(1);
	      auto e2 = face->edge(2);

	      auto p0 = e0.node_ids()[0];
	      auto p1 = e0.node_ids()[1];
	      auto p2 = e1.node_ids()[0];
	      
	      if(p2==p1 || p2==p0){
	        p2 = e1.node_ids()[1];
	      }

	      segment seg_1, seg_2, seg_3, seg_4, seg_5, seg_6;
	      
	      if(p0 == e1.node_ids()[0] || p0 == e1.node_ids()[1]){ // p0 opposite to e2
		seg_1 = seg;
		seg_1.neighbor = edges_to_poly[p0][e1.id()];
		seg_2 = seg;
		seg_2.neighbor = edges_to_poly[p0][e0.id()];
		seg_3 = seg;
		seg_3.neighbor = edges_to_poly[p1][e2.id()];
		seg_4 = seg;
		seg_4.neighbor = edges_to_poly[p1][e0.id()];
		seg_5 = seg;
		seg_5.neighbor = edges_to_poly[p2][e2.id()];
		seg_6 = seg;
		seg_6.neighbor = edges_to_poly[p2][e1.id()];
	      
		segments[p0][edges_to_poly[p0][e0.id()]].push_back(seg_1);

		segments[p0][edges_to_poly[p0][e1.id()]].push_back(seg_2);
	      
		close_cell[p0][edges_to_poly[p0][e0.id()]] = p1;
		close_cell[p0][edges_to_poly[p0][e1.id()]] = p2;
	      
		segments[p1][edges_to_poly[p1][e0.id()]].push_back(seg_3);
	      
		segments[p1][edges_to_poly[p1][e2.id()]].push_back(seg_4);
	      
		close_cell[p1][edges_to_poly[p1][e0.id()]] = p0;
		close_cell[p1][edges_to_poly[p1][e2.id()]] = p2;
	      
		segments[p2][edges_to_poly[p2][e1.id()]].push_back(seg_5);
	      
		segments[p2][edges_to_poly[p2][e2.id()]].push_back(seg_6);
	      
		close_cell[p2][edges_to_poly[p2][e2.id()]] = p0;
		close_cell[p2][edges_to_poly[p2][e1.id()]] = p1;
              }else{ // p0 opposite to e1
		seg_1 = seg;
		seg_1.neighbor = edges_to_poly[p0][e2.id()];
		seg_2 = seg;
		seg_2.neighbor = edges_to_poly[p0][e0.id()];
		seg_3 = seg;
		seg_3.neighbor = edges_to_poly[p1][e1.id()];
		seg_4 = seg;
		seg_4.neighbor = edges_to_poly[p1][e0.id()];
		seg_5 = seg;
		seg_5.neighbor = edges_to_poly[p2][e2.id()];
		seg_6 = seg;
		seg_6.neighbor = edges_to_poly[p2][e1.id()];
	      
		segments[p0][edges_to_poly[p0][e0.id()]].push_back(seg_1);

		segments[p0][edges_to_poly[p0][e2.id()]].push_back(seg_2);
	      
		close_cell[p0][edges_to_poly[p0][e0.id()]] = p1;
		close_cell[p0][edges_to_poly[p0][e2.id()]] = p2;
	      
		segments[p1][edges_to_poly[p1][e0.id()]].push_back(seg_3);
	      
		segments[p1][edges_to_poly[p1][e1.id()]].push_back(seg_4);
	      
		close_cell[p1][edges_to_poly[p1][e0.id()]] = p0;
		close_cell[p1][edges_to_poly[p1][e1.id()]] = p2;
	      
		segments[p2][edges_to_poly[p2][e1.id()]].push_back(seg_5);
	      
		segments[p2][edges_to_poly[p2][e2.id()]].push_back(seg_6);
	      
		close_cell[p2][edges_to_poly[p2][e2.id()]] = p0;
		close_cell[p2][edges_to_poly[p2][e1.id()]] = p1;
              }
	    }
	  }
	}


	
	// Final step requires to build the voronoi cells
	std::vector<CellType> cells;
	cells.resize(n);
	for(auto i = 0; i < n; i++){
	  std::vector<ConvexPolygon> Cell_i;
	  // Cell_i.resize(segments[i].size());
	  for(std::size_t j = 0; j < segments[i].size();j++){ // USE PUSH BACK WITH IF HERE
	    if(segments[i][j].size()>0){
	      ConvexPolygon polygon = build_convex_polygon(segments[i][j],close_cell[i][j], nodes_.row(i)); // Check with Palummo nodes_[i,] for the site, needed for ordering of the normal
	      if(polygon.points().size()!=0){ // The polygon is not degenerate
	        Cell_i.push_back(polygon);
	      }
	    }
	  }
	  cells[i].set_inner_faces(Cell_i);
	  cells[i].set_id(i);
	}

        // Add the cells to internal memory storage
	for(auto i = 0; i< n; i++){
	  this->cells_.insert({i, cells[i]}); 
	}
      };

      // This routine computes the surface polygons associated to each cell. The algorithm mimics the one of Yan, Wang, Levy, Liu (2013).
      // The core idea is to employ a FIFO queue:
      // We start from incident cell-triangle couples. We analyze it, and we can propose a successive incident cell-triangle couple in the queue
      // We employ a boolean matrix to avoid to analyze a couple twice
      void cut_surface_to_cells(){
	// The matrix is needed to avoid double analysis of the same cell-surface_tringle couple
	DMatrix<int> incident_set;
	
	// Note: we can extract the info of adiacency by using the surface method of the mesh.
	// However, we have the problem of treating the map between this surface and the triangles that we use in the normal mesh. However, is this really a problem? Actually no, if we don't need to take care of the underlying tetrahedra ---> but we need to do that!!
	// Discuss with PALUMMO!!!
	auto surface_struct_ = mesh_->surface();
	auto surface_ = surface_struct_.triangulation;
	auto surface_node_map_ = surface_struct_.node_map;
	auto surface_face_map_ = surface_struct_.cell_map;
	auto surface_iterator = surface_.cells_begin(); // CHECK Syntax with Palummo
	
	UInt n_faces = surface_.n_cells();
	UInt n_nodes = mesh_->n_nodes();

	incident_set.resize(n_nodes, n_faces);
	for(UInt i=0; i<n_nodes; i++){
	  for(UInt j=0; j<n_faces; j++){
	    incident_set(i,j)=0;
	  }
	}

	// Create and initialize the FIFO queue.
	struct incident_cell_face_cuples{ // NB this structure contains just the indexes of the 
	  UInt cell;
	  UInt face;
	};

	std::list<incident_cell_face_cuples> FIFO_queue;
	incident_cell_face_cuples current_couple;
	current_couple.cell = surface_node_map_.at(surface_iterator->node_ids()[0]); // NBBBBBBB THIS IS THE NODE IDS OF THE SURFACE; IT IS NOT SAID THAT IS OK ALSO FOR THE GLOBAL MESH
	current_couple.face = surface_iterator->id(); // Syntax ceck with Palummo
	FIFO_queue.push_back(current_couple);
	incident_set(current_couple.cell, current_couple.face)=1;

	// Prepare the structures that will be used to build actually the surfaces, eventually substituted by apposite methods in voronoi cell (push_back in voronoi cell? Little space).
	std::vector<std::vector<ConvexPolygon>> Surfaces_polygons_collection;
	Surfaces_polygons_collection.resize(mesh_->n_nodes());
	
	// Main construction cycle
	while(FIFO_queue.size() != 0){
	  incident_cell_face_cuples actual = FIFO_queue.front();
	  FIFO_queue.pop_front();
	  auto actual_face_iterator = surface_.cells_begin() + actual.face; // Possible problems with iterators plus numbers?
	  
	  // Create the tethrahedron mean
	  SVector<3> tet_mean;
	  tet_mean(0)=0;
	  tet_mean(1)=0;
	  tet_mean(2)=0;
	  auto close_tet = mesh_->cells_begin() + surface_face_map_.at(actual_face_iterator->id()); // Surface_face_map returns the index of the neighboring tet i nthe original triangulation
	  for(auto i=0; i<4;i++){
	    SVector<3> node = nodes_.row(close_tet->node_ids()[i]);
	    tet_mean = tet_mean + node;
	  }
	  tet_mean = tet_mean/4;
	  
	  ConvexPolygon triangle_surface(actual_face_iterator, this->nodes_, surface_node_map_, tet_mean); // NBBB NOT SAID THAT NODE IDS ARE OK WITH OUR NODES_ !!!!

	  // Cut against all the planes of the Voronoi cell (use specific routine inside voronoi cell?)
	  const std::vector<ConvexPolygon> & inner_faces = this->cells_[actual.cell].get_inner_faces();
	  for(std::size_t i = 0; i<inner_faces.size(); i++){
	    bool cut = triangle_surface.cut_by_plane(inner_faces[i].supporting_plane(),i);
	    if(cut){
	      if(incident_set(inner_faces[i].adiacent_cell(), actual.face)==0){
		incident_set(inner_faces[i].adiacent_cell(), actual.face) = 1;
		current_couple.face = actual.face;
		current_couple.cell = inner_faces[i].adiacent_cell();
		FIFO_queue.push_back(current_couple);
	      }
	    }
	  }

	  // Finally, if we have still some adi acency with the surface polygons->add the couples actual cell and the close surface
	  std::vector<UInt> ad_poly = triangle_surface.adiacent_polygons();
	  for(std::size_t j = 0; j < ad_poly.size(); j++){
	    if(ad_poly[j]<0){
	      if(incident_set(actual.cell, -1-ad_poly[j])==0){
		incident_set(actual.cell, -1-ad_poly[j]) = 1;
		current_couple.face = -1-ad_poly[j];
		current_couple.cell = actual.cell;
		FIFO_queue.push_back(current_couple);
	      }
	    }
	  }

	  // Finally add the cut surface to the apporpriate zone
	  if(triangle_surface.points().size()>2){ // Only if the polygon is not degenerate!!!
	    Surfaces_polygons_collection[actual.cell].push_back(triangle_surface);
	  }
	}
	
	// To conclude I need to add the set_surface_faces to each vornoi cell
	for(auto i = 0; i < this->nodes_.rows(); i++){
	  this->cells_[i].set_surface_faces(Surfaces_polygons_collection[i], this->nodes_.row(i));
	}
      }
	      
      // cell data structure
      class VoronoiCell {
      private:
	const Voronoi* v_;
	UInt id_ = 0;
	std::vector<ConvexPolygon> InnerFaces_; // Set when unconstrained Voronoi is built   // New
	std::vector<ConvexPolygon> SurfaceFaces_; // Set with Yan's algorithm
		
	void clip_to_surface(SVector<3> site){ // This methods clips the Inner Polygons (that may go out of the surface) to surface limits. This implementation relies on the adiacent polygons of the surface polygons.
	  // STEP 1: find segments of interest
	  std::vector<std::list<SMatrix<3,2>>> cut_segments; // This vector contains the segments to cut each Polygon, obtained by the surfaces cells;
	  std::vector<std::list<UInt>> cut_indexes; // This vector contains the indexes of the cells corresponding to the cut we are prforming on the surface (to provide the final adiacency)
		  
	  cut_segments.resize(this->InnerFaces_.size());
	  cut_indexes.resize(this->InnerFaces_.size());	   
			
	  // Now cycle over surface and look for segments we need to cut the inner faces
	  for(UInt j = 0; j < this->SurfaceFaces_.size();j++ ){
	    std::vector<SVector<3>> points = SurfaceFaces_[j].points();
	    std::vector<UInt> adiacent_polygons = SurfaceFaces_[j].adiacent_polygons();
	    for(UInt k = 0; k < points.size(); k++){
	      UInt adiacent = adiacent_polygons[k];
	      if(adiacent>=0){ // NOTE: if < 0 be careful, need to add 1 before abs!!!
		SMatrix<3,2> segment;
		segment.col(0) = points[k];
		if(k==points.size()-1){
		  segment.col(1) = points[0];
		}else{
		  segment.col(1) = points[k+1];
		}
		cut_segments[adiacent].push_back(segment);
		cut_indexes[adiacent].push_back(-j -1); // NOTE: if < 0 be careful, need to add 1 before abs!!!
	      }
	    }
	  }
		  
	  // Now we need to sort the segments for each cell cut by the border!!
	  // Reamark: this implementation can becaome very inefficient for a large number of Internal faces - sruface faces couple, due to push_front. However,in practice  I expect 2-3 intsertions at most
	  for(UInt i = 0; i < this->InnerFaces_.size(); i++){ // Each possible internal face
	    UInt final_size = cut_segments[i].size() + 1;
	    if(cut_segments[i].size()  > 0){// at least one segment present
		      
	      // Tolerance for equivalence // Check with Palummo
	      double tol = 0.0000001;
		
	      UInt number_of_borders = 0;
	      bool finished_local_border = false;
	      
	      while(cut_segments[i].size()!=0 && cut_indexes[i].size()!=0 ){	// It may happen that the same polygon may be subject to two distinct borders, related t two different surfaces
	      
	      number_of_borders++; // Increment the number of borders analyzed      
	       
	      // Create a vector of sorted points:
	      std::list<SVector<3>> cell_limit;
	      std::list<UInt> limit_neighbors;
	      
	      // Fill with the first segment
	      SMatrix<3,2> segment = (cut_segments[i]).front();
	      (cut_segments[i]).pop_front();
	      
	      UInt adiacent = (cut_indexes[i]).front();
	      (cut_indexes[i]).pop_front();
		      
	      // Check it is counterclockwise ordered
	      bool counterclock = InnerFaces_[i].check_Direction(segment); // Check direction missing, we need to implement this!!
		      
	      cell_limit.push_front(segment.col(0));
	      if(counterclock){
		cell_limit.push_back(segment.col(1));
	      }else{
		cell_limit.push_front(segment.col(1));
	      }
	      
	      limit_neighbors.push_back(adiacent);
	      
	      finished_local_border = false;

	      while(cut_segments[i].size()!=0 && cut_indexes[i].size()!=0 && !finished_local_border){

		auto iterator_segments = (cut_segments[i]).begin();
		auto iterator_adiacent = (cut_indexes[i]).begin();
		bool found = false;
		while(iterator_segments != cut_segments[i].end() && found == false){
		    
		  if((iterator_segments->col(0) - cell_limit.front()).norm() < tol){
		    cell_limit.push_front(iterator_segments->col(1)); // add the other point
		    limit_neighbors.push_front(*iterator_adiacent);
		    cut_segments[i].erase(iterator_segments);
		    cut_indexes[i].erase(iterator_adiacent);
		    found=true;
		  }else if((iterator_segments->col(1) - cell_limit.front()).norm() < tol){
		    cell_limit.push_front(iterator_segments->col(0)); // add the other point
		    limit_neighbors.push_front(*iterator_adiacent);
		    cut_segments[i].erase(iterator_segments);
		    cut_indexes[i].erase(iterator_adiacent);
		    found=true;
		  }else if((iterator_segments->col(0) - cell_limit.back()).norm() < tol){
		    cell_limit.push_back(iterator_segments->col(1)); // add the other point
		    limit_neighbors.push_back(*iterator_adiacent);
		    cut_segments[i].erase(iterator_segments);
		    cut_indexes[i].erase(iterator_adiacent);
		    found=true;
		  }else if((iterator_segments->col(1) - cell_limit.back()).norm() < tol){
		    cell_limit.push_back(iterator_segments->col(0)); // add the other point
		    limit_neighbors.push_back(*iterator_adiacent);
		    cut_segments[i].erase(iterator_segments);
		    cut_indexes[i].erase(iterator_adiacent);
		    found=true;
		  }
		    
		  if(!found){
		    iterator_segments++;
		    iterator_adiacent++;
		  }
		}

              if(found ==false ){ // We didn't find anything, so the actual border is finished. We need to close the cycle or analyze another border
                finished_local_border = true;
              }

	      }
	      
	      if((cell_limit.front() - cell_limit.back()).norm()<tol){
	        // remove the last point and the last adiacent polygon due to repetition;
	        cell_limit.pop_back(); // Note: now cell limit and limit neighbors have the same size which needs to be handled by the update to surface method
	      }

	      // Find the planes associated with extreme segments
	      const ConvexPolygon * poly_begin = & SurfaceFaces_[-1-limit_neighbors.front()]; // 1 - front to transofrm from negative to positive
	      const ConvexPolygon * poly_end = & SurfaceFaces_[-1-limit_neighbors.back()];
	      
	      //SVector<3> site = v_->site(this->id_);
		
	      // Note: the line of segmetns is already orderd in the right direction. Now we need to cut the faces and finally to throw away the faces that are still unlimited
	      InnerFaces_[i].update_to_surface(cell_limit, limit_neighbors, poly_begin, poly_end, site);//, v_);
	    }
	    }
	    // MISSING (and probably will not be implemented)
	    // InnerFaces_[i].erase_out_of_bound(); // Note: we need to erase the polygons that are still infinite (i.e. case with two or one circumcenter out of surface and ramaining missing)
	  }
            
	}
	  
      public:
	VoronoiCell() = default;
	VoronoiCell(UInt id, const Voronoi* v) : v_(v), id_(id), InnerFaces_(v_->cells_.at(id_).get_inner_faces()), SurfaceFaces_(v_->cells_.at(id_).get_surface_faces()){} 
	void set_id(UInt id){id_=id;}
	void set_inner_faces(const std::vector<ConvexPolygon> & InnerFaces){InnerFaces_ = InnerFaces;}
	void set_surface_faces(const std::vector<ConvexPolygon> & SurfaceFaces, SVector<3> site){
	  SurfaceFaces_ = SurfaceFaces;
	  this->clip_to_surface(site); // Clip the inner faces against surface

	  std::vector<UInt> indices_out; // There may exist a Polygon of the cell that needs to be cut out
	  //indices_out.resize(0);

	  for(std::size_t i = 0; i<InnerFaces_.size(); i++){
	    if(!InnerFaces_[i].closed()){
	      indices_out.push_back(i);
	    }
	  }

	  if(indices_out.size()>0){
	    for(int i = int(indices_out.size()-1); i>-1 ;i--){ // mmhhh... this loop seems strange
	      InnerFaces_.erase(InnerFaces_.begin()+indices_out[i]); // Throw out the extreme faces (very rare situation, in principle)
	    }
	  }
	  return;
	}
	
	// Getters
	const std::vector<ConvexPolygon> & get_inner_faces() const {return this->InnerFaces_;}
	const std::vector<ConvexPolygon> & get_surface_faces() const {return this->SurfaceFaces_;}
	  
	double measure() const { // This part is still critical, we will think in futura how to solve it
	  double volume = 0;
	  SVector<3> site = v_->nodes_.row(this->id_);
	  for (std::size_t i = 0; i < InnerFaces_.size(); ++i) {
	    double area = (InnerFaces_[i]).measure();
	    //SVector<3> site = v_->nodes_[i]; // Check
	    HyperPlane<2,3> supp_plane = ((InnerFaces_[i]).supporting_plane());
	    double sign_dist = (site - supp_plane.project(site)).dot(supp_plane.normal());
	    volume = volume + area * sign_dist;
	  }
            
	  for (std::size_t j = 0; j < SurfaceFaces_.size(); ++j) {
	    double area = (SurfaceFaces_[j]).measure();
	    //SVector<3> site = v_->nodes_[j]; // Check
	    HyperPlane<2,3> supp_plane = ((SurfaceFaces_[j]).supporting_plane());
	    double sign_dist = (site - supp_plane.project(site)).dot(supp_plane.normal());
	    volume = volume + area * sign_dist ; // If negative, correctly subtracts to the total volume (concavity)
	  }
            
	  return  volume/3;
	}
	ConvexPolygon InnerFace(std::size_t i) const {
	  fdapde_assert(i < InnerFaces_.size());
	  return InnerFaces_[i];
	}
	ConvexPolygon SurfaceFace(std::size_t i) const {
	  fdapde_assert(i < SurfaceFaces_.size());
	  return SurfaceFaces_[i];
	}
	bool on_boundary() const {
	  if(SurfaceFaces_.size()>0){
	    return true;
	  }
	  return false;   // no edge on boundary
	}
	bool contains(const SVector<3>& p) const { return v_->locate(p)[0] == id_; }
      };
      
      // getters
      const DMatrix<double>& sites() const { return mesh_->nodes(); }
      SVector<embed_dim> site(UInt id) const { return mesh_->node(id); }
      const BinaryVector<fdapde::Dynamic>& boundary_vertices() const { return nodes_markers_; } // NOTE: set in the constructor, when clipped to surface called, still missing
      const Triangulation<3, 3>& dual() const { return *mesh_; }
      UInt n_nodes() const { return nodes_.rows(); }
      UInt n_cells() const { return mesh_->n_nodes(); }
      using CellType = VoronoiCell;
      CellType cell(UInt id) const { return VoronoiCell(id, this); } // Needs constructor from the Voronoi pointer!!! Look below
      // iterators
      class cell_iterator : public index_based_iterator<cell_iterator, CellType> {
	using Base = index_based_iterator<cell_iterator, CellType>;
	using Base::index_;
	friend Base;
	const Voronoi* voronoi_;
	cell_iterator& operator()(UInt i) {
	  Base::val_ = voronoi_->cell(i);
	  return *this;
	}
      public:
	cell_iterator(UInt index, const Voronoi* voronoi) : Base(index, 0, voronoi->n_cells()), voronoi_(voronoi) {
	  if (index_ < voronoi_->n_cells()) this->val_ = voronoi_->cell(index_);
	}
      };
      cell_iterator cells_begin() const { return cell_iterator(0, this); }
      cell_iterator cells_end() const { return cell_iterator(n_cells(), this); }
      // perform point location for set of points p_1, p_2, \ldots, p_n
      DVector<UInt> locate(const DMatrix<double>& locs) const {
	fdapde_assert(locs.cols() == embed_dim);
	// find delanuay cells containing locs
	DVector<UInt> dual_locs = mesh_->locate(locs);
	for (UInt i = 0; i < locs.rows(); ++i) {
	  if (dual_locs[i] == -1) continue;   // location outside domain
	  // find nearest cell to i-th location
	  typename Triangulation<3, 3>::CellType f = mesh_->cell(dual_locs[i]);
	  SMatrix<1, Triangulation<3, 3>::n_nodes_per_cell> dist =
	    (f.nodes().colwise() - locs.row(i).transpose()).colwise().squaredNorm();
	  UInt min_index;
	  double min = dist.minCoeff(&min_index);
	  DVector<int> neighbors = f.neighbors();
	  for(UInt j = 0; j<neighbors.size(); ++j){ // First level neighbors
	    if(neighbors[j] >= 0){
	      typename Triangulation<3, 3>::CellType f_1 = mesh_->cell(neighbors[j]);
	      SMatrix<1, Triangulation<3, 3>::n_nodes_per_cell> dist_1 =
	        (f_1.nodes().colwise() - locs.row(i).transpose()).colwise().squaredNorm();
	      UInt min_index_1;
	      double min_1 = dist_1.minCoeff(&min_index_1);
	      if(min_1<min){
	        min_index = min_index_1;
	      }
	    }
	  }
	  dual_locs[i] = f.node_ids()[min_index];
	}
	return dual_locs;
      }  
    private:
      const Triangulation<3, 3>* mesh_;
      DMatrix<double> nodes_;                             // voronoi vertices
      BinaryVector<fdapde::Dynamic> nodes_markers_;       // i-th element true if i-th vertex is on boundary
      std::unordered_map<UInt, CellType> cells_;   // for each cell id, the associated voronoiCell
    };


template <> class Voronoi<Triangulation<2, 2>> {
   public:
    static constexpr int local_dim = Triangulation<2,2>::local_dim;
    static constexpr int embed_dim = Triangulation<2,2>::embed_dim;

    Voronoi() = default;
    Voronoi(const Triangulation<2, 2>& mesh) : mesh_(&mesh) {   // constructs voronoi diagram from Delanoy triangulation
        int n_delaunay_faces = mesh_->n_cells();
        int n_delaunay_boundary_edges = mesh_->n_boundary_edges();
        nodes_.resize(n_delaunay_faces + n_delaunay_boundary_edges + mesh_->n_boundary_nodes(), embed_dim);
        nodes_markers_.resize(nodes_.rows());
        int k = n_delaunay_faces;
        for (typename Triangulation<2, 2>::cell_iterator it = mesh_->cells_begin(); it != mesh_->cells_end(); ++it) {
            nodes_.row(it->id()) = it->circumcenter();
            for (int v : it->node_ids()) { cells_[v].push_back(it->id()); }
            if (it->on_boundary()) {
                for (typename Triangulation<2, 2>::CellType::edge_iterator jt = it->edges_begin();
                     jt != it->edges_end(); ++jt) {
                    if (jt->on_boundary()) {
                        nodes_.row(k) = jt->supporting_plane().project(nodes_.row(it->id()));
                        nodes_markers_.set(k);
                        for (int v : jt->node_ids()) { cells_[v].push_back(k); }
                        k++;
                    }
                }
            }
        }
        // augment node set with boundary vertices, sort each cell clockwise (around its mean point)
        for (auto& [key, value] : cells_) {
            if (mesh_->is_node_on_boundary(key)) {
                nodes_.row(k) = mesh_->node(key);
                nodes_markers_.set(k);
                value.push_back(k);
                k++;
            }
            SVector<embed_dim> mean = SVector<embed_dim>::Zero();
            auto compare = clockwise_order<SVector<embed_dim>>(
              std::accumulate(
                value.begin(), value.end(), mean, [&](const auto& c, int a) { return c + nodes_.row(a).transpose(); }) /
              value.size());
            std::sort(value.begin(), value.end(), [&](int i, int j) { return compare(nodes_.row(i), nodes_.row(j)); });
        }
    }

    // cell data structure
    class VoronoiCell {
       private:
        const Voronoi* v_;
        int id_ = 0;
        int n_edges_ = 0;
       public:
        VoronoiCell() = default;
        VoronoiCell(int id, const Voronoi* v) : v_(v), id_(id), n_edges_(v_->cells_.at(id_).size()) { }
        // matrix of edge identifiers
        DMatrix<int> edges() const {
            DMatrix<int> result;
            result.resize(n_edges_, local_dim);
            for (int j = 0; j < n_edges_; ++j) {
                for (int k = 0; k < local_dim; ++k) result(j, k) = v_->cells_.at(id_)[(j + k) % n_edges_];
            }
            return result;
        }
        double measure() const {
            double area = 0;
            for (int j = 0; j < n_edges_; ++j) {
                // compute doubled area of triangle connecting the j-th edge and the center (use cross product)
                SVector<embed_dim> x = v_->vertex(v_->cells_.at(id_)[j]);
                SVector<embed_dim> y = v_->vertex(v_->cells_.at(id_)[(j + 1) % n_edges_]);
                area += x[0] * y[1] - x[1] * y[0];
            }
            return 0.5 * std::abs(area);
        }
        Simplex<1, 2> edge(int i) const {
            fdapde_assert(i < n_edges_);
            SMatrix<embed_dim, 2> coords;
            for (int k = 0; k < embed_dim; ++k) { coords.col(k) = v_->vertex(v_->cells_.at(id_)[(i + k) % n_edges_]); }
            return Simplex<1, 2>(coords);
        }
        bool on_boundary() const {
            for (int j = 0; j < n_edges_; ++j) {
                bool boundary = true;
                for (int k = 0; k < local_dim; ++k)
                    boundary &= v_->nodes_markers_[v_->cells_.at(id_)[(j + k) % n_edges_]];
                if (boundary == true) return true;
            }
            return false;   // no edge on boundary
        }
        bool contains(const SVector<embed_dim>& p) const { return v_->locate(p)[0] == id_; }
    };
    // getters
    const DMatrix<double>& sites() const { return mesh_->nodes(); }
    SVector<embed_dim> vertex(int id) const { return nodes_.row(id); }
    SVector<embed_dim> site(int id) const { return mesh_->node(id); }
    const BinaryVector<fdapde::Dynamic>& boundary_vertices() const { return nodes_markers_; }
    const DMatrix<double>& vertices() const { return nodes_; }
    const Triangulation<2, 2>& dual() const { return *mesh_; }
    int n_nodes() const { return nodes_.rows(); }
    int n_cells() const { return mesh_->n_nodes(); }
    // compute matrix of edges
    DMatrix<int> edges() const {
        std::unordered_set<std::array<int, local_dim>, std_array_hash<int, local_dim>> visited;
        std::array<int, local_dim> edge;
        for (const auto& [key, value] : cells_) {
            int n_edges = value.size();
            for (int j = 0; j < n_edges; ++j) {
                for (int k = 0; k < local_dim; ++k) { edge[k] = value[(j + k) % n_edges]; }
                std::sort(edge.begin(), edge.end());
                if (visited.find(edge) == visited.end()) { visited.insert(edge); }
            }
        }
        DMatrix<int> result;
        result.resize(visited.size(), local_dim);
        int i = 0;
        for (const auto& e : visited) {
            for (int k = 0; k < local_dim; ++k) result(i, k) = e[k];
            i++;
        }
        return result;
    }
    using CellType = VoronoiCell;
    CellType cell(int id) const { return VoronoiCell(id, this); }
    // iterators
    class cell_iterator : public index_based_iterator<cell_iterator, CellType> {
        using Base = index_based_iterator<cell_iterator, CellType>;
        using Base::index_;
        friend Base;
        const Voronoi* voronoi_;
        cell_iterator& operator()(int i) {
            Base::val_ = voronoi_->cell(i);
            return *this;
        }
       public:
        cell_iterator(int index, const Voronoi* voronoi) : Base(index, 0, voronoi->n_cells()), voronoi_(voronoi) {
            if (index_ < voronoi_->n_cells()) this->val_ = voronoi_->cell(index_);
        }
    };
    cell_iterator cells_begin() const { return cell_iterator(0, this); }
    cell_iterator cells_end() const { return cell_iterator(n_cells(), this); }
    // perform point location for set of points p_1, p_2, \ldots, p_n
    DVector<int> locate(const DMatrix<double>& locs) const {
        fdapde_assert(locs.cols() == embed_dim);
        // find delanuay cells containing locs
        DVector<int> dual_locs = mesh_->locate(locs);
        for (int i = 0; i < locs.rows(); ++i) {
            if (dual_locs[i] == -1) continue;   // location outside domain
            // find nearest cell to i-th location
            typename Triangulation<2, 2>::CellType f = mesh_->cell(dual_locs[i]);
            SMatrix<1, Triangulation<2, 2>::n_nodes_per_cell> dist =
              (f.nodes().colwise() - locs.row(i).transpose()).colwise().squaredNorm();
            int min_index;
            double min = dist.minCoeff(&min_index);
            DVector<int> neighbors = f.neighbors();
	    for(UInt j = 0; j<neighbors.size(); ++j){ // First level neighbors
	      typename Triangulation<2, 2>::CellType f_1 = mesh_->cell(neighbors[j]);
	      SMatrix<1, Triangulation<2, 2>::n_nodes_per_cell> dist_1 =
	        (f_1.nodes().colwise() - locs.row(i).transpose()).colwise().squaredNorm();
	      UInt min_index_1;
	      double min_1 = dist_1.minCoeff(&min_index_1);
	      if(min_1<min){
	        min_index = min_index_1;
	      }
	    }
            dual_locs[i] = f.node_ids()[min_index];
        }
        return dual_locs;
    }  
   private:
    const Triangulation<2, 2>* mesh_;
    DMatrix<double> nodes_;                             // voronoi vertices
    BinaryVector<fdapde::Dynamic> nodes_markers_;       // i-th element true if i-th vertex is on boundary
    std::unordered_map<int, std::vector<int>> cells_;   // for each cell id, the ids of the vertices composing it
};

template <> class Voronoi<Triangulation<1, 1>> {
   public:
    static constexpr int local_dim = Triangulation<1, 1>::local_dim;
    static constexpr int embed_dim = Triangulation<1, 1>::embed_dim;

    Voronoi() = default;
    Voronoi(const Triangulation<1, 1>& mesh) : mesh_(&mesh) {   // constructs voronoi diagram from Delanoy triangulation
        int n_delaunay_faces = mesh_->n_cells();
        nodes_.resize(n_delaunay_faces + 2, embed_dim);
        nodes_markers_.resize(nodes_.rows());
        int k = n_delaunay_faces;
        for (typename Triangulation<1, 1>::cell_iterator it = mesh_->cells_begin(); it != mesh_->cells_end(); ++it) {
            nodes_.row(it->id()) = it->circumcenter();
            for (int v : it->node_ids()) { cells_[v].push_back(it->id()); }
            if (it->on_boundary()) {
                for (int i = 0; i < Triangulation<1, 1>::n_nodes_per_cell; ++i) {
                    if (mesh_->is_node_on_boundary(it->node_ids()[i])) {
                        nodes_.row(k) = mesh_->node(it->node_ids()[i]);
                        nodes_markers_.set(k);
                        cells_[it->node_ids()[i]].push_back(k);
                        k++;
                    }
                }
            }
        }
        // sort each cell clockwise (around its mean point)
        for (auto& [key, value] : cells_) {
            if (value[1] < value[0]) std::swap(value[0], value[1]);
        }
    }
    // cell data structure
    class VoronoiCell {
       private:
        const Voronoi* v_;
        int id_ = 0;
       public:
        VoronoiCell() = default;
        VoronoiCell(int id, const Voronoi* v) : v_(v), id_(id) { }
      double measure() const { return (v_->vertex(v_->cells_.at(id_)[1]) - v_->vertex(v_->cells_.at(id_)[0])).norm(); }
        bool on_boundary() const {
            return v_->nodes_markers_[v_->cells_.at(id_)[0]] || v_->nodes_markers_[v_->cells_.at(id_)[1]];
        }
        bool contains(const SVector<embed_dim>& p) const { return v_->locate(p)[0] == id_; }
    };
    // getters
    const DVector<double>& sites() const { return mesh_->nodes(); }
    SVector<embed_dim> vertex(int id) const { return nodes_.row(id); }
    SVector<embed_dim> site(int id) const { return mesh_->node(id); }
    const BinaryVector<fdapde::Dynamic>& boundary_vertices() const { return nodes_markers_; }
    const DMatrix<double>& vertices() const { return nodes_; }
    const Triangulation<1, 1>& dual() const { return *mesh_; }
    int n_nodes() const { return nodes_.rows(); }
    int n_cells() const { return mesh_->n_nodes(); }
    using CellType = VoronoiCell;
    CellType cell(int id) const { return VoronoiCell(id, this); }
    // iterators
    class cell_iterator : public index_based_iterator<cell_iterator, CellType> {
        using Base = index_based_iterator<cell_iterator, CellType>;
        using Base::index_;
        friend Base;
        const Voronoi* voronoi_;
        cell_iterator& operator()(int i) {
            Base::val_ = voronoi_->cell(i);
            return *this;
        }
       public:
        cell_iterator(int index, const Voronoi* voronoi) : Base(index, 0, voronoi->n_cells()), voronoi_(voronoi) {
            if (index_ < voronoi_->n_cells()) this->val_ = voronoi_->cell(index_);
        }
    };
    cell_iterator cells_begin() const { return cell_iterator(0, this); }
    cell_iterator cells_end() const { return cell_iterator(n_cells(), this); }
    // perform point location for set of points p_1, p_2, \ldots, p_n
    DVector<int> locate(const DMatrix<double>& locs) const {
        fdapde_assert(locs.cols() == embed_dim);
        // find delanuay cells containing locs
        DVector<int> dual_locs = mesh_->locate(locs);
        for (int i = 0; i < locs.rows(); ++i) {
            if (dual_locs[i] == -1) continue;   // location outside domain
            // find nearest cell to i-th location
            typename Triangulation<1, 1>::CellType f = mesh_->cell(dual_locs[i]);
            SMatrix<1, Triangulation<1, 1>::n_nodes_per_cell> dist =
              (f.nodes().colwise() - locs.row(i).transpose()).colwise().squaredNorm();
            int min_index;
            dist.minCoeff(&min_index);
            dual_locs[i] = f.node_ids()[min_index];
        }
        return dual_locs;
    }
   private:
    const Triangulation<1, 1>* mesh_;
    DMatrix<double> nodes_;                             // voronoi vertices
    BinaryVector<fdapde::Dynamic> nodes_markers_;       // i-th element true if i-th vertex is on boundary
    std::unordered_map<int, std::vector<int>> cells_;   // for each cell id, the ids of the vertices composing it
};
  
}   // namespace core
}   // namespace fdapde

#endif   // __VORONOI_H__
