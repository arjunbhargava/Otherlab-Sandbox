#include <ompl/base/spaces/SO2StateSpace.h>
#include <ompl/geometric/planners/rrt/RRT.h>
#include <ompl/geometric/planners/rrt/RRTConnect.h>
#include <ompl/geometric/planners/rrt/RRTstar.h>
#include <ompl/geometric/planners/kpiece/KPIECE1.h>
#include <ompl/geometric/planners/est/EST.h>
#include <ompl/geometric/planners/prm/PRM.h>
#include <ompl/geometric/planners/stride/STRIDE.h>
#include <ompl/tools/benchmark/Benchmark.h>
#include <ompl/base/spaces/SO3StateSpace.h>
#include <ompl/base/goals/GoalRegion.h>
#include <boost/math/constants/constants.hpp>
#include <boost/format.hpp>
#include <geode/vector/Rotation.h>
#include <geode/vector/Frame.h>
#include <geode/python/wrap.h>
#include <geode/vector/convert.h>
#include <geode/array/Nested.h>
#include <geode/openmesh/TriMesh.h>
#include <geode/geometry/SimplexTree.h>
#include <geode/geometry/ParticleTree.h>
#include <geode/geometry/traverse.h>
#include <geode/mesh/TriangleSoup.h>
#include <geode/exact/collision.h>
#include <geode/utility/curry.h>

#include <math.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <limits>


//using namespace std;
namespace ob = ompl::base;
namespace og = ompl::geometric;
using namespace geode;

using std::vector;
using std::cout;
using std::endl;
using std::cin;
using std::string;

namespace other {
namespace{
const double pi = boost::math::constants::pi<double>();
int sample_counter = 0;
//double upper_kr16_bounds [] = {185, 35, 154, 350, 130, 350};
//double lower_kr16_bounds [] = {-185, -155,-130, -350, -130, -350};
double upper_kr16_bounds [] = {185, 125, 64, 165, 130, 350};
double lower_kr16_bounds [] = {-185, -65, -210, -165, -130, -350};
//Creating a struct of each linkage, used for rendering/ forward kinematics
struct link_t {
	Vector<real,3> axis;
	double rotation_max;
	double rotation_min;
	Vector<real,3> offsets;
	bool negate_rotation;
};

struct System {
	unsigned int links;
	unsigned int robot_number;
	unsigned int stateDimension;
	vector<link_t> axis_information;
	vector<Vector<real,3>> initial_location;
	Vector<real,3> target_position;
	double initial_angle;
	double tolerance;
	vector<vector<Ref<TriMesh>>> robotMeshes;
	vector<Ref<TriMesh>> obstacleMeshes;
	vector<Ref<SimplexTree<Vector<real,3>,2>>> face_trees;
	vector<Array<Vector<real, 3>>> positions;
}; 


vector<Vector<real,3>> samples; 
Array<real> anglesFromState(System sys, const ob::State * state); 

vector<Frame<Vector<real, 3>>> frame_from_state(System sys, Array<real> &state_angles) {

	vector<Frame<Vector<real,3>>> frames;

	for(unsigned int i = 0; i < sys.stateDimension; ++i) {	
		int index = i%(sys.links);		
		auto f = Frame<Vector<real,3>>(); 
		int factor = sys.axis_information[i].negate_rotation ? -1 : 1; //For Kuka robots to move how we want them to

		if(index == 0) {  
			f.t = sys.axis_information[i].offsets + sys.initial_location[0]; //initial_location of robot 1
			auto rotation_object = Rotation<Vector<real,3>>(factor * state_angles[i], sys.axis_information[i].axis); 
			f.r = rotation_object;
			if(i == 6) {
				auto axis_rotation_object = Rotation<Vector<real,3>>(sys.initial_angle, sys.axis_information[i].axis);
				frames.push_back(Frame<Vector<real, 3>>(axis_rotation_object * sys.axis_information[i].offsets 
				+ sys.initial_location[1], rotation_object * axis_rotation_object)); //initial locatin of robot 2
			} else frames.push_back(f);

		} else {
			f = frames.back();
			auto rotation_object = Rotation<Vector<real,3>>(factor * state_angles[i], f.r * sys.axis_information[i].axis);
			frames.push_back(Frame<Vector<real,3>>(f.t + f.r * sys.axis_information[i].offsets, rotation_object * f.r ));
		}		
	}		
	return frames;
}

class ChainSpace : public ob::CompoundStateSpace {

	public: 
		//Constructor for new space called chainspace
		ChainSpace(System &sys): ob::CompoundStateSpace(), sys_(sys)//Inherits from compound state space
		{
			for(unsigned int i = 0; i < sys.stateDimension; ++i) {
				ob::StateSpacePtr subspace = ob::StateSpacePtr(new ob::SO2StateSpace());
				addSubspace(subspace, 1.0);
			}
			lock(); //Doesn't allow future changes to the subspace
		}

		//Define a distance between 2 states
		double distance(const ob::State *state1, const ob::State *state2) const
    	{	
	    	vector<Vector<real,3>> end_positions1 = effectorPositions(state1,1);
	    	vector<Vector<real,3>> end_positions2 = effectorPositions(state2,1);
	    	double distance = magnitude(end_positions1[0] - end_positions2[0]);
	    	double distance2 = magnitude(end_positions1[1] - end_positions2[1]);
	    	return sqrt(distance*distance + distance2*distance2);
   		}

   		Vector<real,2> distanceToGoal(const ob::State * state) const 
   		{
	        vector<Vector<real,3>> end_positions = effectorPositions(state,1);
		    double distance1 = magnitude(end_positions[0] -  sys_.target_position);
		  	double distance2 = magnitude(end_positions[1] - sys_.target_position);
		    return Vector<real,2>(distance1, distance2);
   		}

   		vector<Vector<real,3>> effectorPositions(const ob::State * state, int effector_flag) const
   		{
   			Array<real> joint_angles = anglesFromState(sys_, state);
   			vector<Frame<Vector<real,3>>> link_frames =	frame_from_state(sys_, joint_angles);
   			if(effector_flag == 1) {
   				vector<Vector<real,3>> end_positions;
	   			for(unsigned int i = 0; i < sys_.robot_number; i++) {
			    		int index = sys_.links* (i+1) -1;
			    		Frame<Vector<real,3>> effector = link_frames[index];
			    		end_positions.push_back(effector.t + effector.r * (100 * sys_.axis_information[index].axis));
			    }
		    return end_positions;
		   
		    } else {
		    	vector<Vector<real,3>> link_positions;
		    	for(unsigned int i = 0; i < sys_.stateDimension; i++) {
		    		Frame<Vector<real,3>> effector = link_frames[i];
		    		link_positions.push_back(effector.t);
		    	}
		    return link_positions;
		    }
   		}

	protected:
		System sys_;

};

Array<real> anglesFromState(System sys, const ob::State * state) 
{
	Array<real> joint_angles(sys.stateDimension);
		const ChainSpace::StateType *cstate1 = static_cast<const ChainSpace::StateType*>(state); 
			for(unsigned int i = 0; i < sys.stateDimension; ++i) {
				joint_angles[i] = cstate1->components[i]->as<ob::SO2StateSpace::StateType>()->value;
			}

		return joint_angles;
}

class SO26ValidityChecker : public ob::StateValidityChecker
{
	public:
	    
	    SO26ValidityChecker(const ob::SpaceInformationPtr &si, System &sys): ob::StateValidityChecker(si), sys_(sys)
	    {
			for(unsigned int i = 0; i < sys.obstacleMeshes.size(); i++) {
				obstacle_trees.push_back(sys.obstacleMeshes[i]->face_tree());
				obstacle_positions.push_back(obstacle_trees[i]->X.copy());
			}
 		 }

	    bool isValid(const ob::State *state) const
	    {	
	    	Array<real> joint_angles(sys_.stateDimension);
			const ChainSpace::StateType *cstate1 = static_cast<const ChainSpace::StateType*>(state); 
			for(unsigned int i = 0; i < sys_.stateDimension; ++i) {
				double angle = cstate1->components[i]->as<ob::SO2StateSpace::StateType>()->value * 180/pi;
				int index = i%sys_.links;
				if(angle > upper_kr16_bounds[index] || angle < lower_kr16_bounds[index]) {
					//cout << "Angle failure on linkage " << index << endl;
					return false;
				}
				joint_angles[i] = angle * pi/180;
			}
			
			generate_mesh(joint_angles);
			if(mesh_self_collisions_wrapper()) return false;
			if(mesh_other_collisions_wrapper()) return false;
			if(sys_.obstacleMeshes.size() > 0)
				return !robot_other_collisions_wrapper();
			return true;
	    }

	private:

		bool mesh_collisions(Ref<SimplexTree<Vector<real,3>,2>> face_tree1, Ref<SimplexTree<Vector<real,3>,2>> face_tree2) const {
			auto X = face_tree1->X;
			auto X2 = face_tree2->X;
			const TriangleSoup& faces = face_tree1->mesh;
			const TriangleSoup& face2 = face_tree2->mesh;

			//Find intersectino between faces
		    struct {
		      const Ref<const SimplexTree<Vector<real,3>,2>> face_tree1;
		      const Ref<const SimplexTree<Vector<real,3>,2>> face_tree2;
		      const RawArray<const Vector<real,3>> X;
		      const RawArray<const Vector<real,3>> X2;
		      int counter;

		      bool cull(const int ne, const int nf) const { return false; }

		      void leaf(const int ne, const int nf) {
		      	auto leaf_size1 = face_tree1->leaf_size;
		      	auto leaf_size2 = face_tree2->leaf_size;

		      	for(auto i = 0; i < leaf_size1; i++) {
		      		for(auto j = 0; j < leaf_size2; j++) {
				        const int face1 = face_tree1->prims(ne)[i],
				                  face2 = face_tree2->prims(nf)[j];
				        const auto fv1 = face_tree1->mesh->elements[face1];
				        const auto fv2 = face_tree2->mesh->elements[face2];
				           
				        if (triangle_triangle_intersection(X[fv1.x], X[fv1.y], X[fv1.z], X2[fv2.x], X2[fv2.y], X2[fv2.z])) 		
							throw -1;
					}
				}
		      
		      }
		    } helper({face_tree1, face_tree2, X, X2, 0});

		    try {
			    double_traverse(*face_tree1,*face_tree2,helper);
			} catch(int e) {
				return true;
			}		
		    return false;
		}
		
		bool robot_other_collisions_wrapper() const {
			for(unsigned int i = 0; i < sys_.stateDimension; i++) {
				for(unsigned int j = 0; j < obstacle_trees.size(); j++) {
					if(mesh_collisions(sys_.face_trees[i], obstacle_trees[j])) {
				//		cout << "Hit a bunny" << endl;
						return true;
					}
				}
			}
			return false;
		}

		bool mesh_other_collisions_wrapper() const {
			for(unsigned int i = 0; i < sys_.links; i++) {
				for(unsigned int j = i; j < sys_.links; j++) {
					int index = i;
					int index2 = sys_.links + j;
					if(mesh_collisions(sys_.face_trees[index], sys_.face_trees[index2])) {
				//		cout << "Got collision on  meshes:" << i << ", " << j << endl;
				//		cout << "Mesh 1 has: " << sys_.robotMeshes[0][i]->n_faces() << " faces. Mesh 2 has: " << sys_.robotMeshes[0][j]->n_faces() << " faces." << endl;
						return true;
					}
				}
			}	
			return false;
		}

		bool mesh_self_collisions_wrapper() const {
			for(unsigned int i = 0; i < sys_.robot_number; i++) {
				for(unsigned int j = 0; j < sys_.links-1; j++) {
					for(unsigned int k = j+1; k < sys_.links; k++) {
						if(k-j > 1) {
							int index = i * sys_.links + j;
							int index2 = i * sys_.links + k;
							if(mesh_collisions(sys_.face_trees[index], sys_.face_trees[index2])) {
				//				cout << "Got self collision comparing robot" << i << "And meshes:" << j << ", " << k << endl;
				//				cout << "Mesh 1 has: " << sys_.robotMeshes[i][j]->n_faces() << " faces. Mesh 2 has: " << sys_.robotMeshes[i][k]->n_faces() << " faces." << endl;
								return true;
							}
						}
					}	
				}
			}
			//cout << "This time, there were " << counter << " self collisions" << endl;
			return false;
		}
  		
		void generate_mesh(Array<real> joint_angles) const { 

			vector<Frame<Vector<real,3>>> frames = frame_from_state(sys_, joint_angles);
			for(unsigned int i = 0; i < sys_.robot_number; i++) {
				for(unsigned int j = 0; j < sys_.robotMeshes[i].size(); j++) {
					int index = i * sys_.links  + j;
					Frame<Vector<real,3>> this_frame = frames[index];
					sys_.face_trees[index]->X.const_cast_().copy(this_frame * sys_.positions[index]);
					sys_.face_trees[index]->update();
				}
			}
		}

	protected:
		System sys_;
		vector<Ref<SimplexTree<Vector<real,3>,2>>> obstacle_trees;
		vector<Array<Vector<real, 3>>> obstacle_positions;
};

vector<Vector<real,3>> computeJacobian(System &sys, ob::State* state) {
	
	ChainSpace*  space_ = new ChainSpace(sys);
	vector<Vector<real,3>> Jacobian_matrix;
	
	Array<real> joint_angles = anglesFromState(sys, state);
	vector<Frame<Vector<real,3>>> link_frames =	frame_from_state(sys, joint_angles);
	vector<Vector<real,3>> end_positions = space_->effectorPositions(state,1);
	vector<Vector<real,3>> link_positions = space_->effectorPositions(state, 0);
	Vector<real,3> delta;

	for(unsigned int i =0; i < sys.stateDimension; i++) {
		
		Vector<real,3> world_axis = link_frames[i].r * sys.axis_information[i].axis;
		world_axis = world_axis/magnitude(world_axis);

		if(i < sys.stateDimension/sys.robot_number) 
			delta = end_positions[0] - link_positions[i];
		else
			delta = end_positions[1] - link_positions[i];

		Vector<real,3> jacobian_entry = cross(world_axis, delta);
		if(sys.axis_information[i].negate_rotation) jacobian_entry = jacobian_entry * -1;

		Jacobian_matrix.push_back(jacobian_entry);
	}
	return Jacobian_matrix;
}

class nearestWeightedSampler : public ob::StateSampler
{
public:
    nearestWeightedSampler(const ob::StateSpace *space, System *sys) : ob::StateSampler(space), sys_(*sys)
    {
        closest1 = numeric_limits<double>::infinity();
        closest2 = numeric_limits<double>::infinity();
        closest_state1 = NULL;
        closest_state2 = NULL;
        this_space = static_cast<const ChainSpace *>(space_);

    }
  
    virtual void sampleUniform(ob::State *state)
    {
		Vector<real,2> distances = this_space->distanceToGoal(state);
		vector<Vector<real,3>> ends =  this_space->effectorPositions(state, 1);
		double distance = sqrt(distances[0] * distances[0] + distances[1]*distances[1]);

	   	cout << ends[0] << ", " << ends[1] << endl;
	    cout << distances[0] << ", " << distances[1] << ":" << distance << endl;
	    cout << closest1 << ", " << closest2 << endl;
	    cout << "-----" << endl;
	    string a;
		std::getline(cin, a);

	    ChainSpace::StateType* cstate = static_cast<ChainSpace::StateType*>(state);

	    if(closest1 > distances[0] || closest_state1 == NULL) {
	    	closest1 = distances[0];
	    	closest_state1 = state;
	    } 

	     if(closest2 > distances[1] || closest_state2 == NULL) {
	    	closest2 = distances[1];
	    	closest_state2 = state;
	    } 

	    double fork = rng_.uniformReal(0, 1);
	    vector<Vector<real,3>> j;
	    Array<real> joint_angles;
	    computeJacobian(sys_,state);
	  
	   	ChainSpace::StateType* best_state = static_cast<ChainSpace::StateType*>(state);
	   	ChainSpace::StateType* cstate1 = static_cast<ChainSpace::StateType*>(closest_state1);
        ChainSpace::StateType* cstate2 = static_cast<ChainSpace::StateType*>(closest_state2);

	   	for(unsigned int i = 0; i < sys_.stateDimension;i++) {
	   		if(i < sys_.stateDimension/sys_.robot_number){
	   			best_state->as<ob::SO2StateSpace::StateType>(i)->value = cstate1->as<ob::SO2StateSpace::StateType>(i)->value;
        	} else {
        		best_state->as<ob::SO2StateSpace::StateType>(i)->value = cstate2->as<ob::SO2StateSpace::StateType>(i)->value;
        	}
	   	}

	    j = computeJacobian(sys_, best_state);
		joint_angles = anglesFromState(sys_, best_state);	    

	  	vector<double> changed_distances = getRotatedDistances(frame_from_state(sys_, joint_angles), j, ends);

    	for(unsigned int i = 0; i < sys_.stateDimension; i++ ) {
    		if(fork > 0) {
	    		if(i < sys_.stateDimension/sys_.robot_number){
	    			if(distances[0] > sys_.tolerance) {
	    			cstate->as<ob::SO2StateSpace::StateType>(i)->value = 
	    				best_state->as<ob::SO2StateSpace::StateType>(i)->value + rng_.uniformReal(0, .1);// * changed_distances[i];
	    			} else {
	    				cstate->as<ob::SO2StateSpace::StateType>(i)->value = cstate->as<ob::SO2StateSpace::StateType>(i)->value;
	    			}
	    		} else {
	    			if(distances[1] > sys_.tolerance) {
	    			cstate->as<ob::SO2StateSpace::StateType>(i)->value = 
	    				best_state->as<ob::SO2StateSpace::StateType>(i)->value + rng_.uniformReal(0, .1);//* changed_distances[i];
	    			} else {
	    				cstate->as<ob::SO2StateSpace::StateType>(i)->value = cstate->as<ob::SO2StateSpace::StateType>(i)->value;
	    			}
	    		}
    		} else {
    			if(rng_.uniformReal(0,1)>.9)
	    			cstate->as<ob::SO2StateSpace::StateType>(i)->value = cstate->as<ob::SO2StateSpace::StateType>(i)->value + rng_.uniformReal(-0.5 ,0.5);
    		}
		}
    }

	virtual void sampleUniformNear(ob::State *state, const ob::State *near, const double distance) {
		cout << "Should never be in sampleUniformNear" << endl;
	 	exit(0);
	}

	virtual void sampleGaussian(ob::State *state, const ob::State *mean, const double stdDev) {
		cout << "Should never be in sampleGaussian" << endl;
	 	exit(0);
	}

private:
	vector<double> getRotatedDistances(vector<Frame<Vector<real,3>>> link_frames, vector<Vector<real,3>> j, vector<Vector<real,3>> ends) {

		double theta_eps = 1;
		vector<double> j_distances;
	    for(unsigned int i = 0; i < sys_.stateDimension; i++) {
	    	double diff;
	    	double diff2; 
	    	if(i < sys_.stateDimension/sys_.robot_number) {
		    	diff = magnitude(sys_.target_position - (ends[0] + j[i]*theta_eps ));
		    	diff2 = magnitude(sys_.target_position - (ends[0] - j[i]*theta_eps));
		    	if(diff2 < diff) diff = -diff2;
		    } else {
		    	diff = magnitude(sys_.target_position - (ends[1] + j[i]*theta_eps));
		    	diff2 = magnitude(sys_.target_position - (ends[1] - j[i]*theta_eps));
		    	if(diff2 <  diff) diff = -diff2;
		    }

	    	j_distances.push_back(diff);
	    }

	   	double mag1 = 0;
	   	double mag2 = 0;
	  	
	  	for(unsigned int i = 0; i < j_distances.size(); i++) {
	  		if(i < sys_.stateDimension/sys_.robot_number) 
		  		mag1 += abs(j_distances[i]);
		  	else mag2 += abs(j_distances[i]);
	  	}

	  	for(unsigned int i = 0; i < j_distances.size(); i++) {
	  		if(i < sys_.stateDimension/sys_.robot_number)
			  	j_distances[i] *= 1/mag1;
			else j_distances[i] *= 1/mag2;
	 	}

	 	int index1, index2 = 0;
	 	double current_max1, current_max2 = 1;
	 	for(unsigned int i = 0; i < j_distances.size(); i++) {
	 		if(i < sys_.stateDimension/sys_.robot_number) {
		 		if(abs(j_distances[i]) < current_max1) {
		 			index1 = i;
	 				current_max1 = j_distances[i];
	 			}
	 		} else {
	 			if(abs(j_distances[i]) < current_max2) {
		 			index2 = i;
	 				current_max2 = j_distances[i];
	 			}
	 		}
	 	}

	// 	cout << index1 << ", " << index2 << endl;
	 	for(unsigned int i = 0; i < j_distances.size(); i++) {
	  		j_distances[i] = 0;
	 	}

	 	if(current_max1 > 0)
		 	j_distances[index1] = 1;
		 else j_distances[index1] = -1;

		if(current_max2 > 0) 
		 	j_distances[index2] = 1;
		 else j_distances[index2] = -1;

		return j_distances;
	}

protected:
	System sys_;
    double closest1;
    double closest2;
    ob::State* closest_state1;
    ob::State* closest_state2;
    const ChainSpace* this_space;

};


ob::StateSamplerPtr allocMyValidStateSampler(System *sys, const ob::StateSpace * space)
{
    return ob::StateSamplerPtr(new nearestWeightedSampler(space, sys));
}


class  pointNormalGoal : public ob::GoalRegion
{
	public:
	    pointNormalGoal(const ob::SpaceInformationPtr &si, System &sys) : ob::GoalRegion(si), sys_(sys), space_(new ChainSpace(sys))
	    {
	    	setThreshold(sys.tolerance);
	    }
	  
	    virtual double distanceGoal(const ob::State *st) const 
	    {	
	    	vector<Vector<real,3>> end_positions1 = space_->effectorPositions(st,1);
	    	samples.push_back(end_positions1[0]);
	    	samples.push_back(end_positions1[1]);
	    	Vector<real,2> distances = space_->distanceToGoal(st);
//	    	cout << distances[0] << ", " << distances[1];

	        return sqrt(distances[0] * distances[0] + distances[1]*distances[1]);
	    }

	    	  
	protected: 
	   	System sys_;
	   	ChainSpace* space_;

};

void initializeAxes(System &sys, Array<Vector<real,3>,2> &offsets) 
{	
	vector<link_t> axis_information;
	Vector<real, 3> x_axis(1,0,0);
	Vector<real, 3> y_axis(0,1,0);
	Vector<real, 3> z_axis(0,0,1);

	for(int o = 0; o < (int) sys.robot_number; ++o) {
		auto robot_offsets = offsets[o];

		for(unsigned int i = 0; i < sys.links; ++i) {
			link_t this_node;
			this_node.rotation_min = lower_kr16_bounds[i];
			this_node.rotation_max = upper_kr16_bounds[i];
			this_node.offsets = i==0 ? robot_offsets[i] : robot_offsets[i] - robot_offsets[i-1]; //nodes.back().offsets;//;
			
			if(i == 3 || i == 5)
				this_node.axis = x_axis;
			else if(i == 0)
				this_node.axis = z_axis;
			else
				this_node.axis = y_axis;

			if(this_node.axis == z_axis || this_node.axis == x_axis)
				this_node.negate_rotation = true;
			else 
				this_node.negate_rotation = false;

			axis_information.push_back(this_node);
		}
	}
	sys.axis_information = axis_information;
}


static Nested<real> plan(unsigned int links, double robot_number, vector<Vector<real,3>> goalState, Array<Vector<real,3>,2> parsed_offsets, 
	vector<vector<Ref<TriMesh>>> robotMeshes, vector<Ref<TriMesh>> obstacleMeshes, double resolution, double range, double solve_time, 
	double initial_angle, vector<Vector<real, 3>> initial_location, double tolerance) 
{
	/*We start by setting up a system so that everything can be passed around without a giant fucking method every time*/
	System sys;
	sys.links = links;
	sys.robot_number = robot_number;
	sys.robotMeshes = robotMeshes;
	sys.obstacleMeshes = obstacleMeshes;
	sys.initial_angle = initial_angle;
	sys.initial_location = initial_location;
	sys.tolerance = tolerance;
	sys.stateDimension = links * robot_number;
	sys.target_position = Vector<real,3>(1300, 0, 1700);

	//Set up the node structure from the given files. 
	initializeAxes(sys, parsed_offsets);

	for(unsigned int i = 0; i < robot_number; i++) {
		for(unsigned int j = 0; j < robotMeshes[i].size(); j++) {
			int index = i * links  + j;
			sys.face_trees.push_back(robotMeshes[i][j]->face_tree());
			sys.positions.push_back(sys.face_trees[index]->X.copy());
		}
	}
	//Create a subspace with the given number of links
	//Create the space information pointer that everything else takes...
	ChainSpace* space_p = new ChainSpace(sys);

	ob::StateSpacePtr space(space_p);
	space->setStateSamplerAllocator(curry(allocMyValidStateSampler, &sys));

	ob::SpaceInformationPtr si(new ob::SpaceInformation(space));
	si->setStateValidityChecker(ob::StateValidityCheckerPtr(new SO26ValidityChecker(si, sys)));
	si->setStateValidityCheckingResolution(resolution);

	ob::ScopedState<ob::CompoundStateSpace> start(space); 
	ob::ScopedState<ob::CompoundStateSpace> end_state(space);	
   	end_state.random();
   	ob::ScopedState<ob::CompoundStateSpace> goal(space);
	og::PathSimplifier simplifier(si);

   	
	//Path planning for multiple goals	
   	Nested<real,false> path;
	for(size_t k = 0; k < goalState.size(); k++) {
		

		og::RRT* solver  = new og::RRT(si);

   		ob::ProblemDefinitionPtr pdef(new ob::ProblemDefinition(si));
   		solver->setRange(range);	
		ob::PlannerPtr planner(solver);

	//	sys.target_position = goalState[k];
		for(unsigned int index = 0; index < si->getStateDimension(); ++index) {
			if(k == 0) {
				start->as<ob::SO2StateSpace::StateType>(index)->value = 0;
		//		goal->as<ob::SO2StateSpace::StateType>(index)->value = goalState[k][index];
			} else { 
				start = end_state;
			}
		}
	//	cout << start << endl;
	///	cout << goal << endl;


		//pdef->setStartAndGoalStates(start, goal);
		pdef->addStartState(start);
		pdef->setGoal(ob::GoalPtr(new pointNormalGoal(si, sys)));
	    planner->setProblemDefinition(pdef);
	    planner->setup();

	  // cout << "We're now going to print the space settings" << endl;
	    si->printSettings(cout);

	   // cout << "We're now going to print the problem settings" << endl;
	   	pdef->print(cout);

	    //cout << "We've finished printing the problem settings, printing solution..." << endl;
	
	    // attempt to solve the problem within ten seconds of planning time
	    ob::PlannerStatus solved = planner->solve(solve_time);
		vector<real> step_vector;

	    if (solved)
	    {
	        ob::PathPtr solution_path = pdef->getSolutionPath();
	        og::PathGeometric geopath = dynamic_cast<og::PathGeometric &>(*solution_path);
	      //  if(simplifier.shortcutPath(geopath, 0, 0, .1)) {
	        //	simplifier.smoothBSpline(geopath);
	       // 	cout << "Simplified the path? " << endl;
	        //}
	       	vector<ob::State *> path_states = geopath.getStates();

	        for(unsigned int i = 0; i < path_states.size(); ++i) {
	        	ChainSpace::StateType* cstate1 = static_cast<ChainSpace::StateType*>(path_states[i]);
	        	for(unsigned int j = 0; j < si->getStateDimension(); ++j) {
	        		step_vector.push_back(cstate1->components[j]->as<ob::SO2StateSpace::StateType>()->value);
	        		if(i == path_states.size() - 1) {
	        			end_state->as<ob::SO2StateSpace::StateType>(j)->value = 
	        				cstate1->components[j]->as<ob::SO2StateSpace::StateType>()->value;
	        		}
	        	}
	      		path.append(asarray(step_vector).copy());
	      		step_vector.clear();
	        }
	        cout << end_state << endl;
	        cout << "Found solution" << endl;
			//solution_patth->print(cout);

	    } else {
	        cout << "No solution found" << endl;
	        step_vector.push_back(0);
	        path.append(asarray(step_vector).copy());
	        exit(0);
		}
	}
	return path;
}

vector<Vector<real,3>> sample_function(unsigned int links, double robot_number, vector<Vector<real,3>> goalState, Array<Vector<real,3>,2> parsed_offsets, 
	vector<vector<Ref<TriMesh>>> robotMeshes, vector<Ref<TriMesh>> obstacleMeshes, double resolution, double range, double solve_time, 
	double initial_angle, vector<Vector<real, 3>> initial_location, double tolerance) 
{
	plan( links, robot_number, goalState, parsed_offsets, robotMeshes, obstacleMeshes,  resolution,  range,  solve_time, initial_angle,  initial_location,  tolerance);
	return samples;
}

}
}
using namespace other;

void wrap_restructured_goal_planner(){
	geode::python::function("plan_goals",plan);
	geode::python::function("sample_path", sample_function);
}
