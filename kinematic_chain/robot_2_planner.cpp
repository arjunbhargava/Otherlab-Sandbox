#include <ompl/base/spaces/SO2StateSpace.h>
#include <ompl/geometric/planners/rrt/RRT.h>
#include <ompl/geometric/planners/rrt/RRTConnect.h>
#include <ompl/geometric/planners/rrt/RRTstar.h>
#include <ompl/geometric/planners/kpiece/KPIECE1.h>
#include <ompl/geometric/planners/est/EST.h>
#include <ompl/geometric/planners/prm/PRM.h>
#include <ompl/geometric/planners/stride/STRIDE.h>
#include <ompl/tools/benchmark/Benchmark.h>
#include <ompl/base/goals/GoalSampleableRegion.h>
#include <ompl/base/spaces/SO3StateSpace.h>
#include <ompl/base/goals/GoalRegion.h>
#include <ompl/base/goals/GoalStates.h>
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
#include <ompl/util/Exception.h>
#include <geode/utility/curry.h>
#include <ompl/util/RandomNumbers.h>
#include "RobotSystem.h"

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
namespace {

const double pi = boost::math::constants::pi<double>();
int sample_counter = 0;
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
	unsigned int stateDimension;
	vector<link_t> axis_information;
	Vector<real,3> initial_location;
	Vector<real,3> target_position;
	double initial_angle;
	double tolerance;
	vector<Ref<TriMesh>> robotMesh;
	vector<Ref<TriMesh>> obstacleMeshes;
	vector<Ref<SimplexTree<Vector<real,3>,2>>> face_trees;
	vector<Array<Vector<real, 3>>> positions;
	double effector_offset;
}; 

vector<Vector<real,3>> samples; 
Array<real> anglesFromState(const System &sys, const ob::State * state); 
class ChainSpace;

vector<Frame<Vector<real, 3>>> frame_from_state(System sys, Array<real> &state_angles) {

	vector<Frame<Vector<real,3>>> frames;

	for(unsigned int i = 0; i < sys.stateDimension; ++i) {	
		auto f = Frame<Vector<real,3>>(); 
		int factor = sys.axis_information[i].negate_rotation ? -1 : 1; //For Kuka robots to move how we want them to

		if(i == 0) {  
			auto rotation_object = Rotation<Vector<real,3>>(factor * state_angles[i], sys.axis_information[i].axis); 
			auto axis_rotation_object = Rotation<Vector<real,3>>(sys.initial_angle, sys.axis_information[i].axis);
			frames.push_back(Frame<Vector<real, 3>>(axis_rotation_object * sys.axis_information[i].offsets 
				+ sys.initial_location, rotation_object * axis_rotation_object)); 
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
	    	Vector<real,3> end_position1 = effectorPositions(state1,1);
	    	Vector<real,3> end_position2 = effectorPositions(state2,1);
	    	double distance = magnitude(end_position1 - end_position2);
	    	return distance;
   		}

   		double distanceToGoal(const ob::State * state) const 
   		{
	        Vector<real,3> end_positions = effectorPositions(state,1);
		    double distance1 = magnitude(end_positions -  sys_.target_position);
		    return distance1;
   		}

   		Vector<real,3> effectorPositions(const ob::State * state, int effector_flag) const
   		{
   			Array<real> joint_angles = anglesFromState(sys_, state);
   			vector<Frame<Vector<real,3>>> link_frames =	frame_from_state(sys_, joint_angles);
   			if(effector_flag == 1) {
   				Vector<real,3> end_position;
	    		int index = sys_.stateDimension - 1;
	    		Frame<Vector<real,3>> effector = link_frames[index];
	    		end_position = effector.t + effector.r * (sys_.effector_offset * sys_.axis_information[index].axis);
		    	return end_position;
			}
			return Vector<real,3>(0, 0,0);
   		}

	protected:
		System sys_;

};


Array<real> anglesFromState(const System &sys, const ob::State * state) 
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
	    
	    SO26ValidityChecker(const ob::SpaceInformationPtr &si, System &sys): ob::StateValidityChecker(si), sys_(sys), space_(new ChainSpace(sys))
	    {
			for(unsigned int i = 0; i < sys.obstacleMeshes.size(); i++) {
				obstacle_trees.push_back(sys.obstacleMeshes[i]->face_tree());
				obstacle_positions.push_back(obstacle_trees[i]->X.copy());
			}
 		 }

	    bool isValid(const ob::State *state) const
	    {	
	    	Vector<real,3> end_position = space_->effectorPositions(state,1);
	    	samples.push_back(end_position);
	    	cout << "Something's happening? " << endl;

	    	Array<real> joint_angles = anglesFromState(sys_, state);
	    	for(int i = 0; i < joint_angles.size(); i++) {
	    		double angle = joint_angles[i] * 180/pi;
	    		if(angle > upper_kr16_bounds[i] || angle < lower_kr16_bounds[i]) {
					//cout << "Angle failure on linkage " << index << endl;
					return false;
				}
			}
			
			generate_mesh(joint_angles);
			if(mesh_self_collisions_wrapper()) return false;
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

		      bool cull(const int ne, const int nf) const { return false; }

		      void leaf(const int ne, const int nf) {
		      	auto leaf_size1 = face_tree1->leaf_size;
		      	auto leaf_size2 = face_tree2->leaf_size;

		      	for(auto face1 : face_tree1->prims(ne)) {
		      		for(auto face2 : face_tree2->prims(nf)) {
				        const auto fv1 = face_tree1->mesh->elements[face1];
				        const auto fv2 = face_tree2->mesh->elements[face2];
				           
				        if (triangle_triangle_intersection(X[fv1.x], X[fv1.y], X[fv1.z], X2[fv2.x], X2[fv2.y], X2[fv2.z])) 		
							throw -1;
					}
				}
		      
		      }

		    } helper({face_tree1, face_tree2, X, X2});

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
						cout << "Hit a bunny" << endl;
						return true;
					}
				}
			}
			return false;
		}

		bool mesh_self_collisions_wrapper() const {
			for(unsigned int i = 0; i < sys_.stateDimension-1; i++) {
				for(unsigned int j = i+1; j < sys_.stateDimension; j++) {
					if(j - i > 1) {
						if(mesh_collisions(sys_.face_trees[i], sys_.face_trees[j])) {
						cout << "Got self collision comparing robot" << i << "And meshes:" << j << endl;
					//		cout << "Mesh 1 has: " << sys_.robotMeshes[i][j]->n_faces() << " faces. Mesh 2 has: " << sys_.robotMeshes[i][k]->n_faces() << " faces." << endl;
							return true;
						}
					}
				}	
			}
			//cout << "This time, there were " << counter << " self collisions" << endl;
			return false;
		}
  		
		void generate_mesh(Array<real> joint_angles) const { 

			vector<Frame<Vector<real,3>>> frames = frame_from_state(sys_, joint_angles);
			for(unsigned int i = 0; i < sys_.robotMesh.size(); i++) {
				Frame<Vector<real,3>> this_frame = frames[i];
				sys_.face_trees[i]->X.const_cast_().copy(this_frame * sys_.positions[i]);
				sys_.face_trees[i]->update();
			}

			for(unsigned int i = 0; i < obstacle_trees.size(); i++) {
				Frame<Vector<real,3>> end_frame = frames[sys_.stateDimension-1];
				Vector<real,3> additional_offset = end_frame.r * sys_.axis_information[sys_.stateDimension-1].axis * sys_.effector_offset;
				obstacle_trees[i]->X.const_cast_().copy(end_frame * obstacle_positions[i] + additional_offset);
				obstacle_trees[i]->update();
			}
		}

	protected:
		System sys_;
		vector<Ref<SimplexTree<Vector<real,3>,2>>> obstacle_trees;
		vector<Array<Vector<real, 3>>> obstacle_positions;
		ChainSpace* space_;
};

class  pointNormalGoal : public ob::GoalStates
{
	public:
	    pointNormalGoal(const ob::SpaceInformationPtr &si, System &sys) : ob::GoalStates(si), sys_(sys), space_(new ChainSpace(sys))
	    {
	    	setThreshold(sys.tolerance);
	    }

	    virtual double distanceGoal(const ob::State *st) const 
	    {	
	    	Vector<real,3>  end_positions1 = space_->effectorPositions(st,1);
	    	samples.push_back(end_positions1);
	    	double distance = space_->distanceToGoal(st);
	        return distance;
	    }

	    	  
	protected: 
	   	System sys_;
	   	ChainSpace* space_;

};


void initializeAxes(System &sys, Array<Vector<real,3>> &offsets) 
{	
	vector<link_t> axis_information;
	Vector<real, 3> x_axis(1,0,0);
	Vector<real, 3> y_axis(0,1,0);
	Vector<real, 3> z_axis(0,0,1);

	auto robot_offsets = offsets;

	for(unsigned int i = 0; i < sys.stateDimension; ++i) {
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

	sys.axis_information = axis_information;
}


static vector<vector<vector<real>>> plan(unsigned int links, vector<Vector<real,3>> goalState, Array<Vector<real,3>> parsed_offsets, 
	vector<Ref<TriMesh>> robotMesh, vector<Ref<TriMesh>> obstacleMeshes, double resolution, double range, double solve_time, 
	double initial_angle, Vector<real, 3> initial_location, double tolerance, vector<vector<Array<real>>> goalAngles, bool smoothingFlag, double effector_offset) 
{
	/*We start by setting up a system so that everything can be passed around without a giant fucking method every time*/
	System sys;
	sys.robotMesh = robotMesh;
	sys.obstacleMeshes = obstacleMeshes;  
	sys.initial_angle = initial_angle;
	sys.initial_location = initial_location;
	sys.tolerance = tolerance;
	sys.stateDimension = links;
	sys.effector_offset = effector_offset;

	RobotSystem * sys2 = new RobotSystem(links, parsed_offsets, robotMesh, obstacleMeshes, initial_angle, initial_location, effector_offset);
	cout << sys2->getStateDimension() << endl;
	exit(0);
	//Set up the node structure from the given files. 
	initializeAxes(sys, parsed_offsets);
	//Initialize the facetrees from the default positions of the robots. 
	//Mostly just want to have this original configuration, everything else is in terms of angles
	for(unsigned int i = 0; i < robotMesh.size(); i++) {
		sys.face_trees.push_back(robotMesh[i]->face_tree());
		sys.positions.push_back(sys.face_trees[i]->X.copy());
	}

	ob::StateSpacePtr temp_space(new ChainSpace(sys));
	ob::ScopedState<ob::CompoundStateSpace> start(temp_space);
	ob::ScopedState<ob::CompoundStateSpace> end_state(temp_space);	
	end_state.random();
//	ob::ScopedState<ob::CompoundStateSpace> goal(space);
	//Path planning for multiple goals	
	vector<vector<vector<real>>> path;
	for(size_t k = 0; k < goalState.size(); k++) 
	{
		sys.target_position = goalState[k];			
		ob::StateSpacePtr space(new ChainSpace(sys));
		ob::SpaceInformationPtr si(new ob::SpaceInformation(space));
		si->setStateValidityChecker(ob::StateValidityCheckerPtr(new SO26ValidityChecker(si, sys)));
		si->setStateValidityCheckingResolution(resolution);

		for(unsigned int index = 0; index < si->getStateDimension(); ++index) 
		{
			if(k == 0) {
				start->as<ob::SO2StateSpace::StateType>(index)->value = 0;
			} else { 
				start = end_state;
				cout << start << endl;
				break;
			}
		}

		ob::GoalStates* goalSet = new pointNormalGoal(si,sys);

		for(unsigned int i = 0 ; i < goalAngles[k].size(); i++) 
		{	
			ob::State *state = space->allocState();
			ChainSpace::StateType *cstate = static_cast<ChainSpace::StateType*>(state); 
			for(int j = 0; j < goalAngles[k][i].size(); j++) {		
				cstate->components[j]->as<ob::SO2StateSpace::StateType>()->value = goalAngles[k][i][j];
			}
			goalSet->addState(state);
		}
		
		og::PathSimplifier simplifier(si);
		og::RRTConnect* solver  = new og::RRTConnect(si);
   		solver->setRange(range);	
   		ob::PlannerPtr planner(solver);

   		ob::ProblemDefinitionPtr pdef(new ob::ProblemDefinition(si));
		pdef->addStartState(start);
	   	pdef->setGoal(ob::GoalPtr(goalSet));
	    planner->setProblemDefinition(pdef);
	    planner->setup();

	  // cout << "We're now going to print the space settings" << endl;
	    si->printSettings(cout);

	   // cout << "We're now going to print the problem settings" << endl;
	   	pdef->print(cout);
	   	cout << "I know everything got set up" << endl;
	    //cout << "We've finished printing the problem settings, printing solution..." << endl;
	    // attempt to solve the problem within ten seconds of planning time
	    ob::PlannerStatus solved = planner->solve(solve_time);
	    cout << "Did we get here??" << endl;
		vector<real> step_vector;
		vector<vector<real>> all_steps;
	    if (solved)
	    {
	        ob::PathPtr solution_path = pdef->getSolutionPath();
	        og::PathGeometric geopath = dynamic_cast<og::PathGeometric &>(*solution_path);
	       	if(smoothingFlag) {
		        if(simplifier.shortcutPath(geopath, 0, 0, .1)) {
		        	simplifier.smoothBSpline(geopath);
		        	cout << "Simplified the path? " << endl;
		        }
	    	}
	       	vector<ob::State *> path_states = geopath.getStates();

	        for(unsigned int i = 0; i < path_states.size(); ++i) 
	        {
	        	ChainSpace::StateType* cstate1 = static_cast<ChainSpace::StateType*>(path_states[i]);
	        	for(unsigned int j = 0; j < si->getStateDimension(); ++j)
	        	{
	        		step_vector.push_back(cstate1->components[j]->as<ob::SO2StateSpace::StateType>()->value);
	        		if(i == path_states.size() - 1) {
	        			end_state->as<ob::SO2StateSpace::StateType>(j)->value = 
	        				cstate1->components[j]->as<ob::SO2StateSpace::StateType>()->value;
	        		}
	        	}
	        	all_steps.push_back(step_vector);
	      		step_vector.clear();
	        }
	        path.push_back(all_steps);
	        cout << end_state << endl;
	        cout << "Found solution" << endl;

	    } else {
	        cout << "No solution found" << endl;
	        step_vector.push_back(0);
	        all_steps.push_back(step_vector);
	        path.push_back(all_steps);
		}
	}
	return path;
}

vector<Vector<real,3>> sample_function(unsigned int links, vector<Vector<real,3>> goalState, Array<Vector<real,3>> parsed_offsets, 
	vector<Ref<TriMesh>> robotMesh, vector<Ref<TriMesh>> obstacleMeshes, double resolution, double range, double solve_time, 
	double initial_angle, Vector<real, 3> initial_location, double tolerance, vector<vector<Array<real>>> goalAngles, bool smoothingFlag, double effector_offset) 
{
	plan(links, goalState, parsed_offsets, robotMesh, obstacleMeshes,  resolution,  range,  solve_time, initial_angle,  initial_location,  tolerance, goalAngles, smoothingFlag, effector_offset);
	return samples;
}

}
}
using namespace other;

void wrap_robot_1_planner(){
	geode::python::function("plan_1",plan);
	geode::python::function("sample_path_1", sample_function);
}
