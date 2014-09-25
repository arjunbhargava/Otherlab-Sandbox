/*This is basically the same as effector_robot_1.cpp but there's a little more collision checking, look at 
the validity checker. If there's a function you can't find it's probably in RobotSystem.h */ 

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
#include <ompl/base/spaces/RealVectorStateSpace.h>
#include <math.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <limits>

#include "RobotSystem.h"


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

vector<Vector<real,3>> samples; 
Array<real> anglesFromState(const ob::State * state, RobotSystem* sys); 


//The state space is changed here a little bit. In addition to the 6 SO2 spaces, we had a 7th dimension for time. The
//RRT then tries to sample in time as well as configuration space. 
class ChainSpace : public ob::CompoundStateSpace {

	public: 
		//Constructor for new space called chainspace
		ChainSpace(RobotSystem* sys): ob::CompoundStateSpace(), sys_(sys)//Inherits from compound state space
		{
			for(unsigned int i = 0; i < sys->getStateDimension(); ++i) {
				ob::StateSpacePtr subspace = ob::StateSpacePtr(new ob::SO2StateSpace());
				addSubspace(subspace, 1.0);
			}
			ob::RealVectorStateSpace* time_space = new ob::RealVectorStateSpace();
			time_space->addDimension(0, 1);
			addSubspace(ob::StateSpacePtr(time_space), 1.0);
			lock(); //Doesn't allow future changes to the subspace
		}

		//Define a distance between 2 states
		double distance(const ob::State *state1, const ob::State *state2) const
    	{	
	    	Vector<real,3> end_position1 = sys_->effectorPositions(anglesFromState(state1, sys_));
	    	Vector<real,3> end_position2 = sys_->effectorPositions(anglesFromState(state2, sys_));
	    	double distance = magnitude(end_position1 - end_position2);
	    	return distance;
   		}

   		// double distanceToGoal(const ob::State * state) const 
   		// {
	    //     Vector<real,3> end_positions = sys_->effectorPositions(anglesFromState(state, sys_));
		   //  double distance1 = magnitude(end_positions -  sys_->getTargetPosition());
		   //  return distance1;
   		// }

	protected:
		RobotSystem* sys_;

};

Array<real> anglesFromState(const ob::State * state, RobotSystem* sys) 
{
	Array<real> joint_angles(sys->getStateDimension());
	const ChainSpace::StateType *cstate1 = static_cast<const ChainSpace::StateType*>(state); 
	for(int i = 0; i < joint_angles.size(); ++i) {
		joint_angles[i] = cstate1->components[i]->as<ob::SO2StateSpace::StateType>()->value;
	}
	return joint_angles;
}

class SO26ValidityCheckerWP : public ob::StateValidityChecker
{
	public:
	    
	    SO26ValidityCheckerWP(const ob::SpaceInformationPtr &si, RobotSystem* sys, RobotSystem* sys2, vector<Array<real>> first_motion): 
	    ob::StateValidityChecker(si), sys_(sys), sys2_(sys2), first_motion_(first_motion), space_(new ChainSpace(sys))
	    {

			for(unsigned int i = 0; i < sys2_->getObstacleMeshes().size(); i++) {
				obstacle_trees.push_back(sys2_->getObstacleMeshes()[i]->face_tree());
				obstacle_positions.push_back(obstacle_trees[i]->X.copy());
			}
 		 }

	    bool isValid(const ob::State *state) const
	    {	
	    	Array<real> joint_angles = anglesFromState(state, sys_);
	    	Vector<real,3> end_position = sys_->effectorPositions(joint_angles);
	    	samples.push_back(end_position);
	    	for(int i = 0; i < joint_angles.size(); i++) {
	    		double angle = joint_angles[i] * 180/pi;
	    		if(angle > upper_kr16_bounds[i] || angle < lower_kr16_bounds[i]) {
					cout << "Angle failure on linkage " << i << endl;
					return false;
				}
			}

			generate_mesh(joint_angles, 0, 0);
			if(mesh_self_collisions_wrapper()) return false;


//This is the only spot in the code thats' different. We need to get the timing of the sampling right, so this 
//part of the code just makes sure that we're somewhere between two valid states so that when we line them up later
//we don't introduce collisions. 
			
			const ChainSpace::StateType* cstate = static_cast<const ChainSpace::StateType*>(state);
			double timestep = cstate->components[6]->as<ob::RealVectorStateSpace::StateType>()->values[0];
			int upper = 0;
			int lower = 0;
			for(size_t i = 0; i < first_motion_.size(); i++) {
				if(timestep >= first_motion_[i][6]) 
					lower = i;

				if(timestep <= first_motion_[i][6]) {
					upper = i;
					break;
				}
			}

			// cout << timestep << ", " << first_motion_[lower][6] << ", " << first_motion_[upper][6] << endl;
			// cout << "testing lower bounds: " << lower << endl;
			generate_mesh(joint_angles, 1, lower);
			if(!robot_other_collisions_wrapper()) {
				// cout << "testing upper bounds: " << upper << endl;
				generate_mesh(joint_angles, 1, upper);
				return !robot_other_collisions_wrapper();
			} else {
			 	return false;
			}
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
		
		bool mesh_self_collisions_wrapper() const {
			for(unsigned int i = 0; i < sys_->getStateDimension()-1; i++) {
				for(unsigned int j = i+1; j < sys_->getStateDimension(); j++) {
					if(j - i > 1) {
						if(mesh_collisions(sys_->getFaceTree()[i], sys_->getFaceTree()[j])) {
						cout << "Got self collision comparing robot" << i << "And meshes:" << j << endl;
							return true;
						}
					}
				}	
			}
			//cout << "This time, there were " << counter << " self collisions" << endl;
			return false;
		}

		bool robot_other_collisions_wrapper() const {
			for(unsigned int i = 0; i < sys_->getStateDimension(); i++) {
				for(unsigned int j = 0; j < sys2_->getObstacleMeshes().size(); j++) {
					if(mesh_collisions(sys_->getFaceTree()[i], obstacle_trees[j])) {
						cout << "Hit a bike" << endl;
						return true;
					}
				}

				for(unsigned int k = 0; k < sys2_->getStateDimension(); k++) {
					if(mesh_collisions(sys_->getFaceTree()[i], sys2_->getFaceTree()[k])) {
						cout << "Hit robot 1 mesh: " << k << " with mesh " << i << endl;
						return true;
					}
				}
			}
			return false;
		}

		void generate_mesh(Array<real> joint_angles, int flag, int index) const { 

			if(!flag)
				sys_->update_mesh(joint_angles);
			else
			{
				sys_->update_mesh(joint_angles);
				sys2_->update_mesh(first_motion_[index]);

				for(unsigned int i = 0; i < obstacle_trees.size(); i++) {
					int end = sys2_->getStateDimension() -1;
				 	Frame<Vector<real,3>> end_frame = sys2_->frame_from_state(first_motion_[index])[end];
				 	Vector<real,3> additional_offset = end_frame.r * sys2_->getAxisInformation()[end].axis * sys2_->getEffectorOffset();
				 	obstacle_trees[i]->X.const_cast_().copy(end_frame * obstacle_positions[i] + additional_offset);
				 	obstacle_trees[i]->update();
				}
			}			
		}

	protected:
		RobotSystem* sys_;
		RobotSystem* sys2_;
		vector<Array<real>> first_motion_;
		vector<Ref<SimplexTree<Vector<real,3>,2>>> obstacle_trees;
		vector<Array<Vector<real, 3>>> obstacle_positions;
		ChainSpace* space_;
};

// class  pointNormalGoal : public ob::GoalStates
// {
// 	public:
// 	    pointNormalGoal(const ob::SpaceInformationPtr &si) : ob::GoalStates(si), space_(new ChainSpace(sys))
// 	    {
// 	    	setThreshold(sys.tolerance);
// 	    }

// 	    virtual double distanceGoal(const ob::State *st) const 
// 	    {	
// 	    	double distance = space_->distanceToGoal(st);
// 	        return distance;
// 	    }

// 	protected: 
// 	   	ChainSpace* space_;
// };


static vector<vector<vector<real>>> plan(unsigned int links, vector<Vector<real,3>> goalState, Array<Vector<real,3>> parsed_offsets, 
	vector<Ref<TriMesh>> robotMesh, vector<Ref<TriMesh>> obstacleMeshes, double resolution, double range, double solve_time, 
	double initial_angle, Vector<real, 3> initial_location, double tolerance, vector<vector<Array<real>>> goalAngles, bool smoothingFlag, 
	double effector_offset, vector<vector<Array<real>>> first_path) 
{
	/*We start by setting up a system so that everything can be passed around without a giant fucking method every time*/
	RobotSystem* sys = new RobotSystem(links, parsed_offsets, robotMesh, vector<Ref<TriMesh>>(), initial_angle, initial_location, 0);
	RobotSystem* sys2 = new RobotSystem(links, parsed_offsets, robotMesh, obstacleMeshes, 0, Vector<real,3>(0, 0, 0), effector_offset);

	ob::StateSpacePtr temp_space(new ChainSpace(sys));
	ob::ScopedState<ob::CompoundStateSpace> start(temp_space);
	ob::ScopedState<ob::CompoundStateSpace> end_state(temp_space);	
	end_state.random();

	//Path planning for multiple goals	
	vector<vector<vector<real>>> path;
	for(size_t k = 0; k < goalState.size(); k++) 
	{
		sys->setTargetPosition(goalState[k]);			
		ob::StateSpacePtr space(new ChainSpace(sys));
		ob::SpaceInformationPtr si(new ob::SpaceInformation(space));
		SO26ValidityCheckerWP* checker = new SO26ValidityCheckerWP(si, sys, sys2, first_path[k]);
		si->setStateValidityChecker(ob::StateValidityCheckerPtr(checker));
		si->setStateValidityCheckingResolution(resolution);

		if(k == 0) {
			for(unsigned int index = 0; index < si->getStateDimension(); ++index) 
			{
				if(index < sys->getStateDimension()) 
					start->as<ob::SO2StateSpace::StateType>(index)->value = 0;
			}
		} else { 
			start = end_state;
		}
		start->as<ob::RealVectorStateSpace::StateType>(6)->values[0] = 0;


		ob::GoalStates* goalSet = new ob::GoalStates(si);//pointNormalGoal(si,sys);

		for(unsigned int i = 0 ; i < goalAngles[k].size(); i++) 
		{	
			ob::State *state = space->allocState();
			ChainSpace::StateType *cstate = static_cast<ChainSpace::StateType*>(state); 
			for(unsigned int j = 0; j < sys->getStateDimension(); j++) {		
				cstate->components[j]->as<ob::SO2StateSpace::StateType>()->value = goalAngles[k][i][j];
			}
			cstate->components[si->getStateDimension()-1]->as<ob::RealVectorStateSpace::StateType>()->values[0] = 
				goalAngles[k][i][si->getStateDimension()-1];

			if(checker->isValid(state)) 
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

	    cout << "Did we get through the setup?" << endl;
	  // cout << "We're now going to print the space settings" << endl;
//	    si->printSettings(cout);

	   // cout << "We're now going to print the problem settings" << endl;
//	   	pdef->print(cout);
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
	        		if(j == si->getStateDimension()-1)
	        			step_vector.push_back(cstate1->components[j]->as<ob::RealVectorStateSpace::StateType>()->values[0]);
	        		else step_vector.push_back(cstate1->components[j]->as<ob::SO2StateSpace::StateType>()->value);
	        		
	        		if(i == path_states.size() - 1) {
	        			if(j == si->getStateDimension()-1) {
	        				end_state->as<ob::RealVectorStateSpace::StateType>(j)->values[0] = 
	        					cstate1->components[j]->as<ob::RealVectorStateSpace::StateType>()->values[0];
	        			} else {
	        				end_state->as<ob::SO2StateSpace::StateType>(j)->value = 
	        					cstate1->components[j]->as<ob::SO2StateSpace::StateType>()->value;
	        			}
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
	double initial_angle, Vector<real, 3> initial_location, double tolerance, vector<vector<Array<real>>> goalAngles, bool smoothingFlag, double effector_offset, vector<vector<Array<real>>> first_path) 
{
	plan(links, goalState, parsed_offsets, robotMesh, obstacleMeshes,  resolution,  range, 
	 solve_time, initial_angle,  initial_location,  tolerance, goalAngles, smoothingFlag, effector_offset, first_path );
	return samples;
}

}
}
using namespace other;

void wrap_robot_2_planner() {
	geode::python::function("plan_2",plan);
	geode::python::function("sample_path_2", sample_function);
}
