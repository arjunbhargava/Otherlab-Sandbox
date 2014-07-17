#include <ompl/base/spaces/SO2StateSpace.h>
#include <ompl/geometric/planners/rrt/RRT.h>
#include <ompl/geometric/planners/rrt/RRTConnect.h>
#include <ompl/geometric/planners/kpiece/KPIECE1.h>
#include <ompl/geometric/planners/est/EST.h>
#include <ompl/geometric/planners/prm/PRM.h>
#include <ompl/geometric/planners/stride/STRIDE.h>
#include <ompl/tools/benchmark/Benchmark.h>
#include <ompl/base/spaces/SO3StateSpace.h>
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

#include <fstream>
#include <iostream>
#include <cmath>

//using namespace std;
namespace ob = ompl::base;
namespace og = ompl::geometric;
using namespace geode;

using std::vector;
using std::cout;
using std::endl;

namespace other {

const double pi = boost::math::constants::pi<double>();

double upper_kr16_bounds [] = {185, 35, 154, 350, 130, 350};
double lower_kr16_bounds [] = {-185, -155,-130, -350, -130, -350};

//Creating a struct of each linkage, used for rendering/ forward kinematics
struct link_t {
	Vector<real,3> axis;
	double rotation_max;
	double rotation_min;
	Vector<real,3> offsets;
	bool negate_rotation;
};


/*This class creates a compound space with 6 SO2 subspaces. Each of the subsapces 
represents a joint angle in the KUKA R16. The space is created here without bounding, 
and can be any dimension/size.  
*/

class ChainSpace : public ob::CompoundStateSpace {

	public: 
		//Constructor for new space called chainspace
		ChainSpace(unsigned int links, double robot_number): 
			ob::CompoundStateSpace(), links_(links), robot_number_(robot_number) //Inherits from compound state space
		{
			for(unsigned int i = 0; i < links * robot_number; ++i) {
				ob::StateSpacePtr subspace = ob::StateSpacePtr(new ob::SO2StateSpace());
				addSubspace(subspace, 1.0);
			}
			lock(); //Doesn't allow future changes to the subspace
		}

		//Define a distance between 2 states
		double distance(const ompl::base::State *state1, const ompl::base::State *state2) const
    	{
	        const StateType *cstate1 = static_cast<const StateType*>(state1);
	        const StateType *cstate2 = static_cast<const StateType*>(state2);
	        double theta1 = 0., theta2 = 0., dx = 0., dy = 0., dist = 0.;

	        for (unsigned int i = 0; i < getSubspaceCount(); ++i)
	        {
	            theta1 = cstate1->components[i]->as<ob::SO2StateSpace::StateType>()->value;
	            theta2 = cstate2->components[i]->as<ob::SO2StateSpace::StateType>()->value;
	            dx = cos(theta1) - cos(theta2);
	            dy = sin(theta1) - sin(theta2);
	            dist += sqrt(dx * dx + dy * dy);
	        }
	        return dist;
   		}

	    double robotNumber() const {
	    	return robot_number_;
	    }

	protected: 
		unsigned int links_;
		double robot_number_;		
};

/*We define a validity checker on SO2 subspaces. Basically we're going to bound the angles based on the 
kuka r16 bounds for each of the joints. Currently working on mesh collision with objects in real space, 
with the current object in real space being the other Kuka robot. */ 

class SO26ValidityChecker : public ob::StateValidityChecker
{
	public:
	    SO26ValidityChecker(const ob::SpaceInformationPtr &si,vector<link_t> &nodes,  Ref<TriMesh> mesh):
			ob::StateValidityChecker(si),  mesh_(mesh), mesh_tree(mesh->face_tree()), nodes_(nodes), 
				stateDimension(si->getStateDimension())
	    {
	    }
	     
	    bool isValid(const ob::State *state) const
	    {
	    	return true;
	    	Array<real> joint_angles(stateDimension);
	        //const ChainSpace* space = static_cast<const ChainSpace*>(si_->getStateSpace().get());
			const ChainSpace::StateType *cstate1 = static_cast<const ChainSpace::StateType*>(state); 
			for(unsigned int i = 0; i < stateDimension; ++i) {
				double angle = cstate1->components[i]->as<ob::SO2StateSpace::StateType>()->value * 180/pi;
				if(angle > upper_kr16_bounds[i] || angle < lower_kr16_bounds[i]) {
					cout << "Failed here " << endl;
					return false;
				}
				joint_angles[i] = angle * pi/180;
			 }

			return true;// computeDistance(effector_from_state(joint_angles));

	    }

	private:

		Vector<real, 3> effector_from_state(const Array<real> &state_angles) const {
	
			vector<Frame<Vector<real,3>>> frames;
			//Rotation<Vector<real, 3>> rotation_object = Rotation<Vector<real, 3>>(0, nodes[0].axis);

			for(unsigned int i = 0; i < stateDimension; ++i) {			
				auto f = Frame<Vector<real,3>>(); //Let's create a frame
				int factor = nodes_[i].negate_rotation ? -1 : 1; //For Kuka robots to move how we want them to

				if(i == 0) {
					f.t = nodes_[0].offsets; //Offset of the base shouldn't ever change..
					//The axis doesn't either, but the rotation angle itself does.  
					auto rotation_object = Rotation<Vector<real,3>>(factor * state_angles[i], nodes_[i].axis); 
					f.r = rotation_object;
					frames.push_back(f);
					
				} else {
					//If it's not the first link, then get the previous link 
					f = frames.back();
					//Your rotation axis should be previous guy's rotation applied to your axis.
					auto rotation_object = Rotation<Vector<real,3>>(factor * state_angles[i], f.r * nodes_[i].axis);

					//Your position is the offset of the other guy + ABSOLUTE rotation * your offset. Your absolute rotation 
					//is the previous absolute rotation multiplied by your relative rotation. Note the order of matrix multiplication
					frames.push_back(Frame<Vector<real,3>>(f.t + f.r * nodes_[i].offsets, rotation_object * f.r ));
				}		
			}		
			
			return frames[stateDimension-1].t;
		}

		bool computeDistance(Vector<real,3> effector_position) const {
			auto difference = mesh_tree->distance(effector_position);
			return difference > 10;
		}

	protected: 
		Ref<TriMesh> mesh_;
		Ref<SimplexTree<Vector<real,3>,2>> mesh_tree;
		vector<link_t> nodes_;
		unsigned int stateDimension;
		int max_distance; 
};

/* Initializes the axes for each linkage (sets up frames) on which the forward kinematics
are eventually computed. These are then used for validity checking. */

void initializeAxes(ob::SpaceInformationPtr &si, vector<link_t> &nodes, Array<Vector<real,3>,2> &offsets) {
	
	Vector<real, 3> x_axis(1,0,0);
	Vector<real, 3> y_axis(0,1,0);
	Vector<real, 3> z_axis(0,0,1);

	for(unsigned int i = 0; i < si->getStateDimension(); ++i) {
		link_t this_node;
		this_node.rotation_min = lower_kr16_bounds[i];
		this_node.rotation_max = upper_kr16_bounds[i];
		this_node.offsets = i==0 ? offsets[0][i] : offsets[0][i] - offsets[0][i-1]; //nodes.back().offsets;//;
		
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

		nodes.push_back(this_node);
	}
}

static Nested<real> plan(unsigned int links, double robot_number, Array<real> goalState, Array<Vector<real,3>,2> parsed_offsets, 
	Ref<TriMesh> obstacleMesh) {

	//Create a subspace with the given number of links
	//Create the space information pointer that everything else takes...
	ob::StateSpacePtr space(new ChainSpace(links, robot_number)); 
	ob::SpaceInformationPtr si(new ob::SpaceInformation(space));

	si->setStateValidityCheckingResolution(.0000001);

	//Set up the node structure from the given files. 
	vector<link_t> nodes;
	initializeAxes(si, nodes, parsed_offsets);

	//We initialize as a 000000 state for the start configuration
	ob::ScopedState<ob::CompoundStateSpace> start(space); 
	ob::ScopedState<ob::CompoundStateSpace> goal(space);

	for(unsigned int index = 0; index < si->getStateDimension(); ++index) {
		start->as<ob::SO2StateSpace::StateType>(index)->value = 0;
		if(index == 5) // I don't really like this, I'll start the robots looking the same way in later iterations.
			start->as<ob::SO2StateSpace::StateType>(index)->value = pi;
		goal->as<ob::SO2StateSpace::StateType>(index)->value = goalState[index];
	}

	//Eventually this will be where I do validity checking 
	si->setStateValidityChecker(ob::StateValidityCheckerPtr(new SO26ValidityChecker(si, nodes, obstacleMesh)));

	//Code from the box: set up a problem, solve it. Resolution parameters are changed here. 
	ob::ProblemDefinitionPtr pdef(new ob::ProblemDefinition(si));
	pdef->setStartAndGoalStates(start, goal);

    // create a planner for the defined space
   	og::RRTConnect* solver  = new og::RRTConnect(si);
   	solver->setRange(.1);
    ob::PlannerPtr planner(solver);
    planner->setProblemDefinition(pdef);
    planner->setup();

    cout << "We're now going to print the space settings" << endl;
    si->printSettings(cout);

    cout << "We're now going to print the problem settings" << endl;
    pdef->print(cout);

    cout << "We've finished printing the problem settings, printing solution..." << endl;

    // attempt to solve the problem within ten seconds of planning time
    ob::PlannerStatus solved = planner->solve(10.0);
	Nested<real,false> path;
	vector<real> step_vector;


    if (solved)
    {
        ob::PathPtr solution_path = pdef->getSolutionPath();
        og::PathGeometric geopath = dynamic_cast<og::PathGeometric &>(*solution_path);
       	vector<ob::State *> path_states = geopath.getStates();

        for(unsigned int i = 0; i < path_states.size(); ++i) {
        	ChainSpace::StateType* cstate1 = static_cast<ChainSpace::StateType*>(path_states[i]);
        	for(unsigned int j = 0; j < si->getStateDimension(); ++j) {
        		step_vector.push_back(cstate1->components[j]->as<ob::SO2StateSpace::StateType>()->value);
        	}
      		path.append(asarray(step_vector).copy());
      		step_vector.clear();
        }

        cout << "Found solution" << endl;
        //solution_path->print(cout);

    } else {
        cout << "No solution found" << endl;
        step_vector.push_back(0);
        path.append(asarray(step_vector).copy());
	}

	return path;
}

}

using namespace other;

void wrap_12d_kinematics_with_collision(){
	geode::python::function("plan",plan);
}


/*
void wrap_chain_space(){

typedef ChainSpace Self;
Class<Self>("ChainSpace")
.GEODE_INIT() //constructor
.GEODE_FIELD(distance)
//.GEODE_METHOD(whatever)

}
*/


