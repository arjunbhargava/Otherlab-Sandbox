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

#include <geode/python/wrap.h>
#include <geode/vector/convert.h>
#include <geode/array/Nested.h>
#include <geode/openmesh/TriMesh.h>

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

/* Order of things to do: 

1.Generate a few compound spaces and set them to be random. See if you can even do that DONE
2.Using a true condition on isValid, see if you can solve any of the paths DONE
3. Set limits on the space 
4. Try a few different planners/visualizing?? 
5. Try building an actual isValid function 

*/
/*This class creates a compound space with 6 SO2 subspaces. Each of the subsapces 
represents a joint angle in the KUKA R16. The space is created here without bounding, 
and can be any dimension/size.  
*/

const double pi = boost::math::constants::pi<double>();

double upper_kr16_bounds [] = {185, 35, 154, 350, 130, 350};
double lower_kr16_bounds [] = {-185, -155,-130, -350, -130, -350};

//Creating a struct of each linkage so we can refer to all that information later. 
struct link_t {

	Vector<real,3> axis;
	double rotation_max;
	double rotation_min;
	Array<Vector<real,3>> offsets;
	bool negate_rotation;

} nodes[6];

class ChainSpace : public ob::CompoundStateSpace {

	public: 
		//Constructor for new space called chainspace
		ChainSpace(unsigned int links, double linkLength): 
			ob::CompoundStateSpace(), links_(links), linkLength_(linkLength) //Inherits from compound state space
		{
			for(unsigned int i = 0; i < links; ++i) {
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
	        return dist * linkLength_;
   		}

	    double linkLength() const {
	    	return linkLength_;
	    }

	protected: 
		unsigned int links_;
		double linkLength_;
		

};

/*We define a validity checker on SO2 subspaces. Basically we're going to bound the angles based on the 
kuka r16 bounds for each of the joints. */ 

class SO26ValidityChecker : public ob::StateValidityChecker
{
	public:
	    SO26ValidityChecker(const ob::SpaceInformationPtr &si, Ref<TriMesh> mesh):
			ob::StateValidityChecker(si), mesh_(mesh), mesh_tree(mesh->face_tree())
	    {
	    }
	     
	    bool isValid(const ob::State *state) const
	    {
	        const ChainSpace* space = static_cast<const ChainSpace*>(si_->getStateSpace().get());
			const ChainSpace::StateType *cstate1 = static_cast<const ChainSpace::StateType*>(state); 

			for(unsigned int i = 0; i < si_->getStateDimension(); ++i) {
				double angle = cstate1->components[i]->as<ob::SO2StateSpace::StateType>()->value * 180/pi;
				if(angle > upper_kr16_bounds[i] || angle < lower_kr16_bounds[i])
					return false;
			 }



			return true;
	    }
	private:

/* Things to do: 
	1. Define the offsets properly 
	2. Figure out how the self-axis rotates
	3. Struct defining the axes of rotation, etc... 
	4. What is this whole compound thing? 
*/
		/*Vector<real, 3> effector_from_state(const vector<real>& joint_angles) {
			vector<Frame<Vector<real,3>>> frames;

			for(int i = 0; i < si_->getStateDimension(); ++i) {
				Frame<Vector<real,3>> f = frames.size() > 0 ? frames.back() : Frame<Vector<real, 3>>();
				frames.push_back(Frame<Vector<real,3>>(offsets[i], )
			}		

			Frame<TV> compound;
  for(int i =0; i < 6; ++i){
    compound = frames[i]*compound;
  }	
		}


	*/
	protected: 
		Ref<TriMesh> mesh_ ;
		Ref<SimplexTree<Vector<real,3>,2>> mesh_tree;

};

bool this_isValid(ob::ScopedState<ob::CompoundStateSpace> state, const ob::SpaceInformationPtr &si) {

	for(unsigned int i = 0; i < si->getStateDimension(); ++i) {
		double angle = state[i] * 180/pi;
		if(angle > upper_kr16_bounds[i] || angle < lower_kr16_bounds[i])
			return false;
	 }
	return true;
}

static Nested<real> plan(unsigned int links, double linkLength, Array<real> goalState, 
	/*Ref<TriMesh> obstacleMesh,*/ Array<Vector<real, 3>> parsed_offsets) {
	offsets = parsed_offsets;
//	Ref<TriMesh> obstacleMesh = linkLength;
	//auto mesh_tree = mesh->face_tree();

	//auto face_point = mesh_tree->closest_point("asdf");

	//Create a subspace with the given number of links
	//Create the space information pointer that everything else takes...
	ob::StateSpacePtr space(new ChainSpace(links, linkLength)); 
	ob::SpaceInformationPtr si(new ob::SpaceInformation(space));

	si->setStateValidityCheckingResolution(.0000001);
	//We initialize as a 000000 state for the start configuration
	ob::ScopedState<ob::CompoundStateSpace> start(space); 
	ob::ScopedState<ob::CompoundStateSpace> goal(space);

	for(unsigned int index = 0; index < si->getStateDimension(); ++index) {
		start->as<ob::SO2StateSpace::StateType>(index)->value = 0;
		goal->as<ob::SO2StateSpace::StateType>(index)->value = goalState[index];
	}

	//Eventually this will be where I do validity checking 
	si->setStateValidityChecker(ob::StateValidityCheckerPtr(new SO26ValidityChecker(si), obstacleMesh)));

	//Code from the box: set up a problem, solve it. Resolution parameters are changed here. 
	ob::ProblemDefinitionPtr pdef(new ob::ProblemDefinition(si));
	pdef->setStartAndGoalStates(start, goal);

    // create a planner for the defined space
   	og::RRTConnect* solver  = new og::RRTConnect(si);
   	solver->setRange(.001);
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

        cout << "Found solution:" << endl;
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

void wrap_6d_kinematics_with_collision(){
	geode::python::function("plan_2",plan);
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


