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
#include <geode/geometry/traverse.h>
#include <geode/mesh/TriangleSoup.h>
#include <geode/exact/collision.h>

#include <math.h>
#include <cmath>
#include <fstream>
#include <iostream>


//using namespace std;
namespace ob = ompl::base;
namespace og = ompl::geometric;
using namespace geode;

using std::vector;
using std::cout;
using std::endl;


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
	    SO26ValidityChecker(const ob::SpaceInformationPtr &si, unsigned int robot_number, vector<link_t> &axis_information, 
	     vector<vector<Ref<TriMesh>>> meshes, vector<Ref<TriMesh>> obstacleMeshes, double initial_angle, vector<Vector<real,3>> initial_location):
			ob::StateValidityChecker(si),  meshes_(meshes), obstacleMeshes_(obstacleMeshes),axis_information_(axis_information),
			 robot_number_(robot_number), stateDimension(si->getStateDimension()), 
			 initial_angle_(initial_angle), initial_location_(initial_location)
	    {
	    	for(unsigned int i = 0; i < robot_number; i++) {
	    		for(unsigned int j = 0; j < meshes[i].size(); j++) {
					int index = i * (stateDimension/robot_number_)  + j;
					face_trees.push_back(meshes[i][j]->face_tree());
					positions.push_back(face_trees[index]->X.copy());
				}
			}

			//cout << positions[0][0] << " , " << positions[6][0] << endl;
			//cout << positions[1][0] << ", " << positions[7][0] << endl;

			for(unsigned int i = 0; i < obstacleMeshes.size(); i++) {
				obstacle_trees.push_back(obstacleMeshes[i]->face_tree());
				obstacle_positions.push_back(obstacle_trees[i]->X.copy());
			}
 		 }



	    bool isValid(const ob::State *state) const
	    {	
	    	Array<real> joint_angles(stateDimension);
			const ChainSpace::StateType *cstate1 = static_cast<const ChainSpace::StateType*>(state); 
			for(unsigned int i = 0; i < stateDimension; ++i) {
				double angle = cstate1->components[i]->as<ob::SO2StateSpace::StateType>()->value * 180/pi;
				int index = i%(stateDimension/robot_number_);
				if(angle > upper_kr16_bounds[index] || angle < lower_kr16_bounds[index]) {
					cout << "Angle failure on linkage " << index << endl;
					return false;
				}
				joint_angles[i] = angle * pi/180;
			}
			
			generate_mesh(joint_angles);
		//	cout << face_trees[0]->X[0] << " , " << face_trees[6]->X[0] << endl;
		//	cout << face_trees[1]->X[0] << " , " << face_trees[7]->X[0] << endl;
			//for(int i = 0; i < face_trees[11]->X.size(); i++) {
			//	cout << face_trees[11]->X[i] << endl;
			//}
			//exit(0);
			//Modify obstacle mesh if so desired... No function currently written for that. Pretty simple tho
			if(mesh_self_collisions_wrapper()) return false;
			if(mesh_other_collisions_wrapper()) return false;
			if(obstacleMeshes_.size() > 0)
				return !robot_other_collisions_wrapper();
			return true;
	    }

	private:

		bool mesh_collisions(Ref<SimplexTree<Vector<real,3>,2>> face_tree1, Ref<SimplexTree<Vector<real,3>,2>> face_tree2) const {
			//GEODE_ASSERT(face_tree1->leaf_size==4);
			//GEODE_ASSERT(face_tree2->leaf_size==4);
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
			for(unsigned int i = 0; i < stateDimension; i++) {
				for(unsigned int j = 0; j < obstacle_trees.size(); j++) {
					if(mesh_collisions(face_trees[i], obstacle_trees[j])) {
						cout << "Hit a bunny" << endl;
						return true;
					}
				}
			}
			return false;
		}

		bool mesh_other_collisions_wrapper() const {
			for(unsigned int i = 0; i < stateDimension/robot_number_; i++) {
				for(unsigned int j = i; j < stateDimension/robot_number_; j++) {
					int index = i;
					int index2 = (stateDimension/robot_number_) + j;
					//cout << index << " , " << index2 << endl;
					if(mesh_collisions(face_trees[index], face_trees[index2])) {
						cout << "Got collision on  meshes:" << i << ", " << j << endl;
						cout << "Mesh 1 has: " << meshes_[0][i]->n_faces() << " faces. Mesh 2 has: " << meshes_[0][j]->n_faces() << " faces." << endl;
						return true;
					}
				}
			}	

			return false;
		}

		bool mesh_self_collisions_wrapper() const {
			
			//First let's just see if we can resolve self-collisions
			int counter = 0;
			for(unsigned int i = 0; i < robot_number_; i++) {
				for(unsigned int j = 0; j < stateDimension/robot_number_-1; j++) {
					for(unsigned int k = j+1; k < stateDimension/robot_number_; k++) {
						if(k-j > 1) {
							int index = i * (stateDimension/robot_number_) + j;
							int index2 = i * (stateDimension/robot_number_) + k;
							if(mesh_collisions(face_trees[index], face_trees[index2])) {
								cout << "Got self collision comparing robot" << i << "And meshes:" << j << ", " << k << endl;
								cout << "Mesh 1 has: " << meshes_[i][j]->n_faces() << " faces. Mesh 2 has: " << meshes_[i][k]->n_faces() << " faces." << endl;
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

			vector<Frame<Vector<real,3>>> frames = frame_from_state(joint_angles);

			for(unsigned int i = 0; i < robot_number_; i++) {
				for(unsigned int j = 0; j < meshes_[i].size(); j++) {
					int index = i * (stateDimension/robot_number_)  + j;
					Frame<Vector<real,3>> this_frame = frames[index];
					face_trees[index]->X.const_cast_().copy(this_frame * positions[index]);
					face_trees[index]->update();
				}
			}
		}

		vector<Frame<Vector<real, 3>>> frame_from_state(const Array<real> &state_angles) const {
	
			vector<Frame<Vector<real,3>>> frames;

			for(unsigned int i = 0; i < stateDimension; ++i) {	
				int index = i%(stateDimension/robot_number_);		
				
				auto f = Frame<Vector<real,3>>(); //Let's create a frame
				int factor = axis_information_[i].negate_rotation ? -1 : 1; //For Kuka robots to move how we want them to

				if(index == 0) { //Comes back in here when it's the nTh base. 
					f.t = axis_information_[i].offsets + initial_location_[0]; //Offset of the base shouldn't ever change..
					//The axis doesn't either, but the rotation angle itself does.  We're gonna change that now. 
					auto rotation_object = Rotation<Vector<real,3>>(factor * state_angles[i], axis_information_[i].axis); 
					f.r = rotation_object;
					if(i == 6) {
						auto axis_rotation_object = Rotation<Vector<real,3>>(initial_angle_, axis_information_[i].axis);
					//	printVector(axis_information_[i].offsets);
					//	printVector(initial_location_[1]);
						frames.push_back(Frame<Vector<real, 3>>(axis_rotation_object * axis_information_[i].offsets + initial_location_[1], rotation_object * axis_rotation_object));
					} else 
						frames.push_back(f);
					
				} else {
					//If it's not the first link, then get the previous link 
					f = frames.back();
					//Your rotation axis should be previous guy's rotation applied to your axis.
					auto rotation_object = Rotation<Vector<real,3>>(factor * state_angles[i], f.r * axis_information_[i].axis);

					//Your position is the offset of the other guy + ABSOLUTE rotation * your offset. Your absolute rotation 
					//is the previous absolute rotation multiplied by your relative rotation. Note the order of matrix multiplication
					frames.push_back(Frame<Vector<real,3>>(f.t + f.r * axis_information_[i].offsets, rotation_object * f.r ));
					//for(int m = 0; m < 3; m++) {
					//cout << f.t[m]<< endl;
					//}
				}		
			}		
			return frames;
		}

	protected: 
		vector<vector<Ref<TriMesh>>> meshes_;
		vector<Ref<TriMesh>> obstacleMeshes_;
		vector<link_t> axis_information_;
		unsigned int robot_number_;
		unsigned int stateDimension;
		vector<Ref<SimplexTree<Vector<real,3>,2>>> face_trees;
		vector<Ref<SimplexTree<Vector<real,3>,2>>> obstacle_trees;
		vector<Array<Vector<real, 3>>> positions;
		vector<Array<Vector<real, 3>>> obstacle_positions;
		double initial_angle_;
		vector<Vector<real,3>> initial_location_;

};

/* Initializes the axes for each linkage (sets up frames) on which the forward kinematics
are eventually computed. These are then used for validity checking. */



void initializeAxes(ob::SpaceInformationPtr &si, vector<link_t> &axis_information, Array<Vector<real,3>,2> &offsets, 
	double robot_number) {
	
	Vector<real, 3> x_axis(1,0,0);
	Vector<real, 3> y_axis(0,1,0);
	Vector<real, 3> z_axis(0,0,1);

	for(int o = 0; o < (int) robot_number; ++o) {
		auto robot_offsets = offsets[o];
		//if((initial_angle != 0 || location.size() > 0) && o == 1) {
		//	Rotation<Vector<real, 3>> axis_rotation(initial_angle, z_axis);
		//	Frame<Vector<real, 3>> new_location(location, axis_rotation);
		//	robot_offsets = new_location * robot_offsets;	
		//xw}
		//for(int m = 0; m < 6; m++) {
		//		cout << robot_offsets[m] << endl;
	//		}	
		for(unsigned int i = 0; i < si->getStateDimension()/robot_number; ++i) {
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
}

static Nested<real> plan(unsigned int links, double robot_number, vector<Array<real>> goalState, Array<Vector<real,3>,2> parsed_offsets, 
	vector<vector<Ref<TriMesh>>> robotMeshes, vector<Ref<TriMesh>> obstacleMeshes, double resolution, double range, double solve_time, 
	double initial_angle, vector<Vector<real, 3>> initial_location) {

	//Create a subspace with the given number of links
	//Create the space information pointer that everything else takes...
	ob::StateSpacePtr space(new ChainSpace(links, robot_number)); 
	ob::SpaceInformationPtr si(new ob::SpaceInformation(space));
	
	si->setStateValidityCheckingResolution(resolution);

	//Set up the node structure from the given files. 
	vector<link_t> axis_information;
	initializeAxes(si, axis_information, parsed_offsets, robot_number);
	//We initialize as a 000000 state for the start configuration
	ob::ScopedState<ob::CompoundStateSpace> start(space); 
	ob::ScopedState<ob::CompoundStateSpace> goal(space);

	og::PathSimplifier simplifier(si);

	//Eventually this will be where I do validity checking 
	si->setStateValidityChecker(ob::StateValidityCheckerPtr(new SO26ValidityChecker(si, 2, axis_information, robotMeshes, obstacleMeshes, initial_angle, initial_location)));
   
	Nested<real,false> path;
		
	for(size_t k = 0; k < goalState.size(); k++) {
		og::RRTConnect* solver  = new og::RRTConnect(si);
   		solver->setRange(range);	
		ob::PlannerPtr planner(solver);
	//Code from the box: set up a problem, solve it. Resolution parameters are changed here. 
		ob::ProblemDefinitionPtr pdef(new ob::ProblemDefinition(si));

		Array<real> substate = goalState[k];
		ob::ScopedState<ob::CompoundStateSpace> new_goal(space);

		for(unsigned int index = 0; index < si->getStateDimension(); ++index) {
			if(substate[index] > pi) substate[index] = substate[index] - 2*pi;
			else if(substate[index] < -pi) substate[index] = substate[index] + 2*pi;

			if(k == 0) {
					start->as<ob::SO2StateSpace::StateType>(index)->value = 0;
			}

			else start->as<ob::SO2StateSpace::StateType>(index)->value = goal->as<ob::SO2StateSpace::StateType>(index)->value;
			new_goal->as<ob::SO2StateSpace::StateType>(index)->value = substate[index];
		}

		goal = new_goal;

		cout << start << endl;
		cout << goal << endl;

		pdef->setStartAndGoalStates(start, goal);
	    // create a planner for the defined space
	    planner->setProblemDefinition(pdef);
	    planner->setup();

	  // cout << "We're now going to print the space settings" << endl;
	    //si->printSettings(cout);

	   // cout << "We're now going to print the problem settings" << endl;
	   // pdef->print(cout);

	    //cout << "We've finished printing the problem settings, printing solution..." << endl;
	
	    // attempt to solve the problem within ten seconds of planning time
	    ob::PlannerStatus solved = planner->solve(solve_time);
		vector<real> step_vector;

	    if (solved)
	    {
	        ob::PathPtr solution_path = pdef->getSolutionPath();
	        og::PathGeometric geopath = dynamic_cast<og::PathGeometric &>(*solution_path);
	        if(simplifier.shortcutPath(geopath, 0, 0, .1)) {
	        	simplifier.smoothBSpline(geopath);
	        	//simplifier.simplify(geopath, 10);
	        	cout << "Simplified the path? " << endl;
	        }
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
	        exit(0);
		}
	}
//	}
	return path;
}

// void printVector(Vector<real,3> vectortoprint) {
// 	for(int i = 0; i < vectortoprint.size(); i++) {
// 		cout << vectortoprint[i] << endl;
// 	}
// }

}
}

using namespace other;

void wrap_multi_planner(){
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


