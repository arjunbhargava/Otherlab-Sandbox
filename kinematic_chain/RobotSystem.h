#pragma once
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


#include <math.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <limits>

namespace ob = ompl::base;
namespace og = ompl::geometric;
using namespace geode;
using std::vector;
using std::cout;
using std::endl;
using std::cin;
using std::string;

namespace other {



/*
This class defines the meshes for the robot arms. Useful helper functions are defined in here as well. 
*/
class RobotSystem {


public:

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

	RobotSystem(unsigned int links, Array<Vector<real,3>> parsed_offsets, vector<Ref<TriMesh>> robotMesh, 
		vector<Ref<TriMesh>> obstacleMeshes, double initial_angle, Vector<real, 3> initial_location, double effector_offset);

	vector<Frame<Vector<real, 3>>> frame_from_state(Array<real> state_angles);


	void update_mesh(Array<real> joint_angles);

	vector<Ref<SimplexTree<Vector<real,3>,2>>> getFaceTree() 
	{
		return sys.face_trees;
	}

	vector<Ref<TriMesh>> getObstacleMeshes() {
		return sys.obstacleMeshes;
	}

	unsigned int getStateDimension() 
	{
		return sys.stateDimension;
	}

	Vector<real,3> getTargetPosition() {
		return sys.target_position;
	}

	double getEffectorOffset() {
		return sys.effector_offset;
	}

	vector<link_t> getAxisInformation() {
		return sys.axis_information;
	}



	void setTargetPosition(Vector<real,3> target) 
	{
		sys.target_position = target;
	}

	Vector<real,3> effectorPositions(Array<real> joint_angles);


private:

	void initializeAxes(Array<Vector<real,3>> parsed_offsets);
	System sys;

};
}



