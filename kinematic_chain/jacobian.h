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
#include <limits>

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


class Jacobian { 

public: 
	Jacobian(System sys, Array<real> current_angles) {


	}

private: 
	System sys_;
	Array<real> current_angles_;

};