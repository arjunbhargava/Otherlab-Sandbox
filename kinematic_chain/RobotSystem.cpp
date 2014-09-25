/* 
This class handles the robot meshes for bike welding. The input parameters are specified in the constructor. 
This makes it easier to handle/update the meshes/get the positions of the end effectors. 
Creting 2 robot meshes can be done much more easily, as the state of the robots is rememebred in this class. 
*/

#include "RobotSystem.h"

//Constructor: Intialize the mesh for the robot, and remember if there are any associated 

other::RobotSystem::RobotSystem(unsigned int links, Array<Vector<real,3>> parsed_offsets, vector<Ref<TriMesh>> robotMesh, 
		vector<Ref<TriMesh>> obstacleMeshes, double initial_angle, Vector<real, 3> initial_location, double effector_offset) 
{	
	sys.robotMesh = robotMesh;
	sys.obstacleMeshes = obstacleMeshes;  
	sys.initial_angle = initial_angle;
	sys.initial_location = initial_location;
	sys.stateDimension = links;
	sys.effector_offset = effector_offset;

	other::RobotSystem::initializeAxes(parsed_offsets);

	for(unsigned int i = 0; i < robotMesh.size(); i++) {
		sys.face_trees.push_back(robotMesh[i]->face_tree());
		sys.positions.push_back(sys.face_trees[i]->X.copy());
	}
}

//Calculate the frame from any given set of angles.

vector<Frame<Vector<real, 3>>> other::RobotSystem::frame_from_state(Array<real> state_angles) 
{
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


//Initializes the axes, as in the constructor. Only useful once, this is what really differentiates the 2 meshes. 
void other::RobotSystem::initializeAxes(Array<Vector<real,3>> offsets) 
{
	vector<link_t> axis_information;
	Vector<real, 3> x_axis(1,0,0);
	Vector<real, 3> y_axis(0,1,0);
	Vector<real, 3> z_axis(0,0,1);

	double upper_kr16_bounds [] = {185, 125, 64, 165, 130, 350};
	double lower_kr16_bounds [] = {-185, -65, -210, -165, -130, -350};

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


//Updates the mesh to the given state by updating the face tree associated with that robot. 
void other::RobotSystem::update_mesh(Array<real> joint_angles) 
{
	vector<Frame<Vector<real,3>>> frames = frame_from_state(joint_angles);
	for(unsigned int i = 0; i < sys.robotMesh.size(); i++) {
		Frame<Vector<real,3>> this_frame = frames[i];
		sys.face_trees[i]->X.const_cast_().copy(this_frame * sys.positions[i]);
		sys.face_trees[i]->update();
	}
}

//Finds the effector position for the given robot. 
Vector<real,3> other::RobotSystem::effectorPositions(Array<real> joint_angles)
{
	vector<Frame<Vector<real,3>>> link_frames =	frame_from_state(joint_angles);
	Vector<real,3> end_position;
	int index = sys.stateDimension - 1;
	Frame<Vector<real,3>> effector = link_frames[index];
	end_position = effector.t + effector.r * (sys.effector_offset * sys.axis_information[index].axis);
	return end_position;
}

