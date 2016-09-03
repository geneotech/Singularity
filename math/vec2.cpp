#include "vec2.h"

std::vector<vec2> generate_circle_points(float radius, float last_angle_in_degrees,float starting_angle_in_degrees, unsigned int number_of_points)
{
	std::vector<vec2> result;

	float step = (last_angle_in_degrees - starting_angle_in_degrees) / number_of_points;

	for (float i = starting_angle_in_degrees; i < last_angle_in_degrees; i += step) {
		result.push_back(vec2().set_from_degrees(i).set_length(radius) );
	}

	return result;
}