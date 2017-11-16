//
//  geometry.cpp
//  ms2_anotator
//
//  Created by Aaron Maurais on 11/14/17.
//  Copyright © 2017 Aaron Maurais. All rights reserved.
//

#include "geometry.hpp"

inline bool geometry::valueInRange(double value, double min, double max){
	return (value >= min) && (value <= max);
}

inline double geometry::dist(const geometry::Point& p1, const geometry::Point& p2){
	return sqrt(pow((p2.x - p1.x), 2) + pow((p2.y - p1.y), 2));
}

void geometry::Point::move(const geometry::Vector2D& mv){
	x += mv.getH(); y += mv.getV();
}

bool geometry::Rect::intersects(const geometry::Rect& comp) const
{
	bool xOverlap = valueInRange(x, comp.x, comp.x + comp.width) ||
		valueInRange(comp.x, x, x + width);
	
	bool yOverlap = valueInRange(y, comp.y, comp.y + comp.height) ||
		valueInRange(comp.y, y, y + height);
	
	return xOverlap && yOverlap;
}
