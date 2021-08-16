#pragma once
#include <string>
#include "PolygonCenterLine.h"

inline std::string calculate_partition(std::string geojson, double R){
    CenterLine::PolygonCenterLine solver;
    if(!solver.calcCenterLine(geojson)) return "{}";
	if(!solver.calcPartition(R)) return "{}";
	return solver.parts_geojson();
}
