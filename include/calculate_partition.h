#pragma once
#include <string>
#include "PolygonCenterLine.h"

// 弧度值差距不超过eps被认为是同一个ucs
inline std::string calculate_ucs_partition(std::string geojson, double eps=1e-2, const CenterLine::Context &context=CenterLine::Context()){
    CenterLine::PolygonCenterLine solver(context);
    if(!solver.calcCenterLine(geojson)) return "{}";
    if(!solver.calcUcsPartition(eps)) return "{}";
    return solver.ucs_parts_geojson();
}

inline std::string calculate_partition(std::string geojson, double R, const CenterLine::Context &context=CenterLine::Context()){
    CenterLine::PolygonCenterLine solver(context);
    if(!solver.calcCenterLine(geojson)) return "{}";
    if(!solver.calcPartition(R)) return "{}";
    return solver.parts_geojson();
}
