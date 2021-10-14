#pragma once
#include <string>
#include "PolygonCenterLine.h"

inline std::string Generate_Center_Line(std::string geojson){
    CenterLine::PolygonCenterLine solver;
    if (!solver.calcCenterLine(geojson)) return "{}";
    return solver.centerline_geojson();
}
