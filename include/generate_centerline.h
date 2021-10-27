#pragma once
#include <string>
#include "PolygonCenterLine.h"

inline std::string Generate_Center_Line(std::string geojson, const CenterLine::Context &context){
    CenterLine::PolygonCenterLine solver(context);
    if (!solver.calcCenterLine(geojson)) return "{}";
    return solver.centerline_geojson();
}
