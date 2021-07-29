#ifndef GENERATE_CENTERLINE
#define GENERATE_CENTERLINE
#include <string>

#include "PolygonCenterLine.h"

inline std::string Generate_Center_Line(std::string geojson){
    CenterLine::PolygonCenterLine solver;
    if(solver.calcCenterLine(geojson)){
        return solver.centerline_geojson();
    }
    else return "{}";
}

#endif // GENERATE_CENTERLINE
