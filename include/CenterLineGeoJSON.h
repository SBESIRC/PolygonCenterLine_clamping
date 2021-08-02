#pragma once
#include "stdafx.h"
#include "valdata.h"

namespace CenterLine
{
	std::string out2str(std::vector<std::pair<Point, Point>> out);
	void parse_geojson(const std::string& datastr, parseout& out);
}
