#pragma once
#include "stdafx.h"
#include "valdata.h"

namespace CenterLine
{
	void parse_geojson(const std::string& datastr, parseout& out);
	std::string out2str(std::vector<std::pair<Point, Point>> out);
	std::string out2str(const std::vector<Block> &ucs_blocks, const std::vector<Point> &corr_ucs);
	std::string out2str(const std::vector<Block> &rect_blocks, const std::vector<Block> &centerline_blocks);
}
