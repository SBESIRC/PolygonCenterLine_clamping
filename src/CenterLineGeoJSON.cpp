#include "CenterLineGeoJSON.h"
#include "json/json.h"

namespace CenterLine
{
	void parse_geojson(const std::string& datastr, parseout& out)
	{
		bool res;
		std::string errs;
		Json::Value root;
		Json::CharReaderBuilder readerBuilder;
		std::unique_ptr<Json::CharReader> const jsonReader(readerBuilder.newCharReader());

		res = jsonReader->parse(datastr.c_str(), datastr.c_str() + datastr.length(), &root, &errs);
		if (!res || !errs.empty()) {
			std::cout << "parseJson err. " << errs << std::endl;
		}
		int size = root["features"].size();
		for (auto i = 0; i < size; ++i)
		{
			auto g = root["features"][i]["geometry"];
			auto gtype = g["type"].asString();
			auto coor = g["coordinates"];
			auto p = root["features"][i]["properties"];
			auto cat = p["Category"].asString();
			auto privacy = p["Privacy"];
			auto sw = p["Switch"];
			auto n = p["Name"];
			bool flag = false;
			std::string name = cat;
			std::string ssw, sprivacy;
			if (!sw.isNull())
			{
				ssw = sw.asString();
			}
			if (!privacy.isNull())
			{
				sprivacy = privacy.asString();
			}
			if (!n.isNull())
			{
				name = n.asString();
			}
			if (gtype == "Polygon")
			{
				Block tb;
				std::vector<std::vector<Point>> tpwh;
				for (auto j = 0; j < coor.size(); ++j)
				{
					std::vector<Point> tpoly;
					Point last_p;
					for (auto k = 0; k < coor[j].size(); ++k)
					{
						Point p{ coor[j][k][0].asDouble(),coor[j][k][1].asDouble() };
						if(k && p.x == last_p.x && p.y == last_p.y) continue;
						if (p.x < out.minx) out.minx = p.x;
						if (p.x > out.maxx) out.maxx = p.x;
						if (p.y < out.miny) out.miny = p.y;
						if (p.y > out.maxy) out.maxy = p.y;
						tpoly.push_back(p);
					}
					if(tpoly.size() > 2 || j == 0) tpwh.push_back(tpoly);
				}
				tb.coords = tpwh;
				tb.category = cat;
				tb.sw = ssw;
				tb.privacy = sprivacy;
				out.data.push_back(tb);
			}
			//if (gtype == "Point")
			//{
			//	Point p{ coor[0].asDouble(),coor[1].asDouble() };
			//	if (p.x < out.minx) out.minx = p.x;
			//	if (p.x > out.maxx) out.maxx = p.x;
			//	if (p.y < out.miny) out.miny = p.y;
			//	if (p.y > out.maxy) out.maxy = p.y;
			//	FireExt fe;
			//	fe.coords = p;
			//	fe.category = name;
			//	out.dataF.push_back(fe);
			//}
		}
	}

	std::string out2str(std::vector<std::pair<Point, Point>> out)
	{
		//auto pt = out.p;
		//auto type = out.type;
		//auto oS = out.oS;
		//double sample_degree = pi / sample_num;
		Json::Value root, rf, lst;
		Json::Value oc, og, op, coords, point, pc;
		for (auto pr : out) {
			coords.clear();
			point.clear();
			point.append(pr.first.x);
			point.append(pr.first.y);
			coords.append(point);
			point.clear();
			point.append(pr.second.x);
			point.append(pr.second.y);
			coords.append(point);
			og["type"] = "LineString";
			og["coordinates"] = coords;
			oc["type"] = "Feature";
			oc["geometry"] = og;
			rf.append(oc);
		}

		root["features"] = rf;
		root["type"] = "FeatureCollection";
		return root.toStyledString();
	}

	Json::Value dump_block(const Block &block){ // polygon_with_holes: [polygon: [point: [x, y] ] ]
		Json::Value polys, poly, point;
		for(auto lst : block.coords){
			poly.clear();
			for(auto p : lst){
				point.clear();
				point.append(p.x); point.append(p.y);
				poly.append(point);
			}
			polys.append(poly);
		}
		return polys;
	}
	std::string out2str(const std::vector<Block> &rect_blocks, const std::vector<Block> &centerline_blocks) {
		//auto pt = out.p;
		//auto type = out.type;
		//auto oS = out.oS;
		//double sample_degree = pi / sample_num;
		Json::Value root, features;
		Json::Value feature, geometry, poly, coords, point, pc;
		for (const auto &block : rect_blocks) {
			geometry["type"] = "Polygon";
			//coords = dump_block(block);
			coords.clear();
			for(auto lst : block.coords){
				poly.clear();
				for(auto p : lst){
					point.clear();
					point.append(p.x); point.append(p.y);
					poly.append(point);
				}
				coords.append(poly);
			}
			geometry["coordinates"] = coords;
			feature["type"] = "Feature";
			feature["geometry"] = geometry;
			feature["properties"]["is_centerline_covered"] = false;
			features.append(feature);
		}
		for (const auto &block : centerline_blocks) {
			geometry["type"] = "Polygon";
			//coords = dump_block(block);
			coords.clear();
			for(auto lst : block.coords){
				poly.clear();
				for(auto p : lst){
					point.clear();
					point.append(p.x); point.append(p.y);
					poly.append(point);
				}
				coords.append(poly);
			}
			geometry["coordinates"] = coords;
			feature["type"] = "Feature";
			feature["geometry"] = geometry;
			feature["properties"]["is_centerline_covered"] = true;
			features.append(feature);
		}
		root["features"] = features;
		root["type"] = "FeatureCollection";
		return root.toStyledString();
	}
}// namespace CenterLine
