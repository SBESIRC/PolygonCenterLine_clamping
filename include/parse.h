#ifndef PARSE_H
#define PARSE_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cassert>
#include <regex>
#include <vector>

namespace Parse{
    using namespace std;
struct Point{
    double x;
    double y;
};

struct Block{
    string type;
    string cate;
    string name;
    vector<vector<Point> > data;
};

static void splits(const string& s, vector<string>& tokens, const string& delimiters = " ") {
	string::size_type lastPos = s.find(delimiters, 0);
	if (lastPos == string::npos)
		return;
	lastPos += delimiters.length();
	string::size_type pos = s.find(delimiters, lastPos);
	while (tokens.push_back(s.substr(lastPos, pos - lastPos)), pos != string::npos) {
		lastPos = pos + delimiters.length();
		pos = s.find(delimiters, lastPos);
	}
}

static vector<Block> parse_geojson(string filename) {
	ifstream f(filename);
    
	stringstream ss;
	ss << f.rdbuf();
	string datastr = ss.str();
    vector<Block> blocks;

	vector<string> s;
	splits(datastr, s, "{\n");
	assert(s.size() % 3 == 0);

	for (int i = 0; i < s.size() / 3; ++i) {
		string type, cate, name = "";
		string s1 = s[3 * i + 1], s2 = s[3 * i + 2];
		vector<string> vs1, vs2;
		vector<vector<Point> > coords;
		splits(s1, vs1, "\"");
		type = vs1[2];
		auto start = s1.find_first_of('[');
		auto end = s1.find_last_of(']');
		string scoords = s1.substr(start, end - start + 1);
		regex re("(-[0-9]+(.[0-9]+)?)|([0-9]+(.[0-9]+)?)");
		sregex_iterator pos(scoords.begin(), scoords.end(), re);
		sregex_iterator end_it;
		coords.push_back(vector<Point>());
		if (type == "LineString") {
			for (; pos != end_it; ++pos) {
				double x = stod(pos->str(0));
				double y = stod((++pos)->str(0));
				Point p{x, y};
				coords.back().push_back(p);
			}
		}
		else if (type == "Polygon") {
			for (; pos != end_it; ++pos) {
				if (coords.back().size() > 1 && 
                    coords.back().front().x == coords.back().back().x &&
                    coords.back().front().y == coords.back().back().y)
					coords.push_back(vector<Point>());
				double x = stod(pos->str(0));
				double y = stod((++pos)->str(0));
				Point p{x, y};
				coords.back().push_back(p);
			}
		}

		splits(s2, vs2, "\"");
		cate = vs2[2];	//Catagory

		if (cate == "DrainageFacility")
			name = vs2[6];

		if (cate == "Space") {
			if (vs2[6] == "")	name = "其他区域";
			else				name = vs2[6];
		}

        Block b{type, cate, name, coords};
        blocks.push_back(b);

        cout << type << "\t" << cate << "\t" << name << endl;
		if (cate == "Space") {
			if (name == "其他区域"){
                //其他区域
            }
			else if (name == "停车区域"){
                //停车区域
            }
			else{
                //生活泵房等区域
            }
		}
		else if (cate == "DrainageFacility") {
			if (name == "排水沟"){
                //排水设施之排水沟
            }
			else{
                //排水设施之集水坑与地漏
            }
		}
		else if (cate == "ShearWall"){
            //墙
        }
        else if(cate == "Column") {
			//柱
		}
		else if (cate == "Obstacle") {
			//障碍
		}
	}

	return blocks;
}

} // namespace std
#endif