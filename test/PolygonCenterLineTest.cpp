#include <iostream>

#define PolygonCenterLine_Implementation
#include "PolygonCenterLine.h"
#include "PolygonCenterLineTest.h"

using namespace std;

int main(int argc, char *argv[])
{
    cout << "argc" << argc << endl;
    if (argc < 2) {
        cout << argc << endl;
        return 1;
    }
    //CenterLine::PolygonCenterLine centerline(argv[1]);
    PolygonCenterLineTest centerline(argv[1]);
    //double interval = 1;
    //if(argc >= 3){
    //	stringstream s_in(argv[2]);
    //	s_in >> interval;
    //}
    //cout << "interval = " << interval << endl;
    centerline.showCenterLine();
    //centerline.showCenterLine(interval);
    cout << "Success." << endl;
    return 0;
}
