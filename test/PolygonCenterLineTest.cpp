#include <iostream>
#include <string>

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
    char filename[1024];
    //for(int i = 173;i <= 210;++i){
    //for(int i = 12;i <= 20;++i){
    for(int i = 1;i <= 212;++i){
        sprintf(filename, argv[1], i);
        PolygonCenterLineTest centerline(filename, i);
        centerline.showCenterLine();
        //centerline.showPartition(0);
    }
    //double interval = 1;
    //if(argc >= 3){
    //	stringstream s_in(argv[2]);
    //	s_in >> interval;
    //}
    //cout << "interval = " << interval << endl;
    //centerline.showCenterLine(interval);
    cout << "Success." << endl;
    return 0;
}
