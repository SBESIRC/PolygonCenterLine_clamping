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
    double R;
    //for(int i = 173;i <= 210;++i){
    //for(int i = 12;i <= 20;++i){
    if(argc == 4){
        int from, to;
        sscanf(argv[2], "%d", &from);
        sscanf(argv[3], "%d", &to);
        for(int i = from;i <= to;++i){
            sprintf(filename, argv[1], i);
            CenterLine::Context context;
            context.prefer_ortho = true;
            PolygonCenterLineTest centerline(filename, i, context);
            centerline.showCenterLine();
            //centerline.showUcsPartition();
            //centerline.showPartition(0);
        }
    }
    else { // if(argc == 3) {
        sscanf(argv[2], "%lf", &R);
        CenterLine::Context context;
        context.prefer_ortho = true;
        PolygonCenterLineTest centerline(argv[1], 0, context);
        centerline.showCenterLine();
        centerline.showUcsPartition();
        //centerline.showPartition(R);
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
