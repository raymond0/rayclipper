//
//  rayclipper.hpp
//  navit
//
//  Created by Ray Hunter on 18/11/2016.
//
//

#ifndef rayclipper_hpp
#define rayclipper_hpp

#include <stdio.h>
#include <vector>
    
namespace rayclipper {
    
typedef struct coord {
    int x;
    int y;
} coord;

struct rect {struct coord l,h;};

typedef std::vector<struct coord> Contour;
class Polygon : public Contour
{
public:
    std::vector<Polygon> holes;
    Polygon() {}
    Polygon( const Contour &contour ) { insert( end(), contour.begin(), contour.end() ); }
};

void RayClipPolygon( const Polygon &inputPolygon, struct rect rect, std::vector<Polygon> &outputPolygons );
void CleanPolygon( const Polygon &inputPolygon, Polygon &outputPolygon );
long long PolygonArea(const Polygon &coords);
}


#endif /* rayclipper_hpp */
