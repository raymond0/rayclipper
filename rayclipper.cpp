//
//  rayclipper.cpp
//  navit
//
//  Created by Ray Hunter on 18/11/2016.
//
//

#include "rayclipper.h"
#include <math.h>
#include <assert.h>

#define coord_is_equal(a,b) ((a).x==(b).x && (a).y==(b).y)
#define sq(x) ((double)(x)*(x))


using namespace std;

namespace rayclipper
{
    
int IndexOfPreceedingPointOnEdge( const Polygon &polygon, const struct coord coordinate, struct rect rect );
    
    
typedef enum
{
    EdgeTop = 0,
    EdgeRight = 1,
    EdgeBottom = 2,
    EdgeLeft = 3,
    EdgeNotAnEdge = 4
} EdgeType;
    
    
double CoordDistance(const struct coord &from, const struct coord &to)
{
    double dx = to.x - from.x;
    double dy = to.y - from.y;
    
    return sqrt( ( dx * dx ) + ( dy * dy ) );
}

    
bool PointIsCompletelyWithinRect(const struct coord &coord, const struct rect &r)
{
    if ( coord.x <= r.l.x ) return false;
    if ( coord.x >= r.h.x ) return false;
    if ( coord.y <= r.l.y ) return false;
    if ( coord.y >= r.h.y ) return false;
    
    return true;
}
    
    
bool PointIsInsideRect(const struct coord &coord, const struct rect &r)
{
    if ( coord.x < r.l.x ) return false;
    if ( coord.x > r.h.x ) return false;
    if ( coord.y < r.l.y ) return false;
    if ( coord.y > r.h.y ) return false;
    
    return true;
}

    
void LineRectIntersection( const struct coord &p1, const struct coord &p2, const struct rect &r,
                           const EdgeType edge, struct coord &result)
{
    double dx=p2.x-p1.x;
    double dy=p2.y-p1.y;
    switch(edge) {
        case EdgeLeft:
            result.y=p1.y + ((double)(r.l.x-p1.x)) * dy/dx;
            result.x=r.l.x;
            break;
        case EdgeRight:
            result.y=p1.y + ((double)(r.h.x-p1.x)) * dy/dx;
            result.x=r.h.x;
            break;
        case EdgeTop:
            result.x=p1.x + ((double)(r.h.y-p1.y)) * dx/dy;
            result.y=r.h.y;
            break;
        case EdgeBottom:
            result.x=p1.x + ((double)(r.l.y-p1.y)) * dx/dy;
            result.y=r.l.y;
            break;
        case EdgeNotAnEdge:
            assert( false );
    }
}

    
struct coord InterestionOfRect( const struct coord &p1, const struct coord &p2, const struct rect &r )
{
    assert( PointIsInsideRect(p1, r) != PointIsInsideRect(p2, r) );
    
    struct coord ret;
    
    bool leftPossible =   ( p1.x < r.l.x && p2.x >= r.l.x ) || ( p2.x < r.l.x && p1.x >= r.l.x );
    bool rightPossible =  ( p1.x <= r.h.x && p2.x > r.h.x ) || ( p2.x <= r.h.x && p1.x > r.h.x );
    bool bottomPossible = ( p1.y < r.l.y && p2.y >= r.l.y ) || ( p2.y < r.l.y && p1.y >= r.l.y );
    bool topPossible =    ( p1.y <= r.h.y && p2.y > r.h.y ) || ( p2.y <= r.h.y && p1.y > r.h.y );
    
    // Left + Right
    if ( p1.x != p2.x )
    {
        if ( leftPossible )
        {
            LineRectIntersection(p1, p2, r, EdgeLeft, ret);
            if ( r.l.y <= ret.y && ret.y <= r.h.y ) return ret;
        }
        
        if ( rightPossible )
        {
            LineRectIntersection(p1, p2, r, EdgeRight, ret);
            if ( r.l.y <= ret.y && ret.y <= r.h.y ) return ret;
        }
    }
    
    // Top + bottom
    if ( p1.y != p2.y )
    {
        if ( topPossible )
        {
            LineRectIntersection(p1, p2, r, EdgeTop, ret);
            if ( r.l.x <= ret.x && ret.x <= r.h.x ) return ret;
        }
        
        if ( bottomPossible )
        {
            LineRectIntersection(p1, p2, r, EdgeBottom, ret);
            if ( r.l.x <= ret.x && ret.x <= r.h.x ) return ret;
        }
    }
    
    assert( false );
    struct coord failed;
    return failed;  // Some compiler configs shout error otherwise
}
    
    
vector<struct coord> AllIntersectionsOfRect( const struct coord &p1, const struct coord &p2, const struct rect &r )
{
    bool leftPossible =   ( p1.x < r.l.x && p2.x >= r.l.x ) || ( p2.x < r.l.x && p1.x >= r.l.x );
    bool rightPossible =  ( p1.x <= r.h.x && p2.x > r.h.x ) || ( p2.x <= r.h.x && p1.x > r.h.x );
    bool bottomPossible = ( p1.y < r.l.y && p2.y >= r.l.y ) || ( p2.y < r.l.y && p1.y >= r.l.y );
    bool topPossible =    ( p1.y <= r.h.y && p2.y > r.h.y ) || ( p2.y <= r.h.y && p1.y > r.h.y );
    
    vector<struct coord> intersections;
    struct coord ret;
    
    // Left + Right
    if ( p1.x != p2.x )
    {
        if ( leftPossible )
        {
            LineRectIntersection(p1, p2, r, EdgeLeft, ret);
            if ( r.l.y <= ret.y && ret.y <= r.h.y )
            {
                intersections.emplace_back( ret );
            }
        }
        
        if ( rightPossible )
        {
            LineRectIntersection(p1, p2, r, EdgeRight, ret);
            if ( r.l.y <= ret.y && ret.y <= r.h.y )
            {
                intersections.emplace_back( ret );
            }
        }
    }
    
    // Top + bottom
    if ( p1.y != p2.y )
    {
        if ( topPossible )
        {
            LineRectIntersection(p1, p2, r, EdgeTop, ret);
            if ( r.l.x <= ret.x && ret.x <= r.h.x )
            {
                intersections.emplace_back( ret );
            }
        }
        
        if ( bottomPossible )
        {
            LineRectIntersection(p1, p2, r, EdgeBottom, ret);
            if ( r.l.x <= ret.x && ret.x <= r.h.x )
            {
                intersections.emplace_back( ret );
            }
        }
    }
    
    assert( intersections.size() <= 2 );
    
    return intersections;
}
    
    
//
//  LineClippedToRect - Clip a line to a rect when both points are outside
//
vector<struct coord> LineIntersectionOfRect( const struct coord &p1, const struct coord &p2, const struct rect &rect )
{
    vector<struct coord> intersections = AllIntersectionsOfRect(p1, p2, rect);
    
    vector<struct coord> output;
    
    if ( intersections.size() != 2 )
    {
        return output;
    }
    
    assert ( intersections[0].y == rect.l.y || intersections[0].y == rect.h.y ||
             intersections[0].x == rect.l.x || intersections[0].x == rect.h.x );
    assert ( intersections[1].y == rect.l.y || intersections[1].y == rect.h.y ||
             intersections[1].x == rect.l.x || intersections[1].x == rect.h.x );
    
    if ( coord_is_equal(intersections[0], intersections[1]) )
    {
        assert( ( intersections[0].x == rect.l.x || intersections[0].x == rect.h.x ) &&
                ( intersections[0].y == rect.l.y || intersections[0].y == rect.h.y ) );

        return output;
    }
    
    double d0 = CoordDistance(p1, intersections[0]);
    double d1 = CoordDistance(p1, intersections[1]);
    
    if ( d0 <= d1 )
    {
        output.emplace_back(intersections[0]);
        output.emplace_back(intersections[1]);
    }
    else
    {
        output.emplace_back(intersections[1]);
        output.emplace_back(intersections[0]);
    }
    
    return output;
}
    

//
//  LineClippedToRect - Clip a line to a rect under all circumstances
//
std::vector<struct coord> LineClippedToRect( const struct coord &p1, const struct coord &p2, const struct rect &rect )
{
    std::vector<struct coord> output;
    
    bool p1InRect = rect.l.x <= p1.x && p1.x <= rect.h.x && rect.l.y <= p1.y && p1.y <= rect.h.y;
    bool p2InRect = rect.l.x <= p2.x && p2.x <= rect.h.x && rect.l.y <= p2.y && p2.y <= rect.h.y;
    
    if ( p1InRect && p2InRect )
    {
        output.emplace_back(p1);
        output.emplace_back(p2);
        return output;
    }
    
    if ( !p1InRect && !p2InRect )
    {
        return LineIntersectionOfRect(p1, p2, rect);
    }
    
    if ( p1InRect )
    {
        vector<struct coord> intersections = AllIntersectionsOfRect(p1, p2, rect);
        assert( intersections.size() == 1 );
        output.emplace_back(p1);
        output.emplace_back(intersections[0]);
        return output;
    }
    
    assert( p2InRect );
    vector<struct coord> intersections = AllIntersectionsOfRect(p1, p2, rect);
    assert( intersections.size() == 1 );
    output.emplace_back(intersections[0]);
    output.emplace_back(p2);
    return output;
}
    
    
bool LineIntersetsRect(const struct coord &p1, const struct coord &p2, const struct rect r)
{
    if ( p1.x < r.l.x && p2.x < r.l.x ) return false;
    if ( p1.x > r.h.x && p2.x > r.h.x ) return false;
    if ( p1.y < r.l.y && p2.y < r.l.y ) return false;
    if ( p1.y > r.h.y && p2.y > r.h.y ) return false;
    
    struct coord ret;
    
    // Left + right
    if ( p1.x != p2.x )
    {
        LineRectIntersection(p1, p2, r, EdgeLeft, ret);
        if ( r.l.y <= ret.y && ret.y <= r.h.y ) return true;
        LineRectIntersection(p1, p2, r, EdgeRight, ret);
        if ( r.l.y <= ret.y && ret.y <= r.h.y ) return true;
    }
    
    // Top + bottom
    if ( p1.y != p2.y )
    {
        LineRectIntersection(p1, p2, r, EdgeTop, ret);
        if ( r.l.x <= ret.x && ret.x <= r.h.x ) return true;
        LineRectIntersection(p1, p2, r, EdgeBottom, ret);
        if ( r.l.x <= ret.x && ret.x <= r.h.x ) return true;
    }
    
    return false;
}
    
    
bool LineIntersetsLine(const struct coord &line1Start, const struct coord &line1End,
                       const struct coord &line2Start, const struct coord &line2End,
                       struct coord &intersection)
{
    struct coord line1Min = { min(line1Start.x, line1End.x), min( line1Start.y, line1End.y ) };
    struct coord line1Max = { max(line1Start.x, line1End.x), max( line1Start.y, line1End.y ) };
    struct coord line2Min = { min(line2Start.x, line2End.x), min( line2Start.y, line2End.y ) };
    struct coord line2Max = { max(line2Start.x, line2End.x), max( line2Start.y, line2End.y ) };

    if ( line1Max.x <= line2Min.x ) return false;
    if ( line1Min.x >= line2Max.x ) return false;
    if ( line1Max.y <= line2Min.y ) return false;
    if ( line1Min.y >= line2Max.y ) return false;
    
    double x1 = line1Start.x, x2 = line1End.x, x3 = line2Start.x, x4 = line2End.x;
    double y1 = line1Start.y, y2 = line1End.y, y3 = line2Start.y, y4 = line2End.y;
    
    double d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    // If d is zero, there is no intersection
    if (d <= __DBL_EPSILON__) return false;
    
    // Get the x and y
    double pre = (x1*y2 - y1*x2), post = (x3*y4 - y3*x4);
    double x = ( pre * (x3 - x4) - (x1 - x2) * post ) / d;
    double y = ( pre * (y3 - y4) - (y1 - y2) * post ) / d;
    
    // Check if the x and y coordinates are within both lines
    if ( x < min(x1, x2) - __DBL_EPSILON__ || x > max(x1, x2) + __DBL_EPSILON__ ||
         x < min(x3, x4) - __DBL_EPSILON__ || x > max(x3, x4) + __DBL_EPSILON__)
        return false;
    if ( y < min(y1, y2) - __DBL_EPSILON__ || y > max(y1, y2) + __DBL_EPSILON__ ||
         y < min(y3, y4) - __DBL_EPSILON__ || y > max(y3, y4) + __DBL_EPSILON__ )
        return false;
    
    intersection.x = x;
    intersection.y = y;

    return true;
}
    
    
bool PointIsInsideContour(const Contour &coords, const struct coord &point)
{
    bool inside = false;
    
    for (size_t i = 0; i < coords.size() - 1; i++ )
    {
        if ( (coords[i].y > point.y ) != ( coords[i+1].y > point.y ) &&
            point.x < ( (long long) coords[i+1].x - coords[i].x ) * ( point.y - coords[i].y ) / ( coords[i+1].y - coords[i].y ) + coords[i].x )
        {
            inside = ! inside;
        }
    }
    
    if ( coords.front().x != coords.back().x || coords.front().y != coords.back().y )
    {
        if ( (coords.back().y > point.y ) != ( coords.front().y > point.y ) &&
            point.x < ( (long long) coords.front().x - coords.back().x ) * ( point.y - coords.back().y ) / ( coords.front().y - coords.back().y ) + coords.back().x )
        {
            inside = ! inside;
        }
    }
    
    return inside;
}


void GetAllSubContours( const Contour &inputPolygon, struct rect rect, size_t lastPointOutside, vector<Contour> &allContours )
{
    bool inside = false;
    shared_ptr<Polygon> currentPolygon;
    
    size_t inputSize = inputPolygon.size();
    
    for ( size_t i = 0; i < inputPolygon.size(); i++ )
    {
        struct coord firstPoint = inputPolygon[(i + lastPointOutside) % inputSize];
        struct coord secondPoint = inputPolygon[(i + 1 + lastPointOutside) % inputSize];

        if ( PointIsInsideRect(secondPoint, rect ) )
        {
            if ( inside )
            {
                if ( currentPolygon->size() == 0 || !coord_is_equal( currentPolygon->back(), secondPoint ) )
                {
                    currentPolygon->emplace_back(secondPoint);
                }
            }
            else
            {
                currentPolygon = shared_ptr<Polygon>(new Polygon());
                inside = true;
                struct coord intersection = InterestionOfRect(firstPoint, secondPoint, rect);
                
                assert ( intersection.y == rect.l.y || intersection.y == rect.h.y ||
                         intersection.x == rect.l.x || intersection.x == rect.h.x );
                
                if ( currentPolygon->size() == 0 || !coord_is_equal( currentPolygon->back(), intersection ) )
                {
                    currentPolygon->emplace_back(intersection);
                }
                
                if ( currentPolygon->size() == 0 || !coord_is_equal( currentPolygon->back(), secondPoint ) )
                {
                    currentPolygon->emplace_back(secondPoint);
                }
            }
        }
        else
        {
            if ( inside )
            {
                struct coord intersection = InterestionOfRect(firstPoint, secondPoint, rect);
                
                assert ( intersection.y == rect.l.y || intersection.y == rect.h.y ||
                         intersection.x == rect.l.x || intersection.x == rect.h.x );
                
                if ( currentPolygon->size() == 0 || !coord_is_equal( currentPolygon->back(), intersection ) )
                {
                    currentPolygon->emplace_back( intersection );
                }
                
                allContours.emplace_back( *currentPolygon );
                
                inside = false;
            }
            else
            {
                //
                // Both outside and currently outside - could be intersection
                //
                vector<struct coord> intersections = LineIntersectionOfRect(firstPoint, secondPoint, rect);
            
                if ( intersections.size() == 2 )
                {
                    currentPolygon = shared_ptr<Polygon>(new Polygon());
                    currentPolygon->emplace_back(intersections[0]);
                    currentPolygon->emplace_back(intersections[1]);
                    allContours.emplace_back( *currentPolygon );
                }
            }
        }
    }
    
    assert( ! inside );
}


EdgeType EdgeForCoord( struct coord coord, struct rect rect )
{
    struct coord topLeft = {rect.l.x, rect.h.y};
    struct coord bottomRight = {rect.h.x, rect.l.y};
    
    if ( coord_is_equal(coord, rect.l) )
    {
        return EdgeLeft;
    }
    
    if ( coord_is_equal(coord, topLeft) )
    {
        return EdgeTop;
    }
    
    if ( coord_is_equal(coord, rect.h) )
    {
        return EdgeRight;
    }
    
    if ( coord_is_equal(coord, bottomRight) )
    {
        return EdgeBottom;
    }
    
    if ( coord.y == rect.h.y ) return EdgeTop;
    if ( coord.x == rect.h.x ) return EdgeRight;
    if ( coord.y == rect.l.y ) return EdgeBottom;
    if ( coord.x == rect.l.x ) return EdgeLeft;
    
    return EdgeNotAnEdge;
}


int DistanceAlongEdge( EdgeType edge, struct coord from, struct coord to )
{
    switch ( edge )
    {
        case EdgeTop:
            assert( from.y == to.y );
            return to.x - from.x;
        case EdgeRight:
            assert( from.x == to.x );
            return from.y - to.y;
        case EdgeBottom:
            assert( from.y == to.y );
            return from.x - to.x;
        case EdgeLeft:
            assert( from.x == to.x );
            return to.y - from.y;
        case EdgeNotAnEdge:
            assert( false );
            return -1;
    }
}


struct coord NextCorner( struct coord coord, struct rect rect )
{
    EdgeType edge = EdgeForCoord( coord, rect );
    
    struct coord topLeft = {rect.l.x, rect.h.y};
    struct coord bottomRight = {rect.h.x, rect.l.y};
    
    switch( edge )
    {
        case EdgeTop:
            return rect.h;
        case EdgeRight:
            return bottomRight;
        case EdgeBottom:
            return rect.l;
        case EdgeLeft:
            return topLeft;
        case EdgeNotAnEdge:
            assert( false );
            return rect.l;
    }
    
    assert(false);
}


int DistanceToNextCorner( struct coord from, struct rect rect, struct coord *nextCorner )
{
    EdgeType fromEdge = EdgeForCoord( from, rect );
    struct coord corner = NextCorner( from, rect );
    
    if ( nextCorner != NULL )
    {
        *nextCorner = corner;
    }
    
    return DistanceAlongEdge( fromEdge, from, corner );
}


int ClockwiseDistance( struct coord from, struct coord to, struct rect rect )
{
    int distance = 0;
    
    EdgeType fromEdge = EdgeForCoord( from, rect );
    EdgeType toEdge = EdgeForCoord( to, rect );
    
    if ( fromEdge == toEdge )
    {
        int d = DistanceAlongEdge( fromEdge, from, to );
        if ( d >= 0 )
        {
            return d;
        }
    }
    
    struct coord currentPosition = from;
    EdgeType currentEdge = fromEdge;
    
    while ( currentEdge != toEdge || distance == 0 )  // Check - 0 check forces us all the way round. Is 0 vaild?
    {
        distance += DistanceToNextCorner( currentPosition, rect, &currentPosition );
        currentEdge = (EdgeType) ((currentEdge + 1) % 4);
    }
    
    assert( currentEdge == toEdge );
    distance += DistanceAlongEdge(currentEdge, currentPosition, to );
    
    return distance;
}


void ClosePolygon( Polygon &polygon, vector<Contour> &otherContours, struct rect rect )
{
    while ( true )
    {
        size_t closestPolygon = -1;
        int closestPolygonDistance = INT_MAX;
        
        for ( size_t i = 0; i < otherContours.size(); i++ )
        {
            auto &other = otherContours[i];
            
            int d = ClockwiseDistance(polygon.back(), other.front(), rect );
            if ( d < closestPolygonDistance )
            {
                closestPolygonDistance = d;
                closestPolygon = i;
            }
        }
        
        struct coord nextCorner;
        int distanceToNextCorner = DistanceToNextCorner(polygon.back(), rect, &nextCorner);
        int distanceToEndOfCurrent = ClockwiseDistance(polygon.back(), polygon.front(), rect);

        if ( distanceToEndOfCurrent <= distanceToNextCorner && distanceToEndOfCurrent <= closestPolygonDistance )
        {
            //
            //  No corners or other polygons in the way. We're done.
            //
            
            //
            //  Trim adjacent retracing of points at beginning/end of polygon
            //
            coord back = polygon.back();
            coord front = polygon.front();
            
            if ( polygon.size() < 4 || ! coord_is_equal(back, front) )
            {
                return;
            }

            coord second = polygon[1];
            coord penultimate = polygon[ polygon.size() - 2 ];
            
            //
            //  Prevent adjacent retracing at join
            //
            if ( penultimate.x == back.x && back.x == second.x )     // All on same x
            {
                polygon.pop_back();
                polygon.erase(polygon.begin());
            }
            else if ( penultimate.y == back.y && back.y == second.y )     // All on same x
            {
                polygon.pop_back();
                polygon.erase(polygon.begin());
            }

            return;
        }
        
        if ( closestPolygonDistance <= distanceToNextCorner )
        {
            auto &other = otherContours[closestPolygon];
            coord back = polygon.back();
            coord newFront = other.front();
            
            size_t toSkipAtStart = 0;
            
            if ( coord_is_equal( back, newFront ) )
            {
                toSkipAtStart = 1;
                
                coord penultimate = polygon[ polygon.size() - 2 ];
                coord newSecond = other[1];
                
                if ( penultimate.x == back.x && back.x == newSecond.x )     // All on same x
                {
                    polygon.pop_back();
                }
                else if ( penultimate.y == back.y && back.y == newSecond.y )     // All on same x
                {
                    polygon.pop_back();
                }
            }
            
            
            assert( other.begin() + toSkipAtStart <= other.end() );

            polygon.insert( polygon.end(), other.begin() + toSkipAtStart, other.end() );
            
            otherContours.erase( otherContours.begin() + closestPolygon );
            
            continue;
        }
        
        assert( distanceToNextCorner < distanceToEndOfCurrent && distanceToNextCorner < closestPolygonDistance );
        polygon.emplace_back( nextCorner );
    }
}


void ClosePolygons( vector<Contour> &contours, struct rect rect, vector<Polygon> &closedPolygons )
{
    while ( contours.size() > 0 )
    {
        auto contour = contours.back();
        contours.pop_back();
        Polygon polygon;
        polygon.insert( polygon.end(), contour.begin(), contour.end() );
        ClosePolygon(polygon, contours, rect);
        
        while ( polygon.size() > 1 && coord_is_equal(polygon.front(), polygon.back() ) )
        {
            //
            //  Could be out -> in -> out on the same intersection
            //
            polygon.pop_back();
        }
        
        if ( polygon.size() < 3 )
        {
            //
            // We need polygons so actually be polygon's completed by this stage
            //
            continue;
        }
        
        closedPolygons.emplace_back( polygon );
    }
}
    
    
bool PointIsOnEdge( const struct coord &coord, const struct rect &rect )
{
    return coord.x == rect.l.x || coord.x == rect.h.x ||
           coord.y == rect.l.y || coord.y == rect.h.y;
}
    
    
vector<Polygon> SplitEdgeTouchingPolygons( vector<Polygon> &polygons, struct rect rect )
{
    vector<Polygon> output;
    for ( auto &polygon : polygons )
    {
        assert( PointIsOnEdge ( polygon.front(), rect ) );
        assert( PointIsOnEdge ( polygon.back(), rect ) );
        
        shared_ptr<Polygon> currentPolygon;
        bool inPolygon = false;
        bool anyPointInside = false;

        for ( size_t i = 1; i < polygon.size(); i++ )
        {
            coord first = polygon[ i - 1 ];
            coord second = polygon[ i ];
            
            if ( PointIsCompletelyWithinRect( first, rect ) ||
                 PointIsCompletelyWithinRect( second, rect ) )
            {
                if ( ! inPolygon )
                {
                    anyPointInside = true;
                    inPolygon = true;
                    currentPolygon = shared_ptr<Polygon>(new Polygon());
                    currentPolygon->emplace_back(first);
                }
                currentPolygon->emplace_back(second);
            }
            else
            {
                bool clockwiseOnEdge = false;
                
                EdgeType firstEdge = EdgeForCoord( first, rect );
                EdgeType secondEdge = EdgeForCoord( second, rect );
                if ( PointIsOnEdge( first, rect ) && PointIsOnEdge( second, rect ) )
                {
                    int firstToSecond = ClockwiseDistance(first, second, rect);
                    int secondToFirst = ClockwiseDistance(second, first, rect);
                    clockwiseOnEdge = firstToSecond < secondToFirst;
                }

                if ( firstEdge != secondEdge || clockwiseOnEdge )
                {
                    if ( ! inPolygon )
                    {
                        anyPointInside = true;
                        inPolygon = true;
                        currentPolygon = shared_ptr<Polygon>(new Polygon());
                        currentPolygon->emplace_back(first);
                    }
                    currentPolygon->emplace_back(second);
                }
                else
                {
                    if ( inPolygon )
                    {
                        inPolygon = false;
                        currentPolygon->emplace_back( second );
                        output.emplace_back( *currentPolygon );
                    }
                }
            }
        }
        
        if ( inPolygon )
        {
            inPolygon = false;
            output.emplace_back( *currentPolygon );
        }
        
        if ( ! anyPointInside )
        {
            //
            // All points lie on the boundary - could be intersection or clockwise
            //
            output.emplace_back( polygon );
        }
        
        assert ( !inPolygon );
    }
    
    return output;
}

typedef enum
{
    ContourResultNothingInside,
    ContourResultEntirelyInside,
    ContourResultCompletelyCovers,
    ContourResultHaveContours
} ContourResult;
    
    
ContourResult GetClippedContours( const Contour &inputPolygon, struct rect rect, vector<Contour> &clippedContours )
{
    if ( inputPolygon.size() < 3 )
    {
        return ContourResultNothingInside;
    }
    
    size_t firstPointOutside = SIZE_MAX;
    bool allPointsAbove = true;
    bool allPointsRight = true;
    bool allPointsBelow = true;
    bool allPointsLeft = true;
    
    for ( size_t i = 0; i < inputPolygon.size(); i++ )
    {
        struct coord c = inputPolygon[i];
        if ( ! PointIsInsideRect(c, rect ) )
        {
            if ( firstPointOutside == SIZE_MAX )
            {
                firstPointOutside = i;
            }
        }

        if ( c.y < rect.h.y ) { allPointsAbove = false; }
        if ( c.y > rect.l.y ) { allPointsBelow = false; }
        if ( c.x < rect.h.x ) { allPointsRight = false; }
        if ( c.x > rect.l.x ) { allPointsLeft = false; }
    }
    
    //
    //  Simple cases - all inside, or no possible overlap
    //
    if ( firstPointOutside == SIZE_MAX )
    {
        return ContourResultEntirelyInside;
    }
    
    if ( allPointsAbove || allPointsLeft || allPointsRight || allPointsBelow )
    {
        return ContourResultNothingInside;
    }
    
    //
    //  Search for incident edges
    //
    size_t inputSize = inputPolygon.size();
    size_t firstPointBackInside = -1;
    
    for ( size_t i = 0; i < inputPolygon.size(); i++ )
    {
        size_t outsidePointIndex = (i + firstPointOutside) % inputSize;
        size_t insidePointIndex = (i + firstPointOutside + 1) % inputSize;
        struct coord outsidePoint = inputPolygon[outsidePointIndex];
        struct coord insidePoint = inputPolygon[insidePointIndex];
        
        if ( PointIsInsideRect(insidePoint, rect )  ||
             LineIntersetsRect(outsidePoint, insidePoint, rect) )
        {
            firstPointBackInside = insidePointIndex;
            break;
        }
    }
    
    if ( firstPointBackInside == SIZE_MAX )
    {
        //
        //  No intersections or points inside. Either completely surrounded and included, or completely surrounded with no overlap
        //
        struct coord rectCenter = { ( rect.l.x + rect.h.x ) / 2, ( rect.l.y + rect.h.y ) / 2 };
        if ( ! PointIsInsideContour( inputPolygon, rectCenter ) )
        {
            //
            // We are not covered
            //
            return ContourResultNothingInside;
        }
        
        //
        // We are covered
        //
        return ContourResultCompletelyCovers;
    }
    
    assert ( firstPointBackInside != SIZE_MAX );
    
    size_t lastPointOutside = firstPointBackInside -1;
    if ( lastPointOutside == SIZE_MAX )
    {
        lastPointOutside = inputSize - 1;
    }
    
    //auto area = geom_poly_area(&inputPolygon[0], inputPolygon.size());
    
    GetAllSubContours(inputPolygon, rect, lastPointOutside, clippedContours );
    
    return ContourResultHaveContours;
}
    
    
int PointInPolygon( const Polygon &polygon, const coord coordinate )
{
    //returns 0 if false, +1 if true, -1 if coordinate ON polygon boundary
    int result = 0;
    size_t cnt = polygon.size();
    if (cnt < 3) return 0;
    coord ip = polygon[0];
    
    for(size_t i = 1; i <= cnt; ++i)
    {
        coord ipNext = (i == cnt ? polygon[0] : polygon[i]);
        if (ipNext.y == coordinate.y)
        {
            if ((ipNext.x == coordinate.x) || (ip.y == coordinate.y &&
               ((ipNext.x > coordinate.x) == (ip.x < coordinate.x)))) return -1;
        }
        if ((ip.y < coordinate.y) != (ipNext.y < coordinate.y))
        {
            if (ip.x >= coordinate.x)
            {
                if (ipNext.x > coordinate.x) result = 1 - result;
                else
                {
                    double d = (double)(ip.x - coordinate.x) * (ipNext.y - coordinate.y) -
                    (double)(ipNext.x - coordinate.x) * (ip.y - coordinate.y);
                    if (!d) return -1;
                    if ((d > 0) == (ipNext.y > ip.y)) result = 1 - result;
                }
            }
            else
            {
                if (ipNext.x > coordinate.x)
                {
                    double d = (double)(ip.x - coordinate.x) * (ipNext.y - coordinate.y) -
                    (double)(ipNext.x - coordinate.x) * (ip.y - coordinate.y);
                    if (!d) return -1;
                    if ((d > 0) == (ipNext.y > ip.y)) result = 1 - result;
                }
            }
        }
        ip = ipNext;
    } 
    return result;
}
    
    
bool PointPreceedsOnEdge( EdgeType edge, struct coord first, struct coord second )
{
    switch ( edge )
    {
        case EdgeTop:
            assert( first.y == second.y );
            return first.x <= second.x;
        case EdgeRight:
            assert( first.x == second.x );
            return first.y >= second.y;
        case EdgeBottom:
            assert( first.y == second.y );
            return second.x <= first.x;
        case EdgeLeft:
            assert( first.x == second.x );
            return second.y >= first.y;
        case EdgeNotAnEdge:
            assert( false );
            return false;
    }
}
    
    
int IndexOfPreceedingPointOnEdge( const Polygon &polygon, const struct coord coordinate, struct rect rect )
{
    assert ( EdgeForCoord( coordinate, rect ) != EdgeNotAnEdge );
    
    EdgeType targetEdge = EdgeForCoord( coordinate, rect );
    
    for ( size_t i = 0; i < polygon.size(); i++ )
    {
        const auto &first = polygon[i];
        const auto &second = polygon[(i == polygon.size() - 1) ? 0 : (i + 1)];
        
        // Both points must be on an edge, first point must be on targetEdge
        
        if ( EdgeForCoord( first, rect ) == EdgeNotAnEdge )
        {
            continue;
        }
        
        if ( EdgeForCoord( second, rect ) == EdgeNotAnEdge )
        {
            continue;
        }
        
        EdgeType firstEdge = EdgeForCoord( first, rect );
        
        if ( firstEdge != targetEdge )
        {
            continue;
        }
        
        assert( firstEdge == targetEdge );
        
        if ( ! PointPreceedsOnEdge( targetEdge, first, coordinate) )
        {
            continue;
        }
        
        EdgeType secondEdge = EdgeForCoord( second, rect );
        
        if ( secondEdge != targetEdge || PointPreceedsOnEdge(targetEdge, coordinate, second ) )
        {
            return (int) i;
        }
    }
    
    return -1;
}
    

void PolygonsWithIncidentEdge( Polygon &polygon, struct rect rect, vector<Polygon> &output )
{
    if ( polygon.size() < 2 )
    {
        return;
    }
    
    for ( size_t firstIndex = 0; firstIndex < polygon.size(); firstIndex++ )
    {
        EdgeType firstEdge = EdgeForCoord( polygon[firstIndex], rect );
        
        if ( firstEdge == EdgeNotAnEdge )
        {
            continue;
        }
        
        while ( EdgeForCoord( polygon[ ( firstIndex + 1 ) % polygon.size() ], rect ) != EdgeNotAnEdge )
        {
            firstIndex = ( firstIndex + 1 ) % polygon.size();
            if ( firstIndex == 0 )
            {
                return;     // Back to start again
            }
        }
        
        assert( EdgeForCoord( polygon[ firstIndex ], rect ) != EdgeNotAnEdge );
        
        size_t secondIndex;
        for ( secondIndex = ( firstIndex + 1 ) % polygon.size(); secondIndex != firstIndex; secondIndex = ( secondIndex + 1 ) % polygon.size() )
        {
            if ( EdgeForCoord( polygon[secondIndex], rect ) != EdgeNotAnEdge )
            {
                break;
            }
        }
        
        assert( EdgeForCoord( polygon[ firstIndex ], rect ) != EdgeNotAnEdge );
        assert( EdgeForCoord( polygon[ secondIndex ], rect ) != EdgeNotAnEdge );
        
        Polygon outp;
        if ( firstIndex < secondIndex )
        {
            outp.insert( outp.end(), polygon.begin() + firstIndex, polygon.begin() + secondIndex + 1 );
            output.emplace_back( outp );
        }
        else
        {
            outp.insert( outp.end(), polygon.begin() + firstIndex, polygon.end() );
            outp.insert( outp.end(), polygon.begin(), polygon.begin() + secondIndex + 1 );
            output.emplace_back( outp );
            
            //  Can only look over end -> start once
            return;
        }
        
        firstIndex = secondIndex;
    }
}
  
    
void EmbedIncidentEdgeHoleIntoParent( Polygon &polygon, struct rect rect, Polygon &hole, vector<Polygon> &cutouts )
{
    int insertionIndex = IndexOfPreceedingPointOnEdge( polygon, hole[0], rect );
    int endInsertionIndex = IndexOfPreceedingPointOnEdge( polygon, hole.back(), rect );
    
    if ( insertionIndex == -1 && endInsertionIndex == -1 )
    {
        cutouts.emplace_back( hole );
        return;
    }
    
    //
    //  Slice off all points between start and end. These might become new polygons.
    //  Anyhow, as we're embedding a hole, we need to cut them out.
    //
    assert( insertionIndex != -1 );
    
    bool cutCorners = true;
    if ( endInsertionIndex == - 1 || insertionIndex == endInsertionIndex )
    {
        cutCorners = false;
    }
    else if ( insertionIndex < endInsertionIndex )
    {
        Polygon newPolygon;
        newPolygon.insert(newPolygon.end(), polygon.begin() + insertionIndex + 1, polygon.begin() + endInsertionIndex + 1);
        cutouts.emplace_back(newPolygon);

        polygon.erase( polygon.begin() + insertionIndex + 1, polygon.begin() + endInsertionIndex + 1 );
    }
    else if ( endInsertionIndex < insertionIndex )
    {
        // We're inserting over the start/end boundary of the parent while cutting corners
        size_t nrCutCornersEnd = polygon.size() - insertionIndex - 1;
        size_t nrCutCornersStart = endInsertionIndex + 1;
        
        Polygon newPolygon;
        newPolygon.insert(newPolygon.end(), polygon.begin() + insertionIndex + 1, polygon.end());
        newPolygon.insert(newPolygon.end(), polygon.begin(), polygon.begin() + nrCutCornersStart);
        cutouts.emplace_back(newPolygon);

        if ( nrCutCornersEnd > 0 )
        {
            polygon.erase( polygon.begin() + insertionIndex + 1, polygon.end() );
        }
        
        polygon.erase( polygon.begin(), polygon.begin() + nrCutCornersStart );
    }
    
    //
    //  Recalc insertion point
    //
    if ( cutCorners )
    {
        assert ( EdgeForCoord( hole[0], rect ) != EdgeNotAnEdge && EdgeForCoord( hole.back(), rect ) != EdgeNotAnEdge );
        
        insertionIndex = IndexOfPreceedingPointOnEdge( polygon, hole[0], rect );
        
        if ( insertionIndex == -1 )
        {
            // ToDo
            cutouts.emplace_back(hole);
            return;
        }
        
        assert ( insertionIndex != -1 );
        
        endInsertionIndex = IndexOfPreceedingPointOnEdge( polygon, hole.back(), rect );
    }
    
    polygon.insert( polygon.begin() + insertionIndex + 1, hole.begin(), hole.end() );
    return;
}
 
    
void EmbedEdgeHolesIntoParent( Polygon &polygon, struct rect rect, vector<Polygon> &cutouts )
{
    for ( size_t holeIdx = 0; holeIdx < polygon.holes.size(); holeIdx++ )
    {
        Polygon &hole = polygon.holes[holeIdx];
        
        vector<Polygon> incidentEdges;
        PolygonsWithIncidentEdge( hole, rect, incidentEdges );
        
        if ( incidentEdges.size() == 0 )
        {
            continue;
        }
        
        for ( auto &incidentEdge : incidentEdges )
        {
            EmbedIncidentEdgeHoleIntoParent( polygon, rect, incidentEdge, cutouts );
        }
        
        polygon.holes.erase( polygon.holes.begin() + holeIdx );
        holeIdx--;
    }
}

    
void AssignHolesToOuterPolygons( std::vector<Polygon> &outerPolygons, const std::vector<Polygon> &holes )
{
    for ( auto &hole : holes )
    {
        for ( auto &outerPoly : outerPolygons )
        {
            int pointInPoly = PointInPolygon( outerPoly, hole[0] );
            
            if ( pointInPoly == -1 )
            {
                pointInPoly = PointInPolygon( outerPoly, hole[hole.size() / 2] );
            }
            
            if ( pointInPoly != 0 )
            {
                outerPoly.holes.emplace_back( hole );
                break;
            }
        }
    }
}

    
void RayClipPolygon( const Polygon &inputPolygon, struct rect rect, vector<Polygon> &outputPolygons )
{
    vector<Contour> contours;
    ContourResult outsideResult = GetClippedContours(inputPolygon, rect, contours);
    
    if ( outsideResult == ContourResultNothingInside )
    {
        return;
    }
    
    if ( outsideResult == ContourResultEntirelyInside )
    {
        outputPolygons.emplace_back( inputPolygon );
        
        vector<Polygon> cutouts;
        
        //
        //  Still embed the holes, if there are any
        //
        for ( auto &polygon : outputPolygons )
        {
            EmbedEdgeHolesIntoParent( polygon, rect, cutouts );
        }
        
        for ( auto &cutout : cutouts )
        {
            vector<Polygon> incidentCutouts;
            PolygonsWithIncidentEdge(cutout, rect, incidentCutouts);
            
            for ( auto &incidentCutout : incidentCutouts )
            {
                if ( incidentCutout.size() > 1 )
                {
                    vector<Contour> otherContours;
                    ClosePolygon(incidentCutout, otherContours, rect);
                    outputPolygons.emplace_back( incidentCutout );
                }
            }
        }

        return;
    }
    
    assert ( outsideResult == ContourResultCompletelyCovers || outsideResult == ContourResultHaveContours );
    
    vector<Polygon> internalHoles;
    
    for ( auto &hole : inputPolygon.holes )
    {
        ContourResult insideResult = GetClippedContours( hole, rect, contours );
        
        switch ( insideResult )
        {
            case ContourResultNothingInside:
                continue;
            case ContourResultCompletelyCovers:
                assert ( outsideResult == insideResult );   // If the hole totally covers, the parent must also
                return;
            case ContourResultEntirelyInside:
                internalHoles.emplace_back( hole );
                continue;
            case ContourResultHaveContours:
                //  Result already in contours
                break;
        }
    }
    
    if ( contours.size() == 0 )
    {
        if ( outsideResult == ContourResultCompletelyCovers )
        {
            Polygon rectAsPolygon;
            rectAsPolygon.emplace_back(rect.l);
            rectAsPolygon.emplace_back(coord{rect.l.x, rect.h.y});
            rectAsPolygon.emplace_back(rect.h);
            rectAsPolygon.emplace_back(coord{rect.h.x, rect.l.y});
            outputPolygons.emplace_back( rectAsPolygon );
            
            for ( auto &internalHole : internalHoles )
            {
                rectAsPolygon.holes.emplace_back( internalHole );
            }
        }
        
        return;
    }
    
    ClosePolygons(contours, rect, outputPolygons);
    AssignHolesToOuterPolygons(outputPolygons, internalHoles);
    
    vector<Polygon> cutouts;

    for ( auto &polygon : outputPolygons )
    {
        EmbedEdgeHolesIntoParent( polygon, rect, cutouts );
    }
    
    for ( auto &cutout : cutouts )
    {
        vector<Polygon> incidentCutouts;
        PolygonsWithIncidentEdge(cutout, rect, incidentCutouts);
        
        for ( auto &incidentCutout : incidentCutouts )
        {
            if ( incidentCutout.size() > 1 )
            {
                vector<Contour> otherContours;
                ClosePolygon(incidentCutout, otherContours, rect);
                outputPolygons.emplace_back( incidentCutout );
            }
        }
    }
}
    
    
bool LastSelfIntersection( const Polygon &inputPolygon, const size_t startIndex, size_t &continuationIndex, struct coord &intersectionPoint )
{
    struct coord start = inputPolygon[startIndex];
    struct coord end = inputPolygon[startIndex + 1];
    
    for ( size_t i = inputPolygon.size() - 2; i > startIndex; i-- )
    {
        struct coord testStart = inputPolygon[i];
        struct coord testEnd = inputPolygon[i + 1];
        
        if ( LineIntersetsLine(start, end, testStart, testEnd, intersectionPoint ) )
        {
            continuationIndex = i + 1;
            return true;
        }
    }
    
    return false;
}
    
    
void CleanPolygon( const Polygon &inputPolygon, Polygon &outputPolygon )
{
    for ( size_t i = 0; i < inputPolygon.size() - 1; i++ )
    {
        auto first = inputPolygon[i];
        outputPolygon.emplace_back(first);
        
        size_t continuationIndex = 0;
        struct coord intersectionPoint;
        if ( LastSelfIntersection( inputPolygon, i, continuationIndex, intersectionPoint ) )
        {
            assert( continuationIndex > i + 1 );
            outputPolygon.emplace_back(intersectionPoint);
            i = continuationIndex - 1;
        }
    }
    
    outputPolygon.emplace_back( inputPolygon.back() );
}
    
    
long long PolygonArea(const Polygon &coords)
{
    long long area=0;
    size_t i,j=0;
    
    for (i = 0; i < coords.size(); i++)
    {
        if (++j == coords.size())
            j = 0;
        area += ( long long )(coords[i].x + coords[j].x ) * (coords[i].y - coords[j].y);
    }
    return area / 2;
}

}
