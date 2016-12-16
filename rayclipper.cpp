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
    
    
typedef enum
{
    EdgeTop = 0,
    EdgeRight = 1,
    EdgeBottom = 2,
    EdgeLeft = 3
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
    float dx=p2.x-p1.x;
    float dy=p2.y-p1.y;
    switch(edge) {
        case EdgeLeft:
            result.y=p1.y + ((float)(r.l.x-p1.x)) * dy/dx;
            result.x=r.l.x;
            break;
        case EdgeRight:
            result.y=p1.y + ((float)(r.h.x-p1.x)) * dy/dx;
            result.x=r.h.x;
            break;
        case EdgeTop:
            result.x=p1.x + ((float)(r.h.y-p1.y)) * dx/dy;
            result.y=r.h.y;
            break;
        case EdgeBottom:
            result.x=p1.x + ((float)(r.l.y-p1.y)) * dx/dy;
            result.y=r.l.y;
            break;
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
    
    
bool PointIsInsidePolygon(const Polygon &coords, const struct coord &point)
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


vector<Polygon> GetAllSubPolygons( const Polygon &inputPolygon, struct rect rect, size_t lastPointOutside )
{
    vector<Polygon> allPolygons;
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
                
                //assert ( currentPolygon->size() >= 3 );
                allPolygons.emplace_back( *currentPolygon );
                
                inside = false;
            }
            else
            {
                //
                // Both outside and currently outside - could be intersection
                //
                vector<struct coord> intersections = AllIntersectionsOfRect(firstPoint, secondPoint, rect);
                
                if ( intersections.size() == 2 )
                {
                    currentPolygon = shared_ptr<Polygon>(new Polygon());
                    
                    assert ( intersections[0].y == rect.l.y || intersections[0].y == rect.h.y ||
                             intersections[0].x == rect.l.x || intersections[0].x == rect.h.x );
                    assert ( intersections[1].y == rect.l.y || intersections[1].y == rect.h.y ||
                             intersections[1].x == rect.l.x || intersections[1].x == rect.h.x );
                    
                    assert( ! coord_is_equal(intersections[0], intersections[1]) );
                    
                    double d0 = CoordDistance(firstPoint, intersections[0]);
                    double d1 = CoordDistance(firstPoint, intersections[1]);
                    
                    if ( d0 <= d1 )
                    {
                        currentPolygon->emplace_back(intersections[0]);
                        currentPolygon->emplace_back(intersections[1]);
                    }
                    else
                    {
                        currentPolygon->emplace_back(intersections[1]);
                        currentPolygon->emplace_back(intersections[0]);
                    }

                    allPolygons.emplace_back( *currentPolygon );
                }
            }
        }
    }
    
    assert( ! inside );
    //assert( allPolygons.size() > 0 );
    
    return allPolygons;
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
    
    assert ( false );
    return EdgeLeft;    // Some compiler configs shout error otherwise
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


void ClosePolygon( Polygon &polygon, vector<Polygon> &otherPolygons, struct rect rect )
{
    while ( true )
    {
        size_t closestPolygon = -1;
        int closestPolygonDistance = INT_MAX;
        
        for ( size_t i = 0; i < otherPolygons.size(); i++ )
        {
            auto &other = otherPolygons[i];
            
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
            //  No corners or other polygons in the way
            return;
        }
        
        if ( closestPolygonDistance <= distanceToNextCorner )
        {
            auto &other = otherPolygons[closestPolygon];
            
            size_t toSkipAtStart = 0;
            size_t toSkipAtEnd = 0;
            
            if ( coord_is_equal( polygon.back(), other.front() ) )
            {
                toSkipAtStart = 1;
            }
            
            //if ( coord_is_equal( other.back(), polygon.front() ) )
            //{
            //    toSkipAtEnd = 1;
            //}
            
            assert( other.begin() + toSkipAtStart <= other.end() - toSkipAtEnd );

            polygon.insert( polygon.end(), other.begin() + toSkipAtStart, other.end() - toSkipAtEnd );
            
            //assert ( !coord_is_equal(polygon.front(), polygon.back()) );
            
            otherPolygons.erase( otherPolygons.begin() + closestPolygon );
            
            continue;
        }
        
        assert( distanceToNextCorner < distanceToEndOfCurrent && distanceToNextCorner < closestPolygonDistance );
        polygon.emplace_back( nextCorner );
    }
}


vector<Polygon> ClosePolygons( vector<Polygon> &polygons, struct rect rect )
{
    vector<Polygon> closedPolygons;
    
    while ( polygons.size() > 0 )
    {
        auto polygon = polygons.back();
        polygons.pop_back();
        ClosePolygon(polygon, polygons, rect);
        
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
    
    return closedPolygons;
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
        /*if ( inPolygon )
        {
            output.emplace_back( *currentPolygon );
        }*/
    }
    
    return output;
}


vector<Polygon> RayClipPolygon( const Polygon &inputPolygon, struct rect rect )
{
    vector<Polygon> output;

    if ( inputPolygon.size() < 3 )
    {
        return output;
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
        output.emplace_back(inputPolygon);
        return output;
    }
    
    if ( allPointsAbove || allPointsLeft || allPointsRight || allPointsBelow )
    {
        return output;
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
        if ( ! PointIsInsidePolygon( inputPolygon, rectCenter ) )
        {
            //
            // We are not covered
            //
            return output;
        }
        
        //
        // We are covered
        //
        Polygon rectAsPolygon;
        rectAsPolygon.emplace_back(rect.l);
        rectAsPolygon.emplace_back(coord{rect.l.x, rect.h.y});
        rectAsPolygon.emplace_back(rect.h);
        rectAsPolygon.emplace_back(coord{rect.h.x, rect.l.y});
        output.emplace_back( rectAsPolygon );
        return output;
    }
    
    assert ( firstPointBackInside != SIZE_MAX );
    
    size_t lastPointOutside = firstPointBackInside -1;
    if ( lastPointOutside == SIZE_MAX )
    {
        lastPointOutside = inputSize - 1;
    }
    
    //auto area = geom_poly_area(&inputPolygon[0], inputPolygon.size());
    
    vector<Polygon> allPolygons = GetAllSubPolygons(inputPolygon, rect, lastPointOutside );
    //vector<Polygon> splitPolygons = SplitEdgeTouchingPolygons( allPolygons, rect );
    vector<Polygon> closedPolygons = ClosePolygons(allPolygons, rect);
    //assert ( closedPolygons.size() > 0 );
    return closedPolygons;
}


}
