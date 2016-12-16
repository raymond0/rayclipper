# rayclipper
An O(n) polygon to rectangle clipper for simple polygons. This implementation assumes that the polygon has a negative winding (-1 only, as polygon’s must be simple). In other words, the polygon must be defined in a clockwise direction, with the inside of the polygon to the right hand side of each line segment. If anyone asks I’ll add it for positive windings too.

This clipper does not check that your polygons are valid and simple (i.e. have no self intersections). If the integrity of the simplicity is violated (without being corrected) over the boundary of the rectangle, you will end up with a rectangle sized flood-filled output polygon. There are a number of tools for checking your polygon is valid before you pass it in. I would recommend Angus Johnson’s excellent clipper library (http://www.angusj.com/delphi/clipper.php) for this purpose. In fact, I would recommend you always use his library unless you are clipping to a rectangle and you are in a time critical situation.

To compile, make sure you have C++11 enabled, e.g. g++ -std=c++11 rayclipper.cpp

Input polygons do not need to have a closed path (i.e. first and last coordinates do not need to match).

The input is the polygon you want to clip and the rectangle you want to clip it to. The output is a set of one or more polygons, which are clipped to the rectangle.