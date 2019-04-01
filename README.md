# rayclipper
An O(n) polygon to rectangle clipper for simple polygons. This implementation assumes that the polygon has a negative winding (-1 only, as polygon’s must be simple). In other words, the polygon must be defined in a clockwise direction, with the inside of the polygon to the right hand side of each line segment. Holes are represented by a polygon in the opposite winding. If anyone asks I’ll add positive windings too.

This clipper does not check that your polygons are valid and simple (i.e. have no self intersections). If the integrity of the simplicity is violated (without being corrected) over the boundary of the rectangle, you will end up with a rectangle sized flood-filled output polygon. There are a number of tools for checking your polygon is valid before you pass it in. I would recommend Angus Johnson’s excellent clipper library (http://www.angusj.com/delphi/clipper.php) for this purpose. In fact, I would recommend you always use his library unless you are clipping to a rectangle and you are in a time critical situation.

To compile, make sure you have C++11 enabled, e.g. g++ -std=c++11 rayclipper.cpp

Input polygons do not need to have a closed path (i.e. first and last coordinates do not need to match).

The input is the polygon you want to clip and the rectangle you want to clip it to. The output is a set of one or more polygons, which are clipped to the rectangle.

Here are some sample inputs and outputs:

<img width="1091" alt="clipper1" src="https://cloud.githubusercontent.com/assets/1854581/21275928/1a767600-c39d-11e6-9536-75699360fa42.png">
<img width="1012" alt="clipper2" src="https://cloud.githubusercontent.com/assets/1854581/21275930/1a7c9aa8-c39d-11e6-9e75-df0540e36308.png">
<img width="1003" alt="clipper3" src="https://cloud.githubusercontent.com/assets/1854581/21275929/1a7c75e6-c39d-11e6-979f-38c7675d9827.png">
<img width="1004" alt="clipper4" src="https://cloud.githubusercontent.com/assets/1854581/21275931/1a8024b6-c39d-11e6-8485-1f4cc9560baf.png">
<img width="882" alt="clipper5" src="https://user-images.githubusercontent.com/1854581/55344572-f8b1fb80-54ad-11e9-94f1-00c95fa47b1f.png">
<img width="882" alt="clipper6" src="https://user-images.githubusercontent.com/1854581/55344586-02d3fa00-54ae-11e9-9c77-4386462c0e14.png">
<img width="882" alt="clipper7" src="https://user-images.githubusercontent.com/1854581/55344598-0c5d6200-54ae-11e9-88dc-40da41e471a3.png">
<img width="882" alt="clipper8" src="https://user-images.githubusercontent.com/1854581/55344603-11221600-54ae-11e9-91c9-2bc8750904ef.png">
<img width="882" alt="clipper9" src="https://user-images.githubusercontent.com/1854581/55344607-154e3380-54ae-11e9-9a8c-4666c695f331.png">
<img width="882" alt="clipper10" src="https://user-images.githubusercontent.com/1854581/55344614-18492400-54ae-11e9-97fc-0824fd14bff5.png">

