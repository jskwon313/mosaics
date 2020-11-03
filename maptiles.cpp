/**
 * @file maptiles.cpp
 * Code for the maptiles function.
 */

#include <iostream>
#include <map>
#include "maptiles.h"
//#include "cs225/RGB_HSL.h"

using namespace std;


Point<3> convertToXYZ(LUVAPixel pixel) {
    return Point<3>( pixel.l, pixel.u, pixel.v );
}

MosaicCanvas* mapTiles(SourceImage const& theSource,
                       vector<TileImage>& theTiles)
{
    /**
     * @todo Implement this function!
     */
    // Initializes variables
    int width = theSource.getRows();
    int height = theSource.getColumns();
    MosaicCanvas *canvas = new MosaicCanvas(width, height);
    vector<Point<3>> v;
    map<Point<3>, int> map;
    Point<3> nearestPoint;
    int temp = 0;
    auto itr = map.end();

    // Creates the KDTree based of the TileImages
    for (unsigned i = 0; i < theTiles.size(); i++) {
      v.push_back(convertToXYZ(theTiles[i].getAverageColor()));
      map.insert(pair<Point<3>, int> (convertToXYZ(theTiles[i].getAverageColor()), i));
    }
    KDTree<3> tree(v);

    // Assigns images based on the average color
    for (int i = 0; i < width; i++) {
      for (int j = 0; j < height; j++) {
        nearestPoint = tree.findNearestNeighbor(convertToXYZ(theSource.getRegionColor(i, j)));
        if (!map.empty())
          itr = map.find(nearestPoint);
        if (itr != map.end())
          temp = itr->second;
        canvas->setTile(i, j, &(theTiles[temp]));
      }
    }

    return canvas;
}
