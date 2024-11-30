/**
* @file maptiles.cpp
* Code for the maptiles function.
*/

#include <iostream>
#include <map>
#include <vector>
#include <memory>

#include "maptiles.h"

Point<3> convertToXYZ(LUVAPixel pixel) {
    return Point<3>(pixel.l, pixel.u, pixel.v);
}


MosaicCanvas* mapTiles(SourceImage const& theSource,
                       vector<TileImage>& theTiles)
{
    /**
     * @todo Implement this function!
     */
    MosaicCanvas* mosaic = new MosaicCanvas(theSource.getRows(), theSource.getColumns());
    map<Point<3>, int> srcTileMap;
    vector<Point<3>> currPoints;
    for(size_t i = 0; i < theTiles.size(); i++) {
      Point<3> tmpPoint = convertToXYZ(theTiles[i].getAverageColor());
      currPoints.push_back(tmpPoint);
      srcTileMap[tmpPoint] = i;
    }
    KDTree<3>* currTree = new KDTree<3>(currPoints);
    for(int row = 0; row < theSource.getRows(); row++) {
      for(int col = 0; col < theSource.getColumns(); col++) {
        Point<3> tmpPoint = convertToXYZ(theSource.getRegionColor(row, col));
        Point<3> closePoint = currTree->findNearestNeighbor(tmpPoint);
        int nearestIndex = srcTileMap[closePoint];
        mosaic->setTile(row, col, &theTiles[nearestIndex]);
      }
    }
    return mosaic;
}