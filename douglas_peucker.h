#ifndef _DOUGLAS_PEUCKER_H
#define _DOUGLAS_PEUCKER_H
//
//-------------------------------------------------------------------------------
// Name:        Douglas Peucker
// Description : A very simple implementation of the Ramer - Douglas - Peucker
//              line simplification algorithm. This does not check for 
//              topological errors.The implementation here is not particularly
//              fast, but it does support pre - processing to quickly simplify
//              to any given distance tolerance or number of points.
//
//              References:
//              
//              Urs Ramer(1972).An iterative procedure for the polygonal approximation
//              of plane curves.Computer Graphics and Image Processing, vol. 1, pp. 244 - 256.
//
//              David H.Douglas and Thomas K.Peucker(1973).Algorithms for the reduction of
//              the number of points required to represent a digitized line or its caricature.
//              The Canadian Cartographer, vol. 10, no. 2, pp. 112 - 122.
//
// Author:      Barry Kronenfeld
// License : MIT License
//-------------------------------------------------------------------------------
#include "OpencvLibrary.h"

std::vector<cv::Point2f> douglas_peucker(const std::vector<cv::Point2f>& pts, int numpts);

#endif // _DOUGLAS_PEUCKER_H

