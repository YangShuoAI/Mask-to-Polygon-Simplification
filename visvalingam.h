#ifndef _VISVALINGAM_H
#define _VISVALINGAM_H

#include <iostream>
#include <opencv2/opencv.hpp>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>

std::vector<cv::Point2f> visvalingam(const std::vector<cv::Point2f>& points, int numPoints);


#endif // _VISVALINGAM_H


