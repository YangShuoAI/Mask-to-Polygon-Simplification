#ifndef _ULTRALYTICS_INTERP_H
#define _ULTRALYTICS_INTERP_H
#include <iostream>
#include <vector>
#include "OpencvLibrary.h"

std::vector<cv::Point2f> ultralytics_interp(const std::vector<cv::Point2f>& segments, int n = 1000);


#endif // _ULTRALYTICS_INTERP_H




