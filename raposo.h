#ifndef _RAPOSO_H
#define _RAPOSO_H

#include "OpencvLibrary.h"

std::vector<cv::Point2f> raposo(
    const std::vector<cv::Point2f>& points,
    int n,
    int max_tries = 25,
    bool keep_first_last = true
);


#endif // _RAPOSO_H


