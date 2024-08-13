#ifndef _APSC_H
#define _APSC_H

#include <vector>
#include <iostream>
#include <opencv2/opencv.hpp>

typedef struct _VertexRelation
{
    std::vector<int> ID;
    std::vector<cv::Point2f> XY;
    std::vector<double> error;
    std::vector<int> parent;
    std::vector<int> LC;
    std::vector<int> RC;
    std::vector<int> LS;
    std::vector<int> RS;
    std::vector<bool> current;

}VertexRelation;

typedef struct _PriorityList
{
    float displacement;
    int A;
    int B;
    int C;
    int D;
    cv::Point2f pE;
    int overlap_endpt;
}PriorityList;


std::vector<cv::Point2f> apsc(std::vector<cv::Point2f> pts, int min_pts_num);

int initVertexTree(const std::vector<cv::Point2f>& pts, VertexRelation *vertex_tree);

#endif // _APSC_H

