// YangShuo
// 2024-8-13
#include <opencv2/core.hpp>  // Include OpenCV core module
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>

// 计算三角形面积
double triangleArea(const cv::Point2f& a, const cv::Point2f& b, const cv::Point2f& c)
{
    return 0.5 * std::abs(a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y));
}

// Visvalingam算法实现
std::vector<cv::Point2f> visvalingam(const std::vector<cv::Point2f>& points, int numPoints)
{
    std::vector<cv::Point2f> simplifiedPoints = points;

    while (simplifiedPoints.size() > numPoints) {
        double minArea = std::numeric_limits<double>::max();
        int minIndex = -1;

        // 找到面积最小的顶点
        for (size_t i = 1; i < simplifiedPoints.size() - 1; ++i) {
            double area = triangleArea(simplifiedPoints[i - 1], simplifiedPoints[i], simplifiedPoints[i + 1]);
            if (area < minArea) {
                minArea = area;
                minIndex = i;
            }
        }

        // 移除面积最小的顶点
        if (minIndex != -1) {
            simplifiedPoints.erase(simplifiedPoints.begin() + minIndex);
        }
    }

    return simplifiedPoints;
}



