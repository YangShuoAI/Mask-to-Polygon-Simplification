#include <opencv2/opencv.hpp>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <limits>
#include "OpencvLibrary.h"


// cHexGrid 类定义
class cHexGrid {
private:
    static constexpr double sqrt_3 = 1.7320508075688772; // sqrt(3)
    double top, left, bottom, right, radius, colWid, rowHt;

public:
    cHexGrid(double top, double left, double bottom, double right, double radius)
        : top(top), left(left), bottom(bottom), right(right), radius(radius) {
        colWid = 1.5 * radius;
        rowHt = radius * sqrt_3;
    }

    cv::Point hexRowCol(double x, double y) const {
        int C = static_cast<int>((x + colWid / 2) / colWid);
        double nearX, nearY;
        cv::Point2f pt = hexCenter(0, C);
        nearX = pt.x;
        nearY = pt.y;
        int cOff = (x < nearX) ? -1 : 1;

        int R;
        int rOff;
        if (C % 2 == 0) {
            R = static_cast<int>((y + rowHt / 2) / rowHt);
            rOff = -1;
        }
        else {
            R = static_cast<int>(y / rowHt);
            rOff = 1;
        }

        std::vector<cv::Point> candidates = { cv::Point(R, C), cv::Point(R, C + cOff), cv::Point(R + rOff, C + cOff) };

        cv::Point winner = candidates[0];
        double mind2 = radius * radius * 2;

        for (const auto& candidate : candidates) {
            double cX, cY;
            cv::Point2f pt = hexCenter(candidate.y, candidate.x);
            cX = pt.x;
            cY = pt.y;
            double d2 = (x - cX) * (x - cX) + (y - cY) * (y - cY);
            if (d2 <= mind2) {
                mind2 = d2;
                winner = candidate;
            }
        }
        return winner;
    }

    cv::Point2f hexCenter(int row, int col) const {
        float x = col * colWid;
        float y = row * rowHt;
        if (col % 2 == 1) {
            y += rowHt / 2.0;
        }
        return cv::Point2f(x, y);
    }
};

// 计算均值
cv::Point2f quant_avg(const std::vector<cv::Point2f>& points, int startid, int nextstartid) {
    double meanx = 0.0, meany = 0.0;
    int n = nextstartid - startid;
    for (int i = startid; i < nextstartid; ++i) {
        meanx += points[i].x;
        meany += points[i].y;
    }
    meanx /= n;
    meany /= n;
    return cv::Point2f(meanx, meany);
}

// 计算中点
cv::Point2f quant_mid(const std::vector<cv::Point2f>& points, int startid, int nextstartid) {
    const auto& p1 = points[startid];
    const auto& p2 = points[nextstartid - 1];
    return cv::Point2f((p1.x + p2.x) / 2.0, (p1.y + p2.y) / 2.0);
}

// 简化算法实现
std::vector<cv::Point2f> simplify_raposo(
    const std::vector<cv::Point2f>& points,
    double hexradius,
    bool keep_first_last = true
) {
    std::vector<cv::Point2f> simplified_line;

    std::vector<double> x, y;
    for (const auto& pt : points) {
        x.push_back(pt.x);
        y.push_back(pt.y);
    }

    cHexGrid hexgrid(*std::max_element(y.begin(), y.end()), *std::min_element(x.begin(), x.end()),
        *std::min_element(y.begin(), y.end()), *std::max_element(x.begin(), x.end()), hexradius);

    std::vector<cv::Point> hexRC;
    for (const auto& pt : points) {
        hexRC.push_back(hexgrid.hexRowCol(pt.x, pt.y));
    }

    if (keep_first_last && hexRC[0] == hexRC[1]) {
        simplified_line.push_back(points[0]);
    }

    int startid = 0;
    for (size_t i = 0; i < points.size(); ++i) {
        if (hexRC[i] != hexRC[startid]) {
            simplified_line.push_back(quant_avg(points, startid, i));
            startid = i;
        }
    }

    simplified_line.push_back(quant_avg(points, startid, points.size()));

    if (keep_first_last && hexRC[hexRC.size() - 2] == hexRC[hexRC.size() - 1]) {
        simplified_line.push_back(points.back());
    }

    return simplified_line;
}

// 尝试简化到 n 个点
std::vector<cv::Point2f> raposo(
    const std::vector<cv::Point2f>& points,
    int n,
    int max_tries = 250,
    bool keep_first_last = true
)
{
    std::vector<double> x, y;
    for (const auto& pt : points) {
        x.push_back(pt.x);
        y.push_back(pt.y);
    }

    double dx = *std::max_element(x.begin(), x.end()) - *std::min_element(x.begin(), x.end());
    double dy = *std::max_element(y.begin(), y.end()) - *std::min_element(y.begin(), y.end());

    double low_radius = 0;
    double high_radius = std::sqrt(dx * dx + dy * dy);
    double pivot = (low_radius + high_radius) / 2.0;

    std::vector<cv::Point2f> simplified_line = simplify_raposo(points, pivot, keep_first_last);
    int dn = std::abs(static_cast<int>(simplified_line.size()) - n);
    int numtries = 0;

    while (dn > 0 && numtries <= max_tries) {
        if (simplified_line.size() < n) {
            high_radius = pivot;
        }
        else {
            low_radius = pivot;
        }
        pivot = (low_radius + high_radius) / 2.0;
        simplified_line = simplify_raposo(points, pivot, keep_first_last);
        int new_dn = std::abs(static_cast<int>(simplified_line.size()) - n);
        if (new_dn < dn) {
            numtries = 0;
        }
        else {
            numtries += 1;
        }
        dn = new_dn;
    }

    return simplified_line;
}

int main0()
{
    std::vector<cv::Point2f> pts;
    pts.push_back(cv::Point2f(47.2f, 86.4f));
    pts.push_back(cv::Point2f(63.9f, 103.3f));
    pts.push_back(cv::Point2f(63.9f, 120.0f));
    pts.push_back(cv::Point2f(63.9f, 136.8f));
    pts.push_back(cv::Point2f(63.8f, 153.6f));
    pts.push_back(cv::Point2f(79.6f, 170.4f));
    pts.push_back(cv::Point2f(79.6f, 187.2f));
    pts.push_back(cv::Point2f(96.3f, 202.9f));
    pts.push_back(cv::Point2f(96.2f, 219.7f));
    pts.push_back(cv::Point2f(112.9f, 236.5f));
    pts.push_back(cv::Point2f(129.5f, 253.3f));
    pts.push_back(cv::Point2f(145.3f, 270.2f));
    pts.push_back(cv::Point2f(145.3f, 286.9f));
    pts.push_back(cv::Point2f(162.0f, 302.6f));
    pts.push_back(cv::Point2f(178.7f, 319.4f));
    pts.push_back(cv::Point2f(194.5f, 336.3f));
    pts.push_back(cv::Point2f(211.1f, 353.1f));
    pts.push_back(cv::Point2f(227.8f, 369.9f));
    pts.push_back(cv::Point2f(243.6f, 386.7f));
    pts.push_back(cv::Point2f(243.5f, 403.5f));

    int numpts = 10; // Target number of points
    double hexradius = 1.0; // Hexagon radius

    std::vector<cv::Point2f> simplified_pts = raposo(pts, numpts, 25, true);

    std::cout << "Simplified Points:" << std::endl;
    for (const auto& pt : simplified_pts) {
        std::cout << "(" << pt.x << ", " << pt.y << ")" << std::endl;
    }

    return 0;
}

