#include <iostream>
#include <vector>
#include "OpencvLibrary.h"

// 线性插值函数
std::vector<float> interp1d(const std::vector<float>& xp, const std::vector<float>& fp, const std::vector<float>& x)
{
    std::vector<float> result(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        auto it = std::lower_bound(xp.begin(), xp.end(), x[i]);
        if (it == xp.begin()) {
            result[i] = fp.front();
        }
        else if (it == xp.end()) {
            result[i] = fp.back();
        }
        else {
            size_t idx = it - xp.begin();
            float x0 = xp[idx - 1], x1 = xp[idx];
            float y0 = fp[idx - 1], y1 = fp[idx];
            result[i] = y0 + (x[i] - x0) * (y1 - y0) / (x1 - x0);
        }
    }
    return result;
}

// 重新采样函数
std::vector<cv::Point2f> ultralytics_interp(const std::vector<cv::Point2f>& segments, int n = 1000)
{
    std::vector<cv::Point2f> resampled_segments;
    // 将当前段加上第一个点，形成闭环
    std::vector<cv::Point2f> s = segments;
    s.push_back(segments[0]);

    // 生成插值坐标
    std::vector<float> x(n);
    std::vector<float> xp(s.size());
    for (int j = 0; j < n; ++j) {
        x[j] = j * (s.size() - 1) / static_cast<float>(n - 1);
    }
    for (size_t j = 0; j < s.size(); ++j) {
        xp[j] = j;
    }

    std::vector<float> s_x;
    std::vector<float> s_y;

    for (size_t j = 0; j < s.size(); j++)
    {
        s_x.push_back(s[j].x);
        s_y.push_back(s[j].y);
    }

    std::vector<float> interpolated_x = interp1d(xp, s_x, x);
    std::vector<float> interpolated_y = interp1d(xp, s_y, x);

    for (int j = 0; j < n; ++j)
    {
        resampled_segments.push_back(cv::Point2f(interpolated_x[j], interpolated_y[j]));
    }

    return resampled_segments;
}


