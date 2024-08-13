#include "douglas_peucker.h"


// 计算点到直线的垂直距离
float distancePtToLine(const cv::Point2f& pt, const cv::Point2f& line_start, const cv::Point2f& line_end)
{
    cv::Point2f a = line_start;
    cv::Point2f b = line_end;
    cv::Point2f x = pt;

    float distance;

    if (a == b) { // 如果线段的起点和终点重合
        distance = std::sqrt(std::pow(a.x - x.x, 2) + std::pow(a.y - x.y, 2)); // 计算点到线段起点的距离
    }
    else { // 一般情况
        // 计算线段向量 b-a 的长度
        float line_length = std::sqrt(std::pow(b.x - a.x, 2) + std::pow(b.y - a.y, 2));

        // 计算点到直线的垂直距离
        distance = std::abs((b.x - a.x) * (a.y - x.y) - (a.x - x.x) * (b.y - a.y)) / line_length;
    }

    return distance;
}

// 递归计算误差
float calcErrorsRecursive(const std::vector<cv::Point2f>& pts, int s, int e, std::vector<float>& errors)
{
    // 初始化
    float maxd = 0;
    int furthest_pt = -1;

    // 遍历点集
    for (int i = s + 1; i < e; ++i) {
        // 计算点到直线的垂直距离
        float d = distancePtToLine(pts[i], pts[s], pts[e]);

        // 更新最大距离和最远点
        if (d > maxd) {
            maxd = d;
            furthest_pt = i;
        }
    }

    // 检查是否找到了最远点
    if (furthest_pt != -1) {
        // 更新最远点的误差值
        errors[furthest_pt] = maxd;

        // 递归计算左侧和右侧的误差
        float maxleft = calcErrorsRecursive(pts, s, furthest_pt, errors);
        float maxright = calcErrorsRecursive(pts, furthest_pt, e, errors);

        // 确保当前最大距离至少与子点的误差值一样大
        maxd = std::max({ maxd, maxleft, maxright });

        // 更新最远点的误差值
        errors[furthest_pt] = maxd;

        // 返回最大误差值
        return maxd;
    }
    else {
        return -1.0f;
    }
}

// 计算每个点的容忍误差
std::vector<float> calcErrors(const std::vector<cv::Point2f>& pts)
{
    std::vector<float> errors(pts.size(), 0.0f);

    // 计算递归误差
    calcErrorsRecursive(pts, 0, pts.size() - 1, errors);

    // 找到最大误差值
    float maxe = *std::max_element(errors.begin(), errors.end());

    // 设置起点和终点的误差值
    if (!errors.empty()) {
        errors[0] = maxe + 1.0f;
        errors.back() = maxe + 1.0f;
    }

    return errors;
}

// 计算误差并返回排序后的误差
int getErrorsSortedErrors(const std::vector<cv::Point2f>& pts, std::vector<float>& errors, std::vector<float>& sorted_errors)
{
    // 计算误差
    errors = calcErrors(pts);

    // 创建排序后的误差副本
    sorted_errors = errors;

    // 对误差进行降序排序
    std::sort(sorted_errors.begin(), sorted_errors.end(), std::greater<float>());

    // 返回错误和排序后的错误
    return 0;
}

// 根据给定的容差简化点集，返回误差大于等于容差的点
std::vector<cv::Point2f> simplifyByDistanceTolerance(const std::vector<cv::Point2f>& pts, const std::vector<float>& errors, float tolerance)
{
    std::vector<cv::Point2f> result;

    // 遍历点集，选择误差大于或等于容差的点
    for (size_t i = 0; i < pts.size(); ++i)
    {
        if (errors[i] >= tolerance)
        {
            result.push_back(pts[i]);
        }
    }

    return result;
}


// 根据点数筛选点
std::vector<cv::Point2f> simplifyByNumPts(
    const std::vector<cv::Point2f>& pts,
    int numpts,
    const std::vector<float>& errors,
    const std::vector<float>& sorted_errors)
{
    // 获取给定数量的最大误差值作为容忍度
    float tolerance = sorted_errors[numpts - 1];

    // 使用容忍度筛选点
    return simplifyByDistanceTolerance(pts, errors, tolerance);
}

// Douglas-Peucker 简化算法
std::vector<cv::Point2f> douglas_peucker(const std::vector<cv::Point2f>& pts, int numpts)
{
    float tolerance = 0.1f;
    // 获取误差和排序后的误差
    std::vector<float> errors;
    std::vector<float> sorted_errors;
    getErrorsSortedErrors(pts, errors, sorted_errors);

    // 根据给定的点数简化点集
    std::vector<cv::Point2f> simplified_by_numpts = simplifyByNumPts(pts, numpts, errors, sorted_errors);

    //// 根据给定的误差容忍度简化点集
    //std::vector<cv::Point2f> simplified_by_tolerance = simplifyByDistanceTolerance(pts, errors, tolerance);

    //// 可以根据需要选择返回其中一个简化结果
    //// 这里返回两种简化结果的合并示例（可根据需求调整）
    //std::vector<cv::Point2f> result = simplified_by_numpts;
    //result.insert(result.end(), simplified_by_tolerance.begin(), simplified_by_tolerance.end());

    //// 去重（如果需要）
    //std::sort(result.begin(), result.end(), [](const cv::Point2f& a, const cv::Point2f& b) {
    //    return (a.x < b.x) || (a.x == b.x && a.y < b.y);
    //    });
    //result.erase(std::unique(result.begin(), result.end()), result.end());

    return simplified_by_numpts;
}


