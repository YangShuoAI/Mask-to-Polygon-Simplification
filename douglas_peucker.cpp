#include "douglas_peucker.h"


// ����㵽ֱ�ߵĴ�ֱ����
float distancePtToLine(const cv::Point2f& pt, const cv::Point2f& line_start, const cv::Point2f& line_end)
{
    cv::Point2f a = line_start;
    cv::Point2f b = line_end;
    cv::Point2f x = pt;

    float distance;

    if (a == b) { // ����߶ε������յ��غ�
        distance = std::sqrt(std::pow(a.x - x.x, 2) + std::pow(a.y - x.y, 2)); // ����㵽�߶����ľ���
    }
    else { // һ�����
        // �����߶����� b-a �ĳ���
        float line_length = std::sqrt(std::pow(b.x - a.x, 2) + std::pow(b.y - a.y, 2));

        // ����㵽ֱ�ߵĴ�ֱ����
        distance = std::abs((b.x - a.x) * (a.y - x.y) - (a.x - x.x) * (b.y - a.y)) / line_length;
    }

    return distance;
}

// �ݹ�������
float calcErrorsRecursive(const std::vector<cv::Point2f>& pts, int s, int e, std::vector<float>& errors)
{
    // ��ʼ��
    float maxd = 0;
    int furthest_pt = -1;

    // �����㼯
    for (int i = s + 1; i < e; ++i) {
        // ����㵽ֱ�ߵĴ�ֱ����
        float d = distancePtToLine(pts[i], pts[s], pts[e]);

        // �������������Զ��
        if (d > maxd) {
            maxd = d;
            furthest_pt = i;
        }
    }

    // ����Ƿ��ҵ�����Զ��
    if (furthest_pt != -1) {
        // ������Զ������ֵ
        errors[furthest_pt] = maxd;

        // �ݹ���������Ҳ�����
        float maxleft = calcErrorsRecursive(pts, s, furthest_pt, errors);
        float maxright = calcErrorsRecursive(pts, furthest_pt, e, errors);

        // ȷ����ǰ�������������ӵ�����ֵһ����
        maxd = std::max({ maxd, maxleft, maxright });

        // ������Զ������ֵ
        errors[furthest_pt] = maxd;

        // ����������ֵ
        return maxd;
    }
    else {
        return -1.0f;
    }
}

// ����ÿ������������
std::vector<float> calcErrors(const std::vector<cv::Point2f>& pts)
{
    std::vector<float> errors(pts.size(), 0.0f);

    // ����ݹ����
    calcErrorsRecursive(pts, 0, pts.size() - 1, errors);

    // �ҵ�������ֵ
    float maxe = *std::max_element(errors.begin(), errors.end());

    // ���������յ�����ֵ
    if (!errors.empty()) {
        errors[0] = maxe + 1.0f;
        errors.back() = maxe + 1.0f;
    }

    return errors;
}

// �������������������
int getErrorsSortedErrors(const std::vector<cv::Point2f>& pts, std::vector<float>& errors, std::vector<float>& sorted_errors)
{
    // �������
    errors = calcErrors(pts);

    // ��������������
    sorted_errors = errors;

    // �������н�������
    std::sort(sorted_errors.begin(), sorted_errors.end(), std::greater<float>());

    // ���ش���������Ĵ���
    return 0;
}

// ���ݸ������ݲ�򻯵㼯�����������ڵ����ݲ�ĵ�
std::vector<cv::Point2f> simplifyByDistanceTolerance(const std::vector<cv::Point2f>& pts, const std::vector<float>& errors, float tolerance)
{
    std::vector<cv::Point2f> result;

    // �����㼯��ѡ�������ڻ�����ݲ�ĵ�
    for (size_t i = 0; i < pts.size(); ++i)
    {
        if (errors[i] >= tolerance)
        {
            result.push_back(pts[i]);
        }
    }

    return result;
}


// ���ݵ���ɸѡ��
std::vector<cv::Point2f> simplifyByNumPts(
    const std::vector<cv::Point2f>& pts,
    int numpts,
    const std::vector<float>& errors,
    const std::vector<float>& sorted_errors)
{
    // ��ȡ����������������ֵ��Ϊ���̶�
    float tolerance = sorted_errors[numpts - 1];

    // ʹ�����̶�ɸѡ��
    return simplifyByDistanceTolerance(pts, errors, tolerance);
}

// Douglas-Peucker ���㷨
std::vector<cv::Point2f> douglas_peucker(const std::vector<cv::Point2f>& pts, int numpts)
{
    float tolerance = 0.1f;
    // ��ȡ�������������
    std::vector<float> errors;
    std::vector<float> sorted_errors;
    getErrorsSortedErrors(pts, errors, sorted_errors);

    // ���ݸ����ĵ����򻯵㼯
    std::vector<cv::Point2f> simplified_by_numpts = simplifyByNumPts(pts, numpts, errors, sorted_errors);

    //// ���ݸ�����������̶ȼ򻯵㼯
    //std::vector<cv::Point2f> simplified_by_tolerance = simplifyByDistanceTolerance(pts, errors, tolerance);

    //// ���Ը�����Ҫѡ�񷵻�����һ���򻯽��
    //// ���ﷵ�����ּ򻯽���ĺϲ�ʾ�����ɸ������������
    //std::vector<cv::Point2f> result = simplified_by_numpts;
    //result.insert(result.end(), simplified_by_tolerance.begin(), simplified_by_tolerance.end());

    //// ȥ�أ������Ҫ��
    //std::sort(result.begin(), result.end(), [](const cv::Point2f& a, const cv::Point2f& b) {
    //    return (a.x < b.x) || (a.x == b.x && a.y < b.y);
    //    });
    //result.erase(std::unique(result.begin(), result.end()), result.end());

    return simplified_by_numpts;
}


