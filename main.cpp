#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include "OpencvLibrary.h"

#include "apsc.h"
#include "raposo.h"
#include "visvalingam.h"
#include "douglas_peucker.h"
#include "ultralytics_interp.h"

int main()
{
    std::fstream fs;
    fs.open("./sample_pts.txt", std::ios::in);
    
    // 检查文件是否成功打开
    if (!fs.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        return 1; // 返回非零值表示错误
    }
    std::string line;

    std::vector<cv::Point2f> sample_pts;
    // 逐行读取文件内容
    while (std::getline(fs, line))
    {
        //std::cout << line << std::endl;
        std::stringstream ss;
        ss << line;
        int x, y;
        ss >> x >> y;
        std::cout << x << " " << y << std::endl;
        sample_pts.push_back(cv::Point2f(x, y));
    }
    // 关闭文件
    fs.close();
    int unsample_num = 32;
    std::vector<cv::Point2f> upsample_pts1, upsample_pts2, upsample_pts3, upsample_pts4, upsample_pts5;
    upsample_pts1 = apsc(sample_pts, unsample_num);
    upsample_pts2 = douglas_peucker(sample_pts, unsample_num);
    upsample_pts3 = raposo(sample_pts, unsample_num);
    upsample_pts4 = visvalingam(sample_pts, unsample_num);
    upsample_pts5 = ultralytics_interp(sample_pts, unsample_num);

    cv::Mat img_show = cv::Mat::zeros(cv::Size(360, 640), CV_8UC3);

    std::vector<cv::Point> upsample_output_pts1, upsample_output_pts2, upsample_output_pts3, upsample_output_pts4, upsample_output_pts5;
    for (size_t i = 0; i < upsample_pts1.size(); i++)
    {
        upsample_output_pts1.push_back(cv::Point(upsample_pts1[i].x, upsample_pts1[i].y));
        upsample_output_pts2.push_back(cv::Point(upsample_pts2[i].x, upsample_pts2[i].y));
        upsample_output_pts3.push_back(cv::Point(upsample_pts3[i].x, upsample_pts3[i].y));
        upsample_output_pts4.push_back(cv::Point(upsample_pts4[i].x, upsample_pts4[i].y));
        upsample_output_pts5.push_back(cv::Point(upsample_pts5[i].x, upsample_pts5[i].y));
    }

    cv::Scalar color0 = cv::Scalar(0, 255, 0);
    cv::Scalar color1 = cv::Scalar(0, 255, 255);
    cv::Scalar color2 = cv::Scalar(0, 0, 255);
    cv::Scalar color3 = cv::Scalar(255, 0, 255);
    cv::Scalar color4 = cv::Scalar(255, 255, 0);
    cv::Scalar color5 = cv::Scalar(255, 0, 0);

    //cv::drawContours(img, contours, contour_index, cv::Scalar(0, 255, 0), 1);
    //cv::drawContours(img, upsample_output_pts, -1, cv::Scalar(0, 0, 255), 1);

    std::vector<cv::Point> pts;
    std::copy(sample_pts.begin(), sample_pts.end(), std::back_inserter(pts));
    cv::polylines(img_show, pts, true, color0, 1, 8, 0);



    cv::polylines(img_show, upsample_output_pts1, true, color1, 1, 8, 0);
    cv::polylines(img_show, upsample_output_pts2, true, color2, 1, 8, 0);
    cv::polylines(img_show, upsample_output_pts3, true, color3, 1, 8, 0);
    cv::polylines(img_show, upsample_output_pts4, true, color4, 1, 8, 0);
    cv::polylines(img_show, upsample_output_pts5, true, color5, 1, 8, 0);

    // 字体类型
    int fontFace = cv::FONT_HERSHEY_COMPLEX;

    // 字体缩放因子
    double fontScale = 0.5;

    // 文本颜色 (B, G, R)

    // 线条粗细
    int thickness = 1;
    // 绘制文本
    std::string text0 = "gt";
    std::string text1 = "apsc";
    std::string text2 = "douglas_peucker";
    std::string text3 = "raposo";
    std::string text4 = "visvalingam";
    std::string text5 = "ultralytics_interp";

    cv::Point pt0 = cv::Point(img_show.cols * 0.5, 10);
    cv::Point pt1 = cv::Point(img_show.cols * 0.5, 30);
    cv::Point pt2 = cv::Point(img_show.cols * 0.5, 50);
    cv::Point pt3 = cv::Point(img_show.cols * 0.5, 70);
    cv::Point pt4 = cv::Point(img_show.cols * 0.5, 90);
    cv::Point pt5 = cv::Point(img_show.cols * 0.5, 110);

    cv::putText(img_show, text0, pt0, fontFace, fontScale, color0, thickness);
    cv::putText(img_show, text1, pt1, fontFace, fontScale, color1, thickness);
    cv::putText(img_show, text2, pt2, fontFace, fontScale, color2, thickness);
    cv::putText(img_show, text3, pt3, fontFace, fontScale, color3, thickness);
    cv::putText(img_show, text4, pt4, fontFace, fontScale, color4, thickness);
    cv::putText(img_show, text5, pt5, fontFace, fontScale, color5, thickness);


    cv::namedWindow("image", cv::WINDOW_NORMAL);
    cv::imshow("image", img_show);
    //cv::imwrite(save_file_path + filenames[k], img_show);
    cv::waitKey(0);

    return 0;
}






