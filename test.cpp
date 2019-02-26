//@note 利用opencv画图测试曲线拟合是否成功

#include <opencv2/opencv.hpp>
#include "PolyFit.hpp"

double GetY(double x, const int poly_n, std::vector<double> &A) {
	double sum = 0.0;
	for (int i = 0; i <= poly_n; ++i) {
		sum += A[i]*std::pow(x, i);
	}
	return sum;
}

void GetCurvePoints(const std::vector<cv::Point> &points, 
	std::vector<cv::Point> &curvepoints, const int poly_n, std::vector<double> factor) {
	for (int x = 0; x < 480; x++) {
		curvepoints.push_back(cv::Point(x, GetY(x, poly_n, factor)));
	}
}

//opencv中画曲线只支持points也就是按像素点，因此只支持整型
void PrintImage(const std::vector<cv::Point> &points, std::vector<cv::Point> &curvepoints) {
	cv::Mat image = cv::Mat::zeros(480, 720, CV_8UC3);
	image.setTo(cv::Scalar(100, 0, 0));

	for (int i = 0; i < points.size(); ++i) {
		cv::circle(image, points[i], 5, cv::Scalar(0, 0, 255), 2, 8, 0);
	}
	cv::polylines(image, points, false, cv::Scalar(0, 255, 0), 1, 8, 0);
	cv::polylines(image, curvepoints, false, cv::Scalar(0, 255, 255), 1, 8, 0);
	cv::imshow("CurveFitting", image);
	cv::waitKey(0);
}

int main() {
	std::vector<cv::Point> keypoints;
	keypoints.push_back(cv::Point(30, 58));
	keypoints.push_back(cv::Point(50, 70));
	keypoints.push_back(cv::Point(120, 90));
	keypoints.push_back(cv::Point(250, 140));
	keypoints.push_back(cv::Point(330, 220));
	keypoints.push_back(cv::Point(450, 500));

	std::vector<curvefit::Point<double>> points;
	points.push_back(curvefit::Point<double>{30, 58});
	points.push_back(curvefit::Point<double>{50, 70});
	points.push_back(curvefit::Point<double>{120, 90});
	points.push_back(curvefit::Point<double>{250, 140});
	points.push_back(curvefit::Point<double>{330, 220});
	points.push_back(curvefit::Point<double>{450, 500});

	curvefit::PolyFit t;
	t.Fit(points, 5);
	std::vector<double> factor= t.GetFactor();
	for (int i = 0; i <factor.size(); i++) {
		std::cout << factor[i] << std::endl;
	}

	std::vector<cv::Point> curvepoints;
	GetCurvePoints(keypoints, curvepoints, 5, factor);
	PrintImage(keypoints, curvepoints);
	return 0;
}
