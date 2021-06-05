#include <pcl/io/pcd_io.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/segmentation/progressive_morphological_filter.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/visualization/cloud_viewer.h>
#include <iostream>
#include<time.h>
#include<vector>
#include<string>
#include<io.h>

# define r_max_square 50*50
# define r_min_square 0.5*0.5

using namespace std;

void getFiles(string path, vector<string>& files)
{
	//文件句柄
	intptr_t hFile = 0;
	//文件信息
	struct _finddata_t fileinfo;
	string p;
	if ((hFile = _findfirst(p.assign(path).append("\\*").c_str(), &fileinfo)) != -1)
	{
		do
		{
			//如果是目录,迭代之
			//如果不是,加入列表
			if ((fileinfo.attrib &  _A_SUBDIR))
			{
				if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
					getFiles(p.assign(path).append("\\").append(fileinfo.name), files);
			}
			else
			{
				files.push_back(p.assign(path).append("\\").append(fileinfo.name));
			}
		} while (_findnext(hFile, &fileinfo) == 0);
		_findclose(hFile);
	}
}

void main()
{
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_input(new pcl::PointCloud<pcl::PointXYZRGB>);
	pcl::ExtractIndices<pcl::PointXYZRGB> extract;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr object_cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
	pcl::PCDReader reader;
	vector<string> d_dir;
	vector<double>acc;
	vector<double>tpr;
	vector<double>fpr;
	vector<double>iou1;
	vector<double>iou2;
	vector<double>f1;
	vector<double>run_time;
	double tem_acc;
	double tem_tpr;
	double tem_fpr;
	double tem_iou1;
	double tem_iou2;
	double tem_f1;
	double tem_time;
	clock_t t1, t2;
	
	int TP, FN, FP, TN;
	const char * filePath = "E:\\系统浏览器下载\\pc_pcd";
	getFiles(filePath, d_dir);
	for (int t = 0; t < d_dir.size(); t++) {
		reader.read(d_dir[t], *cloud_input);
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_test(new pcl::PointCloud<pcl::PointXYZRGB>);
		for (int i = 0; i < cloud_input->size(); ++i) {
			if ((cloud_input->points[i].x*cloud_input->points[i].x + cloud_input->points[i].y*cloud_input->points[i].y) < r_max_square && (cloud_input->points[i].x*cloud_input->points[i].x + cloud_input->points[i].y*cloud_input->points[i].y) > r_min_square) {
				cloud_test->push_back(cloud_input->points[i]);
			}
		}
		TP = 0, FN = 0, FP = 0, TN = 0;
		t1 = clock();
		pcl::ModelCoefficients::Ptr coefficients(new pcl::ModelCoefficients);
		pcl::PointIndices::Ptr inliers(new pcl::PointIndices);
		pcl::SACSegmentation<pcl::PointXYZRGB> seg;
		seg.setOptimizeCoefficients(true);
		seg.setModelType(pcl::SACMODEL_PLANE);
		seg.setMethodType(pcl::SAC_RANSAC);
		seg.setDistanceThreshold(0.15);
		seg.setInputCloud(cloud_test);
		seg.segment(*inliers, *coefficients);
		t2 = clock();

		for (int j = 0; j < inliers->indices.size(); ++j) {
			if (cloud_test->points.at(inliers->indices[j]).g != 0) {
				TP++;
			}
			else FP++;
		}
		
		extract.setInputCloud(cloud_test);
		extract.setIndices(inliers);
		extract.setNegative(true);
		extract.filter(*object_cloud);
		for (int i = 0; i < object_cloud->size(); ++i)
		{
			if (object_cloud->points[i].r != 0) {
				TN++;
			}
			else FN++;
		}
		tem_tpr = long double(TP) / (long double(TP) + long double(FN));
		tem_fpr = long double(FP) / (long double(FP) + long double(TN));
		tem_acc = (long double(TP) + long double(TN)) / (long double(TP) + long double(TN) + long double(FP) + long double(FN));
		tem_iou1 = double(TP) / double(TP + FP + FN);
		tem_iou2 = double(TN) / double(TN + FP + FN);
		tem_f1 = 2 * double(TP) / (2 * double(TP) + double(FP) + double(FN));
		tem_time = (double)(t2 - t1);
		tpr.push_back(tem_tpr);
		fpr.push_back(tem_fpr);
		acc.push_back(tem_acc);
		iou1.push_back(tem_iou1);
		iou2.push_back(tem_iou2);
		f1.push_back(tem_f1);
		run_time.push_back(tem_time);
		if (t%10==0) {
			cout << "已经完成：" << t << "帧数据" << endl;
			cout << "tpr:" << tem_tpr << endl;
			cout << "fpr:" << tem_fpr << endl;
			cout << "acc:" << tem_acc << endl;
			cout << "iou1:" << tem_iou1 << endl;
			cout << "iou2:" << tem_iou2 << endl;
			cout << "f1:" << tem_f1 << endl;
			cout << "run time:" << tem_time << endl;
		}
	}
	double acc_sum = accumulate(begin(acc), end(acc), 0.0);
	double acc_mean = acc_sum / acc.size();
	double tpr_sum = accumulate(begin(tpr), end(tpr), 0.0);
	double tpr_mean = tpr_sum / tpr.size();
	double fpr_sum = accumulate(begin(fpr), end(fpr), 0.0);
	double fpr_mean = fpr_sum / fpr.size();
	double iou1_sum = accumulate(begin(iou1), end(iou1), 0.0);
	double iou1_mean = iou1_sum / iou1.size();
	double iou2_sum = accumulate(begin(iou2), end(iou2), 0.0);
	double iou2_mean = iou2_sum / iou2.size();
	double f1_sum = accumulate(begin(f1), end(f1), 0.0);
	double f1_mean = f1_sum / f1.size();
	double time_sum = accumulate(begin(run_time), end(run_time), 0.0);
	double time_mean = time_sum / run_time.size();
	cout << "最终结果如下：" << endl;
	cout << "tpr:" << tpr_mean << endl;
	cout << "fpr:" << fpr_mean << endl;
	cout << "acc:" << acc_mean << endl;
	cout << "iou1:" << iou1_mean << endl;
	cout << "iou2:" << iou2_mean << endl;
	cout << "f1:" << f1_mean << endl;
	cout << "run time" << time_mean << endl;
}