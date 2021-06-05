#include<iostream>
#include<pcl/point_cloud.h>
#include<pcl/point_types.h>
#include<pcl/io/pcd_io.h>
#include<pcl/visualization/cloud_viewer.h>
#include<pcl/segmentation/sac_segmentation.h>
#include <pcl/filters/voxel_grid.h>
#include<vector>
#include <pcl/filters/statistical_outlier_removal.h>
#include<time.h>
#include<string>
#include<io.h>

# define cj (numeric_limits<double>::max)()
# define cj_min -20
# define n_bin 160
# define z_diff 0.2
# define n_segment 180
# define r_max_square 50*50
# define r_min_square 0.5*0.5
# define point_to_line 0.15
# define sensor_hight 1.8


using namespace std;


class Bin {

public:
	double max_z;
	double min_z;
	double min_z_range;
	bool is_ground;
	//double z_z;
private:
	bool has_point_;
public:
	Bin() :has_point_(false) { init(); }
	void addPoint(const pcl::PointXYZRGB &point) {
		const double d = sqrt(point.x * point.x + point.y * point.y);
		addPoint(d, point.z);
	}
	void addPoint(const double &d, const double &z) {
		has_point_ = true;
		//min_z = cj;
		//printf("初始min_z的值为：%f\n", min_z);
		if (z < min_z) {
			min_z = z;
			min_z_range = d;
			//printf("z的值为：%f\n", min_z);
		}
		if (z > max_z) {
			max_z = z;
		}
	}
	void diff_bin() {
		if (max_z - min_z < z_diff) { is_ground = true; }
		else { is_ground = false; }
	}
	bool hasPoint() { return has_point_; }
	void init() {
		max_z = cj_min;
		min_z = cj;
		is_ground = false;
	}

};


class segment {
public:
	std::vector<Bin> bins_;
	segment() { bins_.resize(n_bin); }

	void diff_segment() {
		for (int i = 0; i < n_bin; i++) {
			bins_[i].diff_bin();
			//cout<<"高度差："<<bins_[i].z_z<<"最高点："<<bins_[i].max_z<<"最低点："<<bins_[i].min_z<<endl;
		}
	}
};

class GroundSegmentation {

public:
	GroundSegmentation() {
		segments_.resize(n_segment);
	}
	vector<segment> segments_;
	vector<std::pair<int, int> > bin_index_;


	void insertPoints(boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>> &cloud) {
		bin_index_.resize(cloud->size());
		double segment_step = 2 * M_PI / n_segment;
		double bin_step = (sqrt(r_max_square) - sqrt(r_min_square)) / n_bin;
		double r_min = sqrt(r_min_square);

		for (int i = 0; i < cloud->size(); ++i) {
			double range_square = cloud->points[i].x * cloud->points[i].x + cloud->points[i].y * cloud->points[i].y;
			double range = sqrt(range_square);

			if (range_square <r_max_square && range_square > r_min_square) {
				//double angle = std::atan2(cloud->points[i].y, cloud->points[i].x);
				unsigned int bin_index = (range - r_min) / bin_step;
				unsigned int segment_index = floor((std::atan2(cloud->points[i].y, cloud->points[i].x) + M_PI) / (2 * M_PI / n_segment));
				//cout << "bin_index:" << bin_index << "----segment_index:" << segment_index << endl;
				if (bin_index > n_bin - 1) { bin_index = n_bin - 1; }
				if (segment_index > n_segment - 1) { segment_index = n_segment - 1; }
				segments_[segment_index].bins_[bin_index].addPoint(cloud->points[i]);
				bin_index_[i] = std::make_pair(segment_index, bin_index);
			}
			else {
				bin_index_[i] = std::make_pair<int, int>(-1, -1);
			}

		}

	}
	void diff_ground() {
		for (int i = 0; i < n_segment; ++i) {
			segments_[i].diff_segment();
		}
	}

	boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>> object_output(boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>> &cloud) {
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_object(new pcl::PointCloud<pcl::PointXYZRGB>);
		//pcl::PointXYZ q;
		double r_min = sqrt(r_min_square);
		for (int i = 0; i < cloud->size(); ++i) {
			//double angle = std::atan2(cloud->points[i].y, cloud->points[i].x);
			unsigned int segment_index = (std::atan2(cloud->points[i].y, cloud->points[i].x) + M_PI) / (2 * M_PI / n_segment);
			double range_square = cloud->points[i].x * cloud->points[i].x + cloud->points[i].y * cloud->points[i].y;
			double range = sqrt(range_square);
			double bin_step = (sqrt(r_max_square) - sqrt(r_min_square)) / n_bin;
			unsigned int bin_index = (range - r_min) / bin_step;
			if (segment_index > n_segment - 1) { segment_index = n_segment - 1; }
			if (bin_index > n_bin - 1) { bin_index = n_bin - 1; }
			if (!segments_[segment_index].bins_[bin_index].is_ground) {
				cloud_object->push_back(cloud->points[i]);
			}
		}
		return cloud_object;
	}

	boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>> ground_output(boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>> &cloud) {
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_ground(new pcl::PointCloud<pcl::PointXYZRGB>);
		const double segment_step = 2 * M_PI / n_segment;
		//pcl::PointXYZ q;
		double r_min = sqrt(r_min_square);
		for (int i = 0; i < cloud->size(); ++i) {
			unsigned int segment_index = (std::atan2(cloud->points[i].y, cloud->points[i].x) + M_PI) / (2 * M_PI / n_segment);
			double range_square = cloud->points[i].x * cloud->points[i].x + cloud->points[i].y * cloud->points[i].y;
			double range = sqrt(range_square);
			double bin_step = (sqrt(r_max_square) - sqrt(r_min_square)) / n_bin;
			unsigned int bin_index = (range - r_min) / bin_step;
			if (segment_index > n_segment - 1) { segment_index = n_segment - 1; }
			if (bin_index > n_bin - 1) { bin_index = n_bin - 1; }
			if (segments_[segment_index].bins_[bin_index].is_ground) {
				cloud_ground->push_back(cloud->points[i]);
			}
		}
		return cloud_ground;
	}
	void accuracy(boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>> &cloud, double&tem_tpr, double&tem_fpr, double&tem_acc, double&tem_iou1, double&tem_iou2, double&tem_f1) {
		//const double segment_step = 2 * M_PI / n_segment;
		int TP = 0;
		int FP = 0;
		int FN = 0;
		int TN = 0;
		const double segment_step = 2 * M_PI / n_segment;
		//pcl::PointXYZ q;
		double r_min = sqrt(r_min_square);
		for (int i = 0; i < cloud->size(); ++i) {
			unsigned int segment_index = (std::atan2(cloud->points[i].y, cloud->points[i].x) + M_PI) / (2 * M_PI / n_segment);
			double range_square = cloud->points[i].x * cloud->points[i].x + cloud->points[i].y * cloud->points[i].y;
			double range = sqrt(range_square);
			double bin_step = (sqrt(r_max_square) - sqrt(r_min_square)) / n_bin;
			unsigned int bin_index = (range - r_min) / bin_step;
			if (segment_index > n_segment - 1) { segment_index = n_segment - 1; }
			if (bin_index > n_bin - 1) { bin_index = n_bin - 1; }
			if (segments_[segment_index].bins_[bin_index].is_ground) {
				//cloud_ground->push_back(cloud->points[i]);
				if (cloud->points[i].g != 0) {
					TP++;
				}
				else FP++;
			}
			else {
				if (cloud->points[i].r != 0) {
					TN++;
				}
				else FN++;
			}
		}
		//cout << "TP：" << TP << "--FP：" << FP << "--FN：" << FN << "--TN：" << TN << endl;
		tem_tpr = long double(TP) / (long double(TP) + long double(FN));
		tem_fpr = long double(FP) / (long double(FP) + long double(TN));
		tem_acc = (long double(TP) + long double(TN)) / (long double(TP) + long double(TN) + long double(FP) + long double(FN));
		tem_iou1 = double(TP) / double(TP + FP + FN);
		tem_iou2 = double(TN) / double(TN + FP + FN);
		tem_f1 = 2 * double(TP) / (2 * double(TP) + double(FP) + double(FN));
		//cout << "TPR：" << long double(TP) / (long double(TP) + long double(FN)) << endl;
		//cout << "FPR：" << long double(FP) / (long double(FP) + long double(TN)) << endl;
		//cout << "精确度：" << (long double(TP) + long double(TN)) / (long double(TP) + long double(TN) + long double(FP) + long double(FN));
	}
};

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

void main() {
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_input(new pcl::PointCloud<pcl::PointXYZRGB>);
	pcl::PCDReader reader;
	vector<string> d_dir;
	vector<double>acc;
	vector<double>tpr;
	vector<double>fpr;
	vector<double>iou1;
	vector<double>iou2;
	vector<double>f1;
	vector<double>run_time;
	const char * filePath = "E:\\系统浏览器下载\\pc_pcd";
	clock_t t1, t2;
	getFiles(filePath, d_dir);
	for (int t = 0; t < d_dir.size(); t++) {
		reader.read(d_dir[t], *cloud_input);
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_test(new pcl::PointCloud<pcl::PointXYZRGB>);
		for (int i = 0; i < cloud_input->size(); ++i) {
			if ((cloud_input->points[i].x*cloud_input->points[i].x + cloud_input->points[i].y*cloud_input->points[i].y) < r_max_square && (cloud_input->points[i].x*cloud_input->points[i].x + cloud_input->points[i].y*cloud_input->points[i].y) > r_min_square) {
				cloud_test->push_back(cloud_input->points[i]);
			}
		}
		GroundSegmentation cloud_segment;
		double tem_time;
		t1 = clock();
		cloud_segment.insertPoints(cloud_test);
		cloud_segment.diff_ground();
		t2 = clock();
		tem_time= (double)(t2 - t1);
		double tem_acc;
		double tem_tpr;
		double tem_fpr;
		double tem_iou1;
		double tem_iou2;
		double tem_f1;
		cloud_segment.accuracy(cloud_test, tem_tpr, tem_fpr, tem_acc, tem_iou1, tem_iou2, tem_f1);
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
	double acc_sum= accumulate(begin(acc), end(acc), 0.0);
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