#include<pcl/io/pcd_io.h>
#include<iostream>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/point_types.h>
#include<pcl/point_cloud.h>
#include<vector>
#include <mutex>
#include <atomic>
#include <list>
#include <limits>
#include <cmath>
#include <memory>
#include <pcl/segmentation/extract_clusters.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/features/normal_3d.h>
#include <pcl/kdtree/kdtree.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/features/moment_of_inertia_estimation.h>
#include <pcl/ModelCoefficients.h>
#include <time.h>
#include <algorithm>

# define cj (numeric_limits<double>::max)()
using namespace std;


# define n_bin 38
# define z_diff 0.15
# define n_segment 180


# define r_max_square_2 50*50
# define point_to_line 0.15
# define sensor_hight 1.8
# define end_point 0.1
# define slope 0.3

# define min_angle 65.2
# define f_angle 0.6
#define PI 3.14159265358979323846

class Bin {

public:
	double max_z;
	double min_z;
	double min_z_range;
	int id;
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
		if (z < min_z) {
			min_z = z;
			min_z_range = d;
		}
		if (z > max_z) {
			max_z = z;
		}
	}
	bool hasPoint() { return has_point_; }
	void init() {
		max_z = (numeric_limits<double>::min)();
		min_z = cj;
	}

};

class Segment {
public:
	typedef std::pair<double, double> LocalLine;
	std::vector<Bin> bins_;
	std::list<LocalLine> lines_;
	std::vector<int> seed_inx_;
	std::vector<int> break_point_;
public:
	Segment() {
		bins_.resize(n_bin);
	};

	bool point_line(double d, double z, int bin_num) {
		if (break_point_.empty()) {
			if (lines_.empty()) {
				return true;
			}
			else {
				double hight = lines_.begin()->first*d + lines_.begin()->second;
				if (abs(z - hight) < point_to_line) {
					return false;
				}
				else return true;
			}
		}
		else {

			for (int i = 0; i < break_point_.size(); i++) {
				if (seed_inx_[break_point_[i]] >= bin_num) {
					list <LocalLine>::iterator it = lines_.begin();
					int j = 0;
					for (int j = 0; j < i; j++) {
						if (it != lines_.end()) {
							++it;
						}
					}

					double hight = it->first*d + it->second;
					if (abs(z - hight) < point_to_line) {
						return false;
					}
					else return true;
				}

			}

			double hight = lines_.back().first*d + lines_.back().second;
			if (abs(z - hight) < point_to_line) {
				return false;
			}
			else return true;

		}
	}


	LocalLine fitLocalLine(const std::vector<Bin>& points) {
		const unsigned int n_points = points.size();
		//printf("%d\n", n_points);
		if (n_points >= 2) {
			Eigen::MatrixXd X(n_points, 2);
			Eigen::VectorXd Y(n_points);
			unsigned int counter = 0;
			for (int i = 0; i < points.size(); ++i) {
				X(counter, 0) = points[i].min_z_range;
				X(counter, 1) = 1;
				Y(counter) = points[i].min_z;
				++counter;
			}
			Eigen::VectorXd result = X.colPivHouseholderQr().solve(Y);
			LocalLine line_result;
			line_result.first = result(0);
			line_result.second = result(1);
			return line_result;
		}
		else {
			LocalLine line_result;
			line_result.first = 0;
			line_result.second = 0;
			return line_result;
		}
	}

	//����������ֱ�ߺ���

	int subsection(const std::vector<Bin>& line, LocalLine &sub_lines_tem) {
		double b = line[0].min_z - sub_lines_tem.first*line[0].min_z_range;
		for (int i = 0; i < line.size(); ++i) {
			if (abs(line[i].min_z - (sub_lines_tem.first*line[i].min_z_range + b)) > end_point) {
				return i;
			}
		}
		return 0;
	}

	std::vector<Bin> get_seed_points() {
		std::vector<Bin> current_line_points;
		for (int i = 0; i < bins_.size(); ++i) {
			if (current_line_points.empty()) {
				if (abs(bins_[i].min_z + sensor_hight) < z_diff + 0.15) {
					current_line_points.push_back(bins_[i]);
					seed_inx_.push_back(i);
				}
			}
			else {
				if (abs((bins_[i].min_z - current_line_points.back().min_z) / (bins_[i].min_z_range - current_line_points.back().min_z_range)) < slope) {
					current_line_points.push_back(bins_[i]);
					seed_inx_.push_back(i);
				}
			}
		}
		return current_line_points;
	}

	void fitlines(std::vector<Bin> seed) {

		LocalLine lines_tem = fitLocalLine(seed);
		lines_.push_back(lines_tem);
		std::vector<Bin> seed_next1;
		std::vector<Bin> seed_next2;
		int middle_point = subsection(seed, lines_tem);
		if (break_point_.empty()) {
			if (middle_point != 0 && middle_point != n_bin) {
				break_point_.push_back(middle_point);
			}
		}
		else {
			if (middle_point != 0 && middle_point != n_bin && middle_point > 1) {
				break_point_.push_back(middle_point + break_point_.back());
			}
		}

		if (middle_point != 0) {
			lines_.pop_back();
			for (int j = 0; j < middle_point; j++) {
				seed_next1.push_back(seed[j]);
			}
			for (int i = middle_point; i < seed.size(); i++) {
				seed_next2.push_back(seed[i]);
			}
			if (seed_next1.size() >= 2) {
				LocalLine lines_tem_1 = fitLocalLine(seed_next1);
				lines_.push_back(lines_tem_1);
			}
			if (seed_next2.size() > 2) {
				fitlines(seed_next2);
			}
			else {
				break_point_.pop_back();
			}
		}
	}


	void fitSegmentLines() {
		std::vector<Bin> current_line_points = get_seed_points();
		if (current_line_points.size() > 0) {
			fitlines(current_line_points);
		}
	}

};

class GroundSegmentation {

public:
	GroundSegmentation() {
		segments_.resize(n_segment);
	}
	vector<Segment> segments_;
	vector<std::pair<int, int> > bin_index_;

	int range_number(double range) {
		return int((atan(range / sensor_hight) * 180 / PI - min_angle) / f_angle);
	}

	void insertPoints(boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>> &cloud) {
		bin_index_.resize(cloud->size());
		double segment_step = 2 * M_PI / n_segment;

		for (int i = 0; i < cloud->size(); ++i) {
			double range_square = cloud->points[i].x * cloud->points[i].x + cloud->points[i].y * cloud->points[i].y;
			double range = sqrt(range_square);
			unsigned int bin_index = range_number(range);
			unsigned int segment_index = floor((std::atan2(cloud->points[i].y, cloud->points[i].x) + M_PI) / (2 * M_PI / n_segment));
			if (bin_index > n_bin - 1) { bin_index = n_bin - 1; }
			if (segment_index > n_segment - 1) { segment_index = segment_index % n_segment; }
			segments_[segment_index].bins_[bin_index].addPoint(cloud->points[i]);
			bin_index_[i] = std::make_pair(segment_index, bin_index);

		}

	}
	void lineFit() {
		for (int i = 0; i < n_segment; ++i) {
			segments_[i].fitSegmentLines();
		}
	}
	boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>> object_output(boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>> &cloud) {
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_object(new pcl::PointCloud<pcl::PointXYZRGB>);
		//pcl::PointXYZ q;
		const double segment_step = 2 * M_PI / n_segment;
		for (int i = 0; i < cloud->size(); ++i) {
			double angle = std::atan2(cloud->points[i].y, cloud->points[i].x);
			unsigned int segment_index = (angle + M_PI) / segment_step;
			double range_square = cloud->points[i].x * cloud->points[i].x + cloud->points[i].y * cloud->points[i].y;
			double range = sqrt(range_square);
			unsigned int bin_index = range_number(range);
			if (bin_index > n_bin - 1) { bin_index = n_bin - 1; }
			if (segment_index > n_segment - 1) { segment_index = segment_index % n_segment; }


			if (segments_[segment_index].lines_.empty()) {
				cloud_object->push_back(cloud->points[i]);
			}
			else if (segments_[segment_index].lines_.begin()->first == 0 && segments_[segment_index].lines_.begin()->second == 0) {
				cloud_object->push_back(cloud->points[i]);
			}
			else {
				if (segments_[segment_index].point_line(range, cloud->points[i].z, bin_index)) {
					cloud_object->push_back(cloud->points[i]);
				}
				else {
					//cloud_ground->push_back(cloud->points[i]);
				}
			}
		}
		return cloud_object;
	}
	void accuracy(boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>> &cloud, double&tem_tpr, double&tem_fpr, double&tem_acc, double&tem_iou1, double&tem_iou2, double&tem_f1) {
		//pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_ground(new pcl::PointCloud<pcl::PointXYZ>);
		const double segment_step = 2 * M_PI / n_segment;
		int TP = 0;
		int FP = 0;
		int FN = 0;
		int TN = 0;
		for (int i = 0; i < cloud->size(); ++i) {
			double angle = std::atan2(cloud->points[i].y, cloud->points[i].x);
			unsigned int segment_index = (angle + M_PI) / segment_step;
			double range_square = cloud->points[i].x * cloud->points[i].x + cloud->points[i].y * cloud->points[i].y;
			double range = sqrt(range_square);

			unsigned int bin_index = range_number(range);

			if (bin_index > n_bin - 1) { bin_index = n_bin - 1; }
			if (segment_index > n_segment - 1) { segment_index = segment_index % n_segment; }
			if (segments_[segment_index].lines_.empty()) {
				
				if (cloud->points[i].r != 0) {
					TN++;
				}
				else {
					FN++;
				}
				
			}
			else if (segments_[segment_index].lines_.begin()->first == 0 && segments_[segment_index].lines_.begin()->second == 0) {
				
				if (cloud->points[i].r != 0) {
					TN++;
				}
				else {
					FN++;
				}
				
			}
			else {
				if (segments_[segment_index].point_line(range, cloud->points[i].z, bin_index)) {
					if (cloud->points[i].r != 0) {
						TN++;
					}
					else {
						FN++;
					}
				}
				else {
					if (cloud->points[i].g != 0) {
						TP++;
					}
					else {
						FP++;
					}
				}
			}

		}
		tem_tpr = long double(TP) / (long double(TP) + long double(FN));
		tem_fpr = long double(FP) / (long double(FP) + long double(TN));
		tem_acc = (long double(TP) + long double(TN)) / (long double(TP) + long double(TN) + long double(FP) + long double(FN));
		tem_iou1 = double(TP) / double(TP + FP + FN);
		tem_iou2 = double(TN) / double(TN + FP + FN);
		tem_f1 = 2 * double(TP) / (2 * double(TP) + double(FP) + double(FN));
	}

	boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>> cloud_output(boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>> &cloud) {
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_view_test(new pcl::PointCloud<pcl::PointXYZRGB>);
		const double segment_step = 2 * M_PI / n_segment;
		for (int i = 0; i < cloud->size(); ++i) {
			double angle = std::atan2(cloud->points[i].y, cloud->points[i].x);
			unsigned int segment_index = (angle + M_PI) / segment_step;
			double range_square = cloud->points[i].x * cloud->points[i].x + cloud->points[i].y * cloud->points[i].y;
			double range = sqrt(range_square);

			unsigned int bin_index = range_number(range);

			if (bin_index > n_bin - 1) { bin_index = n_bin - 1; }
			if (segment_index > n_segment - 1) { segment_index = segment_index % n_segment; }

			pcl::PointXYZRGB p;
			if (segments_[segment_index].lines_.empty()) {
				p.x = cloud->points[i].x;
				p.y = cloud->points[i].y;
				p.z = cloud->points[i].z;
				p.r = 255;
				p.g = 0;
				p.b = 0;
				cloud_view_test->push_back(p);
			}
			else if (segments_[segment_index].lines_.begin()->first == 0 && segments_[segment_index].lines_.begin()->second == 0) {
				p.x = cloud->points[i].x;
				p.y = cloud->points[i].y;
				p.z = cloud->points[i].z;
				p.r = 255;
				p.g = 0;
				p.b = 0;
				cloud_view_test->push_back(p);
			}
			else {
				if (segments_[segment_index].point_line(range, cloud->points[i].z, bin_index)) {
					p.x = cloud->points[i].x;
					p.y = cloud->points[i].y;
					p.z = cloud->points[i].z;
					p.r = 255;
					p.g = 0;
					p.b = 0;
					cloud_view_test->push_back(p);
				}
				else {
					p.x = cloud->points[i].x;
					p.y = cloud->points[i].y;
					p.z = cloud->points[i].z;
					p.r = 0;
					p.g = 255;
					p.b = 0;
					cloud_view_test->push_back(p);
				}
			}
		}
		return cloud_view_test;
	}


	boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>> ground_output(boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>> &cloud) {
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_ground(new pcl::PointCloud<pcl::PointXYZRGB>);
		const double segment_step = 2 * M_PI / n_segment;
		for (int i = 0; i < cloud->size(); ++i) {
			double angle = std::atan2(cloud->points[i].y, cloud->points[i].x);
			unsigned int segment_index = (angle + M_PI) / segment_step;
			double range_square = cloud->points[i].x * cloud->points[i].x + cloud->points[i].y * cloud->points[i].y;
			double range = sqrt(range_square);

			unsigned int bin_index = range_number(range);

			if (bin_index > n_bin - 1) { bin_index = n_bin - 1; }
			if (segment_index > n_segment - 1) { segment_index = segment_index % n_segment; }


			if (segments_[segment_index].lines_.empty()) {
			}
			else if (segments_[segment_index].lines_.begin()->first == 0 && segments_[segment_index].lines_.begin()->second == 0) {
			}
			else {
				if (segments_[segment_index].point_line(range, cloud->points[i].z, bin_index)) {
				}
				else {
					cloud_ground->push_back(cloud->points[i]);
				}
			}
		}
		return cloud_ground;
	}
};

void getFiles(string path, vector<string>& files)
{
	//�ļ����
	intptr_t hFile = 0;
	//�ļ���Ϣ
	struct _finddata_t fileinfo;
	string p;
	if ((hFile = _findfirst(p.assign(path).append("\\*").c_str(), &fileinfo)) != -1)
	{
		do
		{
			//�����Ŀ¼,����֮
			//�������,�����б�
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
	const char * filePath = "E:\\ϵͳ���������\\pc_pcd";
	clock_t t1, t2;
	getFiles(filePath, d_dir);
	for (int t = 0; t < d_dir.size(); t++) {
		reader.read(d_dir[t], *cloud_input);
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_test(new pcl::PointCloud<pcl::PointXYZRGB>);
		for (int i = 0; i < cloud_input->size(); ++i) {
			if ((cloud_input->points[i].x*cloud_input->points[i].x + cloud_input->points[i].y*cloud_input->points[i].y) < r_max_square_2) {
				cloud_test->push_back(cloud_input->points[i]);
			}
		}
		GroundSegmentation ground_segement;

		clock_t t1, t2;
		double tem_time;
		t1 = clock();
		ground_segement.insertPoints(cloud_test);
		ground_segement.lineFit();
		t2 = clock();
		tem_time = (double)(t2 - t1);
		double tem_acc;
		double tem_tpr;
		double tem_fpr;
		double tem_iou1;
		double tem_iou2;
		double tem_f1;
		ground_segement.accuracy(cloud_test, tem_tpr, tem_fpr, tem_acc, tem_iou1, tem_iou2, tem_f1);
		tpr.push_back(tem_tpr);
		fpr.push_back(tem_fpr);
		acc.push_back(tem_acc);
		iou1.push_back(tem_iou1);
		iou2.push_back(tem_iou2);
		f1.push_back(tem_f1);
		run_time.push_back(tem_time);
		if (tem_acc>0.95) {
			cout << "�Ѿ���ɣ�" << t << "֡����" << endl;
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
	cout << "���ս�����£�" << endl;
	cout << "tpr:" << tpr_mean << endl;
	cout << "fpr:" << fpr_mean << endl;
	cout << "acc:" << acc_mean << endl;
	cout << "iou1:" << iou1_mean << endl;
	cout << "iou2:" << iou2_mean << endl;
	cout << "f1:" << f1_mean << endl;
	cout << "run time" << time_mean << endl;
}