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

# define r_1 0.2
# define r_2 0.5
# define r_max_square_1 20*20
# define r_max_square_2 50*50
# define cj 999
# define n_bin 160
# define z_diff 0.25
# define n_segment 180



using namespace std;
struct GroundSegmentationParams {
	GroundSegmentationParams() :
		visualize(true),
		radius(0.5),
		//r_min_square(0.4 * 0.4),
		//r_max_square(20.0 * 20.0),
		n_bins(160),
		n_segments(180),
		max_dist_to_line(0.15),
		max_slope(1),
		max_error_square(0.01),
		long_threshold(1.0),
		max_long_height(0.1),
		max_start_height(0.2),
		sensor_height(1.8),
		line_search_angle(0.5) {}

	bool visualize;
	float radius;
	int n_bins;
	int n_segments;
	double max_dist_to_line;
	double max_slope;
	double max_error_square;
	double long_threshold;
	double max_long_height;
	double max_start_height;
	double sensor_height;
	double line_search_angle;
};

class Bin {
public:
	struct MinZPoint {
		MinZPoint() : z(cj), d(cj) {}
		MinZPoint(double d_input, double z_input) : z(z_input), d(d_input) {}
		double z;
		double d;
	};

private:

	bool has_point_;
	double min_z;
	double min_z_range;

public:
	Bin() :has_point_(false) { min_z = cj; }
	void addPoint(const pcl::PointXYZ &point) {
		const double d = sqrt(point.x * point.x + point.y * point.y);
		addPoint(d, point.z);
	}

	void addPoint(const double &d, const double &z) {
		has_point_ = true;
		if (z < min_z) {
			min_z = z;
			min_z_range = d;
		}
	}

	MinZPoint getMinZPoint() {
		MinZPoint point;
		if (has_point_) {
			point.z = min_z;
			point.d = min_z_range;
		}
		return point;
	}
	bool hasPoint() { return has_point_; }

};

class Segment {
public:
	typedef std::pair<Bin::MinZPoint, Bin::MinZPoint> Line;
	typedef std::pair<double, double> LocalLine;
	std::vector<Bin> bins_;
	std::list<Line> lines_;

private:
	double max_slope_;
	double max_error_;
	double long_threshold_;
	double max_long_height_;
	double max_start_height_;
	double sensor_height_;

public:
	const GroundSegmentationParams params_;
	Segment() {
		int bins_ = params_.n_bins;
		max_slope_ = params_.max_slope;
		max_error_ = params_.max_error_square;
		long_threshold_ = params_.long_threshold;
		max_long_height_ = params_.max_long_height;
		max_start_height_ = params_.max_start_height;
		sensor_height_ = params_.sensor_height;
	};

	double getMaxError(const std::list<Bin::MinZPoint>& points, const LocalLine& line) {
		double max_error = 0;
		for (auto it = points.begin(); it != points.end(); ++it) {
			double residual = (line.first * it->d + line.second) - it->z;
			double error = residual * residual;
			if (error > max_error) max_error = error;
		}
		return max_error;
	}

	double getMeanError(const std::list<Bin::MinZPoint>& points, const LocalLine& line) {
		double error_sum = 0;
		for (auto it = points.begin(); it != points.end(); ++it) {
			double residual = (line.first * it->d + line.second) - it->z;
			error_sum += residual * residual;
		}
		return error_sum / points.size();
	}

	Line localLineToLine(const LocalLine& local_line, const std::list<Bin::MinZPoint>& line_points) {
		Line line;
		const double first_d = line_points.front().d;
		const double second_d = line_points.back().d;
		const double first_z = local_line.first * first_d + local_line.second;
		const double second_z = local_line.first * second_d + local_line.second;
		line.first.z = first_z;
		line.first.d = first_d;
		line.second.z = second_z;
		line.second.d = second_d;
		return line;
	}

	double verticalDistanceToLine(const double& d, const double &z) {
		static const double kMargin = 0.1;
		double distance = -1;
		for (auto it = lines_.begin(); it != lines_.end(); ++it) {
			if (it->first.d - kMargin < d && it->second.d + kMargin > d) {
				const double delta_z = it->second.z - it->first.z;
				const double delta_d = it->second.d - it->first.d;
				const double expected_z = (d - it->first.d) / delta_d * delta_z + it->first.z;
				distance = std::fabs(z - expected_z);
			}
		}
		return distance;
	}

	LocalLine fitLocalLine(const std::list<Bin::MinZPoint>& points) {
		const unsigned int n_points = points.size();

		Eigen::MatrixXd X(n_points, 2);
		Eigen::VectorXd Y(n_points);
		unsigned int counter = 0;
		for (auto iter = points.begin(); iter != points.end(); ++iter) {
			X(counter, 0) = iter->d;
			X(counter, 1) = 1;
			Y(counter) = iter->z;
			++counter;
		}
		Eigen::VectorXd result = X.colPivHouseholderQr().solve(Y);
		LocalLine line_result;
		line_result.first = result(0);
		line_result.second = result(1);
		return line_result;
	}

	//这个就是拟合直线函数
	void fitSegmentLines() {
		auto line_start = bins_.begin();
		while (!line_start->hasPoint()) {
			++line_start;
			if (line_start == bins_.end()) return;
		}
		double cur_ground_height = sensor_height_;
		std::list<Bin::MinZPoint> current_line_points(1, line_start->getMinZPoint());
		LocalLine cur_line = std::make_pair(double(0), double(0));
		int jishu = 2;
		for (auto line_iter = line_start + 1; line_iter != bins_.end(); ++line_iter) {
			if (line_iter->hasPoint()) {
				Bin::MinZPoint cur_point = line_iter->getMinZPoint();
				jishu++;
				if (cur_point.d - current_line_points.back().d > long_threshold_) {
				}

				if (current_line_points.size() >= 2) {
					double expected_z = std::numeric_limits<double>::max();
					if (current_line_points.size() > 2) {
						expected_z = cur_line.first * cur_point.d + cur_line.second;
					}
					current_line_points.push_back(cur_point);
					cur_line = fitLocalLine(current_line_points);
					const double error = getMaxError(current_line_points, cur_line);
					if (error > max_error_ ||
						std::fabs(cur_line.first) > max_slope_ && std::fabs(expected_z - cur_point.z) > max_long_height_) {

						current_line_points.pop_back();

						if (current_line_points.size() >= 3) {
							const LocalLine new_line = fitLocalLine(current_line_points);
							lines_.push_back(localLineToLine(new_line, current_line_points));
							cur_ground_height = new_line.first * current_line_points.back().d + new_line.second;
						}

						current_line_points.erase(current_line_points.begin(), --current_line_points.end());
						--line_iter;
					}

					else {}
				}
				else {
					if (cur_point.d - current_line_points.back().d < long_threshold_ &&
						std::fabs(current_line_points.back().z + cur_ground_height) < max_start_height_) {
						current_line_points.push_back(cur_point);
					}
					else {
						current_line_points.clear();
						current_line_points.push_back(cur_point);
					}
				}
			}
		}

		if (current_line_points.size() > 2) {
			const LocalLine new_line = fitLocalLine(current_line_points);
			lines_.push_back(localLineToLine(new_line, current_line_points));
		}
	}

	inline std::vector<Bin>::iterator begin() {
		return bins_.begin();
	}

	inline std::vector<Bin>::iterator end() {
		return bins_.end();
	}
	bool getLines(std::list<Line>* lines);
};

class GroundSegmentation {

public:
	GroundSegmentation() {}
	const GroundSegmentationParams params_;
	vector<Segment> segments_;
	vector<std::pair<int, int> > bin_index_;
	vector<Bin::MinZPoint> segment_coordinates_;

	void insertPoints(boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>> &cloud) {
		//这部分代码就是给每个点确定它自己的区域，换句话说，就是划分成n*m个小格。然后每个小格的最低点，我是写在了bin类里面。
		bin_index_.resize(cloud->size());
		segment_coordinates_.resize(cloud->size());
		segments_.resize(params_.n_segments);
		for (int hi = 0; hi < params_.n_segments; ++hi) {
			segments_[hi].bins_.resize(params_.n_bins);
		}
		double segment_step = 2 * M_PI / params_.n_segments;
		for (int i = 0; i < cloud->size(); ++i) {
			double range_square = cloud->points[i].x * cloud->points[i].x + cloud->points[i].y * cloud->points[i].y;
			double range = sqrt(range_square);
			if (range_square < r_max_square_1) {
				unsigned int bin_index = range / r_1;
				unsigned int segment_index = floor((std::atan2(cloud->points[i].y, cloud->points[i].x) + M_PI) / (2 * M_PI / n_segment));
				if (bin_index > params_.n_bins - 1) { bin_index = params_.n_bins - 1; }
				if (segment_index > params_.n_segments - 1) { segment_index = params_.n_segments - 1; }
				segments_[segment_index].bins_[bin_index].addPoint(range, cloud->points[i].z);
				bin_index_[i] = std::make_pair(segment_index, bin_index);
			}
			else {
				unsigned int bin_index = (range - 20) / r_2 + 100;
				unsigned int segment_index = floor((std::atan2(cloud->points[i].y, cloud->points[i].x) + M_PI) / (2 * M_PI / n_segment));
				if (bin_index > n_bin - 1) { bin_index = n_bin - 1; }
				if (segment_index > n_segment - 1) { segment_index = segment_index % n_segment; }
				segments_[segment_index].bins_[bin_index].addPoint(range, cloud->points[i].z);
				bin_index_[i] = std::make_pair(segment_index, bin_index);
			}
			segment_coordinates_[i] = Bin::MinZPoint(range, cloud->points[i].z);
		}
	}

	void accuracy(boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>> &cloud, double&tem_tpr, double&tem_fpr, double&tem_acc, double&tem_iou1, double&tem_iou2, double&tem_f1) {
		int TP = 0;
		int FP = 0;
		int FN = 0;
		int TN = 0;
		double segment_step = 2 * M_PI / params_.n_segments;
		for (unsigned int i = 0; i < bin_index_.size(); ++i) {
			Bin::MinZPoint point_2d = segment_coordinates_[i];
			int segment_index = bin_index_[i].first;
			if (segment_index >= 0) {
				if (segment_index > params_.n_segments - 1) { segment_index = params_.n_segments - 1; }
				double dist = segments_[segment_index].verticalDistanceToLine(point_2d.d, point_2d.z);
				int steps = 1;
				while (dist == -1 && steps * segment_step < params_.line_search_angle) {
					int index_1 = segment_index + steps;
					while (index_1 >= params_.n_segments) index_1 -= params_.n_segments;
					int index_2 = segment_index - steps;
					while (index_2 < 0) index_2 += params_.n_segments;
					const double dist_1 = segments_[index_1].verticalDistanceToLine(point_2d.d, point_2d.z);
					const double dist_2 = segments_[index_2].verticalDistanceToLine(point_2d.d, point_2d.z);
					if (dist_1 > dist) {
						dist = dist_1;
					}
					if (dist_2 > dist) {
						dist = dist_2;
					}
					++steps;
				}
				if (dist < params_.max_dist_to_line && dist != -1) {
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
		}
		tem_tpr = long double(TP) / (long double(TP) + long double(FN));
		tem_fpr = long double(FP) / (long double(FP) + long double(TN));
		tem_acc = (long double(TP) + long double(TN)) / (long double(TP) + long double(TN) + long double(FP) + long double(FN));
		tem_iou1 = double(TP) / double(TP + FP + FN);
		tem_iou2 = double(TN) / double(TN + FP + FN);
		tem_f1 = 2 * double(TP) / (2 * double(TP) + double(FP) + double(FN));
	}
	void lineFit() {
		for (int i = 0; i < params_.n_segments; ++i) {
			//拟合直线
			segments_[i].fitSegmentLines();
		}
	}
	

	boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>> object_output(boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>> &cloud) {
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_object(new pcl::PointCloud<pcl::PointXYZRGB>);

		double segment_step = 2 * M_PI / params_.n_segments;
		for (unsigned int i = 0; i < bin_index_.size(); ++i) {
			Bin::MinZPoint point_2d = segment_coordinates_[i];
			int segment_index = bin_index_[i].first;
			if (segment_index >= 0) {
				if (segment_index > params_.n_segments - 1) { segment_index = params_.n_segments - 1; }
				double dist = segments_[segment_index].verticalDistanceToLine(point_2d.d, point_2d.z);
				int steps = 1;
				while (dist == -1 && steps * segment_step < params_.line_search_angle) {
					int index_1 = segment_index + steps;
					while (index_1 >= params_.n_segments) index_1 -= params_.n_segments;
					int index_2 = segment_index - steps;
					while (index_2 < 0) index_2 += params_.n_segments;
					const double dist_1 = segments_[index_1].verticalDistanceToLine(point_2d.d, point_2d.z);
					const double dist_2 = segments_[index_2].verticalDistanceToLine(point_2d.d, point_2d.z);
					if (dist_1 > dist) {
						dist = dist_1;
					}
					if (dist_2 > dist) {
						dist = dist_2;
					}
					++steps;
				}
				if (dist < params_.max_dist_to_line && dist != -1) {
				}
				else {
					cloud_object->push_back(cloud->points[i]);
				}
			}
		}
		return cloud_object;
	}


	boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>> ground_output(boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>> &cloud) {
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_ground(new pcl::PointCloud<pcl::PointXYZRGB>);

		double segment_step = 2 * M_PI / params_.n_segments;
		for (unsigned int i = 0; i < bin_index_.size(); ++i) {
			Bin::MinZPoint point_2d = segment_coordinates_[i];
			int segment_index = bin_index_[i].first;
			if (segment_index >= 0) {
				if (segment_index > params_.n_segments - 1) { segment_index = params_.n_segments - 1; }
				double dist = segments_[segment_index].verticalDistanceToLine(point_2d.d, point_2d.z);
				int steps = 1;
				while (dist == -1 && steps * segment_step < params_.line_search_angle) {
					int index_1 = segment_index + steps;
					while (index_1 >= params_.n_segments) index_1 -= params_.n_segments;
					int index_2 = segment_index - steps;
					while (index_2 < 0) index_2 += params_.n_segments;
					const double dist_1 = segments_[index_1].verticalDistanceToLine(point_2d.d, point_2d.z);
					const double dist_2 = segments_[index_2].verticalDistanceToLine(point_2d.d, point_2d.z);
					if (dist_1 > dist) {
						dist = dist_1;
					}
					if (dist_2 > dist) {
						dist = dist_2;
					}
					++steps;
				}
				if (dist < params_.max_dist_to_line && dist != -1) {
					cloud_ground->push_back(cloud->points[i]);
				}
				else {
				}
			}
		}
		return cloud_ground;
	}

};

void cluster(boost::shared_ptr<pcl::PointCloud<pcl::PointXYZ>> &cloud_object, boost::shared_ptr<pcl::PointCloud<pcl::PointXYZ>> &cloud, boost::shared_ptr<pcl::PointCloud<pcl::PointXYZ>> &cloud_ground) {
	pcl::VoxelGrid<pcl::PointXYZ> vg;
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered(new pcl::PointCloud<pcl::PointXYZ>);
	vg.setInputCloud(cloud_object);
	vg.setLeafSize(0.05f, 0.05f, 0.05f);
	vg.filter(*cloud_filtered);

	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
	//pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered(new pcl::PointCloud<pcl::PointXYZ>);
	tree->setInputCloud(cloud_filtered);

	std::vector<pcl::PointIndices> cluster_indices;
	pcl::EuclideanClusterExtraction<pcl::PointXYZ> ec;

	ec.setClusterTolerance(0.3);
	ec.setMinClusterSize(100);
	ec.setMaxClusterSize(25000);
	ec.setSearchMethod(tree);
	ec.setInputCloud(cloud_filtered);
	ec.extract(cluster_indices);

	int j = 0;
	std::vector<pcl::PointXYZ>min_point_OBB;
	std::vector<pcl::PointXYZ>max_point_OBB;
	std::vector<pcl::PointXYZ>position_OBB;
	std::vector<Eigen::Matrix3f>rotational_matrix_OBB;
	char namelist[256][256];
	char format[512];

	for (std::vector<pcl::PointIndices>::const_iterator it = cluster_indices.begin(); it != cluster_indices.end(); ++it)
	{
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_cluster(new pcl::PointCloud<pcl::PointXYZ>);

		for (std::vector<int>::const_iterator pit = it->indices.begin(); pit != it->indices.end(); ++pit)
			cloud_cluster->points.push_back(cloud_filtered->points[*pit]); 
		j++;
	}
	pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer("重庆理工大学"));
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_view(new pcl::PointCloud<pcl::PointXYZRGB>);

	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_object_view(new pcl::PointCloud<pcl::PointXYZRGB>);
	pcl::PointXYZRGB l;
	for (int i = 0; i < cloud_object->size(); i++) {
		l.x = cloud_object->points[i].x;
		l.y = cloud_object->points[i].y;
		l.z = cloud_object->points[i].z;
		l.r = 255;
		l.g = 0;
		l.b = 0;
		cloud_object_view->push_back(l);
		cloud_view->push_back(l);
	}

	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_ground_view(new pcl::PointCloud<pcl::PointXYZRGB>);
	pcl::PointXYZRGB l_v;
	for (int i = 0; i < cloud_ground->size(); i++) {
		l_v.x = cloud_ground->points[i].x;
		l_v.y = cloud_ground->points[i].y;
		l_v.z = cloud_ground->points[i].z;
		l_v.r = 0;
		l_v.g = 255;
		l_v.b = 0;
		cloud_view->push_back(l_v);
	}

	viewer->addPointCloud<pcl::PointXYZRGB>(cloud_view);
	cout << j << endl;
	while (!viewer->wasStopped())
	{
		viewer->spinOnce(100);
		boost::this_thread::sleep(boost::posix_time::microseconds(100000));
	}

}

double this_abs(double x) {
	return fabs(x);
}
boost::shared_ptr<pcl::PointCloud<pcl::PointXYZ>> screen(boost::shared_ptr<pcl::PointCloud<pcl::PointXYZ>> &cloud, int n) {
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_out(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointXYZ c;
	for (int i = 0; i < cloud->size(); i++) {
		if (this_abs(atan(double(cloud->points[i].y / cloud->points[i].x)) > (double(n / 360)*M_PI))) {
			c.x = cloud->points[i].x;
			c.y = cloud->points[i].y;
			c.z = cloud->points[i].z;
			cloud_out->push_back(c);
		}
	}
	return cloud_out;
}

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

//主函数
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
		if (t % 10 == 0) {
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