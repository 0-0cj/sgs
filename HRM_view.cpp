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

# define cj (numeric_limits<double>::max)()
using namespace std;

# define diff_RMSE 0.01
# define n_bin 38
# define r_1 0.2
# define r_2 0.5
# define z_diff 0.15
# define n_segment 180
#define r_max_square 50*50
# define point_to_line 0.3
# define sensor_hight 1.8
# define end_point 0.15
# define slope 0.3
# define neighbor_number 10
# define gpr_l 0.1935
# define gpr_af 0.2415
# define gpr_an 0.0396
# define min_angle 65.2
# define f_angle 0.6
#define PI 3.14159265358979323846

int nu = 0;
class Bin {

public:
	double max_z;
	double min_z;
	double min_z_range;
	int id;
	double weight;
	double k;
	double b;
	double residuals;
	double bisquare;
private:
	bool has_point_;
public:
	Bin() :has_point_(false) { init(); }
	void addPoint(const pcl::PointXYZRGB &point) {
		const double d = sqrt(point.x * point.x + point.y * point.y);
		addPoint(d, point.z);
	}
	void calculate_residuals() {
		if (k != 0 && b != 0) {
			residuals = abs(min_z - (k*min_z_range + b));
		}
		else {
			residuals = 0;
		}
	}
	void calculate_bisquare(double median_residuals) {
		double r_star = residuals / (6 * median_residuals);
		if (abs(r_star) < 1) {
			bisquare = pow(1 - pow(r_star, 2), 2);
		}
		else {
			bisquare = 0;
		}
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
	int point_number;
public:
	Segment() {
		bins_.resize(n_bin);
		/*
		has_break_point = false;
		for (int i = 0; i < n_bin; i++) {
			bins_[i].id = i;
		}
		*/
		//lines_.begin = std::make_pair(double(0), double(0));
	};

	bool point_line(double d, double z, int bin_num) {
		/*
		if (!has_break_point) {
			for (auto iter = lines_.begin(); iter != lines_.end(); ++iter) {
				double hight = iter->first*d + iter->second;
				//double distance = abs(lines_.first*d - 1 * z + lines_.second) / (sqrt(lines_.first*lines_.first + 1));
				if (abs(z - hight) < point_to_line) {
					return false;
				}
			}
			return true;
		}
		else {
			auto iter = lines_.begin();
			for (int count = 0; count < break_point.size(); count++) {
				if (id <= break_point[count]) {
					double hight = iter->first*d + iter->second;
					if (abs(z - hight) < point_to_line) {
						return false;
					}
					else return true;
				}
				iter = iter++;
			}
			double hight = iter->first*d + iter->second;
			if (abs(z - hight) < point_to_line) {
				return false;
			}
			else return true;
		}
		*/



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
					//cout << "i:" << i << endl;
					for (int j = 0; j < i; j++) {
						if (it != lines_.end()) {
							++it;
						}
					}
					/*
					if (lines_.size() == 1) {
						for (int i_j = 0; i_j < break_point_.size(); i_j++) {
							cout << "断点："<<i_j <<"------:"<< break_point_[i_j] << endl;
						}
					}
					*/
					//cout << "j:" << j << endl;
					//cout << "长度:" << break_point_.size() << endl;
					//cout << "线段数目：" << lines_.size() << endl;
					//cout << "j:" << j << endl;
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

		//return true;
	}

	double k_f(double x1, double x2) {
		return gpr_af * gpr_af*exp(-pow(x1 - x2, 2) / (2 * pow(gpr_l, 2))) + pow(gpr_an, 2);
	}

	/*
	void ini_gpr() {
		int count_bin = 0;
		for (int i = 0; i < bins_.size(); i++) {
			if (bins_[i].hasPoint()) {
				count_bin++;
				//cout << "k:" << bins_[i].k << "--b:" << bins_[i].b << endl;
			}
		}
		point_number = count_bin;

	}
	*/

	double GPR(double x_range) {
		int count_bin = 0;
		for (int i = 0; i < bins_.size(); i++) {
			if (bins_[i].hasPoint()) {
				count_bin++;
				//cout << "k:" << bins_[i].k << "--b:" << bins_[i].b << endl;
			}
		}
		Eigen::MatrixXd K(count_bin, count_bin);
		Eigen::MatrixXd K_star(count_bin, 1);
		Eigen::MatrixXd f_y(count_bin, 1);
		int gpr_i = 0;
		for (int i = 0; i < bins_.size(); i++) {
			if (bins_[i].hasPoint()) {
				f_y(gpr_i, 0) = bins_[i].k*bins_[i].min_z_range + bins_[i].b;
				K_star(gpr_i, 0) = k_f(bins_[i].min_z_range, x_range);
				int gpr_j = 0;
				for (int j = 0; j < bins_.size() && gpr_j <= gpr_i; j++) {
					if (bins_[j].hasPoint()) {
						K(gpr_i, gpr_j) = k_f(bins_[i].min_z_range, bins_[j].min_z_range);
						if (gpr_i != gpr_j) {
							K(gpr_j, gpr_i) = K(gpr_i, gpr_j);
						}
						gpr_j++;
					}
				}
				gpr_i++;
			}
		}
		Eigen::MatrixXd gpr_u(1, 1);
		Eigen::MatrixXd gpr_o(1, 1);
		Eigen::MatrixXd K_star_star(1, 1);
		K_star_star(0, 0) = k_f(x_range, x_range);
		gpr_u = K_star.transpose()*K.inverse()*f_y;
		gpr_o = K_star_star - K_star.transpose()*K.inverse()*K_star;
		//cout << "k_star_star：" << gpr_o(0, 0) << endl;
		//cout << "pow(x_range - gpr_u(0, 0), 2) / (2 * pow(gpr_o(0, 0), 2))的值：" <<pow(x_range - gpr_u(0, 0), 2) / (2 * pow(gpr_o(0, 0), 2)) << endl;
		//cout << "加上exp:" << exp(-pow(x_range - gpr_u(0, 0), 2) / (2 * pow(gpr_o(0, 0), 2))) << endl;
		//cout << gpr_af * gpr_af*exp(-pow(x_range - gpr_u(0, 0), 2) / (2 * pow(gpr_o(0, 0), 2))) + pow(gpr_an, 2) << endl;
		//return gpr_af * gpr_af*exp(-pow(x_range-gpr_u(0,0), 2) / (2 * pow(gpr_o(0,0), 2))) + pow(gpr_an, 2);
		return gpr_u(0, 0);
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

	//这个就是拟合直线函数

	int subsection(const std::vector<Bin>& line, LocalLine &sub_lines_tem) {
		/*
		int count = 0;
		for (int i = 0; i < line.size(); ++i) {
			if (abs(line[i].min_z - (sub_lines_tem.first*line[i].min_z_range + sub_lines_tem.second)) > end_point) {
				has_break_point = true;
				break_point.push_back(line[i].id);
				return count;
			}
			count++;
		}
		*/
		for (int i = 0; i < line.size(); ++i) {
			if (abs(line[i].min_z - (sub_lines_tem.first*line[i].min_z_range + sub_lines_tem.second)) > end_point) {
				return i;
			}
		}
		return 0;
	}

	std::vector<Bin> get_seed_points() {
		std::vector<Bin> current_line_points;
		for (int i = 0; i < bins_.size(); ++i) {
			if (current_line_points.empty()) {
				//cout << "尝试第一个最低点：" << bins_[i].min_z << endl;
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
		/*
		LocalLine lines_tem = fitLocalLine(seed);
		lines_.push_back(lines_tem);
		//cout << "参数：" << lines_tem.first << "x+" << lines_tem.second << endl;
		std::vector<Bin> seed_next1;
		std::vector<Bin> seed_next2;
		int middle_point = subsection(seed, lines_tem);
		int count = 0;
		if (middle_point != 0) {
			//printf("参加了一次\n");
			lines_.pop_back();
			for (int j = 0; j < middle_point; j++) {
				seed_next1.push_back(seed[j]);
			}
			for (int i = middle_point; i < seed.size(); ++i) {
				seed_next2.push_back(seed[i]);
			}
			if (seed_next1.size() > 2) {
				fitlines(seed_next1);
			}
			if (seed_next2.size() > 2) {
				fitlines(seed_next2);
			}
		}
		*/
		LocalLine lines_tem = fitLocalLine(seed);
		lines_.push_back(lines_tem);
		//cout << "参数：" << lines_tem.first << "x+" << lines_tem.second << endl;
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
				//cout << "最后一个元素是：" << break_point_.back() << endl;
				break_point_.push_back(middle_point + break_point_.back());
				//cout << "最后再一次元素是：" << break_point_.back() << endl;
			}
		}

		if (middle_point != 0) {
			//printf("参加了一次\n");
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
				//lines_.push_back(lines_tem);
				//fitlines(seed_next1);
				//cout << "左边" << endl;
			}
			if (seed_next2.size() > 2) {
				fitlines(seed_next2);
			}
			else {
				break_point_.pop_back();
			}
		}
	}

	std::vector<Bin> neighbor(int i, Bin m) {
		int j_i = i + 1, k_i = i - 1;
		vector<Bin> neighbor_vector;
		bool left_end = false;
		bool right_end = false;
		if (j_i > bins_.size() - 1) right_end = true;
		if (k_i < 0)left_end = true;
		while ((!left_end || !right_end) && neighbor_vector.size() < neighbor_number) {
			if ((!right_end) && bins_[j_i].hasPoint()) {
				neighbor_vector.push_back(bins_[j_i]);
			}
			if ((!right_end)) {
				j_i++;
				if (j_i > bins_.size() - 1) right_end = true;
			}


			if (neighbor_vector.size() >= neighbor_number) {
				return neighbor_vector;
			}
			if ((!left_end) && bins_[k_i].hasPoint()) {
				neighbor_vector.push_back(bins_[k_i]);
			}
			if ((!left_end)) {
				k_i--;
				if (k_i < 0)left_end = true;
			}

		}
		return neighbor_vector;
	}


	double point_weight(Bin a, Bin b, double max_x_distance) {
		return pow(1 - pow(abs(a.min_z_range - b.min_z_range) / max_x_distance, 3), 3);
	}


	LocalLine LWR(const std::vector<Bin>& points, Bin current) {
		const unsigned int n_points = points.size();
		//printf("%d\n", n_points);
		Eigen::MatrixXd X(n_points, 2);
		if (n_points >= 2) {
			Eigen::MatrixXd X(n_points, 2);
			Eigen::VectorXd Y(n_points);
			Eigen::MatrixXd W(n_points, n_points);
			//cout << "这里有问题吗？" << endl;
			W = Eigen::MatrixXd::Zero(n_points, n_points);
			double max_x_distance = 0;
			for (int j = 0; j < points.size(); j++) {
				if (max_x_distance < abs(points[j].min_z_range - current.min_z_range)) {
					max_x_distance = abs(points[j].min_z_range - current.min_z_range);
				}
			}
			unsigned int counter = 0;
			//cout << "这里有问题吗？" << endl;
			for (int i = 0; i < points.size(); ++i) {
				X(counter, 0) = points[i].min_z_range;
				X(counter, 1) = 1;
				Y(counter) = points[i].min_z;
				W(counter, counter) = point_weight(points[i], current, max_x_distance);
				++counter;
			}
			//Eigen::MatrixXd result1 = (X.transpose()*W*X).inverse()*X.transpose()*W;
			//cout << "这里有问题吗？" << endl;
			Eigen::VectorXd result = (X.transpose()*W*X).inverse()*X.transpose()*W*Y;
			//cout << "这里没有问题" << endl;
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
	void fitSegmentLines() {
		double diff_ = 10;
		double RMSE_current = 0;
		vector <Bin> new_bin;
		while (diff_ > diff_RMSE) {
			//cout << "进入迭代" << endl;
			for (int i = 0; i < bins_.size(); i++) {
				if (bins_[i].hasPoint()) {
					new_bin.clear();
					//cout << "没问题" << endl;
					new_bin = neighbor(i, bins_[i]);
					//cout << "找到邻居节点" << endl;
					//cout << "当前节点的邻域个数：" << new_bin.size() << endl;
					/*
					if (i == 5) {
						cout << "---------------------------------------------------------" << endl;
						for (int j_j = 0; j_j < 5; j_j++) {
							cout <<"横坐标："<< new_bin[j_j].min_z_range << endl;
							cout << "纵坐标：" << new_bin[j_j].min_z << endl;
						}
						cout<< "---------------------------------------------------------" << endl;
					}
					*/
					LocalLine current_line = LWR(new_bin, bins_[i]);
					bins_[i].k = current_line.first;
					bins_[i].b = current_line.second;
					bins_[i].calculate_residuals();
					//cout <<"第"<<i <<"点的残差：" << bins_[i].residuals << endl;
				//	cout << "斜率：" << bins_[i].k << "      截距：" << bins_[i].b << endl;
					//cout << "拟合完毕"<<endl;
				}
			}
			vector<double> residuals_vector;
			double median_residuals;
			for (int i = 0; i < bins_.size(); i++) {
				if (bins_[i].hasPoint()) {
					residuals_vector.push_back(bins_[i].residuals);
				}
			}
			//cout << "这里有问题吗？" << endl;
			//cout << residuals_vector.size() / 2 << endl;
			median_residuals = residuals_vector[residuals_vector.size() / 2];
			//cout << median_residuals << endl;
			//cout << "这里有问题吗？" << endl;
			for (int i = 0; i < bins_.size(); i++) {
				if (bins_[i].hasPoint()) {
					bins_[i].calculate_bisquare(median_residuals);
					if (bins_[i].min_z > (bins_[i].k*bins_[i].min_z_range + bins_[i].b)) {
						bins_[i].min_z = bins_[i].bisquare*bins_[i].min_z;
						new_bin.clear();
						new_bin = neighbor(i, bins_[i]);
						if (new_bin.size() > 0) {
							double neighbor_min_z = new_bin[0].min_z;

							for (int i_i = 0; i_i < new_bin.size(); i_i++) {
								if (new_bin[i_i].min_z < neighbor_min_z) {
									neighbor_min_z = new_bin[i_i].min_z;
								}
							}

							if (neighbor_min_z > bins_[i].min_z) {
								bins_[i].min_z = neighbor_min_z;
							}
						}
					}
				}
			}
			double sum_var = 0;
			int var_count = 0;
			for (int i = 0; i < bins_.size(); i++) {
				if (bins_[i].hasPoint()) {
					sum_var = sum_var + pow(bins_[i].residuals, 2);
					var_count++;
				}
			}
			//cout << "个数：" << var_count << endl;
			double RMSE = sqrt(sum_var / double(var_count));
			//	cout << "RMSE" << RMSE << endl;
			diff_ = abs(RMSE - RMSE_current);
			//	cout << "差值：" << diff_ << endl;
			RMSE_current = RMSE;
		}
		//printf("没问题");
		/*
		cout << "第" << nu << "次" << endl;
		nu++;
		*/
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
		//double bin_step = (sqrt(r_max_square) - sqrt(r_min_square)) / n_bin;
		//double r_min = sqrt(r_min_square);

		for (int i = 0; i < cloud->size(); ++i) {
			double range_square = cloud->points[i].x * cloud->points[i].x + cloud->points[i].y * cloud->points[i].y;
			double range = sqrt(range_square);
			unsigned int bin_index = range_number(range);
			unsigned int segment_index = floor((std::atan2(cloud->points[i].y, cloud->points[i].x) + M_PI) / (2 * M_PI / n_segment));
			if (bin_index > n_bin - 1) { bin_index = n_bin - 1; }
			if (segment_index > n_segment - 1) { segment_index = segment_index % n_segment; }
			segments_[segment_index].bins_[bin_index].addPoint(cloud->points[i]);
			bin_index_[i] = std::make_pair(segment_index, bin_index);
			/*
			if (range_square > r_max_square_1) {
				//double angle = std::atan2(cloud->points[i].y, cloud->points[i].x);
				unsigned int bin_index = (range - 20) / r_2 + 100;
				unsigned int segment_index = floor((std::atan2(cloud->points[i].y, cloud->points[i].x) + M_PI) / (2 * M_PI / n_segment));
				//cout << "bin_index:" << bin_index << "----segment_index:" << segment_index << endl;
				if (bin_index > n_bin - 1) { bin_index = n_bin - 1; }
				if (segment_index > n_segment - 1) { segment_index = segment_index % n_segment; }
				segments_[segment_index].bins_[bin_index].addPoint(cloud->points[i]);
				bin_index_[i] = std::make_pair(segment_index, bin_index);
			}
			else {
				unsigned int bin_index = range / r_1;
				unsigned int segment_index = floor((std::atan2(cloud->points[i].y, cloud->points[i].x) + M_PI) / (2 * M_PI / n_segment));
				if (bin_index > n_bin - 1) { bin_index = n_bin - 1; }
				if (segment_index > n_segment - 1) { segment_index = segment_index % n_segment; }
				segments_[segment_index].bins_[bin_index].addPoint(cloud->points[i]);
				bin_index_[i] = std::make_pair(segment_index, bin_index);
			}
			*/
		}
	}
	void lineFit() {
		//segments_[4].fitSegmentLines();
		/*
		int n_sum = 0;
		for (int i = 0; i < 30; i++) {
			if (segments_[174].bins_[i].hasPoint()) {
				n_sum++;
			}
		}

		cout << n_sum << endl;
		segments_[174].fitSegmentLines();
		*/

		for (int i = 0; i < n_segment; ++i) {
			segments_[i].fitSegmentLines();
			//cout << "第" << i << "个切片完成"<<endl;
		}


	}

	void test_gpr() {
		//segments_[1].GPR();
	}



	boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>> object_output(boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>> &cloud) {
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_object(new pcl::PointCloud<pcl::PointXYZRGB>);
		//pcl::PointXYZ q;
		const double segment_step = 2 * M_PI / n_segment;
		for (int i = 0; i < cloud->size(); ++i) {
			/*
			double angle = std::atan2(cloud->points[i].y, cloud->points[i].x);
			unsigned int segment_index = (angle + M_PI) / segment_step;
			double range_square = cloud->points[i].x * cloud->points[i].x + cloud->points[i].y * cloud->points[i].y;
			double range = sqrt(range_square);
			double bin_step = (sqrt(r_max_square) - sqrt(r_min_square)) / n_bin;
			double r_min = sqrt(r_min_square);
			unsigned int bin_index = (range - r_min) / bin_step;
			if (bin_index > n_bin - 1) { bin_index = n_bin - 1; }
			if (segment_index > n_segment - 1) { segment_index = segment_index % n_segment; }

			if (range*segments_[segment_index].bins_[bin_index].k + segments_[segment_index].bins_[bin_index].b > cloud->points[i].z) {

			}
			else {
				cloud_object->push_back(cloud->points[i]);
			}
			*/
			/*
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
			*/
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
			//unsigned int segment_index = (angle + M_PI) / segment_step;
			double range_square = cloud->points[i].x * cloud->points[i].x + cloud->points[i].y * cloud->points[i].y;
			double range = sqrt(range_square);
			//double bin_step = (sqrt(r_max_square) - sqrt(r_min_square)) / n_bin;
			//double r_min = sqrt(r_min_square);
			//unsigned int bin_index;
			unsigned int bin_index = range_number(range);
			unsigned int segment_index = floor((std::atan2(cloud->points[i].y, cloud->points[i].x) + M_PI) / (2 * M_PI / n_segment));

			/*
			if (range > sqrt(r_max_square_1)) {
				bin_index = (range - 20) / r_2 + 100;
			}
			else {
				bin_index = (range) / r_1;
			}
			*/
			if (bin_index > n_bin - 1) { bin_index = n_bin - 1; }
			if (segment_index > n_segment - 1) { segment_index = segment_index % n_segment; }
			double point_pre = segments_[segment_index].GPR(range);

			if (abs(point_pre - cloud->points[i].z) < z_diff) {
				if (cloud->points[i].g != 0) {
					TP++;
				}
				else {
					FP++;
				}
			}
			else {
				if (cloud->points[i].r != 0) {
					TN++;
				}
				else {
					FN++;
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
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_point(new pcl::PointCloud<pcl::PointXYZRGB>);
		const double segment_step = 2 * M_PI / n_segment;
		for (int i = 0; i < cloud->size(); ++i) {
			double angle = std::atan2(cloud->points[i].y, cloud->points[i].x);
			//unsigned int segment_index = (angle + M_PI) / segment_step;
			double range_square = cloud->points[i].x * cloud->points[i].x + cloud->points[i].y * cloud->points[i].y;
			double range = sqrt(range_square);
			//double bin_step = (sqrt(r_max_square) - sqrt(r_min_square)) / n_bin;
			//double r_min = sqrt(r_min_square);
			//unsigned int bin_index;
			unsigned int bin_index = range_number(range);
			unsigned int segment_index = floor((std::atan2(cloud->points[i].y, cloud->points[i].x) + M_PI) / (2 * M_PI / n_segment));
			/*
			if (range > sqrt(r_max_square_1)) {
				bin_index = (range - 20) / r_2 + 100;
			}
			else {
				bin_index = (range) / r_1;
			}*/
			//unsigned int bin_index = (range - r_min) / bin_step;
			if (bin_index > n_bin - 1) { bin_index = n_bin - 1; }
			if (segment_index > n_segment - 1) { segment_index = segment_index % n_segment; }
			pcl::PointXYZRGB p;
			double point_pre = segments_[segment_index].GPR(range);
			//cout << "原高度：" << cloud->points[i].z << endl;
			//cout << "预测高度：" << point_pre << endl;
			if (abs(point_pre - cloud->points[i].z) < z_diff) {
				p.x = cloud->points[i].x;
				p.y = cloud->points[i].y;
				p.z = cloud->points[i].z;
				p.r = 0;
				p.g = 255;
				p.b = 0;
				cloud_point->push_back(p);
			}
			else {
				p.x = cloud->points[i].x;
				p.y = cloud->points[i].y;
				p.z = cloud->points[i].z;
				p.r = 255;
				p.g = 0;
				p.b = 0;
				cloud_point->push_back(p);
			}
		}
		return cloud_point;
	}


	boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>> ground_output(boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>> &cloud) {
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_ground(new pcl::PointCloud<pcl::PointXYZRGB>);
		const double segment_step = 2 * M_PI / n_segment;
		for (int i = 0; i < cloud->size(); ++i) {
			/*
			double angle = std::atan2(cloud->points[i].y, cloud->points[i].x);
			unsigned int segment_index = (angle + M_PI) / segment_step;
			double range_square = cloud->points[i].x * cloud->points[i].x + cloud->points[i].y * cloud->points[i].y;
			double range = sqrt(range_square);
			double bin_step = (sqrt(r_max_square) - sqrt(r_min_square)) / n_bin;
			double r_min = sqrt(r_min_square);
			unsigned int bin_index = (range - r_min) / bin_step;
			if (bin_index > n_bin - 1) { bin_index = n_bin - 1; }
			if (segment_index > n_segment - 1) { segment_index = segment_index % n_segment; }

			if (range*segments_[segment_index].bins_[bin_index].k + segments_[segment_index].bins_[bin_index].b > cloud->points[i].z) {
				cloud_ground->push_back(cloud->points[i]);
			}
			else {

			}
			*/
			/*
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
			*/
		}
		return cloud_ground;
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
	reader.read("E:\\系统浏览器下载\\pc_pcd\\000197.pcd", *cloud_input);
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_test(new pcl::PointCloud<pcl::PointXYZRGB>);
	for (int i = 0; i < cloud_input->size(); ++i) {
		if ((cloud_input->points[i].x*cloud_input->points[i].x + cloud_input->points[i].y*cloud_input->points[i].y) < r_max_square) {
			cloud_input->points[i].z = cloud_input->points[i].z + sensor_hight;
			cloud_test->push_back(cloud_input->points[i]);
		}
	}
	GroundSegmentation ground_segement;
	ground_segement.insertPoints(cloud_test);
	ground_segement.lineFit();
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_view(new pcl::PointCloud<pcl::PointXYZRGB>);
	cloud_view = ground_segement.cloud_output(cloud_test);
	pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer("HRT"));
	viewer->addPointCloud<pcl::PointXYZRGB>(cloud_view);
	while (!viewer->wasStopped())
	{
		viewer->spinOnce(100);
		boost::this_thread::sleep(boost::posix_time::microseconds(100000));
	}
}