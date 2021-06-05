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
#include <pcl/kdtree/kdtree.h>
#include <pcl/features/normal_3d.h>
#include <pcl/segmentation/extract_clusters.h>
#include <pcl/features/moment_of_inertia_estimation.h>

void main() {
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_input(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PCDReader reader;
	reader.read("E:\\ÏµÍ³ä¯ÀÀÆ÷ÏÂÔØ\\000000.pcd", *cloud_input);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
	for (int i = 0; i < cloud_input->size(); i++) {
		if ((cloud_input->points[i].x*cloud_input->points[i].x + cloud_input->points[i].y*cloud_input->points[i].y) < (20*20) && (cloud_input->points[i].x*cloud_input->points[i].x + cloud_input->points[i].y*cloud_input->points[i].y) > (0.5*0.5)) {
			cloud->push_back(cloud_input->points[i]);
		}
	}
	pcl::ModelCoefficients::Ptr coefficients(new pcl::ModelCoefficients);
	pcl::PointIndices::Ptr inliers(new pcl::PointIndices);
	pcl::SACSegmentation<pcl::PointXYZ> seg;
	seg.setOptimizeCoefficients(true);
	seg.setModelType(pcl::SACMODEL_PLANE);
	seg.setMethodType(pcl::SAC_RANSAC);
	seg.setDistanceThreshold(0.15);
	seg.setInputCloud(cloud);
	seg.segment(*inliers, *coefficients);

	pcl::PointCloud<pcl::PointXYZRGB>::Ptr out_cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
	pcl::PointXYZRGB p;
	for (int i = 0; i < inliers->indices.size(); ++i)
	{
		p.x = cloud->points.at(inliers->indices[i]).x;
		p.y = cloud->points.at(inliers->indices[i]).y;
		p.z = cloud->points.at(inliers->indices[i]).z;
		p.r = 0;
		p.g = 255;
		p.b = 0;
		out_cloud->push_back(p);
	}

	pcl::PointCloud<pcl::PointXYZ>::Ptr object_cloud(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::ExtractIndices<pcl::PointXYZ> extract;
	extract.setInputCloud(cloud);
	extract.setIndices(inliers);
	extract.setNegative(true);
	extract.filter(*object_cloud);
	for (int i = 0; i < object_cloud->size(); ++i)
	{
		p.x = object_cloud->points[i].x;
		p.y = object_cloud->points[i].y;
		p.z = object_cloud->points[i].z;
		p.r = 255;
		p.g = 0;
		p.b = 0;
		out_cloud->push_back(p);
	}
	pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer("viewer"));
	/*
	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
	tree->setInputCloud(object_cloud);
	std::vector<pcl::PointIndices> cluster_indices;
	pcl::EuclideanClusterExtraction<pcl::PointXYZ> ec;

	ec.setClusterTolerance(0.3);
	ec.setMinClusterSize(100);
	ec.setMaxClusterSize(25000);
	ec.setSearchMethod(tree);
	ec.setInputCloud(object_cloud);
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
			cloud_cluster->points.push_back(object_cloud->points[*pit]); 

		pcl::MomentOfInertiaEstimation <pcl::PointXYZ> feature_extractor;
		feature_extractor.setInputCloud(cloud_cluster);
		feature_extractor.compute();
		std::vector <float> moment_of_inertia;
		std::vector <float> eccentricity;
		pcl::PointXYZ min_point_OBB_0;
		pcl::PointXYZ max_point_OBB_0;
		pcl::PointXYZ position_OBB_0;
		Eigen::Matrix3f rotational_matrix_OBB_0;
		float major_value, middle_value, minor_value;
		feature_extractor.getOBB(min_point_OBB_0, max_point_OBB_0, position_OBB_0, rotational_matrix_OBB_0);
		min_point_OBB.push_back(min_point_OBB_0);
		max_point_OBB.push_back(max_point_OBB_0);
		position_OBB.push_back(position_OBB_0);
		rotational_matrix_OBB.push_back(rotational_matrix_OBB_0);
		j++;
	}

	Eigen::Vector3f *position = new Eigen::Vector3f[j];
	Eigen::Quaternionf *quat = new Eigen::Quaternionf[j];

	

	for (int k = 0; k < j; k++) {
		sprintf(namelist[k], "OBB%d", k);
		position[k](0) = position_OBB.at(k).x;
		position[k](1) = position_OBB.at(k).y;
		position[k](2) = position_OBB.at(k).z;

		quat[k] = rotational_matrix_OBB.at(k);
		//fp << rotational_matrix_OBB.at(k) << "\n" << max_point_OBB.at(k).x - min_point_OBB.at(k).x << "\t" << max_point_OBB.at(k).y - min_point_OBB.at(k).y << "\t" << max_point_OBB.at(k).z - min_point_OBB.at(k).z << "\n";
		// Eigen::Quaternionf quat_0(rotational_matrix_OBB.at(k));
		// quat.push_back(quat_0);


		viewer->addCube(position[k], quat[k], max_point_OBB.at(k).x - min_point_OBB.at(k).x, max_point_OBB.at(k).y - min_point_OBB.at(k).y, max_point_OBB.at(k).z - min_point_OBB.at(k).z, namelist[k]);
		viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_REPRESENTATION, pcl::visualization::PCL_VISUALIZER_REPRESENTATION_WIREFRAME, namelist[k]);
	}
	*/
	viewer->addPointCloud<pcl::PointXYZRGB>(out_cloud);
	
	while (!viewer->wasStopped())
	{
		viewer->spinOnce(100);
		boost::this_thread::sleep(boost::posix_time::microseconds(100000));
	}
}