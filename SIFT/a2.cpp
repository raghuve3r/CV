// B657 assignment 2 skeleton code
//
// Compile with: "make"
//
// See assignment handout for command line and project specifications.


//Link to the header file
#include "CImg.h"
#include <ctime>
#include <iostream>
#include <list>
#include <map>
#include <stdlib.h>
#include <string>
#include <vector>
#include <time.h>
#include <Sift.h>

//Use the cimg namespace to access the functions easily
using namespace cimg_library;
using namespace std;

int getNewPoints(int i, int j, const float (*trans)[3], int axis){
	float a = trans[2][0] * i - trans[0][0];
	float b = trans[2][1] * i - trans[0][1];
	float c = trans[0][2] - i * trans[2][2];
	float d = trans[2][0] * j - trans[1][0];
	float e = trans[2][1] * j - trans[1][1];
	float f = trans[1][2] - j * trans[2][2];
	int x = (c*e - b*f)/(a*e - b*d);
	int y = (a*f - c*d)/(a*e - b*d);
	if(axis == 0){
		return x;
	}
	else{
		return y;
	}
}

//___________________________________________________________________________________________________________________//

CImg<float> transformation(const CImg<float> &input_image, const float (*trans)[3]){
	CImg<float> output(input_image.width(),input_image.height(), 1, 3);
	for(int i = 0; i < input_image.width(); i++){
		for(int j = 0; j < input_image.height(); j++){
			int x = getNewPoints(i, j, trans, 0);
			int y = getNewPoints(i, j, trans, 1);
			if(x >= 0 and x < input_image.width() and y >= 0 and y < input_image.height()){
				output.set_linear_atXY(input_image.atXY(x,y,0,0),i,j,0,0);
				output.set_linear_atXY(input_image.atXY(x,y,0,1),i,j,0,1);
				output.set_linear_atXY(input_image.atXY(x,y,0,2),i,j,0,2);
			}
			else{
				output.set_linear_atXY(0,i,j,0,0);
				output.set_linear_atXY(0,i,j,0,1);
				output.set_linear_atXY(0,i,j,0,2);
			}
		}	
	}
	return output;
}

//___________________________________________________________________________________________________________________//

list<list<int> > findMatches(const vector<SiftDescriptor> &descriptors,const vector<SiftDescriptor> &descriptorsSecond, float threshold){
	list<list<int> > matches;
	for(int i=0; i<descriptors.size(); i++){
		int min = 0;
		int secondMin = 0;
		float distMin = std::numeric_limits<float>::max();
		float distSecondMin = std::numeric_limits<float>::max();
		for(int j = 0; j < descriptorsSecond.size(); j++){
			float distance = 0;
			for(int l=0; l<128; l++){
				distance = distance + pow(descriptors[i].descriptor[l] - descriptorsSecond[j].descriptor[l], 2);
			}
			distance = sqrt(distance);
			if(distance < distMin){
				distSecondMin = distMin;
				distMin = distance;
				secondMin = min;
				min = j;
			}
			else if(distance < distSecondMin){
				distSecondMin = distance;
				secondMin = j;
			}
		}
		if(distSecondMin!=0){
			if(distMin/distSecondMin < threshold){
				list<int> match;
				match.push_back(i);
				match.push_back(min);
				matches.push_back(match);
			}
		}
	}
	return matches;
}

//___________________________________________________________________________________________________________________//

list<int> getValueFromIndex(const vector<SiftDescriptor> &descriptors, const vector<SiftDescriptor> &descriptorsSecond, list<list<int> > &matches, int index){
	int i = 0;
	list<int> coords; 
	std::list<list<int> >::iterator it = matches.begin();
	while(i<=index){
		if(i == index){
			list<int> current = *it;
			std::list<int>::iterator itint=current.begin();
			int pic1 = *itint;
			++itint;     
			int pic2 = *itint;
			coords.push_back(descriptors[pic1].col);
			coords.push_back(descriptors[pic1].row);
			coords.push_back(descriptorsSecond[pic2].col);
			coords.push_back(descriptorsSecond[pic2].row);
			return coords;
			i++;
		}
		else{
			i++;
			++it;
		}
	}
	return coords;
}

//___________________________________________________________________________________________________________________//

bool compare(const pair<string, int>  &pair1, const pair<string, int> &pair2){
	return pair1.second > pair2.second;
}

//___________________________________________________________________________________________________________________//


int getCount(const vector<SiftDescriptor> &descriptors, const vector<SiftDescriptor> &descriptorsSecond, const float (*trans)[3], list<list<int> > &matches){
	int count = 0;
	for(std::list<list<int> >::iterator it=matches.begin(); it != matches.end(); ++it){
		list<int> current = *it;
		std::list<int>::iterator itint=current.begin();
		int pic1 = *itint;
		++itint;     
		int pic2 = *itint;
		int x = descriptors[pic1].col;
		int y = descriptors[pic1].row;
		int xd = descriptorsSecond[pic2].col;
		int yd = descriptorsSecond[pic2].row;
		if(abs(getNewPoints(xd, yd, trans, 0) - x) <= 0 and abs(getNewPoints(xd, yd, trans, 1) - y) <= 0){
			count++;
		}
	}
	return count;
}

//___________________________________________________________________________________________________________________//

CImg<float> getTransformation(const vector<SiftDescriptor> &descriptors,const vector<SiftDescriptor> &descriptorsSecond, list<list<int> > &matches){
	int times = 0;
	int maxCount = 0;
	CImg<float> maxTrans(1,8,1,1);
	srand ( time(NULL) ); 
	while(times < 100){
		int x[4];
		int y[4];
		int xd[4];
		int yd[4];
		int randIndices[matches.size()];
		for(int i = 0; i < matches.size(); i++){
			randIndices[i] = i;
		}
		int l = 0;
		while(l<4){
			int randIndex = rand() % matches.size();
			if(randIndices[randIndex] != -1){
				list<int> coords = getValueFromIndex(descriptors, descriptorsSecond, matches, randIndex); 
				std::list<int>::iterator itint=coords.begin();
				x[l] = *itint;
				++itint;
				y[l] = *itint;
				++itint;
				xd[l] = *itint;
				++itint;
				yd[l] = *itint;
				l++;
				randIndices[randIndex] = -1;
			}
		}
		CImg<float> A(8,8,1,1);
		CImg<float> B(1,8,1,1);
		for(int i = 0; i < A.height(); i++){
			if(i%2 == 0){
				A(2,i) = 1;
				A(3,i) = 0;
				A(4,i) = 0;
				A(5,i) = 0;
			}
			else{
				A(5,i) = 1;
				A(0,i) = 0;
				A(1,i) = 0;
				A(2,i) = 0;
			}
		}
		A(0,0) = A(3,1) = x[0];
		A(1,0) = A(4,1) = y[0];
		A(0,2) = A(3,3) = x[1];
		A(1,2) = A(4,3) = y[1];
		A(0,4) = A(3,5) = x[2];
		A(1,4) = A(4,5) = y[2];
		A(0,6) = A(3,7) = x[3];
		A(1,6) = A(4,7) = y[3];
		A(6,0) = -(x[0] * xd[0]);  
		A(7,0) = -(y[0] * xd[0]);
		A(6,1) = -(x[0] * yd[0]); 
		A(7,1) = -(y[0] * yd[0]); 
		A(6,2) = -(x[1] * xd[1]); 
		A(7,2) = -(y[1] * xd[1]);
		A(6,3) = -(x[1] * yd[1]);
		A(7,3) = -(y[1] * yd[1]); 
		A(6,4) = -(x[2] * xd[2]); 
		A(7,4) = -(y[2] * xd[2]);
		A(6,5) = -(x[2] * yd[2]);
		A(7,5) = -(y[2] * yd[2]);
		A(6,6) = -(x[3] * xd[3]); 
		A(7,6) = -(y[3] * xd[3]);
		A(6,7) = -(x[3] * yd[3]);
		A(7,7) = -(y[3] * yd[3]);
		B(0,0) = xd[0];
		B(0,1) = yd[0];
		B(0,2) = xd[1];
		B(0,3) = yd[1];
		B(0,4) = xd[2];
		B(0,5) = yd[2];
		B(0,6) = xd[3];
		B(0,7) = yd[3];
		
		CImg<float> C = B.solve(A);	
		
		float trans[3][3] = {{C(0,0),C(0,1),C(0,2)} , {C(0,3), C(0,4), C(0,5)} , {C(0,6), C(0,7), 1}};
		int count = getCount(descriptors, descriptorsSecond, trans, matches);
		if(count > maxCount){
			maxCount = count;
			maxTrans = C;
		}
		else{
			if(times==0){
				maxTrans = C;
			}
		}
		times++;
	}
	return maxTrans;
}


//___________________________________________________________________________________________________________________//


double gaussrand(){
    static double V1, V2, S;
    static int phase = 0;
    double X;
    if ( phase == 0 ) {
        do {
            double U1 = (double)rand() / RAND_MAX;
            double U2 = (double)rand() / RAND_MAX;
            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while(S >= 1 || S == 0);
        X = V1 * sqrt(-2 * log(S) / S);
    } else
        X = V2 * sqrt(-2 * log(S) / S);
    phase = 1 - phase;
    return X;
}

//_____________________________________________________________________________________________________________________//

vector< vector<double> > summaryDesciptor(vector< vector<double> > trans, vector<SiftDescriptor> SiftDesc, int k, int w){
	double newDes;
	vector< vector<double> > newDesc1;
	for(int i = 0; i < SiftDesc.size(); i++) {
		vector<double> rr;
		for(int q = 0; q < k; q++) {
			for(int j = 0; j < 128; j++) {
				newDes += SiftDesc[i].descriptor[j] * trans[q][j];
			}
			rr.push_back(floor(newDes/w));
		}
		newDesc1.push_back(rr);
	}
	return newDesc1;
}

//___________________________________________________________________________________________________________________//

int main(int argc, char **argv){
	try{
		if(argc < 2){
			cout << "Insufficent number of arguments; correct usage:" << endl;
			cout << "    a2-p1 part_id ..." << endl;
			return -1;
		}
		int flag = 0;
		string part = argv[1];
		if(part == "part1"){
			if(argc == 5){
				part = "part1_1";
			}
			else if(argc == 7){
				string foo = argv[5];
				if(foo.substr(foo.find_last_of(".") + 1)=="jpg"){
					part = "part1_2";
				}
				else{
					part = "part1_4";
				}
			}
			else if(argc >= 106){
				part = "part1_4b";
			}
			else if(argc > 5 and argc < 70){
				part = "part1_2";
			}
			else if(argc > 70){
				part = "part1_3";
			}			
		}
		else if(part == "part2"){
			if(argc == 3){
				flag = 1;
			}
		}
		float threshold;
		string inputFile;
		if(flag == 0){
			threshold = atof(argv[2]);
			inputFile = argv[3];
		}
		else{
			inputFile = argv[2];
		}
		map<string, int> list_of_image;
		float truePositive = 0;
		float falsePositive = 0;
		float precision = 0;	
		const unsigned char green[] = { 0,255,0 };
		if(part == "part1_1"){
			// This is just a bit of sample code to get you started, to
			// show how to use the SIFT library.
			CImg<float> input_image(inputFile.c_str());
			// convert image to grayscale
			CImg<float> grayInput = input_image.get_RGBtoHSI().get_channel(2);
			vector<SiftDescriptor> descriptors = Sift::compute_sift(grayInput);
			string file = argv[4];
			CImg<float> second_image(file.c_str());
			CImg<float> graySecond = second_image.get_RGBtoHSI().get_channel(2);
			vector<SiftDescriptor> descriptorsSecond = Sift::compute_sift(graySecond);
			list<list<int> > matches = findMatches(descriptors, descriptorsSecond, threshold);
			CImg<float> third_image;
            third_image.append(input_image,'x');
            third_image.append(second_image,'x');
			for (std::list<list<int> >::iterator it=matches.begin(); it != matches.end(); ++it){
				list<int> current = *it;
				std::list<int>::iterator itint=current.begin();
				int pic1 = *itint;
				++itint;
                int pic2 = *itint;
                int x1 = descriptors[pic1].col;
                int y1 = descriptors[pic1].row;
                int x2 = input_image.width() + descriptorsSecond[pic2].col;
                int y2 = descriptorsSecond[pic2].row;
                third_image.draw_line(x1,y1,x2,y2,green);
			}
			third_image.get_normalize(0,255).save("sift2.jpg");
		}
		else if(part == "part1_2"){
			CImg<float> input_image(inputFile.c_str());
			// convert image to grayscale
			CImg<float> grayInput = input_image.get_RGBtoHSI().get_channel(2);
			vector<SiftDescriptor> descriptors = Sift::compute_sift(grayInput);
			for(int i = 4; i < argc; i++){
				string listImage = argv[i];
				CImg<float> list_image(listImage.c_str());
				CImg<float> grayListImage = list_image.get_RGBtoHSI().get_channel(2);
				vector<SiftDescriptor> descriptorsSecond = Sift::compute_sift(grayListImage);
				list<list<int> > matches = findMatches(descriptors, descriptorsSecond, threshold);
				list_of_image[listImage] = matches.size();
			}
			vector< pair<string,int> > matchedVector;
			copy(list_of_image.begin(),list_of_image.end(),back_inserter(matchedVector));
			sort(matchedVector.begin(),matchedVector.end(),compare);
			for(int j=0;j<(argc-4);j++){
				cout<<matchedVector[j].first<<" : "<<matchedVector[j].second<<endl;
			}
		}
		else if(part == "part1_3"){
			string delimiter = "_"; 
			string fileName = inputFile.substr(0,inputFile.find(delimiter));
			CImg<float> input_image(inputFile.c_str());
			// convert image to grayscale
			CImg<float> grayInput = input_image.get_RGBtoHSI().get_channel(2);
			vector<SiftDescriptor> descriptors = Sift::compute_sift(grayInput);
            for(int i = 4; i < argc; i++){
				string listImage = argv[i];
				CImg<float> list_image(listImage.c_str());
				CImg<float> grayListImage = list_image.get_RGBtoHSI().get_channel(2);
				vector<SiftDescriptor> descriptorsSecond = Sift::compute_sift(grayListImage);
				list<list<int> > matches = findMatches(descriptors, descriptorsSecond, threshold);
				list_of_image[listImage] = matches.size();
			}
			vector< pair<string,int> > matchedVector;
			copy(list_of_image.begin(),list_of_image.end(),back_inserter(matchedVector));
			sort(matchedVector.begin(),matchedVector.end(),compare);
			for(int j = 0; j < 10; j++){
				cout << matchedVector[j].first << " : " << matchedVector[j].second << endl;
				if(matchedVector[j].first.find(fileName) != std::string::npos)
					truePositive++;
				else
					falsePositive++;
			}
			precision = (truePositive/(truePositive+falsePositive) * 100);
			cout << "Precision :"<< precision << endl;
		}
		else if(part == "part1_4"){
			string second_file = argv[4];
			int k = atoi(argv[5]);				
			int w = atoi(argv[6]);
			CImg<double> input_image(inputFile.c_str());
			// convert image to grayscale
			CImg<double> grayInput = input_image.get_RGBtoHSI().get_channel(2);			
			vector<SiftDescriptor> descriptors = Sift::compute_sift(grayInput);  
			CImg<float> second_image(second_file.c_str());
			CImg<float> graySecond = second_image.get_RGBtoHSI().get_channel(2);
			vector<SiftDescriptor> descriptorsSecond = Sift::compute_sift(graySecond);
			list<list<int> > matches = findMatches(descriptors, descriptorsSecond, threshold);
			CImg<float> third_image;
            third_image.append(input_image,'x');
            third_image.append(second_image,'x');			
			vector< vector<double> > transform;
			for(int y=0;y<k;y++){
				vector<double> row;
				for(int j=0;j<128;j++){
					row.push_back(gaussrand());
				}
				transform.push_back(row);
			}
			vector< vector<double> > newD = summaryDesciptor(transform,descriptors, k, w);				
			vector< vector<double> > newD1 = summaryDesciptor(transform,descriptorsSecond, k, w);
			vector< vector<double> > euDistance;
			for(int x=0;x<descriptors.size();x++){
				vector<double> first_1; 
				for(int y=0;y<descriptorsSecond.size();y++){
					double distance = 0;
					for(int m=0;m<k;m++){
						double dis = newD[x][m] - newD1[y][m];
						distance += sqrt(dis * dis);
					}	
					first_1.push_back(distance);
				}
				euDistance.push_back(first_1);
			}
			vector<int> near;
			for(int xy = 0; xy<euDistance.size(); xy++){
				double maxNumber = std::numeric_limits<double>::max();
				double secondMax = std::numeric_limits<double>::max();
				int one,two;
				for(int yz=0;yz<euDistance[0].size();yz++){
					if(euDistance[xy][yz]<secondMax){
						if(euDistance[xy][yz]<maxNumber){
							secondMax=maxNumber;
							maxNumber=euDistance[xy][yz];
							two = one;
							one = yz;
						}
						else{
							secondMax = euDistance[xy][yz];
							two = yz;
						}
					}
				}
				near.push_back(one);
				near.push_back(two);
			}
			for(int i=0;i<descriptors.size();i++){
				double dist1=0,dist2=0;
				for(int j=0;j<128;j++){	
					double a = descriptors[i].descriptor[j] - descriptorsSecond[near[i*2]].descriptor[j];
					dist1 += sqrt(a * a);
					double b = descriptors[i].descriptor[j] - descriptorsSecond[near[i*2+1]].descriptor[j];
					dist2 += sqrt(b * b);	
				}
				double ratio = dist1/dist2;
				int index = near[i*2];
				if(dist2<dist1){
					ratio = dist2/dist1;
					index = near[i*2+1];
				}
				cout<<ratio<<endl;
				if(ratio<0.6){
					int x1 = descriptors[i].col;
					int y1 = descriptors[i].row;
					int x2 = input_image.width() + descriptorsSecond[index].col;
					int y2 = descriptorsSecond[index].row;
					third_image.draw_line(x1,y1,x2,y2,green);					
				}
			}
			third_image.get_normalize(0,255).save("sift12.jpg");
		}
		else if(part == "part1_4b"){
			int k = atoi(argv[4]);				
			int w = atoi(argv[5]);	
			string delimiter = "_";  
			string fileName = inputFile.substr(0,inputFile.find(delimiter));
			CImg<double> input_image(inputFile.c_str());
			// convert image to grayscale
			CImg<double> grayInput = input_image.get_RGBtoHSI().get_channel(2);			
			vector<SiftDescriptor> descriptors = Sift::compute_sift(grayInput);    						
			vector< vector<double> > transform;
			for(int y=0;y<k;y++){
				vector<double> row;
				for(int j=0;j<128;j++){
					row.push_back(gaussrand());
				}
				transform.push_back(row);
			}
			vector< vector<double> > newD = summaryDesciptor(transform,descriptors,k,w);
            for(int i = 6; i < argc; i++){
            			int counter = 0;
				string listImage = argv[i];
				CImg<float> list_image(listImage.c_str());
				CImg<float> grayListImage = list_image.get_RGBtoHSI().get_channel(2);
				vector<SiftDescriptor> descriptorsSecond = Sift::compute_sift(grayListImage);
				vector< vector<double> > newD1 = summaryDesciptor(transform,descriptorsSecond,k,w);
				vector< vector<double> > euDistance;
				for(int x=0;x<descriptors.size();x++){
					vector<double> first_1; 
					for(int y=0;y<descriptorsSecond.size();y++){
						double distance = 0;
						for(int m=0;m<k;m++){
							double dis = newD[x][m] - newD1[y][m];		
							distance += sqrt(dis * dis);
						}
							first_1.push_back(distance);
						}
						euDistance.push_back(first_1);
					}
					vector<int> near;
					for(int xy = 0;xy<euDistance.size();xy++){
						double maxNumber = std::numeric_limits<double>::max();
						double secondMax = std::numeric_limits<double>::max();
						int one,two;
						for(int yz=0;yz<euDistance[0].size();yz++){
							if(euDistance[xy][yz]<secondMax){
								if(euDistance[xy][yz]<maxNumber){
									secondMax=maxNumber;
									maxNumber=euDistance[xy][yz];
									two = one;
									one = yz;
								}
								else{
									secondMax = euDistance[xy][yz];
									two = yz;
								}
							}
						}
						near.push_back(one);
						near.push_back(two);
					}
					for(int i=0;i<descriptors.size();i++){
						double dist1=0,dist2=0;
						for(int j=0;j<128;j++){		
							double a = descriptors[i].descriptor[j] - descriptorsSecond[near[i*2]].descriptor[j];
							dist1 += sqrt(a * a);
							double b = descriptors[i].descriptor[j] - descriptorsSecond[near[i*2+1]].descriptor[j];
							dist2 += sqrt(b * b);	
						}
						double ratio = dist1/dist2;
						int index = near[i*2];
						if(dist2<dist1){
							ratio = dist2/dist1;
							index = near[i*2+1];
						}	
						if(ratio<0.4){
							counter++;
						}
					}
				list_of_image[listImage] = counter;
			}
			vector< pair<string,int> > matchedVector;
			copy(list_of_image.begin(),list_of_image.end(),back_inserter(matchedVector));
			sort(matchedVector.begin(),matchedVector.end(),compare);
			for(int j = 0; j < 10; j++){
				cout << matchedVector[j].first << " : " << matchedVector[j].second << endl;
				if(matchedVector[j].first.find(fileName) != std::string::npos)
					truePositive = truePositive + 1;
				else
					falsePositive = falsePositive + 1;
			}
			precision = (truePositive/(truePositive+falsePositive) * 100);
			cout << "Precision :"<< precision << endl;
		}
		else if(part == "part2"){
			if(flag == 1){
				string file = "lincoln.jpg";
				CImg<float> sample_input_image(file.c_str());
				float sample_trans[3][3] = {{0.907,0.258,-182} , {-0.153, 1.44, 58} , {-0.000306, 0.000731, 1}};
				CImg<float> sample_output = transformation(sample_input_image, sample_trans);
				sample_output.save("warpedsample.jpg");	
			}
			else{
				CImg<float> input_image_source(inputFile.c_str());
				CImg<float> grayInput = input_image_source.get_RGBtoHSI().get_channel(2);
				vector<SiftDescriptor> descriptors = Sift::compute_sift(grayInput);
				for(int i = 4; i < argc; i++){
					string imageToTransform = argv[i];	
					CImg<float> input_image_imageToTransform(imageToTransform.c_str());
					CImg<float> grayInputTrans = input_image_imageToTransform.get_RGBtoHSI().get_channel(2);
					vector<SiftDescriptor> transDescriptors = Sift::compute_sift(grayInputTrans);
					list<list<int> > matches = findMatches(descriptors, transDescriptors, threshold);
					CImg<float> B = getTransformation(descriptors, transDescriptors, matches); 
					float trans[3][3] = {{B(0,0),B(0,1),B(0,2)} , {B(0,3), B(0,4), B(0,5)} , {B(0,6), B(0,7), 1}};
					CImg<float> output = transformation(input_image_imageToTransform, trans);
					string name = "warped" + imageToTransform.substr(imageToTransform.find_last_of("/") + 1);
					output.save(name.c_str());
				}
			}
		}
		else
			throw std::string("unknown part!");
	}
	catch(const string &err) {
		cerr << "Error: " << err << endl;
	}
}
