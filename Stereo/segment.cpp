// Skeleton code for B657 A4 Part 2.
// D. Crandall
//
//
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <math.h>
#include <CImg.h>
#include <assert.h>

using namespace cimg_library;
using namespace std;

class Point
{
public:
  Point() {}
  Point(int _col, int _row) : row(_row), col(_col) {}
  int row, col;
};

double gaussian(CImg<double> x, double denom, CImg<double> cov_invert, CImg<double> mean){
    CImg<double> temp;
    temp = (x - mean) * cov_invert;
    temp = temp * ((x-mean).transpose());

    double numerat = exp((-0.5)*temp(0,0));

    double xy = numerat/denom;
    //cout << xy;
    return xy;

}

bool findFG(const vector<Point> &fg, int i,int j){
    bool val = false;
    for(int x=0; x<fg.size();x++){
        if(fg[x].row == i && fg[x].col == j){
            val = true;
        }
    }
    return val;
}


bool findBG(const vector<Point> &bg, int i,int j){
    bool val = false;
    for(int x=0; x<bg.size();x++){
        if(bg[x].row == i && bg[x].col == j){
            val = true;
        }
    }
    return val;
}

double mrfGaussian(int d,int i,int j, const CImg<double> &img, double denom, CImg<double> cov_invert, CImg<double> mean, const vector<Point> &fg, const vector<Point> &bg){

   if(d == 1 && findFG(fg,i,j)){
        return 0;
   }
   else if(d == 0 && findFG(fg,i,j)){
        return INFINITY;
   }
   else if(d == 1 && findBG(fg,i,j)){
        return INFINITY;
   }
   else if(d == 0 && findBG(fg,i,j)){
        return 0;
   }
    else{
    CImg<double> temp;
    temp = (img(j,i) - mean) * cov_invert;
    temp = temp * ((img(j,i)-mean).transpose());

    double numerat = exp((-0.5)*temp(0,0));

    double xy = numerat/denom;
    //cout << xy;
    return xy;
    }
}

CImg<double> calcMean(const CImg<double> &img, const vector<Point> &fg){
    CImg<double> mean(3,1);
    double R = 0;
    double G = 0;
    double B = 0;

    double aR = 0;
    double aG = 0;
    double aB = 0;
    for(int x=0; x<fg.size();x++){
        R += img(fg[x].row,fg[x].col,0,0);
        G += img(fg[x].row,fg[x].col,0,1);
        B += img(fg[x].row,fg[x].col,0,2);
    }
    aR = R/fg.size();
    aG = G/fg.size();
    aB = B/fg.size();

    mean(0,0) = aR;
    mean(1,0) = aG;
    mean(2,0) = aB;

    return mean;
}

CImg<double> calcVar(const CImg<double> &img, const vector<Point> &fg, CImg<double> meanA){
     double vR = 0;
     double vG = 0;
     double vB = 0;
     CImg<double> cov(3,3);
     for(int x=0; x<fg.size();x++){
        vR += pow((meanA(0,0) - img(fg[x].row,fg[x].col,0,0)),2);
        vG += pow((meanA(1,0) - img(fg[x].row,fg[x].col,0,1)),2);
        vB += pow((meanA(2,0) - img(fg[x].row,fg[x].col,0,2)),2);
     }

     double varR = vR/fg.size();
     double varG = vG/fg.size();
     double varB = vB/fg.size();

    //cout << varR << varG << varB << endl;
     cov(0,0) = varR;
     cov(0,1) = 0;
     cov(0,2) = 0;
     cov(1,1) = varG;
     cov(1,0) = 0;
     cov(1,2) = 0;
     cov(2,2) = varB;
     cov(2,0) = 0;
     cov(2,1) = 0;

    return cov;
}


CImg<double> naive_segment(const CImg<double> &img, const vector<Point> &fg, const vector<Point> &bg)
{
  // implement this in step 2...
  //  this placeholder just returns a random disparity map

 CImg<double> result(img.width(), img.height());
 CImg<double> mean(3,1);
 CImg<double> cov(3,3);
 CImg<double> cov_invert(3,3);

 for(int x=0; x<fg.size();x++){
     result(fg[x].row,fg[x].col) = 1;
 }

 for(int x=0; x<bg.size();x++){
     result(bg[x].row,bg[x].col) = 0;
 }
 mean = calcMean(img,fg);
 cov = calcVar(img,fg,mean);
 cov_invert = cov;
 cov_invert.invert();
 double det = cov_invert.det();
 double den = sqrt((2 * 3.14) * (2 * 3.14) * (2 * 3.14) * det);

  double beta = -7.1;

  for(int i=0; i<img.height(); i++){
    for(int j=0; j<img.width(); j++){
      CImg<double> x(3,1);
      x(0,0) = img(j,i,0,0);
      x(1,0) = img(j,i,0,1);
      x(2,0) = img(j,i,0,2);
      double cost = 0;
      if(result(j, i) != 1 || result(j, i) != 0){
         //for(int m = 0; m < 3; m++){
            cost = -1 * log(gaussian(x, den, cov_invert, mean));

         //}

      }
      cout << cost << endl;
      //cout << cost << endl;
      if(cost < beta){
         result(j,i) = 1;
        }
        else{
            result(j,i) = 0;

        }
    }}
  return result;
}

vector< vector<double> > intializeGrid(const CImg<double> $img, const vector<Point> &fg, const vector<Point> &bg){

 vector< vector<double> > neighborGrid(5, vector<double> (2));
/*
  for(int i=1; i<img.height(); i++){
    for(int j=1; j<img.width(); j++){
        for(int x=0; x<fg.size();x++){



     }
    }*/
}


CImg<double> mrf_segment(const CImg<double> &img, const vector<Point> &fg, const vector<Point> &bg)
{
  // implement this in step 3...
  //  this placeholder just returns a random disparity map by calling naive_segment


CImg<double> result(img.width(), img.height());
int iterations = 5;
vector<int> d1;
vector<int> d2;
d1.push_back(0);
d1.push_back(1);
d2.push_back(0);
d2.push_back(1);
map<int, map<int, map<int, map<int, vector<double> > > > > message;
CImg<double> mean(3,1);
CImg<double> cov(3,3);
CImg<double> cov_invert(3,3);
mean = calcMean(img,fg);
cov = calcVar(img,fg,mean);
cov_invert = cov;
cov_invert.invert();
int alpha = 0;
double det = cov_invert.det();
double den = sqrt((2 * 3.14) * (2 * 3.14) * (2 * 3.14) * det);
double lambda = 100;
return naive_segment(img,fg,bg);
 for(int x=0; x<fg.size();x++){
     result(fg[x].row,fg[x].col) = 1;
 }

 for(int x=0; x<bg.size();x++){
     result(bg[x].row,bg[x].col) = 0;
 }

for(int x=1;x<iterations;x++){
    for(int i=1; i<(img.height())-1; i++){
       for(int j=1; j<(img.width())-1; j++){
          CImg<double> nei(4,1);
          nei(0,0) = img(j,i+1);
          nei(1,0) = img(j+1,i);
          nei(2,0) = img(j,i-1);
          nei(3,0) = img(j-1,i);
                          double cost_mrf = 0;

          for(int f=0;f<nei.size();f++){
            for(int d_1 = 0;d_1<d1.size();d_1++){
                for(int d_2=0;d_2<d2.size();d_2++){
                    if(f == 0)
                        cost_mrf = -1 * log(mrfGaussian(d_2,i+1,j,img,den,cov_invert,mean,fg,bg));
                    else if(f == 1)
                        cost_mrf = -1 * log(mrfGaussian(d_2,i,j+1,img,den,cov_invert,mean,fg,bg));
                    else if(f == 2)
                        cost_mrf = -1 * log(mrfGaussian(d_2,i-1,j,img,den,cov_invert,mean,fg,bg));
                    else
                        cost_mrf = -1 * log(mrfGaussian(d_2,i,j-1,img,den,cov_invert,mean,fg,bg));

                    if(d_1 != d_2){
                        cost_mrf += cost_mrf + lambda;
                   }
                }
                //cout << cost_mrf << endl;
            }
            if (cost_mrf < alpha){
                   result(j,i) = 1;
            }
            else{
                result(j,i) = 0;
            }
          }

       }
    }



}
return result;
}

  //return naive_segment(img, fg, bg);


// Take in an input image and a binary segmentation map. Use the segmentation map to split the 
//  input image into foreground and background portions, and then save each one as a separate image.
//
// You'll just need to modify this to additionally output a disparity map.
//
void output_segmentation(const CImg<double> &img, const CImg<double> &labels, const string &fname)
{
  // sanity checks. If one of these asserts fails, you've given this function invalid arguments!
  assert(img.height() == labels.height());
  assert(img.width() == labels.width());

  CImg<double> img_fg = img, img_bg = img;

  for(int i=0; i<labels.height(); i++)
    for(int j=0; j<labels.width(); j++)
      {
	if(labels(j,i) == 0)
	  img_fg(j,i,0,0) = img_fg(j,i,0,1) = img_fg(j,i,0,2) = 0;
	else if(labels(j,i) == 1)
	  img_bg(j,i,0,0) = img_bg(j,i,0,1) = img_bg(j,i,0,2) = 0;
	else
	  assert(0);
      }

  img_fg.get_normalize(0,255).save((fname + "_fg.png").c_str());
  img_bg.get_normalize(0,255).save((fname + "_bg.png").c_str());
}

int main(int argc, char *argv[])
{
  if(argc != 3)
    {
      cerr << "usage: " << argv[0] << " image_file seeds_file" << endl;
      return 1;
    }

  string input_filename1 = argv[1], input_filename2 = argv[2];

  // read in images and gt
  CImg<double> image_rgb(input_filename1.c_str());
  CImg<double> seeds_rgb(input_filename2.c_str());

  // figure out seed points 
  vector<Point> fg_pixels, bg_pixels;
  for(int i=0; i<seeds_rgb.height(); i++)
    for(int j=0; j<seeds_rgb.width(); j++)
      {
	// blue --> foreground
	if(max(seeds_rgb(j, i, 0, 0), seeds_rgb(j, i, 0, 1)) < 100 && seeds_rgb(j, i, 0, 2) > 100)
	  fg_pixels.push_back(Point(j, i));

	// red --> background
	if(max(seeds_rgb(j, i, 0, 2), seeds_rgb(j, i, 0, 1)) < 100 && seeds_rgb(j, i, 0, 0) > 100)
	  bg_pixels.push_back(Point(j, i));
      }

  // do naive segmentation
  CImg<double> labels = naive_segment(image_rgb, fg_pixels, bg_pixels);
  output_segmentation(image_rgb, labels, input_filename1 + "-naive_segment_result");

  // do mrf segmentation
  labels = mrf_segment(image_rgb, fg_pixels, bg_pixels);
  output_segmentation(image_rgb, labels, input_filename1 + "-mrf_segment_result");

  return 0;
}
