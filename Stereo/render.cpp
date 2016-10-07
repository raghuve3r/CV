// Skeleton code for B657 A4 Part 1.
// D. Crandall
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

int main(int argc, char *argv[])
{
  if(argc != 3)
    {
      cerr << "usage: " << argv[0] << " image_file disp_file" << endl;
      return 1;
    }

  string input_filename1 = argv[1];
  string input_filename2 = argv[2];

  // read in images
  CImg<double> image_rgb(input_filename1.c_str());
  CImg<double> image_disp(input_filename2.c_str());
  
  //Code written to perform the functionality described in the part1 of assignment 4
  int rgb_width = image_rgb.width();
  int rgb_height = image_rgb.height();
  double disparity;
  
  //Temporary Image to store values 
  CImg<double> temporary(rgb_width, rgb_height,1,3,0); 
  //Get the first channel of the RGB image.
  CImg<double> final_image(rgb_width, rgb_height,1,3,0); 
  
  for(int ctr = 0; ctr < rgb_width; ctr++){
  		for(int ctr1 = 0 ; ctr1 < rgb_height; ctr1++){
  			//Calculate Disparity
  			disparity = ceil(image_disp(ctr, ctr1) * 0.05);
  			int diff = ctr - disparity;
  			if(diff >= 0){
  				if(diff < rgb_width){
  					temporary(diff, ctr1, 0, 0) = image_rgb(ctr, ctr1, 0 , 0);
  					temporary(diff, ctr1, 0, 1) = image_rgb(ctr, ctr1, 0 , 1);
  					temporary(diff, ctr1, 0, 2) = image_rgb(ctr, ctr1, 0 , 2);
  				}
  			}
  		}
  }
  for(int ctr = 0; ctr < rgb_width; ctr++){
  		for(int ctr1 = 0; ctr1 < rgb_height ; ctr1++){
  			final_image(ctr, ctr1, 0, 0) = temporary(ctr, ctr1, 0, 0);
  			final_image(ctr, ctr1, 0, 1) = image_rgb(ctr, ctr1, 0 , 1);
  			final_image(ctr, ctr1, 0, 2) = image_rgb(ctr, ctr1, 0 , 2);
  		}
  }

  final_image.get_normalize(0,255).save((input_filename1 + "-stereogram.png").c_str());
  final_image.save("test.png");

  return 0;
}
