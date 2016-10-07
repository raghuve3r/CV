#include <SImage.h>
#include <SImageIO.h>
#include <cmath>
#include <limits>
#include <queue>
#include <algorithm>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <DrawText.h>
#include <QtGui/qimage.h>

using namespace std;


#define MAX_GRAY_VALUE 255
#define MIN_GRAY_VALUE 0
#define MAX(a,b) ((a)>(b)?(a):(b))
#define DEG2RAD(deg) ((M_PI*deg)/180.0)
#define PNG_FOMAT "PNG"

struct hough_data{
	double rho_theta;
	double vote_count;
	public:
		hough_data(double r_t,int count) : rho_theta(r_t),vote_count(count) {}
  //min_heap
		bool operator<(const hough_data &other) const {
			return vote_count > other.vote_count;
		}
		//min_heap
		bool operator()(const hough_data &a, const hough_data &b){
			return a.vote_count > b.vote_count;
		}
};

int binary_threshold = 128;
double fl_threshold = 1;

// The simple image class is called SDoublePlane, with each pixel represented as
// a double (floating point) type. This means that an SDoublePlane can represent
// values outside the range 0-255, and thus can represent squared gradient magnitudes,
// harris corner scores, etc. 
//
// The SImageIO class supports reading and writing PNG files. It will read in
// a color PNG file, convert it to grayscale, and then return it to you in 
// an SDoublePlane. The values in this SDoublePlane will be in the range [0,255].
//
// To write out an image, call write_png_file(). It takes three separate planes,
// one for each primary color (red, green, blue). To write a grayscale image,
// just pass the same SDoublePlane for all 3 planes. In order to get sensible
// results, the values in the SDoublePlane should be in the range [0,255].
//

// Below is a helper functions that overlays rectangles
// on an image plane for visualization purpose. 

// Draws a rectangle on an image plane, using the specified gray level value and line width.
//

void overlay_rectangle(SDoublePlane &input, int _top, int _left, int _bottom, int _right, double graylevel, int width)
{
  for(int w=-width/2; w<=width/2; w++) {
    int top = _top+w, left = _left+w, right=_right+w, bottom=_bottom+w;

    // if any of the coordinates are out-of-bounds, truncate them 
    top = min( max( top, 0 ), input.rows()-1);
    bottom = min( max( bottom, 0 ), input.rows()-1);
    left = min( max( left, 0 ), input.cols()-1);
    right = min( max( right, 0 ), input.cols()-1);
      
    // draw top and bottom lines
    for(int j=left; j<=right; j++)
	  input[top][j] = input[bottom][j] = graylevel;
    // draw left and right lines
    for(int i=top; i<=bottom; i++)
	  input[i][left] = input[i][right] = graylevel;
  }
}

// DetectedSymbol class may be helpful!
//  Feel free to modify.
//
typedef enum {NOTEHEAD=0, QUARTERREST=1, EIGHTHREST=2} Type;
class DetectedSymbol {
public:
  int row, col, width, height;
  Type type;
  char pitch;
  double confidence;
};

SDoublePlane scaled_image(const char* input_file_name,const double scale_x,const double scale_y);


// Function that outputs the ascii detection output file
void  write_detection_txt(const string &filename, const vector<struct DetectedSymbol> &symbols)
{
  ofstream ofs(filename.c_str());

  for(int i=0; i<symbols.size(); i++)
    {
      const DetectedSymbol &s = symbols[i];
      ofs << s.row << " " << s.col << " " << s.width << " " << s.height << " ";
      if(s.type == NOTEHEAD)
	ofs << "filled_note " << s.pitch;
      else if(s.type == EIGHTHREST)
	ofs << "eighth_rest _";
      else 
	ofs << "quarter_rest _";
      ofs << " " << s.confidence << endl;
    }
}

// Function that outputs a visualization of detected symbols
void  write_detection_image(const string &filename, const vector<DetectedSymbol> &symbols, const SDoublePlane &input)
{
  SDoublePlane output_planes[3];
  for(int i=0; i<3; i++)
    output_planes[i] = input;

  for(int i=0; i<symbols.size(); i++)
    {
      const DetectedSymbol &s = symbols[i];

      overlay_rectangle(output_planes[s.type], s.row, s.col, s.row+s.height-1, s.col+s.width-1, 255, 2);
      overlay_rectangle(output_planes[(s.type+1) % 3], s.row, s.col, s.row+s.height-1, s.col+s.width-1, 0, 2);
      overlay_rectangle(output_planes[(s.type+2) % 3], s.row, s.col, s.row+s.height-1, s.col+s.width-1, 0, 2);

      if(s.type == NOTEHEAD)
	{
	  char str[] = {s.pitch, 0};
	  draw_text(output_planes[0], str, s.row, s.col+s.width+1, 0, 2);
	  draw_text(output_planes[1], str, s.row, s.col+s.width+1, 0, 2);
	  draw_text(output_planes[2], str, s.row, s.col+s.width+1, 0, 2);
	}
    }

  SImageIO::write_png_file(filename.c_str(), output_planes[0], output_planes[1], output_planes[2]);
}

void writeImage(const char* name, SDoublePlane &a) {
  vector<DetectedSymbol> symbols;
  write_detection_image(name, symbols,a);
}



// The rest of these functions are incomplete. These are just suggestions to 
// get you started -- feel free to add extra functions, change function
// parameters, etc.

// Convolve an image with a separable convolution kernel
//
SDoublePlane convolve_general(const SDoublePlane &input, const SDoublePlane &filter);
SDoublePlane convolve_separable(const SDoublePlane &input, const SDoublePlane &row_filter, const SDoublePlane &col_filter)
{  
	SDoublePlane output(input.rows(), input.cols());
	int rows = input.rows();
	int cols = input.cols();
		
	//applying row filter
	for(int j=0;j<cols;j++){
		for(int i=0;i<rows;i++){
			for(int f=0;f<row_filter.rows();f++){
				int dx = i-f-1;
				dx = dx>=0?dx:(dx+rows)%rows;//reflected index
				output[i][j]+=input[dx][j]*row_filter[f][0];
			}
		}		
	}
	 SDoublePlane output2(input.rows(), input.cols());
	 for(int i=0;i<rows;i++){
		 for(int j=0;j<cols;j++){
			 for(int f=0;f<col_filter.cols();f++){
				 int dy = j-f-1;
				 dy = dy>=0?dy:(dy+cols)%cols;//reflected index
				 output2[i][j]+=output[i][dy]*col_filter[0][f];
			 }
		}		
	 }
	 
		
  return output2;
}

// Convolve an image with a separable convolution kernel
//
SDoublePlane convolve_general(const SDoublePlane &input, const SDoublePlane &filter)
{
  SDoublePlane output(input.rows(), input.cols());
  // Convolution code here
  for (int i = 0; i < input.rows() ; i++) {
    for (int j = 0; j < input.cols() ; j++) {
      output[i][j] = 0;
      for (int q = 0; q < filter.rows(); q++) {
	for (int w = 0; w < filter.cols(); w++) {
	  if (i-q+1 >= 0 && (i-q+1) < input.rows()
	      && j-w+1 >= 0 && (j-w+1) < input.cols()) {
	    output[i][j] += filter[q][w] * input[i-q+1][j-w+1];
	  }
	}
      }      
    }
  }  
  return output;
}

SDoublePlane get_gradient_filter(int neigh_size, char x_y ){
	SDoublePlane grad_x(neigh_size,neigh_size);
	
	for(int i =0;i<neigh_size;i++){
		for(int j=0;j<neigh_size;j++){
			if((j==0&&x_y=='X')||(i==0 && x_y=='Y')){
				grad_x[i][j] = -1;
			}else if((j==(neigh_size-1)&&x_y=='X')||(i==(neigh_size-1)&&x_y=='Y')){
				grad_x[i][j] = 1;
			}
		}
	}
	return grad_x;		
}


SDoublePlane transpose(SDoublePlane &matrix){
	int r = matrix.rows();
	int c = matrix.cols();
	
	for(int i=0;i<r;i++){
		for(int j=i+1;j<c;j++){
			int temp = matrix[i][j];
			matrix[i][j] = matrix[j][i];
			matrix[j][i] = temp;
		}
	}
	return matrix;
}

SDoublePlane get_gaussian_filter(double sigma){
	double s = 2*sigma*sigma;
	int r = round(sigma*6+1);
	
	SDoublePlane fltr(r,r);
	double pwr = 0;
	int mid = r/2;
	double sum = 0;
	for(int i=-mid;i<=mid;i++){
		for(int j=-mid;j<=mid;j++){
			pwr = sqrt(i*i + j*j);
			fltr[i+mid][j+mid]= exp(-(r*r)/s)/(M_PI*s);
			sum+=fltr[i+mid][j+mid];
		}
	}
	//normalize
	for(int i=-mid;i<=mid;i++){
		for(int j=-mid;j<=mid;j++){
			fltr[i+mid][j+mid]/= sum;
		}
	}
	
	return fltr;
}

SDoublePlane compute_grad_magnitude(const SDoublePlane &grad_I_X,const SDoublePlane &grad_I_Y){
	int r = grad_I_X.rows();
	int c = grad_I_X.cols();
	SDoublePlane mag_grad_I_X_Y(r,c);
	
	for(int i = 0;i<r;i++){
		for(int j=0;j<c;j++){
			mag_grad_I_X_Y[i][j]=sqrt(grad_I_X[i][j]*grad_I_X[i][j]+grad_I_Y[i][j]*grad_I_Y[i][j]);
		}
	}
	return mag_grad_I_X_Y;
}


// Apply a sobel operator to an image, returns the result
// 
SDoublePlane sobel_gradient_filter(const SDoublePlane &input, bool _gx)
{
  // Implement a sobel gradient estimation filter with 1-d filters
  SDoublePlane filterX(0,0);
  SDoublePlane filterY(0,0);
  if(_gx){
	filterY = SDoublePlane(3,1);
	for(int i = 0;i<3;i++){
	  filterY[i][0] = (i%2==0?1:2)/8.0;
   }
	  filterX = SDoublePlane(1,3);
	  filterX[0][0] = -1;
	  filterX[0][1] = 0;
      filterX[0][2] = 1;
	 
  }else {
	  
	  filterX = SDoublePlane(1,3);
	  for(int i = 0;i<3;i++){
		filterX[0][i] = (i%2==0?1:2)/8.0;
	  }
	  filterY = SDoublePlane(3,1);
	  filterY[0][0] = 1;
	  filterY[1][0] = 0;
      filterY[2][0] = -1;
	  
	 
  }
  
   return convolve_separable(input,filterY,filterX);
}


// Apply an edge detector to an image, returns the binary edge map
// 
SDoublePlane find_edges(const SDoublePlane &input, double thresh=0)
{
  
  SDoublePlane grad_I_X = sobel_gradient_filter(input,true);
  SDoublePlane grad_I_Y = sobel_gradient_filter(input,false);
  
  // Implement an edge detector of your choice, e.g.
  // use your sobel gradient operator to compute the gradient magnitude and threshold
  SDoublePlane mag_grad_I_X_Y = compute_grad_magnitude(grad_I_X,grad_I_Y);
  
  //filter the image with respect to threshold
  int rows = mag_grad_I_X_Y.rows();
  int cols = mag_grad_I_X_Y.cols();
    
  for(int i=0;i<rows;i++){
	  for(int j=0;j<cols;j++){
		  if((mag_grad_I_X_Y[i][j])<thresh){
			  mag_grad_I_X_Y[i][j] = MIN_GRAY_VALUE;
		  }else {
			 //mag_grad_I_X_Y[i][j] = MAX_GRAY_VALUE-30;
		  }
	  }
  }
  
  return mag_grad_I_X_Y;
}

float get_direction(const float theta){
	//return (fmod(theta+M_PI,M_PI)/M_PI)*8.0;
	float degree = (360.0*theta)/(2*M_PI);
	float direction = 0; 
	if(degree>=0 && degree<=45){
		direction = 1;
	}else if(degree>45 && degree<=90){
		direction = 5;
	}else if(degree>90 && degree<=135){
		direction = 6;
	}else if(degree>135 && degree<=180){
		direction = 1;
	}else if(degree>=-45 && degree<0){
		direction = 7;
	}else if(degree>=-90 && degree<-45){
		direction = 4;
	}else if(degree>=-135 && degree<-90){
		direction = 3;
	}else {
		direction = 8;
	}
	
	return direction;
}

SDoublePlane perfomr_nms_suppression_using_gradient(const SDoublePlane &mag,const SDoublePlane grad_x,const SDoublePlane grad_y){
	//using the eight neighbours, deciding the values for the current (i,j) values could be, depending upon the gradient directions
	int rows = mag.rows();
	int cols = mag.cols();
		
	SDoublePlane nms_output(rows,cols);
	//non - maximal supression
	//considering eight neighbors
	//strating from (1,1) for the eight neighbors possible
	for(int i=1;i<rows-1;i++){
		for(int j=1;j<cols-1;j++){
			//centre pixel is (i,j)
			//get all the neighboring boundaries,in eight directions
			int r_up = i-1;
			int r_bt = i+1; 
			int c_lt = j-1;
			int c_rt = j+1;
			
			float theta = atan2(grad_y[i][j],grad_x[i][j]);
			//convert into 1 - 8 directions
			float direction = get_direction(theta);
			
			bool is_more_than_nghbr = ((direction<=1 || direction> 7) && mag[i][j]>mag[i][c_lt] && mag[i][j]>mag[i][c_rt])
									 ||(direction>1 && direction<=3 && mag[i][j]>mag[r_up][c_rt] && mag[i][j]>mag[r_bt][c_lt])
									 ||(direction>3 && direction<=5 && mag[i][j]>mag[r_up][j] && mag[i][j]>mag[r_bt][j])
									 || (direction>5 && direction<=7 && mag[i][j]>mag[r_up][c_lt] && mag[i][j]>mag[r_bt][c_rt]);
			
			if(is_more_than_nghbr){
				nms_output[i][j] = mag[i][j];
			}else {
				nms_output[i][j] = MIN_GRAY_VALUE;//local minima
			}			
		}
	}
	
	
	//fill the first and the last row with 1
	for(int j = 0;j<cols;j++){
		nms_output[0][j] = 1;
		nms_output[rows-1][j]= 1;
	}
	
	//vector<DetectedSymbol> symbols;
    //write_detection_image("detect_nms_temp.png", symbols, nms_output);
	return nms_output;
}


SDoublePlane  apply_hysteresis_threshold(const SDoublePlane &nms_output,double low_thresold,double high_thresold,const int high_gray_value){
	int rows = nms_output.rows();
	int cols = nms_output.cols();
	
	SDoublePlane output(rows,cols);
	
	for(int i=1;i<rows-1;i++){
		for(int j=1;j<cols-1;j++){
			if(nms_output[i][j]>=high_thresold && output[i][j]==0.0){
				//not yet traced current pixels ie. output[i][j]==0 and its passing the maximum threshold value
					output[i][j] = MAX_GRAY_VALUE;				
					int ngr_index = i*cols+j;
					//now see eight neighbors in each directions
					int indx_ngb[8];
					indx_ngb[0] = ngr_index+1;//right
					indx_ngb[1] = ngr_index-1;//left
					indx_ngb[2] = ngr_index+cols;//bottom
					indx_ngb[3] = ngr_index-cols;//up
					indx_ngb[4] = indx_ngb[2]+1;//bottom-right
					indx_ngb[5] = indx_ngb[2]-1;//bottom-left
					indx_ngb[6] = indx_ngb[3]+1;//up-right
					indx_ngb[7] = indx_ngb[3]-1;//up-left
					
					for(int n_i =0;n_i<8;n_i++){
						int r_i = indx_ngb[n_i]/cols;
						int c_i = indx_ngb[n_i]%cols;
						if(r_i<0 ||c_i<0 ||r_i>=rows||c_i>=cols){
							continue;
						}
						
						if(output[r_i][c_i]==0.0 && nms_output[r_i][c_i]>low_thresold){
							output[r_i][c_i] = high_gray_value;//nms_output[i][j];
						}else {
							//output[r_i][c_i] = MIN_GRAY_VALUE;
						}
					}
			}
		}
	}
	return output;
}

SDoublePlane apply_non_maxmal_suppression(const SDoublePlane &mag,const SDoublePlane grad_x,const SDoublePlane grad_y,const double  low_thresold,const double high_thresold,const int high_gray_value){
	SDoublePlane nms_output = mag;
	nms_output = perfomr_nms_suppression_using_gradient(mag,grad_x,grad_y);
	return apply_hysteresis_threshold(nms_output,low_thresold,high_thresold,high_gray_value);
}

SDoublePlane find_edges_non_max_supr(const SDoublePlane &input,double low_thresold,double high_thresold,const int high_gray_value) {
	
	//convolve with the give, gradient filter
	//SDoublePlane grad_x = get_gradient_filter(3,'X');
	//SDoublePlane grad_y = get_gradient_filter(3,'Y');
	//SDoublePlane grad_I_X = convolve_general(input,grad_x);
	//SDoublePlane grad_I_Y = convolve_general(input,grad_y);
	
	//convolve with the give, sobel filter
	SDoublePlane grad_I_X = sobel_gradient_filter(input,true);
	SDoublePlane grad_I_Y = sobel_gradient_filter(input,false);
	
	int r = input.rows();
	int c = input.cols();
	
	SDoublePlane mag_grad_I_X_Y = compute_grad_magnitude(grad_I_X,grad_I_Y);
	
	return apply_non_maxmal_suppression(mag_grad_I_X_Y,grad_I_X,grad_I_Y,low_thresold,high_thresold,high_gray_value);	
}


SDoublePlane hough_transform(const SDoublePlane &input){
	//now apply the polar coordinates transformation, for accumalator
	//construct the dimesnion for the accumalator
	int im_h = input.rows();
	int im_w = input.cols();
	int max_radius = (int)round(sqrt(2.0)*MAX(im_w,im_h));
	
	//polar coordintae accumalator
	SDoublePlane accumalator(max_radius,180);
	
	double c_x = im_w/2.0;
	double c_y = im_h/2.0;
	
	for(int y=0;y<im_h;y++){
		for(int x=0;x<im_w;x++){
			if(input[y][x]==0){
				continue;
			}
			//changed the degree range from 88 to 92, for assignment, otherwise start it from zero to 180
			for(int deg=0;deg<91;deg++){
				//radius -maxRadius/2 to +maxRadius/2
				double rho = (x-c_x)*cos((double)DEG2RAD(deg))+(y-c_y)*sin((double)DEG2RAD(deg));
				accumalator[(int)round(rho+max_radius/2.0)][deg]++;
			}
		}
		
	}
	writeImage("detect_accumulator_hough.png",accumalator);
	return accumalator;
}

pair<SDoublePlane,vector< pair<int,int> > > start_hough(const SDoublePlane &input,double low_thresold,double high_thresold,const int num_lines,const double threshold_vote_count){
	
	double max_pixel_value = 3;
	SDoublePlane edge_map_image=find_edges_non_max_supr(input,low_thresold,high_thresold,max_pixel_value);
	
	SDoublePlane accumalator = hough_transform(edge_map_image);
	
	//use this accumalator to filter out those lines which have below threshold count
	int max_num_lines = 0;
	int max_voting_count = 0;
	int min_voting_count = numeric_limits<int>::max();
	
	double average_voting_count = 0;
	
	int r_acc = accumalator.rows();
	int c_acc = accumalator.cols();
	double halfMaxRadius = r_acc/2.0;
	
	
	std::priority_queue<hough_data> minHeap_h;
	int non_zero_vote_count = 0;
	
	for(int i=0;i<r_acc;i++){
		for(int j=0;j<c_acc;j++){
			if(max_voting_count<accumalator[i][j]){
				max_voting_count = accumalator[i][j];
			}
			if(accumalator[i][j]!=0 && min_voting_count>accumalator[i][j]){
				min_voting_count = accumalator[i][j];
			}
			if(accumalator[i][j]!=0){
				non_zero_vote_count++;
			}
			average_voting_count+=accumalator[i][j];
		}
	}
	average_voting_count/=(non_zero_vote_count!=0?non_zero_vote_count:(r_acc*c_acc));
	
	//its just an huristic to adjust the number of lines
	int range = abs(max_voting_count-min_voting_count+1);
	int new_threshold_voting = min((int)floor(range/4.0+abs(average_voting_count-min_voting_count+1)/2.0), (int)round(threshold_vote_count));
	
	for(int i=0;i<r_acc;i++){
		for(int j=0;j<c_acc;j++){
			if(accumalator[i][j]>=new_threshold_voting){
				max_num_lines++;
			}		
		}
	}
	//get the minimum number of possible lines
	max_num_lines = max_num_lines==0?num_lines:(int)min(max_num_lines,num_lines);
	max_num_lines = min(r_acc*c_acc,max_num_lines );
	
	
	
	//stats
	cout<<"*****Stats for hough:: balance parameter accordingly to see the result*****"<<endl;
	printf("voting-intended: %lf, max-voting-count-in-hough: %d, min_voting_count:%d and Average voting count:%lf \n",threshold_vote_count, max_voting_count,min_voting_count,average_voting_count);
	cout<<"adjusted heuristic votings::"<<new_threshold_voting<<endl;
	
	for(int i=0;i<r_acc;i++){
		for(int j=0;j<c_acc;j++){
			minHeap_h.push(hough_data(i*c_acc+j,accumalator[i][j]));
			if(minHeap_h.size()==max_num_lines+1){
				minHeap_h.pop();
			}
		}
	}
	
	//cout<<"number of lines::"<<minHeap_h.size()<<", num_lines-intended ::"<<num_lines<<"::(maximum_possible_line - with adjusted/intended voting)"<<max_num_lines<<endl;
	
	vector<hough_data> filtered_max_rho_theta;
	while (!minHeap_h.empty()) {
		filtered_max_rho_theta.push_back(minHeap_h.top());
		minHeap_h.pop();
	}
			
	int im_h = edge_map_image.rows();
	int im_w = edge_map_image.cols();
	
	SDoublePlane final_image(im_h,im_w);
	vector<pair<int,int> > x_y_crdnts_hrzntl_lines;
	//apply the logic for the hough transform for each pixel see if it satisfies the voting count
	
	//optimization -testing :: start  (required for assignments, can be removed for generalized code)
	//get all lines only which has slope of particular type
	
	double slope_epsilon = 2.0d;
	int size = filtered_max_rho_theta.size();
	for(int i=0;i<size;i++){
		int deg = ((int)filtered_max_rho_theta.at(i).rho_theta)%c_acc;
		if(abs(deg-90.0)>slope_epsilon){
			//remove those entry
			filtered_max_rho_theta.erase(filtered_max_rho_theta.begin()+i);
			i--;
			size--;
		}
	}	
	
	//optimization ends
	
	bool found_first = false;
	pair<int,int> last;
	
	for(int y=0;y<im_h;y++){
		found_first = false;
		for(int x=0;x<im_w;x++){
			if(edge_map_image[y][x]==0){
				continue;
			}
			//iterate through priority queue and see if it satisfies for current (x,y) -> (rho and theta)
			for(vector<hough_data>::iterator d=filtered_max_rho_theta.begin();d!=filtered_max_rho_theta.end();++d){
				int rho = (*d).rho_theta/c_acc;
				int deg = ((int)(*d).rho_theta)%c_acc;
				int rho_find = (int)round((x-(im_w/2.0))*cos((double)DEG2RAD(deg))+(y-(im_h/2.0))*sin((double)DEG2RAD(deg))+halfMaxRadius);
				if(rho_find==rho){
					final_image[y][x]=MAX_GRAY_VALUE;
					
					if(!found_first){
						x_y_crdnts_hrzntl_lines.push_back(pair<int,int>(y,x)); //(row,column)
						found_first = true;
					}else{
						last.first = y;
						last.second = x;
					}
					//edge_map_image[y][x]=MAX_GRAY_VALUE;
					break;
				}
			}
		}
		if(found_first){
			x_y_crdnts_hrzntl_lines.push_back(pair<int,int>(last.first,last.second));
		}
	}
	
	//return final_image;
	//return edge_map_image;
	return make_pair(final_image,x_y_crdnts_hrzntl_lines);
}

SDoublePlane get_binary(const SDoublePlane &input,int threshold=128) {
  SDoublePlane output(input.rows(), input.cols());
  for (int i = 0; i < input.rows(); i++)
    for (int j = 0; j < input.cols(); j++)
      output[i][j] = input[i][j] >= threshold ? 1 : 0;
  return output;
}

// Return 1 when less than threshold and reverse when above
SDoublePlane negative_binary(const SDoublePlane &input,int threshold=128) {
  SDoublePlane output(input.rows(), input.cols());
  int rows = input.rows() , cols = input.cols();    
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++) 
      output[i][j] = input[i][j] <= threshold ? 1 : 0;
  return output;
}


void printImage(const SDoublePlane &input) {
  static int t = 1000 ;
  int count = 0 ;
  for (int i = 0; i < input.rows();i++) {
    for (int j = 0 ; j < input.cols(); j++) {
      if (count > t)
		return;
      cout << "Image pixel at " << i << ","<<j << " "<<  input[i][j] << "\n";
	  count++;
    }
  }
}

SDoublePlane reverse_and_flip(const SDoublePlane &input) {
  SDoublePlane output(input.rows(), input.cols());
  SDoublePlane output2(input.rows(), input.cols());  
  int rows = input.rows() , cols = input.cols();
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      output[i][j] = input[rows-i-1][j];
    }
  }
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      output2[i][j] = output[i][cols-j-1];
    }
  }
  return output2;
}


SDoublePlane fl(const SDoublePlane &input) {
  SDoublePlane output(input.rows(), input.cols());
  int rows = input.rows() , cols = input.cols();    
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      output[i][j] = input[i][j] <= fl_threshold ? 0 : 255;
  return output;
}

SDoublePlane add_matrix(const SDoublePlane &input,const SDoublePlane &input2) {
  SDoublePlane output(input.rows(), input.cols());
  int rows = input.rows() , cols = input.cols();    
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      output[i][j] = input[i][j] + input2[i][j];
  return output;
}


SDoublePlane get_hamming(const SDoublePlane &a,const SDoublePlane &b) {
  SDoublePlane a1 = get_binary(a);
  SDoublePlane a2 = negative_binary(a);
  
  SDoublePlane b1 = reverse_and_flip(get_binary(b));
  SDoublePlane b2 = reverse_and_flip(negative_binary(b));  
  // Convolve and add them together 
  return add_matrix(convolve_general(a1,b1) , convolve_general(a2,b2));
}

SDoublePlane get_hamming2(const SDoublePlane &a,const SDoublePlane &b) {
  SDoublePlane output(a.rows(), a.cols());
  SDoublePlane a1 = get_binary(a);
  SDoublePlane b1 = get_binary(b);
  // Convolve and add them together
  for (int i = 0; i < a1.rows() ; i++) {
    for (int j = 0; j < a1.cols(); j++) {
      output[i][j] = 0;
      for (int k = 0; k < b1.rows(); k++) {
	for (int l = 0; l < b1.cols(); l++) {
	  if (i+k < a1.rows() && j+l < a1.cols())
	    output[i][j] += a1[i+k][j+l] * b1[k][l] + (1 - a1[i+k][j+l]) * (1 - b1[k][l]);
	}
      }
    }
  }

  // printImage(output);
  return output;
}

void find_match_and_push_to_vector(const SDoublePlane &w, vector<DetectedSymbol> &symbols,const vector<SDoublePlane> tempList,int n,const SDoublePlane &input,double scale=1) {
  double max = 0,min =-1,count = 0 , total = 0,avg;
  for (int i = 0; i < w.rows() ; i++) {
    for (int j = 0; j < w.cols() ; j++) {
      if (w[i][j] > max)
	max = w[i][j];
      if (w[i][j] < min || min == -1)
	min = w[i][j];
      total += w[i][j];
      if (w[i][j] == 1)
	count++;
    }
  }
  avg = total / (w.rows() * w.cols());
  // double th = 
  // cout <<" Ave=rage " << avg <<" , " << "Max " << max <<" , Min " << min << ", COunt"<< count << "\n";
  double thr = 20;//round((max-avg)/4) ;
  for (int i = 0; i < w.rows() ; i++) {
    for (int j = 0; j < w.cols() ; j++) {
      // 25 is the magic number - Need an algo or use some data structure to ignore nearby points
      if (w[i][j] >= max - thr) {
	DetectedSymbol a;
	a.row = i, a.col = j;
	a.width = (tempList[n].cols())* scale,a.height = (tempList[n].rows()) *scale;

	// Giving a padding for quarter_rest - No idea why prob hits 0 for it in tests
	if (n==1) {
	  a.row-=5 ;
	  a.col-=5;
	  a.width += 10;
	  a.height += 10;
	}
	a.type =(Type) n;
	a.confidence = 1;
	symbols.push_back(a);
	for (int a = -tempList[n].rows()/2; a < tempList[n].rows()/2;a++) {
	  for (int b = -tempList[n].cols()/2; b < tempList[n].cols()/2;b++) {
	    if (i+a >= 0 && i+a < input.rows() && j+b >= 0 && j+b < input.cols()) {
	      //input[i+a][j+b] = 255;
	      w[i+a][j+b] = 0;
	    }
	  }
	}
      }
    }
  }
}
double getAvg(const SDoublePlane &w, vector<DetectedSymbol> &symbols,const vector<SDoublePlane> tempList,int n,const SDoublePlane &input) {
  double max = 0,min =-1,count = 0 , total = 0,avg;
  for (int i = 0; i < w.rows() ; i++) {
    for (int j = 0; j < w.cols() ; j++) {
      if (w[i][j] > max)
	max = w[i][j];
      if (w[i][j] < min || min == -1)
	min = w[i][j];
      total += w[i][j];
      if (w[i][j] == 1)
	count++;
    }
  }
  avg = total / (w.rows() * w.cols());
  // double th = 
  // cout <<" Ave=rage " << avg <<" , " << "Max " << max <<" , Min " << min << ", COunt"<< count << "\n";
  return max - avg;
}

void find_match_and_push_to_vector2(const SDoublePlane &w, vector<DetectedSymbol> &symbols,const vector<SDoublePlane> tempList,int n) {
  int min = -1;
  for (int i = 0; i < w.rows() ; i++) {
    for (int j = 0; j < w.cols() ; j++) {
      if (w[i][j] < min || min == -1)
	min = w[i][j];
    }
  }
  
  for (int i = 0; i < w.rows() ; i++) {
    for (int j = 0; j < w.cols() ; j++) {
      if (w[i][j] < min + 15) {
	DetectedSymbol a;
	a.row = i-5 , a.col = j-5;
	a.width = tempList[n].cols() + 10,a.height = tempList[n].rows() + 10;
	a.type =(Type) n;
	symbols.push_back(a);
      }	
    }
  }
}
/**
   Based on question 4 - Uses hamming distance to find image and pushes detected symbols into vector   
 */
void detect_image(const SDoublePlane &input_image,const vector<SDoublePlane> tempList,vector<DetectedSymbol> &symbols, const SDoublePlane &lineTemp) {  
  vector<DetectedSymbol> symbols3;
  SDoublePlane inp(input_image.rows(),input_image.cols());
  for (int i = 0 ; i < input_image.rows(); i++) {
    for (int j = 0; j < input_image.cols(); j++) {
      inp[i][j] = input_image[i][j];
    }
  }
  for (int n = 0; n < tempList.size(); n++) {    
    SDoublePlane w =get_hamming2(inp,tempList[n]);    
    if (n == 0) {
      SDoublePlane sc = find_edges_non_max_supr(w,0,3,255);
      writeImage("scores4.png",sc);
    }
    find_match_and_push_to_vector(w,symbols,tempList,n,inp);
    // const string name = "b" + to_string(i) + ".png";    
  }
}
// Based on Question 7 uses hamming and other methods
double detect_image3(const SDoublePlane &input_image,const vector<SDoublePlane> tempList,vector<DetectedSymbol> &symbols,const SDoublePlane &lineTemp,const vector<string> fname) {
  SDoublePlane la = get_hamming(input_image,lineTemp);
  vector<DetectedSymbol> symbols3;
  SDoublePlane inp(input_image.rows(),input_image.cols());
  for (int i = 0 ; i < input_image.rows(); i++) {
    for (int j = 0; j < input_image.cols(); j++) {
      inp[i][j] = input_image[i][j];
    }
  }
  // Use first notehead to find scale 
  double scale1 = 1;
  double max = -1;
  double scale = 1;
  for (double a = 0.6; a < 1.8; a+=0.2) {
    SDoublePlane tmp = scaled_image(fname[0].c_str(),a,a);
    SDoublePlane w =get_hamming2(inp,tmp);
    double delta =  getAvg(w,symbols,tempList,0,inp)*25/151/(a*a);
    if (delta > max || max == -1) {	
      max = delta;
      scale = a;
    }      
  }
  scale1 = scale;
  cout << "uSing scale "<< scale  << "\n";
  for (int n = 0; n < tempList.size(); n++) {

    SDoublePlane t1 = scaled_image(fname[n].c_str(),scale,scale);
    SDoublePlane w =get_hamming2(inp,t1);
    find_match_and_push_to_vector(w,symbols,tempList,n,inp,scale);
  }
  return scale1;
}

double distance(int i,int j,int a,int b) {
  return sqrt((i-a)*(i-a) + (j-b)*(j-b));
}

SDoublePlane get_distance_to_closestedge2(const SDoublePlane &inp) {
  SDoublePlane output(inp.rows(), inp.cols());
  for (int i = 0; i < inp.rows(); i++) {
    for(int j = 0; j < inp.cols(); j++) {
      double min = -1;
      if (inp[i][j] > 1) {
      	output[i][j] = 0;
      	continue;
      } else {
	for (int a = 0; a < inp.rows();a++) {
	  for (int b = 0;b < inp.cols(); b++) {
	    if (inp[a][b] > 0) {
	      double t = distance(i,j,a,b);
	      if (t < min || min == -1)
		min = t;
	    }
	  }
	}
	output[i][j] = min;
      }
    }
  }
  return output;
}
// TODO - optimize to use only 2 loops once everything works
SDoublePlane get_distance_to_closestedge(const SDoublePlane &inp) {
  SDoublePlane output(inp.rows(), inp.cols());
  
  SDoublePlane x1(inp.rows(), inp.cols());SDoublePlane x1y(inp.rows(), inp.cols());
  SDoublePlane x2(inp.rows(), inp.cols());SDoublePlane x2y(inp.rows(), inp.cols());
  SDoublePlane y1(inp.rows(), inp.cols());SDoublePlane y1y(inp.rows(), inp.cols());
  SDoublePlane y2(inp.rows(), inp.cols());SDoublePlane y2y(inp.rows(), inp.cols());

  SDoublePlane d1(inp.rows(), inp.cols());SDoublePlane d1y(inp.rows(), inp.cols());
  SDoublePlane d2(inp.rows(), inp.cols());SDoublePlane d2y(inp.rows(), inp.cols());
  SDoublePlane d3(inp.rows(), inp.cols());SDoublePlane d3y(inp.rows(), inp.cols());
  SDoublePlane d4(inp.rows(), inp.cols());SDoublePlane d4y(inp.rows(), inp.cols());
  
  int prevX,prevY;
  int prevX2,prevY2;  
  prevX = -1;prevY = -1;
  for (int j = 0; j < inp.cols(); j++) {
    for (int i = 0 ; i < inp.rows(); i++) {      
      if (inp[i][j] > 0) {
	prevX = i;
	prevY = j;
      } else {
	x1[i][j] = prevX;
	x1y[i][j] = prevY;
      }
    }        
  }
  prevX = -1;prevY = -1;
  for (int j = inp.cols(); j >= 0; j--) {
    for (int i = 0 ; i < inp.rows(); i++) {
      if (inp[i][j] > 0) {
	prevX = i;
	prevY = j;
      } else {
	x2[i][j] = prevX;
	x2y[i][j] = prevY;
      }
    }        
  }

  prevX = -1;prevY = -1;
  for (int i = 0 ; i < inp.rows(); i++) {
    for (int j = 0; j < inp.cols(); j++) {
      if (inp[i][j] > 0) {
	prevX = i;
	prevY = j;
      } else {
	y1[i][j] = prevX;
	y1y[i][j] = prevY;
      }
    }
  }
  
  prevX = -1;prevY = -1;
  for (int i = inp.rows()-1 ; i >=0; i--) {
    for (int j = 0; j < inp.cols(); j++) {
      if (inp[i][j] > 0) {
	prevX = i;
	prevY = j;
      } else {
	y2[i][j] = prevX;
	y2y[i][j] = prevY;
      }
    }
  }

  double diag = distance(0,0,inp.rows(),inp.cols());
  prevX = -1;prevY = -1;
  for (int k = 0; k < inp.rows(); k++) {
    for (int i=inp.rows()-1-k,j=0; j < inp.cols(); i--,j++) {
      if (i >= 0 && i < inp.rows() && j >=0 && j < inp.cols())
      if (inp[i][j] > 0) {
	prevX = i;
	prevY = j;
      } else {
	d1[i][j] = prevX;
	d1[i][j] = prevY;
      }
    }
  }
  prevX = -1;prevY = -1;
  for (int k = 0; k < inp.cols(); k++) {
    for (int i=0,j=k; j < inp.cols(); i--,j++) {
      if (i >= 0 && i < inp.rows() && j >=0 && j < inp.cols())
      if (inp[i][j] > 0) {
	prevX = i;
	prevY = j;
      } else {
	d1[i][j] = prevX;
	d1[i][j] = prevY;
      }
    }
  }
  prevX = -1;prevY = -1;
  for (int k = 0; k < inp.cols(); k++) {
    for (int i=0,j=inp.cols()-1-k; j >=0; i++,j--) {
      if (i >= 0 && i < inp.rows() && j >=0 && j < inp.cols())
      if (inp[i][j] > 0) {
	prevX = i;
	prevY = j;
      } else {
	d2[i][j] = prevX;
	d2y[i][j] = prevY;
      }
    }
  }
  prevX = -1;prevY = -1;
  for (int k = 0; k < inp.rows(); k++) {
    for (int i=k,j=0; i < inp.rows(); i++,j--) {
      if (i >= 0 && i < inp.rows() && j >=0 && j < inp.cols())
      if (inp[i][j] > 0) {
	prevX = i;
	prevY = j;
      } else {
	d2[i][j] = prevX;
	d2y[i][j] = prevY;
      }
    }
  }

  prevX = -1;prevY = -1;
  for (int k = 0; k < inp.rows(); k++) {
    for (int i=k,j=0; i < inp.rows(); i++,j++) {
      if (i >= 0 && i < inp.rows() && j >=0 && j < inp.cols())
	if (inp[i][j] > 0) {
	  prevX = i;
	  prevY = j;
	} else {
	  d3[i][j] = prevX;
	  d3y[i][j] = prevY;
	}
    }
  }

  prevX = -1;prevY = -1;
  for (int k = 0; k < inp.cols(); k++) {
    for (int i=0,j=k; j < inp.cols(); i++,j++) {
      if (i >= 0 && i < inp.rows() && j >=0 && j < inp.cols())
	if (inp[i][j] > 0) {
	  prevX = i;
	  prevY = j;
	} else {
	  d3[i][j] = prevX;
	  d3y[i][j] = prevY;
	}
    }
  }

  prevX = -1;prevY = -1;
  for (int k = 0; k < inp.rows(); k++) {
    for (int i=inp.rows()-1-k,j=inp.cols(); i >=0 ; i--,j--) {
      if (i >= 0 && i < inp.rows() && j >=0 && j < inp.cols())
	if (inp[i][j] > 0) {
	  prevX = i;
	  prevY = j;
	} else {
	  d4[i][j] = prevX;
	  d4y[i][j] = prevY;
	}
    }
  }

  prevX = -1;prevY = -1;
  for (int k = 0; k < inp.cols(); k++) {
    for (int i=inp.rows()-1,j=inp.cols()-1-k; j >=0 ; i--,j--) {
      if (i >= 0 && i < inp.rows() && j >=0 && j < inp.cols())
	if (inp[i][j] > 0) {
	  prevX = i;
	  prevY = j;
	} else {
	  d4[i][j] = prevX;
	  d4y[i][j] = prevY;
	}
    }
  }
  for (int i = 0 ; i < inp.rows(); i++) {
    for (int j = 0; j < inp.cols(); j++) {
      double min = -1,d; int a,b;
      if (inp[i][j] > 0) {
	output[i][j] = 0;
	continue;
      }
      a = x1[i][j]; b = x1y[i][j];
      if (a != -1 && b != -1) {
	d = distance(i,j,a,b);
	if (d < min || min == -1)
	  min = d;
      }
      a = x2[i][j]; b = x2y[i][j];
      if (a != -1 && b != -1) {
	d = distance(i,j,a,b);
	if (d < min || min == -1)
	  min = d;
      }
      a = y1[i][j]; b = y1y[i][j];
      if (a != -1 && b != -1) {
	d = distance(i,j,a,b);
	if (d < min || min == -1)
	  min = d;
      }
      a = y2[i][j]; b = y2y[i][j];
      if (a != -1 && b != -1) {
	d = distance(i,j,a,b);
	if (d < min || min == -1)
	  min = d;
      }
      a = d1[i][j]; b = d2y[i][j];
      if (a != -1 && b != -1) {
	d = distance(i,j,a,b);
	if (d < min || min == -1)
	  min = d;
      }

      a = d2[i][j]; b = d2y[i][j];
      if (a != -1 && b != -1) {
	d = distance(i,j,a,b);
	if (d < min || min == -1)
	  min = d;
      }
      a = d3[i][j]; b = d3y[i][j];
      if (a != -1 && b != -1) {
	d = distance(i,j,a,b);
	if (d < min || min == -1)
	  min = d;
      }
      a = d4[i][j]; b = d4y[i][j];
      if (a != -1 && b != -1) {
	d = distance(i,j,a,b);
	if (d < min || min == -1)
	  min = d;
      }
      output[i][j] = min;
    }
  }
  return output;
}

// TODO - Rename and remove inp
SDoublePlane get_template_match(const SDoublePlane &inp,const SDoublePlane &tmp,const SDoublePlane &D) {
  SDoublePlane f(inp.rows(), inp.cols());
  for (int i = 0; i < inp.rows(); i++) {
    for (int j = 0; j < inp.cols(); j++) {
      f[i][j] = 0;
      for (int k=0; k < tmp.rows(); k++) {
	for (int l = 0; l < tmp.cols();l++) {
	  if(i+k < D.rows() && j+l < D.cols())
	    f[i][j] += D[i+k][j+l] * tmp[k][l];
	  else if (i-k > 0 && j-l > 0)
	    f[i][j] += D[i-k][j-l] * tmp[k][l];
	}
      }
    }
  }
  return f;
}




/**
   Based on question 5 - Detection using sobel and individual convolve
 */
void detect_image2(const SDoublePlane &input_image2,const vector<SDoublePlane> tempList,vector<DetectedSymbol> &symbols, const SDoublePlane &lineTemp) {
  SDoublePlane input_image(input_image2.rows(),input_image2.cols());
  for (int i = 0 ; i < input_image2.rows(); i++) {
    for (int j = 0; j < input_image2.cols(); j++) {
      input_image[i][j] = input_image2[i][j];
    }
  }
  SDoublePlane edgemap = find_edges_non_max_supr(input_image,0,3,255);
  writeImage("edges.png",edgemap);  
  SDoublePlane inp = get_binary(edgemap,1);
  SDoublePlane D = get_distance_to_closestedge(inp);
  fl_threshold = 0;
  SDoublePlane a = fl(inp);
  // writeImage("inE.png",a);
  // writeImage("D.png" ,D);
  
  for (int i = 0; i < tempList.size(); i++) {
    SDoublePlane tmp = get_binary(find_edges_non_max_supr(tempList[i],0,3,255));
    SDoublePlane out = get_template_match(inp,tmp,D);
    find_match_and_push_to_vector2(out,symbols,tempList,i);
  }
}
//
// This main file just outputs a few test images. You'll want to change it to do 
//
void test_convolve() ;
void test_hough_lines(SDoublePlane&,vector<pair<int, int> > &);


void print_hough(const char* filename,const SDoublePlane &input_image,const SDoublePlane &hough_output){
	//make a copy of the input_image
	int r = input_image.rows();
	int c = input_image.cols();
	SDoublePlane copy_input_image(r,c);
	
	//create an other image vector for appropriate rgb values
	
	SDoublePlane image_lines_color(r,c);
	
	for(int y=0;y<r;y++){
		for(int x=0;x<c;x++){
			
			if(hough_output[y][x]){
				image_lines_color[y][x] = MAX_GRAY_VALUE;
				copy_input_image[y][x] = MIN_GRAY_VALUE;
			}else if(!input_image[y][x]){
				copy_input_image[y][x] = MIN_GRAY_VALUE;
				image_lines_color[y][x] = MIN_GRAY_VALUE;
			}else{
				copy_input_image[y][x] = MAX_GRAY_VALUE;
				image_lines_color[y][x] = MAX_GRAY_VALUE;
			}			
		}	
	}	
	SImageIO::write_png_file(filename,copy_input_image,copy_input_image,image_lines_color);
}


SDoublePlane scaled_image(const char* input_file_name,const double scale_x,const double scale_y){
  SDoublePlane temp1= SImageIO::read_png_file(input_file_name);
  QImage img(input_file_name);
  QImage img2 = img.scaled(temp1.cols()*scale_x/*width*/,/*height*/temp1.rows()*scale_y, Qt::KeepAspectRatio);
  const char * output_file_name = "scaled_template.png";
  img2.save(output_file_name,PNG_FOMAT);
  SDoublePlane temp1_2= SImageIO::read_png_file(output_file_name);
  //cout<<"unscaled-image:"<<temp1.rows()<<"::::"<<temp1.cols()<<endl;
  //cout<<"scaled-image"<<temp1_2.rows()<<"::::"<<temp1_2.cols()<<endl;
  return temp1_2;  
}

void calculate_pitch(const SDoublePlane &input_image,vector<DetectedSymbol> &symbols,double scale);
int main(int argc, char *argv[])
//  something more interesting!
{
  if(!(argc == 2))
    {
      cerr << "usage: " << argv[0] << " input_image" << endl;
      return 1;
    }

  string input_filename(argv[1]);
  SDoublePlane input_image= SImageIO::read_png_file(input_filename.c_str());

  double variance = 0.3;
  SDoublePlane gaussian_filter = get_gaussian_filter(variance);
  // //we can perfomr 1-d filter for gaussian too, but for time being we are using the general convolve
  SDoublePlane noiseless_image = convolve_general(input_image,gaussian_filter);
  int low_threshold = 0;
  int high_threshold = 3;
  
  // test step 2 by applying mean filters to the input image
  SDoublePlane mean_filter(3,3);
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      mean_filter[i][j] = 1/9.0;  
  // SDoublePlane output_image = reverse_and_flip(input_image);
  
  // return 1;
  // Read all the 3 template images
  vector<string> fname;
  fname.push_back("template1.png");
  fname.push_back("template2.png");
  fname.push_back("template3.png");
  SDoublePlane temp1= SImageIO::read_png_file("template1.png");
  SDoublePlane temp2= SImageIO::read_png_file("template2.png");
  SDoublePlane temp3= SImageIO::read_png_file("template3.png");
  // This is the line template 
  SDoublePlane lineTemp= SImageIO::read_png_file("template4.png");
  
  vector<SDoublePlane> tempList;
  tempList.push_back(temp1);
  tempList.push_back(temp2);
  tempList.push_back(temp3);
  vector<DetectedSymbol> symbols,symbols2,symbols3;
  
  // Has step a and b 
  detect_image(input_image,tempList,symbols,lineTemp);
  calculate_pitch(input_image,symbols,1);
  write_detection_image("detected4.png", symbols, input_image);
  
  // Has step c and d 
  detect_image2(input_image,tempList,symbols2,lineTemp);
  calculate_pitch(input_image,symbols2,1);
  write_detection_image("detected5.png", symbols2, input_image);
  
  // Step e - Staves.png
  SDoublePlane nms_output_img = find_edges_non_max_supr(noiseless_image,low_threshold,high_threshold,MAX_GRAY_VALUE);
  //write_detection_image("detect_nms.png",symbols,nms_output_img);
  int max_numer_of_lines = 5000;
  int vote_threshold_count = 300;//extreme case- will be all possible lines
  pair<SDoublePlane,vector<pair<int,int> > > hough_output = start_hough(noiseless_image,low_threshold,high_threshold,max_numer_of_lines,vote_threshold_count);
  print_hough("staves.png",input_image,hough_output.first);  

  

  // Step f and g
  double scale = detect_image3(input_image,tempList,symbols3,lineTemp,fname);
  calculate_pitch(input_image,symbols3,scale);
  write_detection_image("detected7.png", symbols3, input_image);
  write_detection_txt("detected7.txt", symbols3);
}

void calculate_pitch(const SDoublePlane &input_image2,vector<DetectedSymbol> &symbols,double scale=1) {
  SDoublePlane input_image(input_image2.rows(),input_image2.cols());
  for (int i = 0 ; i < input_image2.rows(); i++) {
    for (int j = 0; j < input_image2.cols(); j++) {
      input_image[i][j] = input_image2[i][j];
    }
  }
  double variance = 0.3;
  SDoublePlane gaussian_filter = get_gaussian_filter(variance);
  SDoublePlane noiseless_image = convolve_general(input_image,gaussian_filter);
  int low_threshold = 0;
  int high_threshold = 3;
  int max_numer_of_lines = 5000;
  int vote_threshold_count = 300;//extreme case- will be all possible lines
  pair<SDoublePlane,vector<pair<int,int> > > hough_output = start_hough(noiseless_image,low_threshold,high_threshold,max_numer_of_lines,vote_threshold_count);
  vector<pair<int,int> > lineList = hough_output.second;
  vector<pair<int,int> > lines;
  vector<pair<int,int> > smallLines;
  int lastLine = -1;
  // Y coordinate of last line
  int lineDist = 8;
  for (int a = 0; a < lineList.size(); a+=2) {
    int x1 =lineList[a].first , y1 = lineList[a].second ;
    int x2 =lineList[a+1].first , y2 = lineList[a+1].second;
    double length = distance(x1,y1,x2,y2);
    // Avoiding consecutive lines
    if (x1 < lastLine+lineDist && length > input_image.cols()/2) {
      lastLine = x2;
      continue;
    }
    
    if (length > ((double)input_image.cols()*3)/4.0) {
      lines.push_back(lineList[a]);
      lines.push_back(lineList[a+1]);
      lastLine = x2;
      // cout << "added Line is "<< x1  <<" ," <<y1<< "\n";
      // cout << "added Line is "<< x2  <<" ," <<y2<< "\n";
    }
  }
  // for (int a = 0 ; a < lines.size() ; a++) {
  //   cout << lines[a].first <<","<<lines[a].second << "\n";
  // }
  vector<pair<int,int> > flines;
  int sz = lines.size();
  // Filter out lines which don't have any other line in a threshold range
  int max_distance = 13;
  for (int a = 0 ; a < lines.size() ; a+=2) {
    int x1 =lines[a].first , y1 = lines[a].second;
    int x2 =lines[a+1].first , y2 = lines[a+1].second;
    // Check if it has something above or below it - else
    if (a >= 2 && a < lines.size()-1 && x1 - lines[a-1].first < max_distance) {
      flines.push_back(lines[a]);
    } else if(a >= 0 && a+1 < lines.size()) {
      if ((lines[a+2].first) - x1 < max_distance)
	flines.push_back(lines[a]);
    }
  }
  //  cout << "Number of lines is " << flines.size()  << "\n";
  lastLine = 0;
  // for (int a = 0 ; a < flines.size() ; a++) {
  //   cout << flines[a].first <<","<<flines[a].second << "\n";
  // }
  /// Iterate over the symbols
  int difThresh = 3 * scale;
  for (int i = 0; i < symbols.size() ;i++) {
    // Binary search
    int r = symbols[i].row ;int h = symbols[i].height ;
    int pos = (r + h)/2;
    int start = 0,end = flines.size();
    int foundPos = 0;
    while(1) {
      int mid = (end+start)/2;
      if (!(mid >= 0 && mid < flines.size())) {
	foundPos = start *2 + 1;
	break;
      }
      // Try mid , mid + .5 , mid - .5     
      int midV = flines[mid].first;
      int midVUp = (flines[mid].first +  flines[mid+1].first)/2;
      int midVDw = (flines[mid].first + flines[mid+1].first)/2;
      // cout <<"R Value is " << r << " : " << mid << ",midv " << midV<< "\n";
      // cout <<"Value at "<< mid << ":" << midV << " ROw us " << r<< "\n";      
      if (mid + 1 < flines.size()) {
	if (r >= midVUp - difThresh && r <= midVUp + difThresh) {
	  // Then found
	  foundPos = mid *2 + 1;
	  break;
	}
      }
      if (end-start < 2) {
	foundPos = mid;
	break;
      }
      if (mid - 1 >=0) {
	if (r >= midVDw - difThresh && r <= midVDw + difThresh) {
	  // Then found
	  foundPos = mid * 2 - 1;
	  break;
	}
      }
      if (r >= midV - difThresh && r <= midV + difThresh) {
	// Then found
	foundPos = mid * 2;
	break;
      }      
      
      if (flines[mid].first > r) {
	end = mid;
      } else {
	start = mid;
      }
    }
    
    // cout << "Got position " << foundPos  << "\n";   
    if (foundPos > 9)
      foundPos+=2;
    int fp = foundPos;
    fp = fp % 7 + 1;
    switch(fp) {
    case 6 : symbols[i].pitch = 'G';break;
    case 7 : symbols[i].pitch = 'F';break;
    case 1 : symbols[i].pitch = 'E';break;
    case 2 : symbols[i].pitch = 'D';break;
    case 3 : symbols[i].pitch = 'C';break;
    case 4 : symbols[i].pitch = 'B';break;
    case 5 : symbols[i].pitch = 'A';break;
    // case 8 : symbols[i].pitch = 'G';break;
    // case 9 : symbols[i].pitch = 'F';break;
    // case 10 : symbols[i].pitch = 'D';break;
    // case 11 : symbols[i].pitch = 'C';break;
    // default :symbols[i].pitch = 'A';break;
    }
    
    // cout << "Found at " << foundPos << "\n";
  }
  write_detection_image("detect_hough.png",symbols,hough_output.first);
}


void test_convolve() {  
  SDoublePlane a(3,3);SDoublePlane b(3,3);
  
  a[0][0] = 1;a[0][1] = 2;a[0][2] = 3;
  a[1][0] = 4;a[1][1] = 5;a[1][2] = 6;
  a[2][0] = 7;a[2][1] = 8;a[2][2] = 9;
  
  b[0][0] = -1;b[0][1] = -2;b[0][2] = -1;
  b[1][0] = 0;b[1][1] = 0;b[1][2] = 0;
  b[2][0] = 1;b[2][1] = 2;b[2][2] = 1;

  printImage(convolve_general(a,b));
}

void test_hough_lines(SDoublePlane& hg_input,vector<pair<int, int> > &hrzntl_lines){
	
	int r = hg_input.rows();
	int c = hg_input.cols();
	int size = hrzntl_lines.size();
	SDoublePlane hg_line_img(r,c);
	
	for(int i=0;i<size-1;i+=2){
	  pair<int,int> y_x = hrzntl_lines.at(i);
		int x1 = y_x.second;
		int y1 = y_x.first;
		
		y_x = hrzntl_lines.at(i+1);
		
		int x2 = y_x.second;
		int y2 = y_x.first;
		while(x1<c){
		  hg_line_img[y1][x1] = MAX_GRAY_VALUE;
		  x1++;
		}
		while(x2>=0){
			hg_line_img[y2][x2] = MAX_GRAY_VALUE;
			x2--;
		}
	}
	writeImage("detect_hrzntl_hough.png",hg_line_img);
}
// test_convolve();
  // detect_image(input_image,tempList,symbols,lineTemp);
  //detect_image2(input_image,tempList,symbols2,lineTemp);  // write_detection_image("ed1.png", symbols, find_edges_non_max_supr(temp1,0,1,200));
  // write_detection_image("ed11.png", symbols, find_edges(temp1,0));
  // write_detection_image("ed2.png", symbols, find_edges_non_max_supr(temp2,0,3,255));
  // write_detection_image("ed3.png", symbols, find_edges_non_max_supr(temp3,0,3,255));
  // for(int i = 0 ; i < 15; i++) {
  //   string name = "det"+to_string(i)+".png";
  //   write_detection_image(name, symbols, find_edges_non_max_supr(temp3,0,i,255));
  // }
  //write_detection_image("detect2.png", symbols2, input_image);
  // write_detection_image("detect3.png", symbols, find_edges_non_max_supr(temp3,0,3,255));

  //
  // detect_image2(input_image,tempList,symbols2,lineTemp);
  // write_detection_image("ed1.png",symbols3,find_edges_non_max_supr(temp2,0,5,255));
  // randomly generate some detected symbols -- you'll want to replace this
  //  with your symbol detection code obviously!
  //detectImage(input_image,tempList,symbols);
  // write_detection_image("con.png", symbols, output_image);
  // for(int i=0; i<10; i++)
  //   {
  //     DetectedSymbol s;
  //     s.row = rand() % input_image.rows();input_image
  //     s.col = rand() % input_image.cols();
  //     s.width = 20;
  //     s.height = 20;
  //     s.type = (Type) (rand() % 3);
  //     s.confidence = rand();
  //     s.pitch = (rand() % 7) + 'A';
  //     symbols.push_back(s);
  //   }

  //write_detection_txt("detected.txt", symbols);
  //write_detection_image("detected.png", symbols, input_image);



  // write_detection_image("detect_sobel.png", symbols, find_edges(input_image,high_threshold));//noiseless_image
  // write_detection_image("detect_nms.png",symbols,find_edges_non_max_supr(noiseless_image,low_threshold,high_threshold,MAX_GRAY_VALUE));
  //test_hough_lines(hough_output.first,hough_output.second);
  //write_detection_image("detect_sobel.png", symbols, find_edges(input_image,high_threshold));//noiseless_image
