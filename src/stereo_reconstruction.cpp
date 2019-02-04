#include <SImage.h>
#include <SImageIO.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <math.h>
#include <ctime>

#include "mrf.h"
#include "GraphCuts.h"
#include "LoopyBP.h"

using namespace std;

//////////////////////////////////////////////////
// Define and Global Variables
//
const int NUM_OF_LABELS = 80;
const int WINDOW_SIZE = 5;

// Input Images (Left and Right)
SDoublePlane image_left;
SDoublePlane image_right;

//////////////////////////////////////////////////
// Cost Functions
unsigned int data_cost(int x, int y, int label)
{
  unsigned int data_cost = 0;

  // window size
  unsigned int window_size = WINDOW_SIZE;

  unsigned int left_val = 0;
  unsigned int right_val = 0;
  unsigned int sum_val = 0;

  for (int dy = -(int)(window_size / 2); dy <= (int)(window_size / 2); dy++) {
    for (int dx = -(int)(window_size / 2); dx <= (int)(window_size / 2); dx++) {
      left_val = image_left[y + dy][x + dx];
      right_val = image_right[y + dy][x + dx - label];
      sum_val += (left_val - right_val) * (left_val - right_val);
    }
  }

  data_cost = sum_val / (window_size * window_size);

  //cout << "x: " << x << ", y: " << y <<" data_cost: " << data_cost << endl;

  return data_cost;
}

SDoublePlane compute_disp_using_block_matching()
{ 
  // image container to store result 
  SDoublePlane result(image_left.rows(), image_left.cols());

  //cout << total << " " << image_left.rows() << " " << image_left.cols() << endl;
  // to assign best label which has minmum cost
  unsigned int min_cost = UINT_MAX;
  // current cost
  unsigned int cost = 0;

  // image width
  unsigned int y = 0;
  // image height
  unsigned int x = 0;

  // initialize
  int dIndex = 0;
  for (int i = 0; i < image_left.rows(); i++) {
    for (int j = 0; j <image_right.cols(); j++) {
      result[i][j] = 0; 
    }
  }

  // window size
  unsigned int window_size = WINDOW_SIZE;
  for (int y = (int)(window_size / 2); y < image_left.rows() - (int)(window_size / 2); y++) {
    for (int x = (int)(window_size / 2); x < image_right.cols() - (int)(window_size / 2); x++) {

      // init min_cost
      min_cost = UINT_MAX;
      for (int j = 0; j < NUM_OF_LABELS; j++) {
        //if (x - j > 0) {
            cost = data_cost(x, y, j);
            // assign new label
            if (cost < min_cost) {
              min_cost = cost;
              //cout << "min_cost: " << min_cost << " ,j = " << j << endl;
              result[y][x] = j * (255/NUM_OF_LABELS);
              //cout << "result[" << y << "][" << x << "]: " << result[y][x] << endl;
            }
        //}
      }
    }
  }

  return result;
}

void compute_data_cost(MRF::CostVal *&data_cost_val)
{
  // width, hegiht, number of labels
  data_cost_val = new MRF::CostVal[image_left.cols() * image_left.rows() * NUM_OF_LABELS];

  // number of pixels without borders
  unsigned int total = 0;
  // current cost
  unsigned int cost = 0;

  unsigned int left_val = 0;
  unsigned int right_val = 0;

  // initialize
  int dIndex = 0;
  for (int y = 0; y < image_left.rows(); y++) {
    for (int x = 0; x <image_left.cols(); x++) {
      for (int k = 0; k < NUM_OF_LABELS; k++) {
        data_cost_val[dIndex++] = UINT_MAX; 
      }
    }
  }

  // window size
  unsigned int window_size = WINDOW_SIZE;
  for (int y = (int)(window_size / 2); y < image_left.rows() - (int)(window_size / 2); y++) {
    for (int x = (int)(window_size / 2); x < image_right.cols() - (int)(window_size / 2); x++) {
      for (int k = 0; k < NUM_OF_LABELS; k++) {
        //if (x - k > 0) {
            cost = data_cost(x, y, k);
            data_cost_val[y * image_left.cols() * NUM_OF_LABELS + x * NUM_OF_LABELS + k] = cost;
        //}
      }
    }
  }

}


SDoublePlane getDisparities(MRF *mrf)
{ 
  // image container to store result 
  SDoublePlane result(image_left.rows(), image_left.cols());

  int n = 0;

  for (int i = 0; i < image_left.rows(); i++) {
    for (int j = 0; j <image_right.cols(); j++) {
      result[i][j] = mrf->getLabel(n++) * (255/NUM_OF_LABELS);
    }
  }

  return result;
}

void mesureErrors(SDoublePlane &out_disp, SDoublePlane &gt_disp)
{
  // Measure error with respect to ground truth
  double err = 0;

  for (int i = 0; i < gt_disp.rows(); i++) {
    for (int j = 0; j < gt_disp.cols(); j++) {
      err += sqrt((out_disp[i][j] - gt_disp[i][j])*(out_disp[i][j] - gt_disp[i][j]));
    }
  }

  cout << "MRF stereo technique mean error = " << err/gt_disp.rows()/gt_disp.cols() << endl;

}

//////////////////////////////////////////////////
// MRF Algorithms
//
// Loopy Belife Propagation
MRF* loopy_belief_propagation(int img_width, int img_height, int num_of_labels, EnergyFunction *energy);
// Graph Cuts: Alpha-Beta Swap
MRF* gt_alpha_beta_swap(int img_width, int img_height, int num_of_labels, EnergyFunction *energy);
// Graph Cuts: Alpha Expansion
MRF* gt_alpha_expansion(int img_width, int img_height, int num_of_labels, EnergyFunction *energy);

// Main Function
int main(int argc, char *argv[])
{
  cout << "MRF based stereo reconstruction start!" << endl;

  if (argc != 5 && argc != 4) {
      cerr << "usage: " << argv[0] << " algorithm image_file1 image_file2 [gt_file]" << endl;
      return 1;
  }

  string alg_name = argv[1];
  string input_filename1 = argv[2], input_filename2 = argv[3];
  string gt_filename;

  // read in images and gt
  image_left = SImageIO::read_png_file(input_filename1.c_str());
  image_right = SImageIO::read_png_file(input_filename2.c_str());
  SDoublePlane gt;

  if (argc == 5)
    gt_filename = argv[4];

  if (gt_filename != "") {
    gt = SImageIO::read_png_file(gt_filename.c_str());
    // gt maps are scaled by a factor of 3, undo this...
    for(int i=0; i<gt.rows(); i++)
      for(int j=0; j<gt.cols(); j++) {
        gt[i][j] = gt[i][j]/3.0 * (255/NUM_OF_LABELS);
        //cout << gt[i][j] << endl;
      }
  }

  MRF *mrf = NULL;
  
  // width, hegiht, number of labels
  // compute data cost
  MRF::CostVal *data_cost_val = NULL;
  compute_data_cost(data_cost_val);

  MRF::CostVal hCue[image_left.cols()*image_left.rows()];
  MRF::CostVal vCue[image_left.cols()*image_left.rows()];

  // generate function
  MRF::CostVal* ptr;
  for (ptr=&hCue[0]; ptr<&hCue[image_left.cols()*image_left.rows()]; ptr++) *ptr = rand() % 3;
  for (ptr=&vCue[0]; ptr<&vCue[image_left.cols()*image_left.rows()]; ptr++) *ptr = rand() % 3;

  MRF::CostVal smoothMax = (MRF::CostVal)25, lambda = (MRF::CostVal)2;

  // Energy Function
  DataCost *data         = new DataCost(data_cost_val);
  SmoothnessCost *smooth = new SmoothnessCost(1,smoothMax,lambda);
  EnergyFunction *energy = new EnergyFunction(data,smooth);

  if (alg_name == "block") {
    cout << "Block matching " << endl;
    int start_s = clock();
    SDoublePlane block_disp = compute_disp_using_block_matching();
    int stop_s = clock();
    cout << "time: " << (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000 << endl;
    SImageIO::write_png_file("block_disp.png", block_disp, block_disp, block_disp);
    //SImageIO::write_png_file("gt_disp.png", gt, gt, gt);
    if (gt_filename != "")
      mesureErrors(block_disp, gt);
  } else if (alg_name == "lbp") {
    int start_s = clock();
    mrf = loopy_belief_propagation(image_left.cols(),image_left.rows(),NUM_OF_LABELS,energy);
    SDoublePlane gt_lbp_disp = getDisparities(mrf);
    int stop_s = clock();
    cout << "time: " << (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000 << endl;
    SImageIO::write_png_file("lbp_disp.png", gt_lbp_disp, gt_lbp_disp, gt_lbp_disp);
    if (gt_filename != "")
      mesureErrors(gt_lbp_disp, gt);
  } else if (alg_name == "gt_swap") {
    int start_s = clock();
    mrf = gt_alpha_beta_swap(image_left.cols(),image_left.rows(),NUM_OF_LABELS,energy);
    SDoublePlane gt_swap_disp = getDisparities(mrf);
    int stop_s = clock();
    cout << "time: " << (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000 << endl;
    SImageIO::write_png_file("gt_swap_disp.png", gt_swap_disp, gt_swap_disp, gt_swap_disp);
    if (gt_filename != "")
      mesureErrors(gt_swap_disp, gt);
  } else if (alg_name == "gt_exp") {
    int start_s = clock();
    mrf = gt_alpha_expansion(image_left.cols(),image_left.rows(),NUM_OF_LABELS,energy);
    SDoublePlane gt_alpha_disp = getDisparities(mrf);
    int stop_s = clock();
    cout << "time: " << (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000 << endl;
    SImageIO::write_png_file("gt_alpha_disp.png", gt_alpha_disp, gt_alpha_disp, gt_alpha_disp);
    if (gt_filename != "")
      mesureErrors(gt_alpha_disp, gt);
  } else {
    cout << "Check the name of algorithm! It should be one of the followings:" << endl;
    cout << ">> block, lbp, gt_swap, gt_exp" << endl;

    // release memory
    delete energy;
    delete smooth;
    delete data;
    delete [] data_cost_val;

    return 0;
  } 

  cout << "MRF Stereo Reconstruction Completed!" << endl;

  // release memory
  if (mrf)
    delete mrf;

  delete energy;
  delete smooth;
  delete data;
  delete [] data_cost_val;

  return 0;
}

MRF* loopy_belief_propagation(int img_width, int img_height, int num_of_labels, EnergyFunction *energy)
{
  MRF *mrf;

  MRF::EnergyVal E;

  float t, tot_t;
  int iter;

  mrf = new LoopyBP(img_width, img_height, num_of_labels, energy);

  cout << "Loopy Belief Propagation." << endl;

  mrf->initialize();
  mrf->clearAnswer();

  E = mrf->totalEnergy();
  printf("Energy at the Start = %g (%g,%g)\n", (float)E,
      (float)mrf->smoothnessEnergy(), (float)mrf->dataEnergy());

  tot_t = 0;

  mrf->optimize(50, t);
  E = mrf->totalEnergy();
  tot_t = tot_t + t ;
  printf("energy = %g (%f secs)\n", (float)E, tot_t);

  return mrf;
}

MRF* gt_alpha_beta_swap(int img_width, int img_height, int num_of_labels, EnergyFunction *energy)
{
  MRF *mrf;

  MRF::EnergyVal E;

  float t, tot_t;
  int iter;
  mrf = new Swap(image_left.cols(),image_left.rows(),NUM_OF_LABELS,energy);

  cout << "Graph Cuts: Alpha-Beta Swap" << endl;

  mrf->initialize();
  mrf->clearAnswer();

  E = mrf->totalEnergy();
  printf("Energy at the Start = %g (%g,%g)\n", (float)E,
      (float)mrf->smoothnessEnergy(), (float)mrf->dataEnergy());

  tot_t = 0;

  mrf->optimize(100, t);
  E = mrf->totalEnergy();
  tot_t = tot_t + t ;
  printf("energy = %g (%f secs)\n", (float)E, tot_t);

  return mrf;
}

MRF* gt_alpha_expansion(int img_width, int img_height, int num_of_labels, EnergyFunction *energy)
{
  MRF *mrf;

  MRF::EnergyVal E;

  float t, tot_t;
  int iter;
  mrf = new Expansion(image_left.cols(),image_left.rows(),NUM_OF_LABELS,energy);

  cout << "Graph Cuts: Alpha Expansion" << endl;

  mrf->initialize();
  mrf->clearAnswer();

  E = mrf->totalEnergy();
  printf("Energy at the Start = %g (%g,%g)\n", (float)E,
      (float)mrf->smoothnessEnergy(), (float)mrf->dataEnergy());

  tot_t = 0;

  mrf->optimize(100, t);
  E = mrf->totalEnergy();
  tot_t = tot_t + t ;
  printf("energy = %g (%f secs)\n", (float)E, tot_t);


  return mrf;
}
