/*
 execute as :
 .x TIPP_code.C++("training_23Nedd_TIPP_Zm2_Shift.root",1, 12, 12 )              for Det 12 strip F12
 .x TIPP_code.C++("training_23Nedd_TIPP_Zm2_Shift.root",2, 10, 30 )              for Det 10 Strip B30
*/

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMarker.h"
#include "TH2F.h"
#include "TMath.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TGraphErrors.h"
#include <vector>
#include <iostream>
#include <numeric>
#include <fstream>

//Global Variables


Int_t gmethod = 1 ; //to specify the selection between a front strip (gmethod=1) and a back strip (gmethod=2)
Int_t nb_pixels ; //for selected front strip we have 48 pixels, for a back strip we have 24 pixels
TString gDataFile ; //The file containing the simulated data
TString Kinematic_file = "datadd.txt"; //the file containing the calculated data
vector<vector<Double_t>> gMean_strip; // vector containing a set of vectors
                                      //each one of this set represents a pixel and contains mean values of x,y,z,Energy and Theta for this pixel
TGraph* gKline; //where the calculated data is stocked
double gLimitXmin = -1 ;
double gLimitXmax = +1 ;
double gLimitYmin = -1 ;
double gLimitYmax = +1 ;
double gLimitZmin = -1 ;
double gLimitZmax = +6 ;

//DECLARATION OF FUNCTIONS

vector<vector<Double_t>> Crunsh_Strip(Int_t Det_selected, Int_t Strip_selected, Int_t nb_pixels);
        //takes the selected detector,strip and pixel's number on the strip as parammeters and returns means values for each pixel as set of vectors
void add_elements(vector<vector<vector<Double_t>>>*Strip_big_vector, Int_t index, Double_t x, Double_t y, Double_t z, Double_t ELab, Double_t Theta);
        //fills the big vector with all events informations for the selected strip
Double_t mean_vector (vector<Double_t> vector_samples);
        //takes a vector of double values (x values per exemple) as parameter and returns the mean value of these values
Double_t std_vector(vector<Double_t> vector_samples, Double_t mean);
        //takes the vector of double values and the mean of this vector as parameters and returns the standard_deviation
Double_t theta_cal (Double_t x, Double_t y, Double_t z);
        //returns the theta calculated for a selected position
Double_t chi2 (const Double_t *param);
        //returns the chi square value
TGraph* GetTheoreticaLine (TString name_file);
        //returns the kinematic graph (theoretical line) from selected file
Double_t GetThetaCritical(vector<vector<Double_t>> vector_final, Int_t nb_pixels);
        //returns the value of theta having the maximum energy
const Double_t* Minimization(std::string name_minimizer, Double_t (*Func)(const Double_t*), Int_t dim);
        //returns the final values after the minimzization and takes the function to follow for the minimization and the number of variables as parameters



//###############  MAIN Function  ###############

        //when excuting the code the name of file containing the simulated data, the method precising if the selected strip is a front strip (1) or
        //back strip(2), the selected detector number (12 for detector on x>0 and 10 for the one on x<0), and the selected strip number should be
        //given as parameters
void TIPP_code(TString name ,const Int_t method , const Int_t Det_selected,const Int_t Strip_selected)
{
  gmethod=method; //put in global variable
  gDataFile=name; //put in global variable

  if (method == 1) //front strip selected
  {
    nb_pixels = 49;
  }
  else if (method == 2) //back strip selected
        {
        nb_pixels = 24;
        }
        else
                {
                 cout << "Please select an available type of method (1 to select a Front strip and 2 to select a Back strip)" << endl;
                }


 //Open simulated data file
  TFile* File = new TFile (gDataFile);
 //Charge the kinematic line (theoritical)
  gKline = GetTheoreticaLine(Kinematic_file);

 //Collect of data (means and std for each pixel) from one strip
 //gMean_strip = < <x_0, y_0, z_0, xstd_0, ystd_0, zstd_0, E_0, Estd_0, Theta_0, Thetastd_0>, <x_2, y_2, ...>, ... >
  gMean_strip = Crunsh_Strip(Det_selected, Strip_selected, nb_pixels);

//for(Int_t i =0; i< (Int_t) gMean_strip.size();i++)
//  {cout <<" x " <<gMean_strip.at(i).at(0) <<" y " <<gMean_strip.at(i).at(1) <<" z " <<gMean_strip.at(i).at(2) <<endl; }

 //Minimization of colllected data from gMean_strip
  const Double_t *xs = Minimization("Minuit2", &chi2, 3);

 //DISPLAY

  //Draw errors
  TCanvas *c1 = new TCanvas("initial and shifted values","initial and shifted values", 250, 100, 800, 500);
  vector<Double_t> ELab_vector;
  vector<Double_t> Theta_vector;
  vector<Double_t> ELab_std_vector;
  vector<Double_t> Theta_std_vector;
  vector<Double_t> Theta_vector_shift;
  //vector<Double_t> Theta_std_vector_shift;

//TO COMMENT AFTER

    for (Int_t i = 0; i < (Int_t) gMean_strip.size(); i ++)
{
        ELab_vector.push_back( (gMean_strip.at(i)).at(6) );
        ELab_std_vector.push_back( (gMean_strip.at(i)).at(7) );
        Theta_vector.push_back( (gMean_strip.at(i)).at(8) );
        Theta_std_vector.push_back( (gMean_strip.at(i)).at(9) );

        Double_t a = theta_cal( (gMean_strip.at(i)).at(0) + xs[0], (gMean_strip.at(i)).at(1)+ xs[1], (gMean_strip.at(i)).at(2) + xs[2] );
        Theta_vector_shift.push_back( a );
      }



 //Initial value graph
    TGraphErrors* Graph_errors_init = new TGraphErrors(gMean_strip.size(), &Theta_vector[0], &ELab_vector[0], &Theta_std_vector[0], &ELab_std_vector[0]);
    Graph_errors_init->SetMarkerColor(2);
    Graph_errors_init->SetMarkerStyle(20);
    Graph_errors_init->SetMarkerSize(0.5);
    Graph_errors_init->SetTitle("initial angles");
 //shifted values graph
    TGraphErrors* Graph_errors_shift = new TGraphErrors(gMean_strip.size(), &Theta_vector_shift[0], &ELab_vector[0], &Theta_std_vector[0], &ELab_std_vector[0]);
    Graph_errors_shift->SetMarkerColor(4);
    Graph_errors_shift->SetMarkerStyle(20);
    Graph_errors_shift->SetMarkerSize(0.5);
    Graph_errors_shift->SetTitle("shifted angles");
 //kinematic line options
    gKline->SetTitle("Theoritical line");
    gKline->SetLineColor(1);

 //Plotting
    c1->cd();
    gKline->Draw("ac");
    gKline->GetXaxis()->SetRangeUser(50,80);
    gKline->GetYaxis()->SetRangeUser(0,22);
    Graph_errors_init->Draw("p");
    Graph_errors_shift->Draw("p");
    c1->BuildLegend();

}







//###############  FUNCTIONS  ###############

vector<vector<Double_t>> Crunsh_Strip(Int_t Det_selected , Int_t Strip_selected , Int_t nb_pixels)
{
 //Crunsh_Strip's variables
  Double_t ELab, x, y, z, Theta, Theta_calc ,Theta_crit;
  Int_t Front, Back, Det;
  TVector3 p2, p3; //the 2 vectors to calculate theta (scalar product)
  vector<vector<vector<Double_t>>> Strip_big_vector (nb_pixels); //Declaration of a global vector which gather a whole strip
  vector<Double_t> one_dim_vector; // vector of double values, these values can be x , or y or z values for example
  vector<vector<Double_t>> vector_final(nb_pixels); //the final vector containing the mean and std values for all the pixels
  vector<vector<Double_t>> vector_final_corrected(nb_pixels); //the corrected final vector after elimination of pixels without events

 //Load the TFile
  TFile* File = new TFile (gDataFile);
  if (!File->IsOpen())
  {
    cout << "file not found! " << endl; //Error
  }

 //Load the TTree
  TTree* Tree = (TTree*)File->Get("SimDataTree");
  if (!Tree)
  {
    cout << "SimDataTree tree not found! " << endl; //Error
  }
 //Number of entries in our Tree
  Int_t nentries = Tree->GetEntries();

 //pointing on branches values
  Tree->SetBranchAddress("ELab", &ELab);
  Tree->SetBranchAddress("SharcX", &x);
  Tree->SetBranchAddress("SharcY", &y);
  Tree->SetBranchAddress("SharcZ", &z);
  Tree->SetBranchAddress("ThetaLab", &Theta);
  Tree->SetBranchAddress("SharcFront", &Front);
  Tree->SetBranchAddress("SharcBack", &Back);
  Tree->SetBranchAddress("SharcDet", &Det);

 //defining the number of vectors for each pixel (subvector) in the Strip big vector
  for(Int_t i = 0; i<nb_pixels; i++)
  {
    for (Int_t j = 0; j<5; j++) //5 values to store : x, y, z, ELab, Theta, for each pixel
    {
      (Strip_big_vector.at(i)).push_back(one_dim_vector);
    }
  }

 //Crunshing
  for (Int_t i = 0; i<nentries; i++)
  {
   //Pointing on the current entry
    Tree->GetEntry(i);

   if (gmethod == 1) //selected front strip
    {
      if (Det==Det_selected && Front==Strip_selected)
      {
       Theta_calc = theta_cal(x,y,z); //Compute the angle
       add_elements(&Strip_big_vector, Back, x, y, z, ELab, Theta_calc); //fill the selected front strip big vector with values
      }
    }
    if (gmethod == 2){
      //Selection of one pixel
      if (Det==Det_selected && Back==Strip_selected && Front>1 &&Front<23)
      {
       Theta_calc = theta_cal(x,y,z); //Compute the angle
        add_elements(&Strip_big_vector, Front, x, y, z, ELab, Theta_calc); //fill the selected back strip big vector with values
      }
    }
   }

 //Means and standard deviations of different variables (x, y, z, ELab and Theta) calculation and storage

  for(Int_t i = 0; i<nb_pixels; i++)
  {
   //Compute the mean with the sorted data of vector_samples
    Double_t x_mean = mean_vector(  (Strip_big_vector.at(i)).at(0)  );
    Double_t y_mean = mean_vector(  (Strip_big_vector.at(i)).at(1)  );
    Double_t z_mean = mean_vector(  (Strip_big_vector.at(i)).at(2)  );
    Double_t ELab_mean = mean_vector(  (Strip_big_vector.at(i)).at(3)  );
    Double_t Theta_mean = mean_vector(  (Strip_big_vector.at(i)).at(4)  );

   //Standard deviation
    Double_t x_std = std_vector( (Strip_big_vector.at(i)).at(0)  , x_mean);
    Double_t y_std = std_vector( (Strip_big_vector.at(i)).at(1)  , y_mean);
    Double_t z_std = std_vector( (Strip_big_vector.at(i)).at(2)  , z_mean);
    Double_t ELab_std = std_vector( (Strip_big_vector.at(i)).at(3)  , ELab_mean);
    Double_t Theta_std = std_vector( (Strip_big_vector.at(i)).at(4)  , Theta_mean);


   //Storage into vector_final
    (vector_final.at(i)).push_back(x_mean);//0
    (vector_final.at(i)).push_back(y_mean);//1
    (vector_final.at(i)).push_back(z_mean);//2
    (vector_final.at(i)).push_back(x_std); //3
    (vector_final.at(i)).push_back(y_std); //4
    (vector_final.at(i)).push_back(z_std); //5
    (vector_final.at(i)).push_back(ELab_mean);//6
    (vector_final.at(i)).push_back(ELab_std); //7
    (vector_final.at(i)).push_back(Theta_mean);//8
    (vector_final.at(i)).push_back(Theta_std); //9
}

  //Compute theta critical
   Theta_crit = GetThetaCritical(vector_final, nb_pixels);

  //filling the corrected final vector (useful one)
  Double_t counter = 0;
  for(Int_t i = 0; i<nb_pixels; i++)
  {
    Double_t Theta_current = vector_final.at(i).at(8);
    if (Theta_current > Theta_crit)
      {
        (vector_final_corrected.at(counter)).push_back( vector_final.at(i).at(0) );//0
        (vector_final_corrected.at(counter)).push_back( vector_final.at(i).at(1) );//1
        (vector_final_corrected.at(counter)).push_back( vector_final.at(i).at(2) );//2
        (vector_final_corrected.at(counter)).push_back( vector_final.at(i).at(3) ); //3
        (vector_final_corrected.at(counter)).push_back( vector_final.at(i).at(4) ); //4
        (vector_final_corrected.at(counter)).push_back( vector_final.at(i).at(5) ); //5
        (vector_final_corrected.at(counter)).push_back( vector_final.at(i).at(6) );//6
        (vector_final_corrected.at(counter)).push_back( vector_final.at(i).at(7) ); //7
        (vector_final_corrected.at(counter)).push_back( vector_final.at(i).at(8) );//8
        (vector_final_corrected.at(counter)).push_back( vector_final.at(i).at(9) ); //9

        counter += 1;
      }
  }
  vector_final_corrected.resize(counter);

  return vector_final_corrected;
}


void add_elements(vector<vector<vector<Double_t>>>*  Strip_big_vector, Int_t index, Double_t x, Double_t y, Double_t z, Double_t ELab, Double_t Theta)
{
  ((Strip_big_vector->at(index)).at(0)).push_back(x);
  ((Strip_big_vector->at(index)).at(1)).push_back(y);
  ((Strip_big_vector->at(index)).at(2)).push_back(z);
  ((Strip_big_vector->at(index)).at(3)).push_back(ELab);
  ((Strip_big_vector->at(index)).at(4)).push_back(Theta);
}


Double_t mean_vector(vector<Double_t> vector_samples)
{
  Double_t mean = accumulate(vector_samples.begin(), vector_samples.end(), 0.0)/(vector_samples.size());
  return mean;
}


Double_t std_vector(vector<Double_t> one_dim_vector, Double_t mean)
{
  Double_t sq_sum = inner_product(one_dim_vector.begin(), one_dim_vector.end(), one_dim_vector.begin(), 0.0);
  Double_t stdev = sqrt(sq_sum / one_dim_vector.size() - mean * mean);
  return stdev;
}


Double_t theta_cal (Double_t x, Double_t y, Double_t z)
{
  TVector3 p2, p3;
  p2.SetXYZ(0, 0, z); p3.SetXYZ(x, y, z); //Value of the 2 vectors
  Double_t Theta_calc = p2.Angle(p3)*180/(TMath::Pi()); //Compute the angle
  return Theta_calc;
}


Double_t chi2 (const Double_t *param)
{
  Double_t chi2_sum = 0;

  for (Int_t i = 0; i < (Int_t) gMean_strip.size(); i++)
  {
   //New parameters with shift
    Double_t x = (gMean_strip.at(i)).at(0) + param[0];
    Double_t y = (gMean_strip.at(i)).at(1) + param[1];
    Double_t z = (gMean_strip.at(i)).at(2) + param[2];
    Double_t ELab = (gMean_strip.at(i)).at(6);
   //Std on ELab
    Double_t ELab_std = (gMean_strip.at(i)).at(7);
   //Compute the new Theta_calc
    Double_t Theta_calc_new = theta_cal(x, y, z);
   //Evaluate the value with theoretical line at position Theta_calc_new
    Double_t ELab_theoric = gKline->Eval(Theta_calc_new);

      Double_t chi2 = ((ELab - ELab_theoric)/ELab_std) * ((ELab - ELab_theoric)/ELab_std); //Little chi2
      chi2_sum += chi2;

  }
  return (chi2_sum);
}


Double_t GetThetaCritical(vector<vector<Double_t>> Mean_Strip, Int_t nb_pixels)
{
  Double_t max_Energy = 0; //starting value
  Double_t current_Energy = 0; //starting value
  Int_t index = 0; //starting value

  for (Int_t i = 0; i<nb_pixels; i++)
  {
    current_Energy = Mean_Strip.at(i).at(6);
    if (current_Energy > max_Energy)
    {
      max_Energy = current_Energy;
      index = i;
    }
  }
  Double_t Theta_crit = Mean_Strip.at(index).at(8);
  return Theta_crit;
}


const Double_t* Minimization(std::string name_minimizer ="Minuit2", Double_t (*Func)(const Double_t*) = &chi2, Int_t dim = 3)
{
  ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer(name_minimizer, "");

 // set tolerance , etc...
  minimum->SetMaxFunctionCalls(10000000);
  minimum->SetTolerance(0.000001);
  minimum->SetPrintLevel(1);


 // create function wrapper for minimizer
  ROOT::Math::Functor f(Func, dim);

 //Step for each variables (x, y, z)
  double step[3] = {0.001, 0.001, 0.001};

  minimum->SetFunction(f);

 // Set the free variables to be minimized !
  minimum->SetVariable(0,"delta_x",0, step[0]);
  minimum->SetVariable(1,"delta_y",0, step[1]);
  minimum->SetVariable(2,"delta_z",0, step[2]);

 //Limits
  minimum->SetVariableLimits(0,gLimitXmin,gLimitXmax);
  minimum->SetVariableLimits(1,gLimitYmin,gLimitYmax);
  minimum->SetVariableLimits(2,gLimitZmin,gLimitZmax);

 // do the minimization
  minimum->Minimize();

  const double *xs = minimum->X();
  cout << "Minimum: f(" << xs[0] << "," << xs[1] << "," << xs[2] << "): "
            << minimum->MinValue()  << endl;

  return xs;
}


TGraph* GetTheoreticaLine (TString name_file = "datadd.txt")
{
  Int_t n = 101;
  Float_t x[n], y[n];
  ifstream fichier(name_file);

  if (fichier)
  {
    for (Int_t i=0; i<n; i++)
    {
      fichier >> x[i];
      fichier >> y[i];
    }
  }

  TGraph *reaction_dd = new TGraph (n, x, y);

  return reaction_dd;
}
