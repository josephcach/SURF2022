#include "Piecewise_trilinear.h"




std::vector<std::vector<std::vector<double> > > Piecewise_trilinear::GetCorners(const int N, 
   const std::function<double(double, double, double)>& fun)
{
  std::vector<std::vector<std::vector<double> > > corners(N, std::vector<std::vector<double> >(N, std::vector<double>(N, 0))); 
  std::vector<double> corner_coords(N);
  double h = 2.0 / (double) (N-1);
  for(int i=0;i<N;i++)
  {
    corner_coords[i] = -1 + i*h;
  }
  
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      for(int k=0;k<N;k++){
        corners[i][j][k] = fun(corner_coords[i],corner_coords[j],corner_coords[k]);
      }
    }
  }
  return corners;
}

double Piecewise_trilinear::Interpolate(const int N, const std::vector<std::vector<std::vector<double> > >&  corners, 
  const double x, const double y, const double z )
  {
    double h = 2.0 / (double) (N-1);
    int xind =  (x+1)/h;
    int yind =  (y+1)/h;
    int zind =  (z+1)/h;
    double xd = (x - (xind*h-1))/h;
    double yd = (y - (yind*h-1))/h;
    double zd = (z - (zind*h-1))/h;

    double c00 =  corners[xind][yind][zind]*(1-xd)+corners[xind+1][yind][zind]*xd;
    double c01 =  corners[xind][yind][zind+1]*(1-xd)+corners[xind+1][yind][zind+1]*xd;
    double c10 =  corners[xind][yind+1][zind]*(1-xd)+corners[xind+1][yind+1][zind]*xd;
    double c11 =  corners[xind][yind+1][zind+1]*(1-xd)+corners[xind+1][yind+1][zind+1]*xd;

    double c0 = c00*(1-yd)+c10*yd;
    double c1 = c01*(1-yd) +c11*yd;
    return c0*(1-zd)+c1*zd;
  }


