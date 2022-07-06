#include <array>
#include <vector>
#include <functional>

class Interpolator_chebyshev
{
private:
    /* data */
public:
    Interpolator_chebyshev();
    ~Interpolator_chebyshev();

    
    //static std::vector<std::vector<double>>  Interpolator_chebyshev::GetPoints(int N1, int N2, int N3);
    
    static double ChebPol(int deg, double x);
    
    static std::vector<double> GetPoints(int N);
    
    static std::vector<std::vector<std::vector<double> > > 
        GetCoefficients(const int N, const std::function<double(double, double, double)>& fun);
    
    static double InterpolateWCoeffs(std::vector<std::vector<std::vector<double> > > coefs, int N, double x, double y, double z);

    static double function(double x, double y, double z);
};