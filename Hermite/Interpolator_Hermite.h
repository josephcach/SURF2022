#ifndef INTERPOLATOR_HERMITE_H
#define INTERPOLATOR_HERMITE_H

#include <array>
#include <vector>
#include <stdexcept>

///class representing a single Hermite interpolation scheme on [-1, 1]^3
///For simplicity we assume the interpolaiton points are the corners of the cube and we 
///only need to define the number of derivatives (>= 0) in addition to the values are stored
class InterpolatorHermite
{
    ///1d hermite interpolation utilized in the 2d scheme
    ///assumes two interpolation points at -1 and 1
    ///the template parameter M denotes the number of
    ///derivative values per point i.e. >= 0
    template <int M>
    class InterpolatorHermite1D
    {
        ///interpolate function that takes the function values at p1 and p2
        ///which are assumed to sit at -1 and 1, repectively,
        ///and returns the interpolated value at -1 <= target <= 1.
        static double Interpolate(const double value_p1, const std::array<double, M>& der_p1, 
            const double value_p2, const std::array<double, M>& der_p2, 
            const double target)
        {
            throw std::logic_error("Not yet impleneted");

        }
    };

    ///2d hermite interpolation based on the 1d scheme. 
    ///The 3d hermite interpolation utilizes 1d and 2d scheme
    ///scheme assume points at at corners of the square [-1,1]^2
    ///template parameter M determines number of derivative values per point i.e. >= 0
    template <int M>
    class InterpolatorHermite2D
    {
        static double Interpolate(/*put a 2D equivalent of the 1D interface here*/)
        {
            ///TODO: think about the interface here
            throw std::logic_error("Not yet implemented");
        }
    };


private:
    ///need a parameter to determine number of derivative 
    ///values per point.
    ///points are given as corners of cube
    ///this parameter needs to be >= 0
    static constexpr int M_ = 3;

public:
    static double Interpolate(/*put a 3D equivalent of the 1D interface here*/)
    {
        ///TODO: think about the interface here and implement function
        throw std::logic_error("Not yet implemented");
    }


};



#endif