#include "hpclib_global.h"

// MathFuncsDll.h

namespace MathFuncs
{
    class MyMathFuncs
    {
    public:
        // Returns a + b
        static HPCLIB_EXPORT double Add(double a, double b);

        // Returns a - b
        static HPCLIB_EXPORT double Subtract(double a, double b);

        // Returns a * b
        static HPCLIB_EXPORT double Multiply(double a, double b);

        // Returns a / b
        // Throws DivideByZeroException if b is 0
        static HPCLIB_EXPORT double Divide(double a, double b);
    };
}
