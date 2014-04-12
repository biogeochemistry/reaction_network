#ifndef RUNGEKUTTA_HPP
#define RUNGEKUTTA_HPP

#include "OdeSolver.hpp"

class RungeKutta : public OdeSolver {
    public:
        double SolveEquation();
};



#endif // RUNGEKUTTA_HPP
