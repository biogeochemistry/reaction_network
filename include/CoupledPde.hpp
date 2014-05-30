#ifndef COUPLED_PDE_HPP
#define COUPLED_PDE_HPP

#include "BvpPde.hpp"

class CoupledPde {
    private:
        BvpPde *Pde1;
        BvpPde *Pde2;
        int counter;
        double mdt, mT;
        void WriteSolution();

    public:
        CoupledPde(BvpPde *FirstPde, BvpPde *SecondPde);
        void SolveSystem();
};

#endif // COUPLED_PDE_HPP
