#include <iostream>
#include <fstream>
#include <cassert>
#include "BvpOde.hpp"


void BvpOde::Solve(){
    mpLinearSolver = new LinearSolver(*mpLhsMat, *mpRhsVec);
    solution = mpLinearSolver->SolveLinearSystem();
    WriteSolutionFile();
}

BvpOde::BvpOde(SecondOrderOde* pOde,BoundaryConditions* pBcs, int numNodes){
    mpOde = pOde; 
    mpBconds = pBcs;
    mNumNodes = numNodes;
    mpGrid = new FiniteDifferenceGrid(mNumNodes, pOde->mXmin, pOde->mXmax);
    mpRhsVec = new VectorXd(mNumNodes);
    mpLhsMat = new SparseMatrix<double> (mNumNodes, mNumNodes);
    mFilename = "ode_output.dat";
    PopulateMatrix6thOrder();
    PopulateVector();
    ApplyBoundaryConditions();

}

BvpOde::~BvpOde() {
    delete mpRhsVec;
    delete mpLhsMat;
    delete mpGrid;
}

void BvpOde::PopulateMatrix() {
    for (int i=1; i<mNumNodes-1; i++) {
        double xm = mpGrid->mNodes[i-1].C.x; 
        double x = mpGrid->mNodes[i].C.x;
        double xp = mpGrid->mNodes[i+1].C.x;
        double diffusion_alpha = 2.0/(xp-xm)/(x-xm);
        double diffusion_beta = -2.0/(xp-x)/(x-xm);
        double diffusion_gamma = 2.0/(xp-xm)/(xp-x);
        (*mpLhsMat).insert(i,i-1) = (mpOde->mCoeffOfUxx)*diffusion_alpha - (mpOde->mCoeffOfUx)/(xp-xm);
        (*mpLhsMat).insert(i,i) = (mpOde->mCoeffOfUxx)*diffusion_beta + mpOde->mCoeffOfU;
        (*mpLhsMat).insert(i,i+1) = (mpOde->mCoeffOfUxx)*diffusion_gamma +(mpOde->mCoeffOfUx)/(xp-xm);
    }
}

void BvpOde::PopulateMatrix6thOrder() {
    // This works only with uniform grid!!!
    double h = mpGrid->mNodes[1].C.x - mpGrid->mNodes[0].C.x; 
    double D = mpOde->mCoeffOfUxx;
    double w = mpOde->mCoeffOfUx;
    double k = mpOde->mCoeffOfU;

    int i1 = 1;
    (*mpLhsMat).insert(i1,i1+6) =   11*D/(180*h*h);
    (*mpLhsMat).insert(i1,i1+5) =  -90*D/(180*h*h) +   (+2*w)/(60*h);
    (*mpLhsMat).insert(i1,i1+4) =  324*D/(180*h*h) +  (-15*w)/(60*h);
    (*mpLhsMat).insert(i1,i1+3) = -670*D/(180*h*h) +  (+50*w)/(60*h);
    (*mpLhsMat).insert(i1,i1+2) =  855*D/(180*h*h) + (-100*w)/(60*h);
    (*mpLhsMat).insert(i1,i1+1) = -486*D/(180*h*h) + (+150*w)/(60*h);
    (*mpLhsMat).insert(i1,i1)   =  -70*D/(180*h*h) +  (-77*w)/(60*h) + k;
    (*mpLhsMat).insert(i1,i1-1) =  126*D/(180*h*h) +  (-10*w)/(60*h);
    
    int i2 = 2;
    (*mpLhsMat).insert(i2,i2+5) =  -2*D/(180*h*h);
    (*mpLhsMat).insert(i2,i2+4) =  +16*D/(180*h*h) +    (-w)/(60*h);
    (*mpLhsMat).insert(i2,i2+3) =  -54*D/(180*h*h) +  (+8*w)/(60*h);
    (*mpLhsMat).insert(i2,i2+2) =  +85*D/(180*h*h) + (-30*w)/(60*h);
    (*mpLhsMat).insert(i2,i2+1) = +130*D/(180*h*h) + (+80*w)/(60*h);
    (*mpLhsMat).insert(i2,i2)   = -378*D/(180*h*h) + (-35*w)/(60*h) + k;
    (*mpLhsMat).insert(i2,i2-1) = +214*D/(180*h*h) + (-24*w)/(60*h);
    (*mpLhsMat).insert(i2,i2-2) =  -11*D/(180*h*h) +  (+2*w)/(60*h);


    int im2 = mNumNodes-3;
    (*mpLhsMat).insert(im2,im2-5) =  -2*D/(180*h*h);
    (*mpLhsMat).insert(im2,im2-4) =  +16*D/(180*h*h) +    (+w)/(60*h);
    (*mpLhsMat).insert(im2,im2-3) =  -54*D/(180*h*h) +  (-8*w)/(60*h);
    (*mpLhsMat).insert(im2,im2-2) =  +85*D/(180*h*h) + (+30*w)/(60*h);
    (*mpLhsMat).insert(im2,im2-1) = +130*D/(180*h*h) + (-80*w)/(60*h);
    (*mpLhsMat).insert(im2,im2)   = -378*D/(180*h*h) + (+35*w)/(60*h) + k;
    (*mpLhsMat).insert(im2,im2+1) = +214*D/(180*h*h) + (+24*w)/(60*h);
    (*mpLhsMat).insert(im2,im2+2) =  -11*D/(180*h*h) +  (-2*w)/(60*h);
    
    int im1 = mNumNodes-2;
    (*mpLhsMat).insert(im1,im1-6) =   11*D/(180*h*h);
    (*mpLhsMat).insert(im1,im1-5) =  -90*D/(180*h*h) +   (-2*w)/(60*h);
    (*mpLhsMat).insert(im1,im1-4) =  324*D/(180*h*h) +  (+15*w)/(60*h);
    (*mpLhsMat).insert(im1,im1-3) = -670*D/(180*h*h) +  (-50*w)/(60*h);
    (*mpLhsMat).insert(im1,im1-2) =  855*D/(180*h*h) +  (100*w)/(60*h);
    (*mpLhsMat).insert(im1,im1-1) = -486*D/(180*h*h) + (-150*w)/(60*h);
    (*mpLhsMat).insert(im1,im1)   =  -70*D/(180*h*h) +  (+77*w)/(60*h) + k;
    (*mpLhsMat).insert(im1,im1+1) =  126*D/(180*h*h) +  (+10*w)/(60*h);

    for (int i=3; i<mNumNodes-3; i++) {
        double xm = mpGrid->mNodes[i-1].C.x; 
        double x = mpGrid->mNodes[i].C.x;
        double xp = mpGrid->mNodes[i+1].C.x;
        
        // Checking that grid is uniform
        assert(h - (x-xm) < 1e-5);
        assert(h - (xp-x) < 1e-5);

        // Constructing 6th order finite difference scheme:
        /*
         diffusion_alpha_m3 = alpha(-3) 
         diffusion_alpha_m2 = alpha(-2) 
         diffusion_alpha_m1 = alpha(-1) 
         diffusion_alpha_0  = alpha(0)
         diffusion_alpha_p1 = alpha(+1) 
         diffusion_alpha_p2 = alpha(+2) 
         diffusion_alpha_p3 = alpha(+3) 
         */
        double diffusion_alpha_m3 =   2*D/(180*h*h) +   (-w)/(60*h);
        double diffusion_alpha_m2 = -27*D/(180*h*h) + (+9*w)/(60*h);
        double diffusion_alpha_m1 = 270*D/(180*h*h) +(-45*w)/(60*h);
        double diffusion_alpha_0  =-490*D/(180*h*h) + k;
        double diffusion_alpha_p1 = 270*D/(180*h*h) +(+45*w)/(60*h);
        double diffusion_alpha_p2 = -27*D/(180*h*h) + (-9*w)/(60*h);
        double diffusion_alpha_p3 =   2*D/(180*h*h) +   (+w)/(60*h);

        (*mpLhsMat).insert(i,i-3) = diffusion_alpha_m3;
        (*mpLhsMat).insert(i,i-2) = diffusion_alpha_m2;
        (*mpLhsMat).insert(i,i-1) = diffusion_alpha_m1;
        (*mpLhsMat).insert(i,i)   = diffusion_alpha_0;
        (*mpLhsMat).insert(i,i+1) = diffusion_alpha_p1;
        (*mpLhsMat).insert(i,i+2) = diffusion_alpha_p2;
        (*mpLhsMat).insert(i,i+3) = diffusion_alpha_p3;
    }
}

void BvpOde::PopulateVector() {
    for (int i=1; i<mNumNodes-1; i++) {
        double x = mpGrid->mNodes[i].C.x;
        (*mpRhsVec)(i) = mpOde->mpRhsFunc(x);
    }
}

void BvpOde::ApplyBoundaryConditions() {
    bool left_bc_applied = false; 
    bool right_bc_applied = false;
    double D = mpOde->mCoeffOfUxx;
    double w = mpOde->mCoeffOfUx;
    double k = mpOde->mCoeffOfU;        
    if (mpBconds->mX0BcIsConst) {
        (*mpLhsMat).insert(0,0) = 1.0;
        (*mpRhsVec)(0) = mpBconds->mX0BcValue; 
        left_bc_applied = true;
    }

    if (mpBconds->mXNBcIsConst) {
        (*mpLhsMat).insert(mNumNodes-1,mNumNodes-1) = 1.0;
        (*mpRhsVec)(mNumNodes-1) = mpBconds->mXNBcValue; 
        right_bc_applied = true;
    }

    if (mpBconds->mX0BcIsNoFlux) {
        assert(left_bc_applied == false);
        double h = mpGrid->mNodes[1].C.x - mpGrid->mNodes[0].C.x;
        (*mpLhsMat).insert(0,0) = -2.0*D/h; 
        (*mpLhsMat).insert(0,1) = 2.0*D/h;
        (*mpRhsVec)(0) = 0;//No Flux =0 
        left_bc_applied = true;
    }

    if (mpBconds->mXNBcIsNoFlux) {
        assert(right_bc_applied == false);
        double h = mpGrid->mNodes[mNumNodes-1].C.x - mpGrid->mNodes[mNumNodes-2].C.x; 
        (*mpLhsMat).insert(mNumNodes-1,mNumNodes-1) = -2.0*D/h; 
        (*mpLhsMat).insert(mNumNodes-1,mNumNodes-2) = 2.0*D/h; 
        (*mpRhsVec)(mNumNodes-1) = 0; //No flux =0 
        right_bc_applied = true;
    }

    if (mpBconds->mX0BcIsFlux) {
        assert(left_bc_applied == false);
        double F = mpBconds->mX0BcValue;
        double h = mpGrid->mNodes[1].C.x - mpGrid->mNodes[0].C.x;
        (*mpLhsMat).insert(0,0) = - 2.0*D/(h*h) - 2.0*w/h - w*w/D + k; 
        (*mpLhsMat).insert(0,1) = 2.0*D/(h*h);
        (*mpRhsVec)(0) = F * (2.0/h + w/D); 
        left_bc_applied = true;
    }

    if (mpBconds->mXNBcIsFlux) {
        assert(right_bc_applied == false);
        double F = mpBconds->mXNBcValue;
        double h = mpGrid->mNodes[mNumNodes-1].C.x - mpGrid->mNodes[mNumNodes-2].C.x; 
        (*mpLhsMat).insert(mNumNodes-1,mNumNodes-2) = 2*D/(h*h);
        (*mpLhsMat).insert(mNumNodes-1,mNumNodes-1) = -2*D/(h*h) + 2*w/h + w*w/D + k; 
        (*mpRhsVec)(mNumNodes-1) = F * (2.0/h + w/D); 
        cout <<endl<<endl<<endl<<(*mpRhsVec)(mNumNodes-1)<<endl<<endl<<endl;
        right_bc_applied = true;
    }
    assert(right_bc_applied & left_bc_applied);
}

void BvpOde::WriteSolutionFile() {
    std::ofstream output_file(mFilename.c_str()); 
    assert(output_file.is_open());
    for (int i=0; i<mNumNodes; i++) {
        double x = mpGrid->mNodes[i].C.x;
          output_file << x << "  " << solution(i) << "\n";
    }
   output_file.flush();
   output_file.close();
   std::cout<<"Solution written to "<<mFilename<<"\n";
}

void BvpOde::SetFilename(const std::string& name){
    mFilename = name;
}

