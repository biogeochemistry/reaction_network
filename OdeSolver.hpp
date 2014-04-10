#ifndef ODESOLVERDEF
#define ODESOLVERDEF

#include <iostream>
#include <iosfwd>
#include <fstream>

class OdeSolver {
public:
  double current_time;
  double current_value;  
  double stepSize;
  double initialTime;
  double finalTime;
  double initialValue;
  void SetStepSize (double h);
  void SetTimeInterval (double t0, double t1);
  void SetInitialValue (double y0);
  double RightHandSide (double t, double y);
  void WritingResultsIntoFile(std::ofstream &file, double t, double y);
  virtual double SolveEquation() = 0;
};

#endif // ODESOLVERDEF