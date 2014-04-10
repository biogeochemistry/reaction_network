#include "Runge-Kutta.hpp"


double RungeKutta::SolveEquation() {
  std::ofstream results;
  results.open("results.txt");
  results << "time\t|\tvalue\n";
  current_value = initialValue; 
  double k1, k2, k3, k4;
  for (current_time = initialTime; current_time < finalTime; current_time += stepSize)
  {
    OdeSolver::WritingResultsIntoFile(results, current_time, current_value);
    k1 = stepSize * OdeSolver::RightHandSide(current_time, current_value);
    k2 = stepSize * OdeSolver::RightHandSide(current_time + 0.5 * stepSize, current_value + 0.5 * k1);
    k3 = stepSize * OdeSolver::RightHandSide(current_time + 0.5 * stepSize, current_value + 0.5 * k2);
    k4 = stepSize * OdeSolver::RightHandSide(current_time + stepSize, current_value + k3);
    current_value += (k1 + 2.f * k2 + 2.f * k3 + k4) / 6.f; 
   
  }
  results.close();
  return 0;
}