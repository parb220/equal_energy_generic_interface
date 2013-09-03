#include<vector>
#include <gsl/gsl_poly.h>
#include <cmath>

using namespace std; 

bool Solve_Polynomial_Equation(vector<double> &t, size_t n, double t0, double tk_1)
{
	t.resize(n); 

	if (n < 2)
		return false; 
	if (n == 2)
		t[1]=t[0] = t0; 
	else 
	{
		double *coefficients = new double[n-1]; 
		coefficients[0] = t0-tk_1; 
		for (int i=1; i<n-1; i++)
			coefficients[i] = 1.0; 
		double *Z = new double[(n-2)*2]; 
		gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc(n-1); 
		gsl_poly_complex_solve(coefficients, n-1, w, Z);
		
		double gamma;
                bool continue_flag = true;
                for (int i=0; i<n-2 && continue_flag; i++)
                {
                        if (Z[2*i]>0 && abs(Z[2*i+1]) <= 1.0e-6)
                        {
                                gamma = Z[2*i];
                                continue_flag = false;
                        }
                }
		delete [] Z; 
		delete [] coefficients; 
		if (continue_flag)
			return false; 
		t[0] = t0; 
		t[n-2] = tk_1; 
		for (int i=1; i<n-2; i++)
			t[i] = t[i-1]+pow(gamma, (double)i); 
		t[n-1] = t[n-2]+pow(gamma, (double)(n-1)); 
	}
	return true; 
}
