
#include "dw_MSStateSpace.h"

#include "dsgelinv1_comfuns.h"

T_MSStateSpace* Create_MSStateSpaceFromStateModel(TStateModel *TZstatemodel, int nlags_encoded)
{
  int i, t;
  T_MSStateSpace *statespace;
  TMarkovStateVariable *sv=TZstatemodel->sv;
  TSdsgemodel *TZdsgemodel=(TSdsgemodel*)(TZstatemodel->theta);

  // Create T_MSStateSpace structure
  //statespace=Create_MSStateSpace(TZdsgemodel->ny,TZdsgemodel->nz,0,TZdsgemodel->nshocks,
  statespace=Create_MSStateSpace(TZdsgemodel->ny,TZdsgemodel->nz,0,TZdsgemodel->nshocks_kalman,
                    sv->nbasestates,nlags_encoded,TZdsgemodel->nSample);

  //printf("basestates: %d  lags encoded: %d  new lags encoded: %d\n",sv->nbasestates,sv->nlags_encoded,nlags_encoded);
  //getchar();

  for (t=statespace->nobs; t > 0; t--)
    for (i=statespace->ny-1; i >= 0; i--)
      ElementV(statespace->y[t],i)=ElementM(TZdsgemodel->SampleData_dm,t-1,i);

  statespace->TZstatemodel=TZstatemodel;

  return statespace;
}

/*
   Initializes MSStateSpace structure.  Returns 1.
*/
int Initialize_MSStateSpace(TStateModel *model)
{
  T_MSStateSpace *statespace=(T_MSStateSpace*)(model->theta);
  TSkalfilmsinputs_1stapp *TZKalman=((TSdsgemodel*)(statespace->TZstatemodel->theta))->kalfilmsinputs_1stapp_ps;
  int i, j;

  // Update Tao's structures
  RefreshEverything(statespace->TZstatemodel);

  // Update T_MSStateSpace structure
  // Parameters
  for (i=statespace->nbasestates-1; i >= 0; i--)
    {
      ColumnVector(statespace->a[i],TZKalman->at_dm,i);
      EquateMatrix(statespace->H[i],TZKalman->Ht_dc->C[i]);
      ColumnVector(statespace->b[i],TZKalman->bt_dm,i);
      EquateMatrix(statespace->F[i],TZKalman->Ft_dc->C[i]);
      EquateMatrix(statespace->R[i],TZKalman->Rt_dc->C[i]);
      EquateMatrix(statespace->V[i],TZKalman->Vt_dc->C[i]);
      EquateMatrix(statespace->G[i],TZKalman->Gt_dc->C[i]);
      if (statespace->nu > 0) EquateMatrix(statespace->Phiy[i],TZKalman->Psiut_dc->C[i]);
      if (statespace->nepsilon > 0) EquateMatrix(statespace->Phiz[i],TZKalman->Psiet_dc->C[i]);

      /* Check - Comment out except when debugging *
      TMatrix X;
      PRECISION norm;
      if (statespace->nu > 0) 
	{
	  X=ProductTransposeMM((TMatrix)NULL,statespace->Phiy[i],statespace->Phiy[i]);
	  norm=MatrixNormEuclidean(SubtractMM(X,X,statespace->R[i]));
	  FreeMatrix(X);
	}
      else
	norm=MatrixNormEuclidean(statespace->R[i]);
      if (norm > 1e-9)
	{
	  printf("R[%d] and Phiy[%d]Phiy[%d]' are not equal (%lg)- press enter to continue...\n",i,i,i,norm);
	  getchar();
	}
      if (statespace->nepsilon > 0) 
	{
	  X=ProductTransposeMM((TMatrix)NULL,statespace->Phiz[i],statespace->Phiz[i]);
	  norm=MatrixNormEuclidean(SubtractMM(X,X,statespace->V[i]));
	  FreeMatrix(X);
	}
      else
	norm=MatrixNormEuclidean(statespace->V[i]);
      if (norm > 1e-9)
	{
	  printf("V[%d] and Phiz[%d]Phiz[%d]' are not equal (%lg) - press enter to continue...\n",i,i,i,norm);
	  getchar();
	}
      /**/
    }

  for (i=statespace->nstates-1; i >= 0; i--)
    {
      j=model->sv->lag_index[i][0];
      ColumnVector(statespace->Ez1[1][i],TZKalman->z0_dm,j);
      EquateMatrix(statespace->Ezz1[1][i],TZKalman->P0_dc->C[j]);
    }

  // Flags
  statespace->NonZeroG=0;
  statespace->t0=0;

  return 1;
}
