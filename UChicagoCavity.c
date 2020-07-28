#include <math.h>
#include <stdlib.h>
#include "quac.h"
#include "operators.h"
#include "error_correction.h"
#include "solver.h"
#include "dm_utilities.h"
#include "quantum_gates.h"
#include "quantum_circuits.h"
#include "petsc.h"
#include "qasm_parser.h"
#include "qsystem.h"

PetscErrorCode ts_monitor(TS,PetscInt,PetscReal,Vec,void*);

qvec wf_in;
qvec rho;

PetscErrorCode ts_monitor(TS ts,PetscInt step,PetscReal time, Vec rho1, void *ctx){
  //Print out things at each time step, if desired
  //PetscPrintf(PETSC_COMM_WORLD,"Step: %d, time: %f\n",step,time);
  PetscReal fid;
  get_superfidelity_qvec(wf_in,rho,&fid);
  printf("Fidelity: %f\n",fid);
  /* print_qvec(rho); */
  PetscFunctionReturn(0);
}

int main(int argc,char **args){
  
  circuit circ;
  qsystem system;
  operator *operators;
  PetscReal gamma_amp,gamma_dep,dt,time_max,fidelity;
  PetscInt num_modes,num_stpes,steps_max,i;
  PetscScalar mat_val;
  PetscReal chi[] = {-7539822.3686155 , -5466371.21724624, -4209734.15581032,
       -3581415.62509236, -2953097.09437441, -2324778.56365645,
       -2261946.71058465, -2010619.29829747};
  PetscReal kappa[] = {56548.66776462, 56674.33147076, 32861.05915655, 26326.54643708,
       20357.52039526,  8796.45943005,  9299.11425463,  4523.89342117};
  PetscReal gamma_c1[] = {495.04950495, 467.28971963, 487.80487805, 531.91489362,
       558.65921788, 543.47826087, 497.51243781, 515.46391753};
  PetscReal gamma_c2[] = {346.02076125, 343.64261168, 396.82539683, 364.96350365,
       436.68122271, 392.15686275, 420.16806723, 396.82539683};

  PetscReal alpha = -892212313.6195012;
  char wf_in_string[PETSC_MAX_PATH_LEN]="wf_in.dat", dm_out_string[PETSC_MAX_PATH_LEN]="dm_out.dat";
  /* Initialize QuaC */
  QuaC_initialize(argc,args);

  PetscOptionsGetInt(NULL,NULL,"-num_modes",&num_modes,NULL);
  PetscOptionsGetString(NULL,NULL,"-wf_in",wf_in_string,sizeof(wf_in_string),NULL);
  PetscOptionsGetString(NULL,NULL,"-dm_out",dm_out_string,sizeof(dm_out_string),NULL);

  operators = malloc((num_modes+1)*sizeof(struct operator));

  gamma_amp = 12658.227848101267;
  gamma_dep = 17241.379310344826;

  initialize_system(&system);

  create_op_sys(system,2,&operators[0]);

  for(i=1;i<=num_modes;i++){
      create_op_sys(system,7,&operators[i]);
  }

  add_ham_term(system,alpha,4,operators[0]->dag,operators[0]->dag,operators[0],operators[0]);

  for(i=1;i<=num_modes;i++){
      add_ham_term(system,kappa[(i-1)],4,operators[i]->dag,operators[i]->dag,operators[i],operators[i]);
      add_ham_term(system,chi[(i-1)],4,operators[0]->dag,operators[0],operators[i]->dag,operators[i]);
  }

  add_lin_term(system,gamma_amp,1,operators[0]);
  add_lin_term(system,gamma_dep,1,operators[0]->n);

  for(i=1;i<=num_modes;i++){
      add_lin_term(system,gamma_c1[(i-1)],1,operators[1]);
      add_lin_term(system,gamma_c2[(i-1)],1,operators[i]->n);
  }

  //Construct the matrix now that we are done adding to it
  construct_matrix(system);

  //Read in the input wavefunction
  read_qvec_wf_binary(&(wf_in),wf_in_string);


  //Time step until time max=1
  time_max  = 0.0001;
  dt        = 0.00000000001;
  steps_max = 100000;

  /* Set the ts_monitor to print results at each time step, if desired */
  set_ts_monitor_sys(system,ts_monitor,NULL);

  create_qvec_sys(system,&rho);
  //Copy the wf into the DM
  copy_qvec(wf_in,rho);
  //print_sparse_mat(system->mat_A);
  //Run the evolution, with the error channels above
  time_step_sys(system,rho,0.0,time_max,dt,steps_max);

  //Get superfidelity - bound for fidelity, but equals fidelity if one state is pure
  //get_superfidelity_qvec(wf_in,rho,&fidelity);
  //printf("Fidelity: %f\n",fidelity);
  //print_qvec(rho);

 destroy_op(&operators[0]);
 for (i=1;i<num_modes;i++){
    destroy_op(&operators[i]);
  }
  destroy_qvec(&rho);
  destroy_qvec(&wf_in);
  //free(chi);
  //free(kappa);
  //free(gamma_c1);
  //free(gamma_c2);
  QuaC_finalize();
  return 0;

}




