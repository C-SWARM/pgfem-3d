/**
 * Authors:
 *  Matt Mosby, [1], <mmosby1@nd.edu>
 *  Sangmin Lee, [1], <slee43@nd.edu>
 *  Alberto Salvadori, [1], <asalvad2@nd.edu>
 *  [1] - University of Notre Dame, Notre Dame, IN
 */

#include "data_structure.h"
#include "state_variables.h"
#include "utils.h"


int State_variables::initialization(const size_t n_Fs,
                                    const size_t n_vars,
                                    const size_t n_flags)
{
  int err = 0;
  this->n_Fs = n_Fs;
  this->Fs = new Matrix<double>[n_Fs];

  for(size_t i = 0; i < n_Fs; i++)
  {
    this->Fs[i].initialization(3,3);
    this->Fs[i](1,1) = this->Fs[i](2,2) = this->Fs[i](3,3) = 1.0;
  }

  /* state variables is column vector */
  this->state_vars = new Matrix<double>[1];
  this->state_vars[0].initialization(n_vars,1);

  this->n_flags = n_flags;
  this->flags = new int[n_flags]();

  return err;
}

size_t State_variables::state_variables_get_packed_size(void)
{
  /* get the length of all Fs. Note that we are computing the length
     of each tensor just in case. */
  size_t len_Fs = 0;
  for (int i = 0, e = this->n_Fs; i < e; i++) {
    len_Fs += this->Fs[i].m_col * this->Fs[i].m_row;
  }
  const size_t len_vars = this->state_vars->m_col * this->state_vars->m_row;
  return ( (len_Fs + len_vars) * sizeof(double) + this->n_flags * sizeof(int) );
}

int State_variables::state_variables_pack(char *buffer,
                                          size_t *pos)
{
  int err = 0;

  /* pack the deformation gradients */
  for (size_t i = 0, e = this->n_Fs; i < e; i++) {
    const size_t len = this->Fs[i].m_row * this->Fs[i].m_col;
    const double *data = this->Fs[i].m_pdata;
    pack_data(data, buffer, pos, len, sizeof(*data));
  }

  /* pack the state variables */
  {
    const size_t len = this->state_vars->m_row * this->state_vars->m_col;
    const double *data = this->state_vars->m_pdata;
    pack_data(data, buffer, pos, len, sizeof(*data));
  }

  /* pack the flags */
  pack_data(this->flags, buffer, pos, this->n_flags, sizeof(*(this->flags)));

  return err;
}

int State_variables::state_variables_unpack(const char *buffer,
                                            size_t *pos)
{
  int err = 0;

  /* pack the deformation gradients */
  for (size_t i = 0, e = this->n_Fs; i < e; i++) {
    const size_t len = this->Fs[i].m_row * this->Fs[i].m_col;
    double *data = this->Fs[i].m_pdata;
    unpack_data(buffer, data, pos, len, sizeof(*data));
  }

  /* pack the state variables */
  {
    const size_t len = this->state_vars->m_row * this->state_vars->m_col;
    double *data = this->state_vars->m_pdata;
    unpack_data(buffer, data, pos, len, sizeof(*data));
  }

  /* pack the flags */
  unpack_data(buffer, this->flags, pos, this->n_flags, sizeof(*(this->flags)));

  return err;
}

/// compute mid-point value of state variables
///
/// \param[in] id_nm1 id of state variable at t(n-1)
/// \param[in] id_n   id of state variable at t(n)
/// \param[in] id_np1 id of state variable at t(n+1)
/// \param[in] npa    if npa = 0: returns (1-alpha)*value(id_nm1) + alpha*value(id_n)
///                      npa = 1: returns (1-alpha)*value(id_n)   + alpha*value(id_np1)
/// \param[in] alpha  mid-point alpha
/// \return    mid-point value of the state variable
double State_variables::compute_state_vars_npa(const int id_nm1,
                                               const int id_n,
                                               const int id_np1,
                                               const int npa,
                                               const double alpha){
  double *vars = this->state_vars[0].m_pdata;
  double value_npa = vars[id_np1];
  switch(npa)
  {
    case 0:
      mid_point_rule(&value_npa, vars + id_nm1, vars + id_n, alpha, 1);
      break;
    case 1:  
      mid_point_rule(&value_npa, vars + id_n, vars + id_np1, alpha, 1);
      break;
  }
  return value_npa;
}

int
Model_var_info::print_variable_info(FILE *f)
{
  int err = 0;
  fprintf(f,"F names: ");
  for(int i = 0, e = this->n_Fs; i < e; i++) fprintf(f,"%s ",this->F_names[i]);
  fprintf(f,"\nVar names: ");
  for(int i = 0, e = this->n_vars; i < e; i++) fprintf(f,"%s ",this->var_names[i]);
  fprintf(f,"\nFlag names: ");
  for(int i = 0, e = this->n_flags; i < e; i++) fprintf(f,"%s ",this->flag_names[i]);
  fprintf(f,"\n");
  return err;
}

Model_var_info::~Model_var_info()
{
  /// deallocate internal memory
  for (size_t ia=0; ia<this->n_Fs; ia++) {
    if (this->F_names[ia]) free(this->F_names[ia]);
  }
  if(this->F_names) free(this->F_names);

  for(size_t ia=0; ia<this->n_vars; ia++) {
    if (this->var_names[ia]) free(this->var_names[ia]);
  }

  if(this->var_names) free(this->var_names);

  for(size_t ia=0; ia<this->n_flags; ia++) {
    if(this->flag_names[ia]) free(this->flag_names[ia]);
  }
  if (this->flag_names) free(this->flag_names);
}

/// Set constitutive model information such as variable names.
/// This function can be used in each constitutive model to set name for variables,
/// tensors, and flags 
/// \param[in] info    Object for storing state variables at an integration point
/// \param[in] Vars_no number of variabls
/// \param[in] Vars    container of variable names
/// \param[in] Tens_no number of tensors
/// \param[in] Tens    container of tensor names
/// \param[in] Flag_no number of flags
/// \param[in] Flag    container of flag names
/// \return non-zero on internal error
int constitutive_model_info(Model_var_info &info,
                            const int Vars_no, const CMVariableNames *Vars,
                            const int Tens_no, const CMVariableNames *Tens,
                            const int Flag_no, const CMVariableNames *Flag)
{
  int err = 0;
  int test_size = 1024;
  // make sure I don't leak memory
  
  // set variable names    
  info.n_vars = Vars_no;
  info.var_names = (char **)malloc(sizeof(char*)*Vars_no);
  if(info.var_names == NULL){    
    ++err;
    return err;
  }
    
  for(int a=0; a<Vars_no; a++){
    info.var_names[a] = (char *)malloc(sizeof(char)*test_size);
    if(info.var_names[a] == NULL)
      ++err;
  }

  if(err>0)
    return err;
  
  // set tensor names  
  for(int a=0; a<Vars_no; a++)
    sprintf(info.var_names[a], Vars[a].name);

  info.n_Fs = Tens_no;
  info.F_names = (char **)malloc(sizeof(char*)*Tens_no);
  if(info.F_names==NULL){
    ++err;
    return err;
  }
  
  for(int a=0; a<Tens_no; a++){
    info.F_names[a] = (char *)malloc(sizeof(char)*test_size);
    if(info.F_names[a] == NULL)
      ++err;
  }
  
  if(err>0)
    return err;
  
  for(int a=0; a<Tens_no; a++)
    sprintf(info.F_names[a], Tens[a].name);
    
    
  // set flag names  
  info.n_flags = Flag_no;
  info.flag_names = (char **)malloc(sizeof(char*)*Flag_no);
  if(info.flag_names == NULL){
    ++err;
    return err;
  }
  
  for(int a=0; a<Flag_no; a++){
    info.flag_names[a] = (char *)malloc(sizeof(char)*test_size);
    if(info.flag_names[a]==NULL)
      ++err;
  }

  if(err>0)
    return err;
    
  for(int a=0; a<Flag_no; a++)
    sprintf(info.flag_names[a], Flag[a].name);

  return err;  
} 