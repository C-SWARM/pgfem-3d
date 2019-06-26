/**
 * Unit test for issue #43 add integer flag array to state variables
 * AUTHORS:
 *  Matt Mosby, UNiversity of Notre Dame, <mmosby1@nd.edu>
 */
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "allocation.h"
#include "constitutive_model.h"
#include <cassert>

int main(int argc, char **argv)
{
  int err = 0;
  int test_all = 1;
  int test_single = -1;
  if (argc > 1) {
    test_all = 0;
    test_single = atoi(argv[1]);
  }

  HOMMAT mat;
  mat.nu = 0.3;
  mat.E = 1.0;
  mat.devPotFlag = 1;
  mat.volPotFlag = 2;

  if (test_all) {
    for (int i = 0; i < NUM_MODELS; i++) {
      Model_parameters *p = new Model_parameters[1];
      err += construct_Model_parameters(&p,0,i);
      err += p->initialization(&mat, i);
      assert(err == 0);

      /* Create and print the model info */
      Model_var_info info;
      err += p->get_var_info(info);
      err += info.print_variable_info(stdout);
      delete [] p;
      printf("\n");
    }
  } else {
      Model_parameters *p = new Model_parameters[1];
      err += construct_Model_parameters(&p,0,test_single);
      err += p->initialization(&mat, test_single);
      assert(err == 0);

      /* Create and print the model info */
      Model_var_info info;
      err += p->get_var_info(info);
      err += info.print_variable_info(stdout);
      delete [] p;
      printf("\n");
  }

  return err;
}
