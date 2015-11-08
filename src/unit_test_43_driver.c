/**
 * Unit test for issue #43 add integer flag array to state variables
 * AUTHORS:
 *  Matt Mosby, UNiversity of Notre Dame, <mmosby1@nd.edu>
 */

#include "constitutive_model.h"
#include <assert.h>

int main(int argc, char **argv)
{
  int err = 0;
  int test_all = 1;
  int test_single = -1;
  if (argc > 1) {
    test_all = 0;
    test_single = atoi(argv[1]);
  }

  Model_parameters *p = malloc(sizeof(*p));
  Model_var_info *info = NULL;

  if (test_all) {
    for (int i = 0; i < NUM_MODELS; i++) {
      err += model_parameters_construct(p);
      err += model_parameters_initialize(p, NULL, NULL, NULL, i);
      assert(err == 0);

      /* Create and print the model info */
      err += p->get_var_info(&info);
      err += model_var_info_print(stdout, info);
      printf("\n");

      /* cleanup */
      err += model_var_info_destroy(&info);
      err += model_parameters_destroy(p);
    }
  } else {
    err += model_parameters_construct(p);
    err += model_parameters_initialize(p, NULL, NULL, NULL, test_single);
    assert(err == 0);

    /* Create and print the model info */
    err += p->get_var_info(&info);
    err += model_var_info_print(stdout, info);
    printf("\n");

    /* cleanup */
    err += model_var_info_destroy(&info);
    err += model_parameters_destroy(p);
  }

  free(p);

  return err;
}
