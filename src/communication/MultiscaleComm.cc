#include "pgfem3d/Communication.hpp"
#include "pgfem3d/MultiscaleCommunication.hpp"
#include "PGFEM_io.h"

using namespace pgfem3d::net;

namespace pgfem3d {

  MultiscaleComm::MultiscaleComm(const PGFem3D_Comm comm_world,
				 Network *n)
  {
    int myrank = 0;
    
    /* duplicate communicator */
    n->comm_dup(comm_world, &world);
    macro = world;
    micro_all = world;
    micro = world;
    mm_inter = world;
    
    /* store rank */
    n->comm_rank(comm_world, &myrank);
    rank_world = myrank;
    rank_macro = myrank;
    rank_micro_all = myrank;
    rank_micro = myrank;
    rank_mm_inter = myrank;
    
    /* set valid flag */
    valid_macro = 1;
    valid_micro_all = 1;
    valid_micro = 1;
    valid_mm_inter = 1;

    /* save network handle */
    net = n;
  }
  
  MultiscaleComm::~MultiscaleComm()
  {
    net->comm_free(&world);
    
    if(valid_macro) {
      net->comm_free(&macro);
    }
    
    if(valid_micro_all) {
      net->comm_free(&micro_all);
    }
    
    if(valid_micro) {
      net->comm_free(&micro);
      net->comm_free(&worker_inter);
    }
    
    if(valid_mm_inter){
      net->comm_free(&mm_inter);
    }
  }
  
  void MultiscaleComm::MM_split(const int macro_nproc,
				const int micro_group_size)
  {
    int nproc_world = 0;
    net->comm_size(world, &nproc_world);
    
    /* error check group sizes */
    const int micro_nproc_all = nproc_world - macro_nproc;
    if (micro_nproc_all < micro_group_size
	|| (micro_nproc_all%micro_group_size != 0)){
      if (rank_world == 0){
	PGFEM_printerr("ERROR: incorrect group sizes in %s!\n",__func__);
	PGFEM_printerr("WORLD = %d || MACRO = %d || GROUP = %d\n",
		       nproc_world,macro_nproc,micro_group_size);
	PGFEM_Abort();
      }
    }
    
    /* split into macro and micro communicators */
    {
      const int macro_color = (rank_world < macro_nproc) ? 1 : NET_UNDEFINED;
      net->comm_split(world,macro_color, rank_world, &macro);
      
      const int micro_color = (rank_world < macro_nproc) ? NET_UNDEFINED : 1;
      net->comm_split(world, micro_color, rank_world, &micro_all);
    }
    
    /* get ranks on new communicators */
    if (macro == NET_COMM_NULL) {
      valid_macro = 0;
      rank_macro = NET_UNDEFINED;
    } else {
      net->comm_rank(macro, &rank_macro);
    }
    
    if (micro_all == NET_COMM_NULL){
      valid_micro_all = 0;
      rank_micro_all = NET_UNDEFINED;
    } else {
      net->comm_rank(micro_all, &rank_micro_all);
    }
    
    /* split micro communicator into work groups */
    if (valid_micro_all && micro_group_size > 0) {
      int color = NET_UNDEFINED;
      /* integer division for group id */
      if (rank_micro_all != NET_UNDEFINED){
	color = rank_micro_all / micro_group_size;
      }
      net->comm_split(micro_all, color, rank_micro_all, &micro);
      net->comm_rank(micro, &rank_micro);
      
      /*Create inter-micro communicators between equivalent procs on
	different microstructures. This allows direct communication
	between workers using the rank as the server id. Split
	micro_all communicator by rank in micro using rank_micro as the
	color and rank_micro_all as the key. */
      net->comm_split(micro_all, rank_micro, rank_micro_all, &worker_inter);
      net->comm_rank(worker_inter, &server_id);
    } else {
      micro = NET_COMM_NULL;
      valid_micro = 0;
      rank_micro = NET_UNDEFINED;
    }
    
    /* create the micro-macro intercommunicator */
    {
      int color = NET_UNDEFINED;
      if(rank_macro != NET_UNDEFINED
	 || rank_micro == 0){
	color = 1;
      }
      net->comm_split(world, color, rank_world, &mm_inter);
    }
    
    /* get rank on new communicator */
    if (mm_inter == NET_COMM_NULL) {
      valid_mm_inter = 0;
      rank_mm_inter = NET_UNDEFINED;
    } else {
      net->comm_rank(mm_inter, &rank_mm_inter);
    }
  }
  
} // end namespace pgfem3d