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
    mm_inter_ROM = world;
    micro_1 = world;
    micro_ROM = world;   
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
    valid_mm_inter_ROM = 1;
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

    if(valid_mm_inter_ROM){
      net->comm_free(&mm_inter_ROM);
    }

  }
  
  void MultiscaleComm::MM_split(const int macro_nproc,
				const int micro_group_size,
        const int micro_ROM_group_size,
        const int full_micro_np)
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
    //NET is typically MPI 
    /* split into macro and micro communicators */
    {
      const int macro_color = (rank_world < macro_nproc) ? 1 : NET_UNDEFINED;
      net->comm_split(world,macro_color, rank_world, &macro);
      
      const int micro_color = (rank_world < macro_nproc) ? NET_UNDEFINED : 1;
      net->comm_split(world, micro_color, rank_world, &micro_all);
    }

    // split micro into micro 1 and 2
    {
      const int micro1_color = (rank_world < macro_nproc + full_micro_np) ? 1 : 2;
      net->comm_split(micro_all,micro1_color, micro_all, &micro_1);

      const int micro_ROM_color = (rank_world < macro_nproc + full_micro_np) ? 2 : 1;
      net->comm_split(micro_all,micro_ROM_color, micro_all, &micro_ROM);
      
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
   
    if (micro_1 == NET_COMM_NULL) {
      valid_micro_1 = 0;
      rank_micro_1 = NET_UNDEFINED;
    } else {
      net->comm_rank(micro_1, &rank_micro_1);
    }

    if (micro_ROM == NET_COMM_NULL) {
      valid_micro_ROM = 0;
      rank_micro_ROM = NET_UNDEFINED;
    } else {
      net->comm_rank(micro_ROM, &rank_micro_ROM);
    }


 
    /* split micro 1 communicator into work groups */
    if (valid_micro_1 && micro_group_size > 0) {
      int color = NET_UNDEFINED;
      /* integer division for group id */
      if (rank_micro_1 != NET_UNDEFINED){
        color = rank_micro_1 / micro_group_size;
      }
      net->comm_split(micro_1, color, rank_micro_1, &micro);
      net->comm_rank(micro, &rank_micro_1);

      /*Create inter-micro communicators between equivalent procs on
	different microstructures. This allows direct communication
	between workers using the rank as the server id. Split
	micro_all communicator by rank in micro using rank_micro as the
	color and rank_micro_all as the key. */
      net->comm_split(micro_1, rank_micro_1, rank_micro_all, &worker_inter);
      net->comm_rank(worker_inter, &server_id);
    } else {
      micro = NET_COMM_NULL;
      valid_micro = 0;
      rank_micro = NET_UNDEFINED;
    }


    /* split micro 2 communicator into work groups */
    if (valid_micro_ROM && micro_ROM_group_size > 0) {
      int color = NET_UNDEFINED;
      /* integer division for group id */
      if (rank_micro_ROM != NET_UNDEFINED){
        color = rank_micro_ROM / micro_ROM_group_size;
      }
      net->comm_split(micro_ROM, color, rank_micro_ROM, &micro_ROM);
      net->comm_rank(micro_ROM, &rank_micro_ROM);

      /*Create inter-micro communicators between equivalent procs on
    different microstructures. This allows direct communication
    between workers using the rank as the server id. Split
    micro_all communicator by rank in micro using rank_micro as the
    color and rank_micro_all as the key. */
      net->comm_split(micro_ROM, rank_micro_ROM, rank_micro_all, &worker_inter_ROM);
      net->comm_rank(worker_inter_ROM, &server_id);
    } else {
      micro = NET_COMM_NULL;
      valid_micro = 0;
      rank_micro = NET_UNDEFINED;
    }

    
    /* create the micro-macro intercommunicator */
    {
      int color = NET_UNDEFINED;
      if(rank_macro != NET_UNDEFINED || rank_micro_1 == 0){
      	color = 1;
      }
      net->comm_split(world, color, rank_world, &mm_inter);
    }
    
    /* create the micro2-macro intercommunicator */
    {
      int color = NET_UNDEFINED;
      if(rank_macro != NET_UNDEFINED || rank_micro_ROM == 0){
        color = 1;
      }
      net->comm_split(world, color, rank_world, &mm_inter_ROM);
    }

    /* get rank on new communicator */
    if (mm_inter == NET_COMM_NULL) {
      valid_mm_inter = 0;
      rank_mm_inter = NET_UNDEFINED;
    } else {
      net->comm_rank(mm_inter, &rank_mm_inter);
    }
    /* get rank on new communicator */
    if (mm_inter_ROM == NET_COMM_NULL) {
      valid_mm_inter_ROM = 0;
      rank_mm_inter_ROM = NET_UNDEFINED;
    } else {
      net->comm_rank(mm_inter_ROM, &rank_mm_inter_ROM);
    }

  }
  
} // end namespace pgfem3d
