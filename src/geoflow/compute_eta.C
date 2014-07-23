/*******************************************************************
 * Copyright (C) 2003 University at Buffalo
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Author: 
 * Description: 
 *
 *******************************************************************
 * $Id: step.C 164 2014-07-10 12:11:22 haghakha $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "../header/hpfem.h"

void calc_map_area(HashTable* El_Table, StatProps* statprops);

double compute_eta(HashTable* El_Table,StatProps* statprops){

  calc_map_area( El_Table, statprops);

  double local_eta=0.,eta=0. ,*dx, phi=0.0;
  HashEntryPtr* buck = El_Table->getbucketptr();

  for(int i=0; i<El_Table->get_no_of_buckets(); i++)
    if(*(buck+i))
    {
      HashEntryPtr currentPtr = *(buck+i);
      while(currentPtr)
      {
	Element* Curr_El=(Element*)(currentPtr->value);
	if(Curr_El->get_adapted_flag()>0) {

	  dx=Curr_El->get_dx();
	  phi=*(Curr_El->get_state_vars());

	  local_eta+=4*dx[0]*dx[1]*phi*(phi*phi-1);
	}
	currentPtr=currentPtr->next;
      }
    }
  local_eta/=statprops->map_area;
  MPI_Allreduce(&local_eta,&eta,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);


  return eta;
}

void calc_map_area(HashTable* El_Table, StatProps* statprops){

  double map_area=0.0,*dx;

  HashEntryPtr* buck = El_Table->getbucketptr();

  for(int i=0; i<El_Table->get_no_of_buckets(); i++)
    if(*(buck+i))
    {             
      HashEntryPtr currentPtr = *(buck+i);
      while(currentPtr)             
      {                                         
	Element* Curr_El=(Element*)(currentPtr->value);       
	if(Curr_El->get_adapted_flag()>0) { 
	  dx=Curr_El->get_dx();

	  map_area+=dx[0]*dx[1];

	}
	currentPtr=currentPtr->next;
      }
    }

  MPI_Allreduce(&map_area,&(statprops->map_area),1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  return;
}


