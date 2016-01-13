/* 
* Copyright (C) 2008, Brian Tanner

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 */

#include <string.h> /*strcmp*/
#include <stdio.h> /*printf*/
#include <stdlib.h>
#include <boost/shared_ptr.hpp>

#include <rlglue/Environment_common.h>/* env_ function prototypes and RL-Glue types */
#include <rlglue/utils/C/RLStruct_util.h> /* helpful functions for allocating structs and cleaning them up */

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <lifev/core/LifeV.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include "Heart.hpp"

#include "RippleMacroParams.hpp"

observation_t this_observation;
reward_observation_terminal_t this_reward_observation;
boost::shared_ptr<LifeV::Heart> ptrHeart = NULL;
bool flg_mpi = false;

void observe(void)
{
   if( ptrHeart ){
      for(int i = 0; i < SIZE_OBS; i++){
         this_observation.doubleArray[i] = 0.0; //*(ptrHeart->M_Uptr)[i];
      }
   }
}


const char* env_init()
{    
   
   int argc = 3;
   char** argv = static_cast<char**>(calloc(sizeof(char*), argc));
   for(int i=0; i < argc; i++){
      argv[i] = static_cast<char*>(calloc(sizeof(char), 20 ));
   }
   sprintf(argv[0], "RippleEnvironment");
   sprintf(argv[1], "-f");
   sprintf(argv[2], "data");
  
   if(!flg_mpi){
      MPI_Init (&argc, &argv);
      Epetra_MpiComm Comm (MPI_COMM_WORLD);
      if ( Comm.MyPID() == 0 ){
          std::cout << "% using MPI" << std::endl;
      }
      flg_mpi = true;
   }

   GetPot command_line (argc, argv);
   const std::string data_file_name = command_line.follow ("data", 2, "-f", "--file");
   GetPot dataFile (data_file_name);
   ptrHeart.reset( new LifeV::Heart(dataFile) );

   for(int i = 0; i < argc; i++){
      free(argv[i]);
   }
   free(argv);

   if(SIZE_OBS != ptrHeart->M_Uptr->size()){
      std::cout << "Observation size mismatch!" << std::endl;
      std::cout << "SIZE_OBS:" << SIZE_OBS << std::endl;
      std::cout << "M_Uptr size:" << ptrHeart->M_Uptr->size() << std::endl;
      return NULL;
   }
	char* task_spec = "VERSION RL-Glue-3.0 PROBLEMTYPE episodic DISCOUNTFACTOR 1.0 OBSERVATIONS DOUBLES^SIZE_OBS ACTIONS INTS (0 1)  REWARDS (-1.0 1.0)  EXTRA RippleEnvironment(C/C++) by Naoki Tomii.";

	/* Allocate the observation variable */
	allocateRLStruct(&this_observation,0,SIZE_OBS,0);
	
   /* Setup the reward_observation variable */
	this_reward_observation.observation=&this_observation;
	this_reward_observation.reward=0;
	this_reward_observation.terminal=0;
   
   std::cout << task_spec << std::endl;
   return task_spec;
}

const observation_t *env_start()
{ 
   observe();
  	return &this_observation;
}

const reward_observation_terminal_t *env_step(const action_t *this_action)
{
  static int cnt = 0;
  int episode_over=0;
  double the_reward=0;

  std::cout << "env_step @" << ++cnt << std::endl;
  if(cnt>=2000){
    episode_over=1;
    the_reward=-1;
  }

  // Make sure the action is valid 
  try{
    assert(this_action->numInts==1);
    if(ptrHeart){
      if(this_action->intArray[0] >= 0){
        std::cout << "Action : " << this_action->intArray[0] << std::endl;
        ptrHeart->stimulate(this_action->intArray[0]);
      }
      ptrHeart->step();
      observe();
    }
  }catch(...){
    episode_over=1;
  }
  this_reward_observation.terminal = episode_over;
  this_reward_observation.reward = the_reward;

  return &this_reward_observation;
}

void env_cleanup()
{
   if(flg_mpi){
      MPI_Finalize();
      flg_mpi = false;
   }
	clearRLStruct(&this_observation);
}

const char* env_message(const char* inMessage) {
	if(strcmp(inMessage,"what is your name?")==0)
		return "my name is RippleEnvironment!";

	return "I don't know how to respond to your message";
}

