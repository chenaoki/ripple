#include <stdio.h>
#include <rlglue/RL_glue.h>
#include <iostream>

int whichEpisode=0;

/* Run One Episode of length maximum cutOff*/
void runEpisode(int stepLimit) {        
    int terminal=RL_episode(stepLimit);
	printf("Episode %d\t %d steps \t%f total reward\t %d natural end \n",whichEpisode,RL_num_steps(),RL_return(), terminal);
	whichEpisode++;
}

int main(int argc, char *argv[]) {

	const char* task_spec;
	const char* responseMessage;
	const reward_observation_action_terminal_t *stepResponse;
	const observation_action_t *startResponse;

	printf("\n\nExperiment starting up!\n");

	printf("\n\n----------Initialization----------\n");
	task_spec=RL_init();
	printf("RL_init called, the environment sent task spec: %s\n",task_spec);

	printf("\n\n----------Response check----------\n");
	responseMessage=RL_agent_message("what is your name?");
	printf("Agent responded to \"what is your name?\" with: %s\n",responseMessage);
	responseMessage=RL_env_message("what is your name?");
	printf("Environment responded to \"what is your name?\" with: %s\n",responseMessage);

	//printf("\n\n----------Running a episode----------\n");
	//runEpisode(0);
	//RL_cleanup();

	printf("\n\n----------Stepping through an episode----------\n");
	startResponse=RL_start();
	stepResponse=RL_step();
	while(stepResponse->terminal!=1){
		stepResponse=RL_step();
      /*
		if(stepResponse->terminal!=1){
			printf("(%d,%d) ",stepResponse->o.intArray[0],stepResponse->a.intArray[0]);
		}*/
	}
   std::cout << "done!" << std::endl;
	printf("\n\n----------Summary----------\n");
	printf("It ran for %d steps, total reward was: %f\n",RL_num_steps(),RL_return());
	RL_cleanup();
	
   return 0;
}
