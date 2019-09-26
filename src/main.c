#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mQuery.h"
#include "mSector.h"
#include "mUtility.h"
#include "mABM.h"
#include "mTest.h"
//#include <Python.h>

// Entry for wrapper.
int simulation(char **args){

	srand(atoi(args[4]));
	
	setenv("TZ", "GMT", 1);
	tzset();
	
	CONF_t config;
	Aircraft_t *Flight,*flight;
	int Nflight;
	SHOCK_t shock;
	SECTOR_t *sectors;
	
	char *input_ABM=args[1];
	char *output_ABM=args[2];
	char *config_file=args[3];
	char output_ABM_nsim[400]; //TODO: change this!
	
	/*Initialization of Variable*/
	init_Sector(&Flight,&Nflight,&config,&shock,input_ABM, config_file,&sectors);
	
	//printf("Finished Initialization\n");
	int i;
	/* run nsim simulation with the same M1 file*/
	for(i=0;i<config.nsim;i++) {
		
		/*Create a backup copy for the flight*/
		copy_flight(Flight,Nflight,&flight);
		
		add_nsim_output(output_ABM_nsim,output_ABM,i);
		collision_counter(output_ABM_nsim,&config);
		
		if( ABM(&flight,Nflight,config,shock,sectors) == 0){
			/*if ABM does not solve the conflicts It run again the simulation*/
			del_flight(&flight, Nflight, Flight);
			i--;
			continue;
		}
		
		//printf("Sim %d\n",i+1);
		/*create and save the ouptput file*/
		save_m3(flight,Nflight,Flight,output_ABM_nsim);
		
		del_flight(&flight, Nflight, Flight);
	}
	//~ 
	/*free the memory*/
	del_flight_pos(&Flight,Nflight,config);
	del_flight(&Flight, Nflight, Flight);
	del_conf(&config);
	del_shock(&shock,&config);
	
	del_bound(&sectors,config.n_sect);
	
	return 0;
}

// Manual entry
int main(int argc,char *argv[]){
	simulation(argv);
	return 0;
}


