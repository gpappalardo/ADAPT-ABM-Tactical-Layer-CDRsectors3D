/*
 *  mQuery.c
 *  ElsaABM
 *
 *  Created by Christian Bongiorno on 06/03/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#include "mQuery.h"
#include "mUtility.h"
#include "mSector.h" 
#include "mABM.h"

#include<string.h>
#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<time.h>

/*Take string as HH:MM:SS and return the number of seconds*/

long double _convert_time(char *c ){
	long double time;
	int H,m,s, i,f;
	
	H=atoi(c);
	for(i=0;c[i]!=':';i++);
	m=atoi(&c[++i]);
	for(;c[i]!=':';i++);
	s=atoi(&c[++i]);
	for(;c[i]!=':';i++);
	f=atoi(&c[++i]);
	
	time= (H*3600.)+(m*60.) + (long double) s + (long double) f/1000000.;
	
	return time;
}

int _calculate_velocity(Aircraft_t *flight,int Nflight){
	int i,j;
	long double t;
	for(i=0;i<Nflight;i++){
		for(j=0;j<(flight[i].n_nvp-1);j++) {
			t=(flight[i].time[j+1]-flight[i].time[j]);

			#ifdef EUCLIDEAN
			flight[i].vel[j]=euclid_dist2d(flight[i].nvp[j], flight[i].nvp[j+1])/t;
			#else
			flight[i].vel[j]=haversine_distance(flight[i].nvp[j], flight[i].nvp[j+1])/t;
			#endif
			
		}	
	}
	return 1;
}

int get_M1(char *m1_file,Aircraft_t **flight,CONF_t *conf){
	
	FILE *rstream;
	
	rstream=fopen(m1_file,"r");
	if(rstream==NULL) BuG("M1 File doesn't exist\n");
	
	char c[R_BUFF];
	if(!fgets(c, R_BUFF, rstream)) BuG("Impossible to read M1\n");
	
	int Nflight=atoi(c);

	int i,j,h;

	(*flight)=( Aircraft_t* ) malloc(Nflight*sizeof(Aircraft_t));
	
	#ifdef BOUND_CONTROL
	int inside,m;
	#endif
	
	#ifdef WORKAROUND_NIGHT
	/*It does not manage date*/
	int over_night;
	#endif
	
	struct tm t = {0};
	
	for(i=0;i<Nflight;i++){
		
		#ifdef WORKAROUND_NIGHT
		over_night=0;
		#endif
	
		#ifdef BOUND_CONTROL
		inside=0;
		#endif
	
		if(!fgets(c, R_BUFF, rstream)) BuG("BUG in M1 File -lx0\n");
		(*flight)[i].ready = 0;
		(*flight)[i].ID=atoi(c);
		//int mytest;
		//mytest = (*flight)[i].ID;
		//if(mytest==11424)
		//	mytest = 0;

		for(j=1;c[j]!='\t';j++);
		
		(*flight)[i].n_nvp=atoi(&c[++j]);

		if((*flight)[i].n_nvp==0) BuG("BUG in M1 File: flight has no navpoint\n");
		
		(*flight)[i].nvp=falloc_matrix((*flight)[i].n_nvp,DPOS);
		(*flight)[i].time=falloc_vec((*flight)[i].n_nvp);
		(*flight)[i].vel=falloc_vec((*flight)[i].n_nvp - 1);

		for(h=0;h<(*flight)[i].n_nvp;h++){
			for(++j;c[j]!='\t'&&c[j]!='\0'&&c[j]!='\n';j++);
			if(c[j]=='\0'||c[j]=='\n') BuG("BUG in M1 File -lx1\n");
			(*flight)[i].nvp[h][0]=atof(&c[++j]);
			
			for(++j;c[j]!=','&&c[j]!='\0'&&c[j]!='\n';j++);
			if(c[j]=='\0') BuG("BUG in M1 File -lx2\n");
			(*flight)[i].nvp[h][1]=atof(&c[++j]);
			
			#ifdef EUCLIDEAN
			project((*flight)[i].nvp[h],(*flight)[i].nvp[h]);
			#endif

			for(++j;c[j]!=','&&c[j]!='\0'&&c[j]!='\n';j++);
			if(c[j]=='\0') BuG("BUG in M1 File -lx3\n");
			(*flight)[i].nvp[h][2]=atof(&c[++j]);
			
			for(++j;c[j]!=','&&c[j]!='\0';j++);
			if(c[j]=='\0') BuG("BUG in M1 File -lx4\n");
			
			
			strptime(&c[++j],"%Y-%m-%d %H:%M:%S",&t);
			(*flight)[i].time[h]=(long double) mktime(&t);
			
			//~ for(++j;c[j]!=' '&&c[j]!='\0';j++);
			//~ if(c[j]=='\0') BuG("BUG in M1 File -lx5\n");
			//~ (*flight)[i].time[h]=_convert_time(&c[++j]);
			
			#ifdef WORKAROUND_NIGHT
			if(h>0 && over_night==0){
				if((*flight)[i].time[h]<(*flight)[i].time[h-1]) over_night = 1;
			}
			if(over_night!=0) ((*flight)[i].time[h])+=DAY;
			#endif
			
			
			#ifdef CAPACITY
			for(++j;c[j]!=','&&c[j]!='\0';j++);
			if(c[j]=='\0') BuG("BUG in M1 File -lx6\n");
			(*flight)[i].nvp[h][4]=atof(&c[++j])+WA_SECTOR_LABEL; // TAKE OUTTTT
			#endif		
			
			#ifdef BOUND_CONTROL
			/*It does not work! */
			if(inside==0) if(point_in_polygon((*flight)[i].nvp[h],(*conf).bound,(*conf).Nbound))inside=1;
			#endif

		}
		#ifdef BOUND_CONTROL
		if(inside==0) {
			printf("Flight %d doesn't cross the boundary\n",(*flight)[i].ID);
			exit(0);
		}
		#endif
	}
	
	/*Evaluate velocity as the mean velocity between two NVPs*/
	_calculate_velocity((*flight),Nflight);
	
	fclose(rstream);
	if((*conf).use_delay_file==1){
		rstream=fopen((*conf).delay_file  ,"r");
		if(rstream==NULL) BuG("Delay File doesn't exist\n");
		if(!fgets(c, R_BUFF, rstream)) BuG("Impossible to read Delay File\n");
	
		long int fid,delay,n;
		for(n=0;fgets(c,500,rstream);n++){
			if(n>=Nflight) break;
			fid =atof(c);
			for(i=0;c[i]!='\t';i++);
			delay=atof(&c[++i]);
			//printf("%d %d\n",fid,delay);
			int tr = -1;
			for(int jj=0;jj<Nflight;jj++){
				if((*flight)[jj].ID == fid){
					tr = jj;
					break;
				}
			}
			if(tr==-1) {
				printf("%ld",fid);
				printf("non trovato!\n");

			}
			else (*flight)[tr].tw = delay;
		}
			fclose(rstream);
	}
	return Nflight;
}

long double _find_value_string(char *config_file,char *label){

	FILE *rstream=fopen(config_file, "r");
	if(rstream==NULL){
		printf("I was looking here %s for the config file\n", config_file);
		BuG("BUG - configuration file doesn't exist\n");
	}
	
	char c[R_BUFF];
	int i,lsize;
	
	for(lsize=0;label[lsize]!='\0';lsize++);
	
	while (fgets(c, R_BUFF, rstream)) {
		if(c[0]=='#'||c[0]=='\n'||c[0]==' ') continue;
		for(i=0;c[i]!='#'&&c[i]!='\0';i++);
		if(c[i]=='\0') BuG("configuration file not standard\n");
		if(!memcmp(&c[++i], label, lsize)) {
			fclose(rstream);
			return atof(c);
		}
	}
	
	fclose(rstream);
	printf("Impossible to find %s in config-file\n",label);
	exit(0);
}
long double _find_value_datetime (char *config_file,char *label){

	FILE *rstream=fopen(config_file, "r");
	if(rstream==NULL){
		printf("I was looking here %s for the config file\n", config_file);
		BuG("BUG - configuration file doesn't exist\n");
	}
	struct tm t = {0};
	char c[R_BUFF];
	int i,lsize;
	
	for(lsize=0;label[lsize]!='\0';lsize++);
	
	while (fgets(c, R_BUFF, rstream)) {
		if(c[0]=='#'||c[0]=='\n'||c[0]==' ') continue;
		for(i=0;c[i]!='#'&&c[i]!='\0';i++);
		if(c[i]=='\0') BuG("configuration file not standard\n");
		if(!memcmp(&c[++i], label, lsize)) {
			fclose(rstream);
			
			strptime(c,"%Y-%m-%d %H:%M:%S",&t);
			return (long double) mktime(&t);	
		}
	}
	
	fclose(rstream);
	printf("Impossible to find %s in config-file\n",label);
	exit(0);
}


char * _find_value_string_char(char *config_file,char *label){
	/*
	Same function as previous for strings.
	*/
	FILE *rstream=fopen(config_file, "r");
	if(rstream==NULL){
		printf("I was looking here %s for the config file\n", config_file);
		BuG("BUG - configuration file doesn't exist\n");
	}
	
	char *c = (char*) malloc((R_BUFF+1)*sizeof(char));
	if(c==NULL) BuG("Memory\n");
	//char c[R_BUFF];
	int i,lsize;
	//char *d;
	for(lsize=0;label[lsize]!='\0';lsize++);
	
	while (fgets(c, R_BUFF, rstream)) {
		if(c[0]=='#'||c[0]=='\n'||c[0]==' ') continue;
		for(i=0;c[i]!='#'&&c[i]!='\0';i++);
		if(c[i]=='\0') BuG("configuration file not standard\n");
		if(!memcmp(&c[++i], label, lsize)) {
			fclose(rstream);
			for(i=0;c[i]!='\t';i++);
			c[i]='\0';
			return c;
		}
	}
	
	fclose(rstream);
	printf("Impossible to find %s in config-file\n",label);
	exit(0);
}

int get_configuration(char *config_file,CONF_t *config){
	
	/*It searches for a #string in the configuration file and it assignes the left value*/
	(*config).conf_ang = _find_value_string(config_file,"max_ang");
	(*config).extr_ang = _find_value_string(config_file,"extr_ang");
	(*config).nsim = (int) _find_value_string(config_file,"nsim");
	(*config).n_polygon = atoi(_find_value_string_char(config_file,"npol"));

	(*config).x_capacity = _find_value_string(config_file,"x_capacity");
	(*config).direct_thr = _find_value_string(config_file,"direct_thr");
	(*config).xdelay = _find_value_string(config_file,"xdelay");
	(*config).pdelay = _find_value_string(config_file,"pdelay");
	(*config).laplacian_vel = (int) _find_value_string(config_file,"laplacian_vel");
	(*config).Nm_shock = (int) _find_value_string(config_file,"Nm_shock");
	(*config).radius = _find_value_string(config_file,"radius");
	(*config).t_w = (int) _find_value_string(config_file,"t_w");
	(*config).t_d = (int) _find_value_string(config_file,"t_d");	
	(*config).t_i = _find_value_string(config_file,"t_i");
	(*config).t_r = _find_value_string(config_file,"t_r");
	(*config).shortest_path = (int) _find_value_string(config_file,"shortest_path");
	(*config).d_thr = (int) _find_value_string(config_file,"d_thr");
	(*config).f_lvl[0] = _find_value_string(config_file, "shock_f_lvl_min");
	(*config).f_lvl[1] = _find_value_string(config_file, "shock_f_lvl_max");
	(*config).geom = _find_value_string(config_file, "geom");
	(*config).sig_V = _find_value_string(config_file, "sig_V");
	(*config).tmp_from_file = _find_value_string(config_file, "tmp_from_file");
	(*config).lifetime = _find_value_string(config_file, "lifetime");
	(*config).ls = _find_value_string(config_file,"ls");
	(*config).as = _find_value_string(config_file,"as");
	(*config).max_t = _find_value_string(config_file,"max_t");
	(*config).rer_active = (int) _find_value_string(config_file,"rer_active");
	//printf("%s\n", config_file);
	//exit(0);
	//(*config).main_dir = "/home/earendil/Documents/ELSA/ABM/ABM_FINAL";
	//(*config).main_dir = _find_value_string_char(config_file, "main_dir");
	(*config).temp_nvp = _find_value_string_char(config_file, "temp_nvp");
	(*config).shock_tmp = _find_value_string_char(config_file, "shock_tmp");
	(*config).bound_file = _find_value_string_char(config_file, "bound_file");
	(*config).capacity_file = _find_value_string_char(config_file, "capacity_file");
	(*config).use_delay_file = (int)_find_value_string(config_file,"use_delay");
	//(*config).delay_file = _find_value_string_char(config_file,"delay_file");
	
	(*config).start_datetime = _find_value_datetime(config_file, "start_datetime");
	(*config).end_datetime = _find_value_datetime(config_file, "end_datetime");
	
	(*config).noise_d_thr = _find_value_string(config_file,"noise_d_thr");
	
	return 1;
}

int get_boundary( CONF_t *config ){
	FILE *rstream=fopen((*config).bound_file, "r");
	
	if(rstream==NULL) BuG("No Bound File found\n");
	
	//Finding the number of Points on the boundaries
	char c[R_BUFF];
	int Nbound,i,j;
	for(Nbound=0;fgets(c, R_BUFF, rstream);Nbound++);
	(*config).Nbound=Nbound;
	fclose(rstream);
	
	//Actually reading the file
	rstream=fopen((*config).bound_file, "r");
	(*config).bound = falloc_matrix(Nbound, 2);
	for(i=0;fgets(c, R_BUFF, rstream)&&i<Nbound;i++){
		(*config).bound[i][0]=atof(c);
		for(j=0;c[j]!='\t'&&c[j]!=' ';j++);
		(*config).bound[i][1]=atof(&c[++j]);
	}
	fclose(rstream);
	
	if((*config).bound[0][0]!=(*config).bound[(*config).Nbound-1][0] || (*config).bound[0][1]!=(*config).bound[(*config).Nbound-1][1] )
		BuG("Last point of the Boundary has to be the same than the first.\n");
	
	return 1;
}

int get_temp_shock(CONF_t *conf){
	/*
	char* rep_tail = "/abm_tactical/config/shock_tmp.dat";
	char * rep = malloc(snprintf(NULL, 0, "%s%s", (*conf).main_dir, rep_tail) + 1);
	sprintf(rep, "%s%s", (*conf).main_dir, rep_tail);
	*/

	FILE *rstream=fopen((*conf).shock_tmp,"r");
	if(rstream==NULL) BuG("Impossible to read shock_tmp.dat\n");
	
	int i,j,N;
	char c[500];
	for(N=0;fgets(c,500,rstream);N++);
	fclose(rstream);
	
	(*conf).n_point_shock = N;
	(*conf).point_shock = falloc_matrix(N, 2);
	
	rstream=fopen((*conf).shock_tmp,"r");
	for(i=0;fgets(c,500,rstream)&&i<N;i++){
		(*conf).point_shock[i][0]=atof(c);
		for(j=0;c[j]!='\t'&&c[j]!=' '&&c[j]!='\0';j++);
		if(c[j]=='\0') BuG("Error in shock_tmp.dat");
		(*conf).point_shock[i][1]=atof(&c[++j]);
		
		#ifdef EUCLIDEAN
		project((*conf).point_shock[i],(*conf).point_shock[i]);
		#endif 
	}
	fclose(rstream);
	//free(rep);
	return 1;
	
}

int add_nsim_output(char *file_out,char *file_in, int n){

	int i;
	for(i=0;file_in[i]!='\0';i++);
	for(;file_in[i]!='.';i--);
	int j;
	char str[50];
	sprintf(str, "%d", n);
	
	for(j=0;j<i;j++) file_out[j]=file_in[j];
	file_out[j]='\0';
	
	strcat(file_out,"_");
	strcat(file_out,str);
	strcat(file_out,".dat");
	
	return 1;
	
}

int collision_counter(char *file_output,CONF_t *c){
	int i,j;
      	for(j=0;file_output[j]!='\0';j++);
        for(;file_output[j]!='.';j--);

	for(i=0;i<j;i++) (*c).output_CDF[i] = file_output[i];
	(*c).output_CDF[i] = '\0';
	
	strcat((*c).output_CDF,"_Counter.dat");
	
	FILE *wstream = fopen((*c).output_CDF,"w");
	if(wstream==NULL) BuG("No Counter File\n");
	fclose(wstream);
	
	return 1;
	
}

int print_operation(long double *point ,long double tp,int idf,char *typ,CONF_t cnf){
	FILE *wstream = fopen(cnf.output_CDF,"a");
	if(wstream==NULL) BuG("Counter file is no more there\n");
	//~ 
	time_t pT;
	struct tm pTm={0};
	char buff[100];
	
				
	pT = (time_t) tp;
	localtime_r(&pT, &pTm);
	
	strftime(buff,100,"%Y-%m-%d %H:%M:%S",&pTm);
	
	
	fprintf(wstream,"%s\t%s\t(%Lf,%Lf)\t%d\n",typ,buff,point[0],point[1],idf);
	
	fclose(wstream);
	
	return 1;
}
		


int get_capacity(char *file_r,CONF_t *conf){
	
	FILE *rstream=fopen(file_r,"r");
	if(rstream==NULL) BuG("Miss Capacity file\n");
	
	char c[R_BUFF];
	for((*conf).n_sect=0;fgets(c, R_BUFF, rstream);((*conf).n_sect)++) if(c[0]=='#') ((*conf).n_sect)--;
	fclose(rstream);
	
	((*conf).n_sect) += 2;
	(*conf).capacy = ialloc_vec((*conf).n_sect);
	
	rstream=fopen(file_r,"r");
	int i,j;
	for(i=1;fgets(c, R_BUFF, rstream);i++){
		if(c[0]=='#'){
			i--;
			continue;
		}
		if(atoi(c)!=i) BuG("Not Regular Capacity file, miss index\n");
		for(j=0;c[j]!='\t';j++);
		(*conf).capacy[i]=atoi(&c[j+1])*(*conf).x_capacity;
		

	}
	(*conf).capacy[0]=SAFE;
	(*conf).capacy[(*conf).n_sect-1] = SAFE;
	fclose(rstream);
	
	return 1;	
}
