/*
 *  mSector.c
 *  ElsaABM_v1
 *
 *  Created by Christian Bongiorno on 07/03/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */
#include "mQuery.h"
#include "mUtility.h"
#include "mSector.h"
#include "mTest.h"
#include "mABM.h"

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>

//Find extream value of a matrix respect to the column c. if Max=1 find max if max!=1 find min
long double _find_extrem(long double **p,int np,int c,int Max){
	int i;
	long double x;
	if(Max==1){
		for(i=1,x=p[0][c];i<np;i++) if(p[i][c]>x) x=p[i][c];
		return x;
	}
	else{
		for(i=1,x=p[0][c];i<np;i++) if(p[i][c]<x) x=p[i][c];
		return x;		
	}
}

int _point_proj(long double *P,long double *A,long double *B,long double *H){
	
	if(B[0]==A[0]){
		H[0]=B[0];
		H[1]=P[1];
		return 1;
	}
	
	long double m=(B[1]-A[1])/(B[0]-A[0]);
	long double c=A[1]-A[0]*m;
	
	H[0]=(P[0]/m + P[1] - c )/(m + 1./m);
	H[1]=m*((P[0]/m + P[1] - c )/(m + 1./m)) +c;
	
	return 1;
}

int _check_tmp_point(long double *p,CONF_t conf){
	int i;
	if(!point_in_polygon(p, conf.bound, conf.Nbound)) return 1;
	for(i=0;i<(conf.Nbound-1);i++) if( distance_point_segment(conf.bound[i],conf.bound[i+1],p)<DT ) return 1;
	
	return 0;
}

int generate_temporary_point(CONF_t *config){
	
	(*config).n_tmp_nvp=NTMP;
	
	(*config).tmp_nvp = falloc_matrix((*config).n_tmp_nvp, 2);
	int n;
//#ifdef TMP_FROM_FILE
	/*
	char* rep_tail = "/abm_tactical/config/temp_nvp.dat";
	char * rep = malloc(snprintf(NULL, 0, "%s%s", (*config).main_dir, rep_tail) + 1);
	sprintf(rep, "%s%s", (*config).main_dir, rep_tail);
	*/

	if((*config).tmp_from_file){
		//printf("Attention! read temporary nvp from file\n");
		int i;
		char c[500];

		FILE *rstream=fopen((*config).temp_nvp,"r");
		if(rstream==NULL) BuG("Impossible to open temp_nvp.dat\n");
		for(n=0;fgets(c,500,rstream);n++){
			if(n>=NTMP) break;
			(*config).tmp_nvp[n][0]=atof(c);
			for(i=0;c[i]!='\t';i++);
			(*config).tmp_nvp[n][1]=atof(&c[++i]);
			
			#ifdef EUCLIDEAN
			project((*config).tmp_nvp[n],(*config).tmp_nvp[n]);
			#endif

#ifdef BOUND__CONTROL
			/* ATTENTION THIS CONTROL DOES NOT WORKS */
			if( !point_in_polygon((*config).tmp_nvp[n],(*config).bound,(*config).Nbound))
				BuG("Temporary Point outside boundary\n");
#endif
		}
		fclose(rstream);
		(*config).n_tmp_nvp = n;
		//free(rep);	
		return 1;
	}
//#endif
	/* ATTENTION !!!! (DEPRECATED)*/
	long double Mx= _find_extrem((*config).bound,(*config).Nbound,0,1);
	long double mx= _find_extrem((*config).bound,(*config).Nbound,0,0);
	long double My= _find_extrem((*config).bound,(*config).Nbound,1,1);
	long double my= _find_extrem((*config).bound,(*config).Nbound,1,0);
	
	for (n=0; n<(*config).n_tmp_nvp ; n++ ) {
		do{
			(*config).tmp_nvp[n][0]=frand(mx, Mx);
			(*config).tmp_nvp[n][1]=frand(my, My);
		}while ( _check_tmp_point( (*config).tmp_nvp[n],(*config) ));
	}
	FILE *wstream=fopen((*config).temp_nvp,"w");
	for(n=0;n<(*config).n_tmp_nvp;n++) fprintf(wstream,"%Lf\t%Lf\n",(*config).tmp_nvp[n][0],(*config).tmp_nvp[n][1]);
	fclose(wstream);

	//free(rep);
	
	return 1;
}

int cheak_traj_intersect_bound(Aircraft_t *flight,int Nflight,CONF_t config){
	int i,j,k,n_inter;
	for(i=0;i<Nflight;i++) {
		for(j=0,n_inter=0;j<flight[i].n_nvp;j++) 
			for(k=0;k<(config.Nbound-1);k++) 
				if(isbetween(config.bound[k], config.bound[k+1], flight[i].nvp[j])){
					n_inter++;
					break;
				}
		if(n_inter<2) BuG("The current version of the ABM require that the trajectories of all flight have two nvp on the boundary of the sector\n");
		
		if(n_inter>2) BuG("There are trajectories with more than two intersections with the boundary of the sector\n");	
	} 
	return 1;
}

int init_traj_intersect_bound(Aircraft_t **flight,int Nfligth,CONF_t config){
	int i,j,k,c;
	
	for(i=0;i<Nfligth;i++) 
		for(j=0,c=0;j<(*flight)[i].n_nvp&&c<2;j++) 
			for(k=0;k<(config.Nbound-1);k++) 
				if(isbetween(config.bound[k],config.bound[k+1],(*flight)[i].nvp[j])){
					(*flight)[i].bound[c++]=j;
					break;
				}
	return 1;
}

int remove_aircraft(Aircraft_t **fligth, int *Nfligth, int sel){
	Aircraft_t *New_fligth=(Aircraft_t*) malloc((*Nfligth-1)*sizeof(Aircraft_t));
	if(New_fligth==NULL) BuG("Memory BUG\n");
	int i,h;
	for(i=0,h=0;i<*Nfligth;i++) if(i!=sel) New_fligth[h++]=(*fligth)[i];
	free((*fligth) );
	
	(*fligth)=New_fligth;
	*Nfligth= *Nfligth -1;
	
	return 1;
}
//~ FIND The sector related to the point p
int _find_sector(long double *p,SECTOR_t **sec,int n_sect){
	int i;
	for(i=1;i<n_sect;i++)
		if(point_in_polygonOK( &((*sec)[i]),p)) return i;
	return 0;
}

int add_nvp_st_pt(Aircraft_t *f){
	
	long double **nvp=falloc_matrix( ((*f).n_nvp+1),DPOS );
	long double *vel=falloc_vec( (*f).n_nvp );
	long double *time = falloc_vec( (*f).n_nvp+1 );
	
	int i,j;
	for(i=0;i<(*f).st_indx;i++) {
		for(j=0;j<DPOS;j++) nvp[i][j]=(*f).nvp[i][j];
		vel[i]=(*f).vel[i];
		time[i]=(*f).time[i];
		printf("%Lf\n",time[i]);
	}

	
	for(j=0;j<DPOS;j++) nvp[i][j] = (*f).st_point[j];
	vel[i] = (*f).vel[i-1];
	time[i] = time[i-1] + haversine_distance( nvp[i-1], nvp[i] ) / vel[i];
	
	for(i++;i<((*f).n_nvp);i++){
		for(j=0;j<DPOS;j++) nvp[i][j]=(*f).nvp[i-1][j];
		vel[i]=(*f).vel[i-1];
		time[i]=(*f).time[i-1];
	}
	
	for(j=0;j<DPOS;j++) nvp[i][j]=(*f).nvp[i-1][j];
	time[i]=(*f).time[i-1];
	
	//ffree_2D( (*f).nvp, (*f).n_nvp);
	for(i=0;i<(*f).n_nvp;i++) free((*f).nvp[i]);
	free( (*f).nvp );
	
	free( (*f).vel );
	free( (*f).time );
		
	(*f).nvp=nvp;
	
	(*f).vel=vel;
	(*f).time=time;
	
	(*f).n_nvp = (*f).n_nvp +1;
	(*f).st_indx = (*f).st_indx +1;
	(*f).bound[1] = (*f).bound[1]  +1;
		
	return 1;
}
//~ ADD n nvp after st_indx in equally spaced
int add_n_nvp(Aircraft_t *f,CONF_t *conf,SECTOR_t **sec, int n,long double dV){
	
	int i,j;
	long double **nvp=falloc_matrix( (*f).n_nvp+n,DPOS );
	long double *vel=falloc_vec( (*f).n_nvp+ n-1 );
	long double *time = falloc_vec(  (*f).n_nvp+n );
	
	for(i=0;i<=(*f).st_indx;i++) {
		for(j=0;j<DPOS;j++) nvp[i][j] = (*f).nvp[i][j];
		vel[i] = (*f).vel[i];
		time[i] = (*f).time[i];		
	}
	for(;i<(*f).n_nvp;i++) {
		for(j=0;j<DPOS;j++) nvp[i+n][j] = (*f).nvp[i][j];
		vel[i+n-1] = (*f).vel[i-1];
		time[i+n] = (*f).time[i];		
	}
	
	
	long double l = haversine_distance((*f).nvp[(*f).st_indx],(*f).nvp[(*f).st_indx+1]);
	long double dl = l/(n+1);
	for(i=0;i<n;i++){
		vel[(*f).st_indx+i] = (*f).vel[(*f).st_indx];
		nvp[(*f).st_indx+i+1][2] = (*f).nvp[(*f).st_indx][2];
		nvp[(*f).st_indx+i+1][3] = 1.;
		coord_NoVel((i+1)*dl,(*f).nvp[(*f).st_indx],(*f).nvp[(*f).st_indx+1],nvp[(*f).st_indx+(i+1)]);
		nvp[(*f).st_indx+i+1][4] =  _find_sector(nvp[(*f).st_indx+(i+1)],sec,(*conf).n_sect);

		
	}
	
	
	ffree_2D((*f).nvp,(*f).n_nvp);
	free((*f).vel);
	free((*f).time);
	
	(*f).nvp = nvp;
	(*f).vel = vel;
	(*f).time = time;
	(*f).n_nvp = (*f).n_nvp+n;
	
	_position(f, (*f).st_point, (*conf).t_d, (*conf).t_i, (*conf).t_r,dV);
	
	
	return 1;
}

int add_fist_nvpInSec(Aircraft_t *f,CONF_t *conf, SECTOR_t **sec,int next_sec){
	int i;
	int prev_sec = (int) (*f).nvp[(*f).st_indx][4];
	//int next_sec = (int) (*f).nvp[(*f).st_indx+1][4];
	
	if(prev_sec==0||next_sec==0) return 1;
	
	int mynull;
	if( next_sec != prev_sec ){
		for(i=1;i<(*conf).t_d;i++) if(point_in_polygonOK( &((*sec)[next_sec]),(*f).pos[i])){
			//~ 
			//~ plot(*f,*conf,"/tmp/F.dat");
			//~ plot_point((*f).pos[i],"/tmp/point.dat");
			//~ plot_pos(*f,*conf,"/tmp/pos.dat");
			//~ plot_sector((*sec)[prev_sec],"/tmp/prev_S.dat");
			//~ plot_sector((*sec)[next_sec],"/tmp/next_sec.dat");
			//~ 
			
			add_nvp_fromPos(f,i,next_sec,conf);
			//~ plot(*f,*conf,"/tmp/OUT.dat");
			//~ printf("DIRECT: %d %d\n",prev_sec,next_sec);
			//~ scanf("%d",&mynull);
			break;
		}
	}
	return 1;
}

int add_nvp_fromPos(Aircraft_t *f,int indx,int next_sect,CONF_t *conf){
	
	long double **nvp=falloc_matrix( ((*f).n_nvp+1),DPOS );
	long double *vel=falloc_vec( (*f).n_nvp );
	long double *time = falloc_vec( (*f).n_nvp+1 );
	
	int i,j;
	for(i=0;i<((*f).st_indx+1);i++) {
		for(j=0;j<DPOS;j++) nvp[i][j]=(*f).nvp[i][j];
		vel[i]=(*f).vel[i];
		time[i]=(*f).time[i];
	}

	
	for(j=0;j<DPOS;j++) nvp[i][j] = (*f).pos[indx][j];
	nvp[i][4] = next_sect;
	
	vel[i] = (*f).vel[i-1];
	time[i] = time[i-1] + haversine_distance( nvp[i-1], nvp[i] ) / vel[i];
	
	for(i++;i<((*f).n_nvp);i++){
		for(j=0;j<DPOS;j++) nvp[i][j]=(*f).nvp[i-1][j];
		vel[i]=(*f).vel[i-1];
		time[i]=(*f).time[i-1];
	}
	
	for(j=0;j<DPOS;j++) nvp[i][j]=(*f).nvp[i-1][j];
	time[i]=(*f).time[i-1];
	
	//ffree_2D( (*f).nvp, (*f).n_nvp);
	for(i=0;i<(*f).n_nvp;i++) free((*f).nvp[i]);
	free( (*f).nvp );
	
	free( (*f).vel );
	free( (*f).time );
		
	(*f).nvp=nvp;
	
	(*f).vel=vel;
	(*f).time=time;
	
	(*f).n_nvp = (*f).n_nvp +1;
	//(*f).st_indx = (*f).st_indx +1;
	(*f).bound[1] = (*f).bound[1]  +1;
	
	_position(f, (*f).st_point, (*conf).t_d, (*conf).t_i, (*conf).t_r, 0.); //
	return 1;
}

/*
int _add_nvp_bound(Aircraft_t **f,int i,int j,long double **bound,int k){
	int st_indx=j+1;
	long double *p=falloc_vec(DPOS);

	find_intersection((*f)[i].nvp[j],(*f)[i].nvp[j+1],bound[k],bound[k+1],p);
	add_nvp( &( (*f)[i] ),&st_indx,p);
	return 1;
}
*/


int project(long double *in, long double *out){
	long double in2[2];
	in2[0]=in[0];
	in2[1]=in[1];
	
	/*Change this to change Projection*/
	gall_peter(in2,out);
	return 1;
}

int _is_to_add(Aircraft_t f,int xp,long double **bound,int k){
	
	if( ( segments_intersect(f.nvp[xp], f.nvp[xp+1], bound[k], bound[k+1]) && !isbetween(bound[k], bound[k+1],f.nvp[xp+1]) ) && !isbetween( bound[k], bound[k+1],f.nvp[xp]) ){
		return 1;
	   }
	return 0;
} 

int modify_traj_intersect_bound(Aircraft_t **flight,int *Nflight,CONF_t config){
	int i,j;
			
	/*Add flag to nvp*/
	for(i=0;i<*Nflight;i++){
		(*flight)[i].nvp[0][3]=0;
		for(j=1;j<((*flight)[i].n_nvp-2);j++) (*flight)[i].nvp[j][3]=1;
		(*flight)[i].nvp[j][3]=0;
	}
	
	return 1;
}

int set_boundary_flag_onFlight(Aircraft_t **f,int *Nflight, CONF_t c){
	int i,j;
	for(i=0;i<*Nflight;i++) {
		(*f)[i].st_indx=1;
		(*f)[i].ready=0;
		for(j=0;j<((*f)[i].n_nvp-1);j++) {
			if(( point_in_polygon((*f)[i].nvp[j], c.bound, c.Nbound)||point_in_polygon((*f)[i].nvp[j+1], c.bound, c.Nbound) )&&(*f)[i].nvp[j][2]>=F_LVL_MIN){
				(*f)[i].nvp[j][3]=1.;
				}
			 else  (*f)[i].nvp[j][3]=0.;
			}
			   
//		if(point_in_polygon((*f)[i].nvp[j], c.bound, c.Nbound)) (*f)[i].nvp[j][3]=1.;
//		else  (*f)[i].nvp[j][3]=0.;
		
		for(j=0;j<((*f)[i].n_nvp-1);j++) if(point_in_polygon((*f)[i].nvp[j], c.bound, c.Nbound)||point_in_polygon((*f)[i].nvp[j+1], c.bound, c.Nbound)){
			(*f)[i].bound[0]=j;
			break;
		}
		if(j==((*f)[i].n_nvp-1)){
			remove_aircraft(f, Nflight, i--);
			printf("Removed %d Aircraft\n",(*f)[i].ID);
			continue;
		}
		//if((*f)[i].bound[0]==0) (*f)[i].bound[0]=1;
		for(++j;j<((*f)[i].n_nvp);j++) if(!point_in_polygon((*f)[i].nvp[j], c.bound, c.Nbound)&&point_in_polygon((*f)[i].nvp[j-1], c.bound, c.Nbound)) {
			(*f)[i].bound[1]=j;
			break;			
		}
		if(j==((*f)[i].n_nvp)) (*f)[i].bound[1]=j-1;
		   
		
		if((*f)[i].bound[0]==(*f)[i].bound[1]){
			printf("SOmetingh wrong\n");
			(*f)[i].bound[1]=(*f)[i].bound[1]+1;
		}
	}
	

	return 1;
}

int is_on_bound(long double *p,long double **bound,int N){
	int i;
	for(i=0;i<(N-1);i++) if(isbetween(bound[i], bound[i+1], p)) return 1;
	return 0;
}

int _alloc_shock( CONF_t conf,SHOCK_t *shock ){
	(*shock).Nshock = (int) (((conf.f_lvl[1]-conf.f_lvl[0])/10.)*conf.Nm_shock*DAY/(conf.t_w*conf.t_r*conf.t_i));
	(*shock).shock=falloc_matrix( (*shock).Nshock , 6);
	return 1;
}

int _set_cross_timeM1(Aircraft_t **f, int N){
	int i;
	for(i=0;i<N;i++) (*f)[i].delta_t[0]=(*f)[i].time[(*f)[i].bound[1]] - (*f)[i].time[(*f)[i].bound[0]];
	return 1;
}

int _alloc_flight_pos(Aircraft_t **f,int N_f,CONF_t *conf){
	int i;
	for(i=0;i<N_f;i++) (*f)[i].pos = falloc_matrix((*conf).t_d, DPOS);
	return 1;
}

int point_in_polygonOK(SECTOR_t *s,long double *p){

	int   i, j=(*s).n_side-1 ;
	int  oddNodes=0;
	

	for (i=0; i<(*s).n_side; i++) {
		if (((*s).bound[i][1]< p[1] && (*s).bound[j][1]>=p[1] ||   (*s).bound[j][1]< p[1] && (*s).bound[i][1]>=p[1])) {
			oddNodes^=(p[1]*(*s).multiple[i]+(*s).constant[i]<p[0]);
		}
		j=i;
	}

  return oddNodes; 
 }



int precalc_values(SECTOR_t *s) {

	int   i, j=(*s).n_side -1 ;

	(*s).constant = falloc_vec((*s).n_side);
	(*s).multiple = falloc_vec((*s).n_side);

	for(i=0; i<(*s).n_side; i++) {
		if((*s).bound[j][1]==(*s).bound[i][1]) {	  
			(*s).constant[i]=(*s).bound[i][0];
			(*s).multiple[i]=0;
		}
		else {
			(*s).constant[i]=(*s).bound[i][0]-((*s).bound[i][1]*(*s).bound[j][0])/((*s).bound[j][1]-(*s).bound[i][1])+((*s).bound[i][1]*(*s).bound[i][0])/((*s).bound[j][1]-(*s).bound[i][1]);
			(*s).multiple[i]=((*s).bound[j][0]-(*s).bound[i][0])/((*s).bound[j][1]-(*s).bound[i][1]); 
		}
		j=i;
	}

	return 1;
}

int del_bound(SECTOR_t **sector,int N_sect){
	int i;
	for(i=1;i<(N_sect);i++){
		ffree_2D( (*sector)[i].bound,(*sector)[i].n_side );
		free((*sector)[i].multiple);
		free((*sector)[i].constant);

	}
	free(*sector);
	return 1;
}

int _get_bound(char *directory,SECTOR_t *sector, int sel_sect){
	
	char file[R_BUFF];
	char num[R_BUFF];
	sprintf(num, "%d", sel_sect);
	
	strcpy(file,directory);
	strcat(file,num);
	strcat(file,"_bound_latlon.dat");
	
	//printf("%s\n",file);
	/* Count the side of the sector */
	FILE *rstream = fopen(file,"r");
	if(rstream==NULL) BuG("Boundary file does not exist\n");
	int Nbound;
	char c[R_BUFF];
	for((*sector).n_side=0;fgets(c, R_BUFF, rstream);((*sector).n_side)++);
	fclose(rstream);
	
	
	//Actually reading the file
	rstream=fopen(file, "r");
	(*sector).bound = falloc_matrix((*sector).n_side, 2);
	int i,j;
	for(i=0;fgets(c, R_BUFF, rstream)&&i<(*sector).n_side ;i++){
		(*sector).bound[i][0]=atof(c);
		for(j=0;c[j]!='\t'&&c[j]!=' ';j++);
		(*sector).bound[i][1]=atof(&c[++j]);
	}
	fclose(rstream);
	
	//Calculate value useful for poin_in_polygon
	precalc_values(sector);
	
	return 1;
}

int _get_sectors_boundary(CONF_t *config,SECTOR_t **sectors){
	
	// Allocation n_sect
	(*sectors) = (SECTOR_t*) malloc((*config).n_sect * sizeof(SECTOR_t));
	if((*sectors)==NULL) BuG("No Memory\n");
	
	int i;
	
	
	for(i=1;i<(*config).n_sect;i++) _get_bound((*config).bound_file,&((*sectors)[i]),i);

	return 1;
		
}

int init_Sector(Aircraft_t **flight,int *Nflight,CONF_t	*config, SHOCK_t *shock,char *input_ABM, char *config_file,SECTOR_t **sectors){
	
	/*Get configuration parameter from file*/
	get_configuration(config_file, config);
	
	/*Get temporary point from file
	 * The generation is deprecated*/
	generate_temporary_point(config);
	
	/*Get the M1 flight plan*/
	*Nflight=get_M1(input_ABM,flight,config);
	
	/*Set 0 on the activation flag for the first and the last nvp*/
	modify_traj_intersect_bound(flight, Nflight, *config);

	/*allocation of memory for the position matrix*/
	_alloc_flight_pos(flight,*Nflight,config);
	
	/*Allocataion and initiazation of the shock*/
	_alloc_shock(*config,shock);
	get_temp_shock(config);
	
	/*Get the capacity constrains from file*/
	#ifdef CAPACITY
	get_capacity((*config).capacity_file, config);
	_get_sectors_boundary(config,sectors);
	#endif
	
	
	
	return 1;
}





