#include <vector>
#include <iostream>
#include <fstream>
#include "core/particle_filter.hpp"
#include "random/normal.h"
#include <core/urand.hpp>
#include <rng/prng.hpp>
#include <unistd.h>
using namespace std;
using namespace mcmc_utilities;
#include <SDL2/SDL.h>

ranlib::Normal<> nrg(0,1);
const int h(480),w(640);

class target
{
private:
  double x;//position
  double vx;//velocity
  double sx;//stddev of x and y
  double svx;//stddev of vx and vy
  double osx;//stddev of observation value of x and y

  double y;//position
  double vy;//velocity
  double sy;//stddev of x and y
  double svy;//stddev of vx and vy
  double osy;//stddev of observation value of x and y

public:
  target()
    :x(0),vx(0),sx(.1),svx(.01),osx(1),y(0),vy(0),sy(.1),svy(.01),osy(1)
  {}

public:
  void evolve()
  {
    x+=(vx+nrg.random()*sx);
    vx+=nrg.random()*svx;

    y+=(vy+nrg.random()*sy);
    vy+=nrg.random()*svy;

  }

  void state(double& x1,double& vx1,double& y1,double& vy1)const
  {
    x1=x;
    vx1=vx;

    y1=y;
    vy1=vy;

  }

  void observe(double& ox,double& oy)const
  {
    ox=x+osx*nrg.random();
    oy=y+osy*nrg.random();
  }
};



class target_model
  :public pf_model<double,std::vector<double>,std::vector<double>,double>
{
public:
  typedef double T_p;
  typedef std::vector<double> T_state;
  typedef std::vector<double> T_obs;
  typedef double T_t;
  std::vector<double> xl,xr;
public:
  target_model()
    :xl(),xr()
  {}

private:
  T_p do_evol_log_prob(const T_state& x,const T_t& t,const T_state& prev_stat,const T_t& prev_t,int n)const
  {
    double sx(.1),svx(.01);
    double sy(.1),svy(.01);

    double kx(.0),ky(.0);
    //double x1=prev_stat.back()[0];
    double x_prev=prev_stat[0];
    double vx_prev=prev_stat[1];

    double y_prev=prev_stat[2];
    double vy_prev=prev_stat[3];
    
    double x2=x_prev+vx_prev;
    double vx2=vx_prev;

    double y2=y_prev+vy_prev;
    double vy2=vy_prev;

    double result=0;
    if(n==0)
      {
	result= -(x[0]-x2)*(x[0]-x2)/(2*sx*sx);
      }
    else if(n==1)
      {
	result= -(x[1]-vx2)*(x[1]-vx2)/(2*svx*svx);
      }
    else if(n==2)
      {
	result= -(x[2]-y2)*(x[2]-y2)/(2*sy*sy);
      }
    else if(n==3)
      {
	result= -(x[3]-vy2)*(x[3]-vy2)/(2*svy*svy);
      }
    else
      {
	result=
	  -(x[0]-x2)*(x[0]-x2)/(2*sx*sx)
	  -(x[1]-vx2)*(x[1]-vx2)/(2*svx*svx)
	  -(x[2]-y2)*(x[2]-y2)/(2*sy*sy)
	  -(x[3]-vy2)*(x[3]-vy2)/(2*svy*svy);
      }
    assert(isfinite(result));
    return result;
  }

  T_p do_obs_log_prob(const T_obs& obs,const T_state& stat,const T_t& t,int n)const
  {
    //result= -(x[0]-y[0])*(x[0]-y[0])/(2*.4*.4);
    double result=0;
    double osx(1),osy(1);
    if(n==0)
      {
	result= -(obs[0]-stat[0])*(obs[0]-stat[0])/(2*osx*osx);
      }
    else if(n==2)
      {
	result= -(obs[1]-stat[2])*(obs[1]-stat[2])/(2*osy*osy);
      }
    else if(n==1||n==3)
      {
	result= 0;
      }
    else
      {
	result= -(obs[0]-stat[0])*(obs[0]-stat[0])/(2*osx*osx)
	  -(obs[1]-stat[2])*(obs[1]-stat[2])/(2*osy*osy);
      }
    assert(isfinite(result));
    return result;
  }

  std::pair<double,double> do_state_var_range(const T_t& t,const T_state& x0,const T_t& prev_t,size_t ndim)const
  {
    
    if(ndim==0||ndim==2)
      {
	return make_pair(x0[ndim]-100,x0[ndim]+100);
      }
    else
      {
	return make_pair(x0[ndim]-20,x0[ndim]+20);
      }
  }
};

void init_display(SDL_Window*& window,SDL_Renderer*& renderer,int h,int w)
{
  SDL_Init(SDL_INIT_VIDEO);              // Initialize SDL2
  window = SDL_CreateWindow(
			"An SDL2 window",                  // window title
			SDL_WINDOWPOS_UNDEFINED,           // initial x position
			SDL_WINDOWPOS_UNDEFINED,           // initial y position
			w,                               // width, in pixels
			h,                               // height, in pixels
			SDL_WINDOW_OPENGL                  // flags - see below
			);
  renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);;
  if (window == NULL) {
    // In the event that the window could not be made...
    printf("Could not create window: %s\n", SDL_GetError());
    assert(0);
  }
}

void clear_sdl(SDL_Window*& window,SDL_Renderer*& renderer)
{
  SDL_DestroyRenderer(renderer);
  SDL_DestroyWindow(window);
}

void regulate(int& ix,int& iy,double x,double y,int w,int h,double scale)
{
  ix=x*scale+w/2;
  iy=y*scale+h/2;

  ix%=w;
  iy%=h;
  while(ix<0)
    {
      ix+=w;
    }
  while(iy<0)
    {
      iy+=h;
    }
}


void draw_vector(SDL_Renderer& renderer,int w,int h,double x0,double y0,double vx,double vy,double scale)
{
  int ix0,iy0;
  int ix1,iy1;

  regulate(ix0,iy0,x0,y0,w,h,scale);
  //regulate(ix1,iy1,x0+vx*10,y0+vy*10,w,h,scale);
  ix1=ix0+vx*10;
  iy1=iy0+vy*10;

  SDL_RenderDrawLine(&renderer,ix0,iy0,ix1,iy1);
}

void draw_target(SDL_Renderer& renderer,int w,int h,double x,double y,double scale)
{
  SDL_Rect rect;
  rect.w=rect.h=2;
  //regulate(rect,w,h,x,y,scale);
  regulate(rect.x,rect.y,x,y,w,h,scale);
  SDL_RenderFillRect(&renderer,&rect);
}


int main()
{

  SDL_Window *window;
  SDL_Renderer *renderer;
  init_display(window,renderer,h,w);
  //SDL_Rect rect;

  double scale=30;
  nrg.seed(time(0));
  int nparticles=500;
  ranlib::Normal<double> ng(0,1);
  target target1;
  target_model tm;
  std::vector<particle<double,std::vector<double> > > particles(nparticles);
  double prev_t;

  for(int i=0;i<particles.size();++i)
    {
      std::vector<double> x(4);
      x[0]=0;
      x[1]=0;
      x[2]=0;
      x[3]=0;
      
      particles[i].state=x;
      particles[i].weight=1;
    }
  std::vector<double> obs(2);
  srand(time(0));
  obs[0]=0;
  obs[1]=0;
  
  //ifstream ifs("noisy_motion.txt");
  double t=0;
  //urand<double> rng;
  prng<double> rng;
#if 0
  std::vector<std::shared_ptr<base_urand<double> > > rng_array;
  for(int i=0;i<nparticles;++i)
    {
      rng_array.push_back(std::shared_ptr<base_urand<double> >{new prng<double>(i)});
    }
  //rng_array.push_back(std::shared_ptr<base_urand<double> >{new urand<double>});
#endif
  
  for(int step=0;step<10000;++step)
    {
      t+=1;
      target1.evolve();
      target1.observe(obs[0],obs[1]);
      
      double x_mean=0;
      double vx_mean=0;
      double y_mean=0;
      double vy_mean=0;

      tm.update_sir(obs,t,particles,prev_t,rng);
     
      double x,vx,y,vy;
      //cout<<t<<" "<<x_mean<<" "<<v_mean<<endl;
      target1.state(x,vx,y,vy);
      
      //cout<<t<<" "<<y<<" "<<y_mean<<" "<<obs[1]<<endl;
      //cout<<t<<" "<<xmean<<" "<<xstd<<endl;


      
      //regulate(rect,w,h,x,y,scale);
      
      SDL_SetRenderDrawColor(renderer,0,0,0,255);
      SDL_RenderClear(renderer);
      SDL_SetRenderDrawColor(renderer,255,255,255,255);
      //SDL_RenderFillRect(renderer,&rect);
      draw_target(*renderer,w,h,x,y,scale);
      draw_vector(*renderer,w,h,x,y,vx,vy,scale);
      SDL_SetRenderDrawColor(renderer,0,0,255,255);
      for(int i=0;i<particles.size();++i)
	{
	  x_mean+=particles[i].state[0];
	  vx_mean+=particles[i].state[1];
	  y_mean+=particles[i].state[2];
	  vy_mean+=particles[i].state[3];
	  
	  //regulate(rect,w,h,particles[i].state[0],particles[i].state[2],scale);
	  //SDL_RenderFillRect(renderer,&rect);
	  draw_target(*renderer,w,h,particles[i].state[0],particles[i].state[2],scale);
	  //cout<<particles[i].state[0]<<" "<<particles[i].weight<<endl;
	}
      //return 0;
      x_mean/=particles.size();
      vx_mean/=particles.size();
      y_mean/=particles.size();
      vy_mean/=particles.size();

      //regulate(rect,w,h,x_mean,y_mean,scale);

      SDL_SetRenderDrawColor(renderer,0,255,255,255);
      //SDL_RenderFillRect(renderer,&rect);
      draw_target(*renderer,w,h,x_mean,y_mean,scale);
      draw_vector(*renderer,w,h,x_mean,y_mean,vx_mean,vy_mean,scale);
      //regulate(rect,w,h,obs[0],obs[1],scale);
      SDL_SetRenderDrawColor(renderer,255,0,255,255);
      //SDL_RenderFillRect(renderer,&rect);
      draw_target(*renderer,w,h,obs[0],obs[1],scale);

      SDL_RenderPresent(renderer);
      SDL_Event event;
      while (SDL_PollEvent(&event))
	{
	  switch(event.type)
	    {
	    case SDL_KEYDOWN:
	      if(event.key.keysym.sym==SDLK_q)
		{
		  exit(0);
		}
	      break;
	    }
	}
      //SDL_Delay(1000);  // Pause execution for 3000 milliseconds, for example
      if(step%10==0)
	{
	  //cout<<step<<endl;
	}
      //cout<<t<<" "<<x<<" "<<x_mean<<" "<<obs[0]<<endl;
      cout<<t<<" "<<y<<" "<<y_mean<<" "<<obs[1]<<endl;
      //cout<<vx_mean<<" "<<vy_mean<<endl;
      //sleep(1);
    }
  clear_sdl(window,renderer);
  SDL_Quit();
}
