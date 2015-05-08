#include <vector>
#include <iostream>
#include <fstream>
using namespace std;
#include "core/particle_filter.hpp"
#include "random/normal.h"
using namespace mcmc_utilities;
#include <SDL2/SDL.h>

ranlib::Normal<> nrg(0,1);
const int h(512),w(512);

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
    :x(0),vx(0),sx(.1),svx(.01),osx(10),y(0),vy(0),sy(.1),svy(.01),osy(10)
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
  typedef std::vector<double> T_stat;
  typedef std::vector<double> T_obs;
  typedef double T_t;
  std::vector<double> xl,xr;
private:
  T_p do_evol_log_prob(const T_stat& x,const T_t& t,const T_stat& prev_stat,const T_t& prev_t)const
  {
    double sx(.1),svx(.01),osx(10);
    double sy(.1),svy(.01),osy(10);

    //double x1=prev_stat.back()[0];
    double x_prev=prev_stat[0];
    double vx_prev=prev_stat[1];

    double y_prev=prev_stat[2];
    double vy_prev=prev_stat[3];
    
    double x2=x_prev+vx_prev;
    double vx2=vx_prev;

    double y2=y_prev+vy_prev;
    double vy2=vy_prev;

    
    return
      -(x[0]-x2)*(x[0]-x2)/(2*sx*sx)
      -(x[1]-vx2)*(x[1]-vx2)/(2*svx*svx)
      -(x[2]-y2)*(x[2]-y2)/(2*sy*sy)
      -(x[3]-vy2)*(x[3]-vy2)/(2*svy*svy);
  }

  T_p do_obs_log_prob(const T_obs& obs,const T_stat& stat,const T_t& t)const
  {
    //return -(x[0]-y[0])*(x[0]-y[0])/(2*.4*.4);
    double osx(10),osy(10);
    return -(obs[0]-stat[0])*(obs[0]-stat[0])/(2*osx*osx)
      -(obs[1]-stat[2])*(obs[1]-stat[2])/(2*osx*osx);
  }

  void do_stat_var_range(const T_stat& x0,T_stat& xl,T_stat& xr)const
  {
    xl.resize(4);
    xr.resize(4);
    //x1[0]=-10;x2[0]=10;
    T_stat x1(x0),x2(x0);
    x1[0]-=100;
    x1[1]-=20;
    x1[2]-=100;
    x1[3]-=20;

    x2[0]+=100;
    x2[1]+=20;
    x2[2]+=100;
    x2[3]+=20;

    xl=x1;
    xr=x2;
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

void regulate(SDL_Rect& rect,int h,int w)
{
  rect.x*=10;
  rect.y*=10;
  rect.x+=w/2;
  rect.y+=h/2;
  rect.x%=w;
  rect.y%=h;
  while(rect.x<0)
    {
      rect.x+=w;
    }
  while(rect.y<0)
    {
      rect.y+=h;
    }
}
	     
int main()
{

  SDL_Window *window;
  SDL_Renderer *renderer;
  init_display(window,renderer,h,w);
  SDL_Rect rect;

  
  nrg.seed(time(0));
  int nparticles=100;
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
  for(int step=0;step<10000;++step)
    {
      t+=1;
      target1.evolve();
      target1.observe(obs[0],obs[1]);
      
      double x_mean=0;
      double vx_mean=0;
      double y_mean=0;
      double vy_mean=0;

      tm.update_sir(obs,t,particles,prev_t);
     
      double x,vx,y,vy;
      //cout<<t<<" "<<x_mean<<" "<<v_mean<<endl;
      target1.state(x,vx,y,vy);
      
      //cout<<t<<" "<<y<<" "<<y_mean<<" "<<obs[1]<<endl;
      //cout<<t<<" "<<xmean<<" "<<xstd<<endl;


      rect.x=x;
      rect.y=y;
      rect.w=2;
      rect.h=2;
      regulate(rect,h,w);
      
      SDL_SetRenderDrawColor(renderer,0,0,0,255);
      SDL_RenderClear(renderer);
      SDL_SetRenderDrawColor(renderer,255,255,255,255);
      SDL_RenderFillRect(renderer,&rect);

      SDL_SetRenderDrawColor(renderer,0,0,255,255);
      for(int i=0;i<particles.size();++i)
	{
	  x_mean+=particles[i].state[0];
	  vx_mean+=particles[i].state[1];
	  y_mean+=particles[i].state[2];
	  vy_mean+=particles[i].state[3];
	  rect.x=particles[i].state[0];
	  rect.y=particles[i].state[2];

	  regulate(rect,h,w);
	  SDL_RenderFillRect(renderer,&rect);
	  //cout<<particles[i].state[0]<<" "<<particles[i].weight<<endl;
	}
      //return 0;
      x_mean/=particles.size();
      vx_mean/=particles.size();
      y_mean/=particles.size();
      vy_mean/=particles.size();

      rect.x=x_mean;
      rect.y=y_mean;
      regulate(rect,h,w);

      SDL_SetRenderDrawColor(renderer,0,255,255,255);
      SDL_RenderFillRect(renderer,&rect);

      rect.x=obs[0];
      rect.y=obs[1];

      regulate(rect,h,w);
      SDL_SetRenderDrawColor(renderer,255,0,255,255);
      SDL_RenderFillRect(renderer,&rect);


      SDL_RenderPresent(renderer);
      //SDL_Delay(10);  // Pause execution for 3000 milliseconds, for example
      if(step%10==0)
	{
	  //cout<<step<<endl;
	}
      cout<<t<<" "<<x<<" "<<x_mean<<" "<<obs[0]<<endl;
    }
  clear_sdl(window,renderer);
  SDL_Quit();
}
