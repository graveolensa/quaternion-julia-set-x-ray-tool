#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>


/*
Quaternion Julia Set X-Ray tool. 

Instead of dipping quaternion julia sets in opengl paint, when compiled, the following
will allow you to construct rotations of quaternion julia sets as rasters. Example
videos which were generated using this code:

http://www.youtube.com/watch?v=zNiYZnd0rn4
http://www.youtube.com/watch?v=au_hTYx_P4Y
http://www.youtube.com/watch?v=zUWNVRidCyE
http://www.youtube.com/watch?v=gCJKJrbcxFo

Owen Maresh
http://owen.maresh.info/

*/


typedef struct{
  double re;
  double i;
  double j;
  double k;
} quaternion;

typedef struct{
  double x;
  double y;
  double z;
} point;

quaternion add(quaternion m, quaternion n)
{
  quaternion out;
  out.re = m.re + n.re;
  out.i=m.i+n.i;
  out.j=m.j+n.j;
  out.k=m.k+n.k;
  return out;
}

double norm2(quaternion m)
{
  double out;
  out = m.re*m.re + m.i*m.i + m.j*m.j + m.k*m.k;
  return out;
}

double norm(quaternion m)
{
  return sqrt(norm2(m));
}


quaternion conjugate(quaternion m)
{
  quaternion out;
  out.re = m.re;
  out.i = -1.0 * m.i;
  out.j = -1.0 * m.j;
  out.k = -1.0 * m.k;
  return out;
}


quaternion multiply(quaternion m, quaternion n)
{
  quaternion out;
  out.re = (m.re*n.re) - (m.i*n.i)  - (m.j*n.j) - (m.k*n.k); 
  out.i  = (m.re*n.i)  + (m.i*n.re) + (m.j*n.k) - (m.k*n.j);
  out.j  = (m.re*n.j)  + (m.j*n.re) + (m.k*n.i) - (m.i*n.k);
  out.k  = (m.re*n.k)  + (m.k*n.re) + (m.i*n.j) - (m.j*n.i);
  return out;
}

quaternion rmult(double s, quaternion m)
{
  quaternion out;
  out.re = m.re * s;
  out.i  = m.i  * s;
  out.j  = m.j  * s;
  out.k  = m.k  * s;
  return out;
}

quaternion reciprocal(quaternion m)
{
  quaternion out;
  out = rmult(1/norm2(m),conjugate(m));
  return out;
}




quaternion divide(quaternion m, quaternion n)
{
  quaternion out;
  if(norm2(n) == 0)
    {
      printf("divide by zero attempted in quaternion divide routine\n"); 
   
    }
  else
    {      
      out = multiply(m,reciprocal(n));

    }
  return out;
} 


quaternion normscale(quaternion m)
{
  /* note: this kind of assumes that the quaternion isn't zero
     but given that it has extremely limited use, I don't have
     to sanity check its input */
  quaternion out;
  out = rmult(1/norm(m),m);
  return out;
}


quaternion negative(quaternion m)
{
  return rmult(-1.0,m);
}

quaternion difference(quaternion m, quaternion n)
{
  quaternion out;
  out.re = m.re - n.re;
  out.i = m.i - n.i;
  out.j = m.j - n.k;
  out.k = m.k - n.k;
  return out;
}

quaternion square(quaternion m)
{
  return multiply(m,m);
}

quaternion cube(quaternion m)
{
  return multiply(square(m),m);
}

quaternion quart(quaternion m)
{
  return square(square(m));
}

quaternion quint(quaternion m)
{
  return multiply(m,quart(m));
}

quaternion threeadd(quaternion m, quaternion n, quaternion o)
{
  quaternion out;
  return add(m,add(n,o));
}

quaternion fouradd(quaternion m, quaternion n, quaternion o, quaternion p)
{
  quaternion out;
  return add(add(m,n),add(o,p));
}

quaternion sixadd(quaternion m, quaternion n, quaternion o, quaternion p, quaternion q, quaternion r)
{
  return add(threeadd(m,n,o),threeadd(p,q,r));
}

quaternion fourmult(quaternion m, quaternion n, quaternion o, quaternion p)
{
  quaternion out;
  out = multiply(multiply(m,n),multiply(o,p));
  return out;
}

quaternion threemult(quaternion m, quaternion n, quaternion o)
{
  quaternion out;
  out = multiply(multiply(m,n),o);
  return out;
}

point interpret(quaternion m)
{
  point out;
  out.x = m.re;
  out.y = m.i;
  out.z = m.j;
  return out;
}


quaternion w1,w2,w3,w4;
quaternion i,j,k,one;
quaternion mi,mj,mk,mone;
quaternion zero;
quaternion current,next;

double tx,ty,tz,rx,ry;

int display[512][512];
int A,B,C;

/* we completely don't care about the order of the points */
point *firechar;
point zull,store;
double tol;
int vba,vbb;
int iter;
int level;
double s;
double costheta,sintheta;
double theta;
int n;
double TX,TY;
int m,l;
FILE *file;
int zullcount;
double vranj,cosvranj,sinvranj;
double ax,az,rz;


int main(int argc, char *argv[])
{
 
  char name_buf[20];
  srand48(time(NULL));

 /* define a non-point as zull */
  zull.x = -23.0;
  zull.y = -33.0;
  zull.z = -10.0;
  /* define zero */
  zero.re = 0.0; zero.i = 0.0 ; zero.j = 0.0 ; zero.k = 0.0 ;

  /* define some units */
  one.re = 1.0 ; one.i = 0.0 ; one.j = 0.0 ; one.k = 0.0;

  i.re   = 0.0 ; i.i = 1.0   ; i.j = 0.0   ; i.k = 0.0;

  j.re   = 0.0 ; j.i = 0.0   ; j.j = 1.0   ; j.k = 0.0;

  k.re   = 0.0 ; k.i = 0.0   ; k.j = 0.0   ; k.k = 1.0;

  /* negate them */
  mone   = negative(one);
  mi     = negative(i);
  mj     = negative(j);
  mk     = negative(k);

  /* now define the unscaled roots */
  w1 = threeadd(i,j,k);
  w2 = threeadd(mi,mj,k);
  w3 = threeadd(i,mj,mk);
  w4 = threeadd(mi,j,mk);
 
  /* rescale the unscaled roots so that they're on the unit sphere */
  w1 = normscale(w1);
  w2 = normscale(w2);
  w3 = normscale(w3);
  w4 = normscale(w4);

  printf("constants defined\n");

  quaternion constant;
  constant.re = -0.450;
  constant.i = -0.477;
  constant.j = 0.181;
  constant.k = 0.306;
  /* define the polynomial */
  quaternion julia(quaternion m)
  {
    return add(square(m),constant);
  }

  printf("newton's method defined\n");

  /* linear point array 256^3 = 16777216 */
  firechar = malloc( 27270901* sizeof(point));

  printf("storage array generated\n");
  /* clear the display array */
  tol = 0.001;


  for(A=0;A<301;A++)
    {
      for(B=0;B<301;B++)
  {
	  for(C=0;C<301;C++)
	    {
	      current.re = ((double) (A - 151))/150.0 + 0.001* (drand48()-1);
	      current.i = ((double) (B - 151))/150.0 + 0.001* (drand48()-1);
	      current.j = ((double) (C - 151))/150.0; + 0.001* (drand48()-1);
	      current.k = 0.001* (drand48()-1);
	      store = interpret(current);
	      firechar[90601*A+301*B+C].x=store.x;
	      firechar[90601*A+301*B+C].y=store.y;
	      firechar[90601*A+301*B+C].z=store.z;
	    }
	}
    }

  /* generate the data array */
  for(A=0;A<301;A++)
    {
      for(B=0;B<301;B++)
	{
	  for(C=0;C<301;C++)
	    {
	      /* printf("I got this far, no really! A=%d B=%d C=%d\n",A,B,C); */

              current.re = ((double) (A - 151))/150.0 + 0.001* (drand48()-1);
	      current.i = ((double) (B - 151))/150.0 + 0.001* (drand48()-1);
	      current.j = ((double) (C - 151))/150.0 + 0.001* (drand48()-1);
	      current.k = 0.001* (drand48()-1);
	      store = interpret(current);
	      iter = 0;
	      do
		{	     
		  next=julia(current);
		  current=next;
		  iter=iter+1;
		} while ((norm(current)<2)&&(iter<30));
	      /* printf("%d, %.10f\n", iter, norm(current)); */
	      if(norm(current)>2)
		{
		  firechar[90601*A+301*B+C] = zull;
		  /* printf("setting a zull\n"); */
		}
	    }
	}
    }

  zullcount=0;
  for(n=0;n<27270901;n++)
    {
      if(firechar[n].x==-23.0)
	{
	  zullcount=zullcount+1;
	}

    }
  printf("total cells occuped %d out of \n", 27270901-zullcount);
  printf("percentage occupied = %f\n", ((double)(27270901-zullcount))/27270901);


  for (level = 0; level < 1024; level++) 
    {
      theta = 6.283185307179586 * ((double) (level))/1024.0;
      costheta = cos(theta);
      sintheta = sin(theta);
      vranj = 2.718281828459045 * theta;
      cosvranj = cos(vranj);
      sinvranj = sin(vranj);
      /* define the zero point for reference purposes */
      
      /* white the display arraw */

        for(vba=0;vba<512;vba++)
	{
	  for(vbb=0;vbb<512;vbb++)
	    {
	      display[vba][vbb]=0;
	    }
	}


      for(n=0;n<27270901;n++)
	{
	  if(firechar[n].x!=-23.0)
	    {
	      tx=firechar[n].x;
	      ty=firechar[n].y;
	      tz=firechar[n].z;
	      ax = tx * cosvranj - tz * sinvranj;
	      az = tx * sinvranj + tz * cosvranj; 
	      rz = az;

	      rx = ax * costheta - ty * sintheta;
	      ry = ax * sintheta + ty * costheta;
	      /* printf("I got this far, no really!\n");     */ 	      
	      TX=1.61*rx+0.00010394;
	      TY=1.61*rz+0.00010394;
	      
	      m = (int) ((TX + 2) / 4 * 511);
	      l = (int) ((TY + 2) / 4 * 511);
	      

	      /* printf("m=%d l=%d\n", m,l); */
	      /* display[m][l]=0; */
	      if(display[l][m]<255)
	      	{
	      	  display[l][m]=display[l][m]+1;
	      	}
	      /* if(display[l][m]>0) */
	      /* 	{ */
	      /* 	  if((m>0)&&(m<512)) */
	      /* 	    { */
	      /* 	      if((l=0)&&(l<512)) */
	      /* 		{ */
	      /* 		  display[l][m]=0; */
	      /* 		} */
	      /* 	    } */
	      /* 	} */
	    }
	}

  
      sprintf(name_buf,"%04d.pnm", level);
      file = fopen(name_buf, "w");
      fprintf(file, "P2\n512 512\n255\n");
      for(vba=0;vba<512;vba++)
	{
	  for(vbb=0;vbb<512;vbb++)
	    {
	      fprintf(file,"%u\n", display[vba][vbb]);
	    }
	}
      fclose(file);
    }






  return 0;
}
