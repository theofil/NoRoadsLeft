#include "TMath.h"
#include <iostream>

#ifndef pnumber_h
#define pnumber_h
class pnumber
{
  public:
  pnumber(){x=0;xE=0;}
  pnumber(float n){x=n;xE=sqrt(n);if(n<40)std::cout << "warning (pnumber) "<< n <<" has non Gaussian error bar" << std::endl;}
  pnumber(float n, float nE){x=n;xE=nE;}
  void set(float n, float nE){x=n;xE=nE;}
   
  pnumber operator + (pnumber c)
  {  
    pnumber tmp; 
    tmp.x = x + c.x;
    tmp.xE = sqrt(c.xE*c.xE + xE*xE);
    return tmp;
  }

  pnumber operator - (pnumber c)
  {  
    pnumber tmp; 
    tmp.x = x - c.x;
    tmp.xE = sqrt(c.xE*c.xE + xE*xE);
    return tmp;
  }

  pnumber operator = (pnumber c)
  {  
    x = c.x; 
    xE = c.xE;
    return *this;
  }

  pnumber operator / (pnumber c)
  {  
    pnumber tmp; 
    tmp.x = x / c.x;
    if( !(x == 0 || c.x ==0) )tmp.xE = ( x / c.x ) * sqrt(  (c.xE*c.xE)/(c.x * c.x) + (xE*xE)/(x*x)  ) ;
    if(  (x == 0 || c.x ==0) )tmp.xE =0;

    return tmp;
  }

  pnumber operator * (pnumber c)
  {  
    pnumber tmp; 
    tmp.x = x * c.x;
    if( !(x == 0 || c.x ==0) )tmp.xE = ( x * c.x ) * sqrt(  (c.xE*c.xE)/(c.x * c.x) + (xE*xE)/(x*x)  ) ;
    if(  (x == 0 || c.x ==0) )tmp.xE =0;
    return tmp;
  }
  
  pnumber pow(float n)
  {
    //this calculates e.g. sqrt(a), where b=(0.5,0)
    pnumber result;
    if(x!=0)
    {
      result.x  = TMath::Power(x,n);
      result.xE = fabs(n*TMath::Power(x,n-1)*xE);
    }else
    {
      result.x  = 0;
      result.xE  = 0;
    }
    return result;
  }

  float x;
  float xE;
};

std::ostream &operator<<(std::ostream &ostr, pnumber v){return ostr <<  v.x  << " +/- " << v.xE;}
#endif
