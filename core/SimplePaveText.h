#ifndef  SimplePaveText_h
#define SimplePaveText_h
#include "TPaveText.h"
#include <string>
class SimplePaveText:public TPaveText
{
  public:
  SimplePaveText(string text, float x1 = 0.68, float y1 = 0.76, float x2 = 0.88, float y2 = 0.96, int textSize = 18, int textFont = 63);
};
SimplePaveText::SimplePaveText(string text, float x1, float y1, float x2, float y2, int textSize, int textFont)
{
  this->SetX1NDC(x1);
  this->SetX2NDC(x2);
  this->SetY1NDC(y1);
  this->SetY2NDC(y2);
  this->SetOption("blNDC");
  this->SetBorderSize(0);
  this->SetFillColor(0);
  this->SetFillStyle(0);
  this->SetTextFont(textFont);
  this->SetTextSize(textSize);
  this->SetTextColor(9);
  this->AddText(text.c_str());
}
#endif

