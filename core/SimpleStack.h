#ifndef SimpleStack_h
#define SimpleStack_h
class SimpleStack:public THStack
{
  public:
  using THStack::THStack;
  std::vector<string> sampleTitles_;
  std::vector<TH1F*>  histoPointers_;
};
#endif
