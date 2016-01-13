#ifndef __STIMELECTRODE__
#define __STIMELECTRODE__

#include <exception>
#include <vector>
#include <map>

namespace ripple{

template <typename Pos, typename Cur>
class StimElectrode
{

  protected:
    std::map<unsigned int, Cur> mapTimeCurrent;

  public:
   StimElectrode(void){};
   StimElectrode( const StimElectrode& obj)
     : mapTimeCurrent( obj.mapTimeCurrent ){};
   virtual ~StimElectrode(void){};

  public:
   void setCurrent(unsigned int time_, const Cur& current_)
   {
     this->mapTimeCurrent[time_] = current_;
   };
   virtual Cur getCurrent(unsigned int time, Pos x, Pos y, Pos z) = 0;
};

template <typename Pos, typename Cur>
class AxisStimElectrode : public StimElectrode<Pos, Cur>
{

  public:
    typedef struct Axis{ 
      unsigned int numAxis; 
      Pos p;
      Pos distance; 
    } Axis;

  private:
    std::vector<Axis> vecAxis;

  public:
   AxisStimElectrode(void)
     : StimElectrode<Pos,Cur>(){};
   AxisStimElectrode(const AxisStimElectrode& obj) 
     : StimElectrode<Pos,Cur>(obj), vecAxis(obj.vecAxis){};
   virtual ~AxisStimElectrode(void){};
  
  public:
   void addAxis(unsigned int numAxis_, Pos p_, Pos distance_)
   {
     if( distance_ < 0 || numAxis_ > 2){
       throw std::string("error");
     }
     Axis ax;
     ax.numAxis = numAxis_; ax.p = p_; ax.distance = distance_;
     this->vecAxis.push_back(ax);
     return;
   };

   Cur getCurrent(unsigned int time, Pos x, Pos y, Pos z){
     bool flag = true;
     Pos dist = 0;
     for(typename std::vector<Axis>::const_iterator it = vecAxis.begin(); it != vecAxis.end(); it++){
       if(it->numAxis == 0 && abs(it->p - x) >= it->distance ) flag = false;
       if(it->numAxis == 1 && abs(it->p - y) >= it->distance ) flag = false;
       if(it->numAxis == 2 && abs(it->p - z) >= it->distance ) flag = false;
     }
     if( vecAxis.size() == 0 || !flag || 
         this->mapTimeCurrent.find(time) == this->mapTimeCurrent.end())
     {
       return 0;
     }else{
       return this->mapTimeCurrent[time];
     }
   };
};

}//namespace
#endif // __STIMELECTRODE__
