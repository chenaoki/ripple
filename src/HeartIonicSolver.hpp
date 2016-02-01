#ifndef __HEARTIONICSOLVER__
#define __HEARTIONICSOLVER__

#include "HeartIonicSolverBase.hpp"
#include "HeartIonicSolverOHR.hpp"

namespace LifeV{

enum EnumHeartIonicSolverType{
      RM = 1, // Rogers & McCulloch
      LR, // Luo & Rudy
      MS, // Mitchell & Schaeffer
      OR // O'Hara & Rudy
};

}

#endif// __HEARTIONICSOLVER__
