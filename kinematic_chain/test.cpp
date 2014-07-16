#include "test.h"
#include <geode/python/Class.h>
    
namespace other{

  GEODE_DEFINE_TYPE(Whatever)

  Whatever::Whatever()
    : mything(10)
  {
  }
}

using namespace other;
using namespace geode;

void wrap_test(){
  typedef Whatever Self;
  Class<Self>("Whatever")
    .GEODE_INIT()
    .GEODE_FIELD(mything);
}
