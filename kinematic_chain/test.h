#pragma once
#include <geode/python/Object.h>

using namespace geode;

namespace other{
  class Whatever : public Object{
  public:
    GEODE_DECLARE_TYPE(GEODE_NO_EXPORT)
      typedef Object Base;
  protected:
    Whatever();
  public:
    int mything;
  };
}
