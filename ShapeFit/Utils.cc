#include "Utils.h"

#include <TString.h>

const char* _Form(const char* fmt, va_list ap)
{
  char str[2048];
  vsnprintf(str, 2048, fmt, ap);
  return Form("%s", str);      // use ROOT's circular buffer
}

__attribute__((format(printf, 1, 2)))
const char* LeakStr(const char* fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  auto s = new TString(_Form(fmt, ap)); // intentional leak
  va_end(ap);
  return s->Data();
}

