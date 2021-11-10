#pragma once

#include <cstdarg>

const char* _Form(const char* fmt, va_list ap);

// "leaky" Form; ensures result won't get clobbered
const char* LeakStr(const char* fmt, ...);
