// To be included by Paths.cc

#include "Utils.h"

#include <TString.h>

#include <cstdarg>
#include <cstdlib>
#include <libgen.h>             // dirname
#include <stdexcept>
#include <string>

// Returns true if $envvar is set to anything except "0".
// Returns false if $envvar is unset or it is set to "0".
static bool checkenv(const char* envvar)
{
  const char* val = getenv(envvar);
  return val && (strcmp(val, "0") != 0);
}

static void ensure_dir(const char* path)
{
  std::string s(path);
  system(Form("mkdir -p %s", dirname(s.data())));
}

// We assume that directories (e.g. LBNL_FIT_OUTDIR) are specified relative to
// the fitter root, but we are running from ShapeFit or toySpectra, so prepend
// ../ if the path is a relative one
static const char* normalize(const char* path)
{
  return (path[0] == '/') ? path
    : LeakStr("../%s", path);
}

static const char* normalized_or(const char* envvar, const char* dflt)
{
  const char* d = getenv(envvar);
  if (!d) d = dflt;
  return normalize(d);
}

static const char* indir()
{
  return normalized_or("LBNL_FIT_INDIR", "input");
}

static const char* outdir()
{
  return normalized_or("LBNL_FIT_OUTDIR", indir());
}

namespace Paths {

__attribute__((format(printf, 1, 2)))
const char* inpath(const char* fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  const char* relpath = _Form(fmt, ap);
  va_end(ap);

  return LeakStr("%s/%s", indir(), relpath);
}

__attribute__((format(printf, 1, 2)))
const char* outpath(const char* fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  const char* relpath = _Form(fmt, ap);
  va_end(ap);

  const char* path = LeakStr("%s/%s", outdir(), relpath);
  ensure_dir(path);
  return path;
}

} // namespace Paths

static int stage_nADs(int istage)
{
  switch (istage) {
    case 0: return 6;
    case 1: return 8;
    case 2: return 7;
  }

  throw std::runtime_error("stage_nADs: Unknown stage");
}

const char* stage_lwc(int istage)
{
  return LeakStr("%dad", stage_nADs(istage));
}

const char* stage_upc(int istage)
{
  return LeakStr("%dAD", stage_nADs(istage));
}
