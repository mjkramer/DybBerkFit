// To be included by Paths.cc

#include "TString.h"

#include <cstdarg>
#include <cstdlib>
#include <libgen.h>             // dirname
#include <stdexcept>
#include <string>

static void ensure_dir(const char* path)
{
  std::string s(path);
  system(Form("mkdir -p %s", dirname(s.data())));
}

static const char* _Form(const char* fmt, va_list ap)
{
  char str[2048];
  vsnprintf(str, 2048, fmt, ap);
  return Form("%s", str);      // use ROOT's circular buffer
}

// We assume that directories (e.g. LBNL_FIT_OUTDIR) are specified relative to
// the fitter root, but we are running from ShapeFit or toySpectra, so prepend
// ../ if the path is a relative one
static const char* normalize(const char* path)
{
  return (path[0] == '/') ? path
    : Form("../%s", path);
}

static const char* normalized_or(const char* envvar, const char* dflt)
{
  const char* d = getenv(envvar);
  if (!d) d = dflt;
  return normalize(d);
}

static const char* outdir()
{
  return normalized_or("LBNL_FIT_OUTDIR", "output");
}

static const char* indir()
{
  return normalized_or("LBNL_FIT_INDIR", "input");
}

__attribute__((format(printf, 1, 2)))
static const char* outpath(const char* fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  const char* relpath = _Form(fmt, ap);
  va_end(ap);

  const char* path = Form("%s/%s", outdir(), relpath);
  ensure_dir(path);
  return path;
}

__attribute__((format(printf, 1, 2)))
static const char* inpath(const char* fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  const char* relpath = _Form(fmt, ap);
  va_end(ap);

  const char* path = Form("%s/%s", indir(), relpath);
  return path;
}

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
  return Form("%dad", stage_nADs(istage));
}

const char* stage_upc(int istage)
{
  return Form("%dAD", stage_nADs(istage));
}
