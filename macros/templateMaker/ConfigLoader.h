#ifndef CONFIGLOADER_H
#define CONFIGLOADER_H

#include "Config.h"
#include <string>

HistogramConfig LoadHistogramConfig(const std::string& filename);

#endif
