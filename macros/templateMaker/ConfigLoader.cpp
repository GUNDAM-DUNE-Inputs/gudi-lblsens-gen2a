#include "ConfigLoader.h"
#include <yaml-cpp/yaml.h>
#include <iostream>

HistogramConfig LoadHistogramConfig(const std::string& filename)
{
    YAML::Node config = YAML::LoadFile(filename);

    if(!config["systematicsInfo"])
    {
        std::cerr << "Missing 'systematicsInfo' section in YAML.\n";
        exit(1);
    }

    YAML::Node systNode = config["systematicsInfo"];

    HistogramConfig cfg;

    cfg.dialNames = systNode["Dials"].as<std::vector<std::string>>();

    cfg.sampleNames = systNode["Samples"].as<std::vector<std::string>>();

    cfg.variationValues = systNode["Variations"].as<std::vector<double>>();

    cfg.binEdges = systNode["binEdges"].as<std::vector<double>>();

    if(cfg.binEdges.size() < 2)
    {
        std::cerr << "Bin edges must contain at least 2 values.\n";
        exit(1);
    }

    return cfg;
}
