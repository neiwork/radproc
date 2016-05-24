#pragma once

#include <boost/property_tree/ptree_fwd.hpp>

/** creates the output folder and returns its location */
std::string prepareOutputfolder();

boost::property_tree::ptree readConfig();

std::string getFileName(boost::property_tree::ptree cfg, std::string outputFolder, std::string fileId);