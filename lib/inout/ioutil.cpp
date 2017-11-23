#include "ioutil.h"

#include <boost/property_tree/json_parser.hpp>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <sstream>
#include <fparameters/parameters.h>
using boost::property_tree::ptree;
namespace fs = boost::filesystem;
using boost::format;
using std::string;

const string& FILE_CONFIG{ "parameters.json" };
const string& FILE_OUTPUT_DIR_PFX{ "output_" };
const string& FS{ "\\" };

std::string time2str(boost::posix_time::ptime now)
{
	using namespace boost::posix_time;
	static std::locale loc(std::cout.getloc(),
		new time_facet("%Y%m%d_%H%M%S"));
	std::ostringstream wss;
	wss.imbue(loc);
	wss << now;
	return wss.str();
}

std::string prepareOutputfolder() {

	using namespace boost::posix_time;
	ptime now = second_clock::universal_time();
	string folder{ FILE_OUTPUT_DIR_PFX + time2str(now) };
	fs::create_directory(folder);
//	fs::create_directory(folder + "/plots");
	fs::copy_file(FILE_CONFIG, folder + FS + FILE_CONFIG);
	return std::move(folder);
}

ptree readConfig(std::string filename) {
	ptree cfg;
	boost::property_tree::read_json(filename, cfg);
	return std::move(cfg);
}

std::string getFileName(std::string outputFolder, std::string fileId) {
	return outputFolder + FS + GlobalConfig.get<string>("file." + fileId);
}