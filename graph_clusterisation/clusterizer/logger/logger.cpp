/* #include "logger.h"

namespace logging = boost::log;
namespace src = boost::log::sources;
namespace sinks = boost::log::sinks;
namespace keywords = boost::log::keywords;


void init(std::string& logFile, int level)
{
    logging::add_file_log(
            keywords::file_name = logFile,
            keywords::format = "%Message%");

    logging::core::get()->set_filter(logging::trivial::severity >= boost::log::trivial::severity_level(level));
}*/