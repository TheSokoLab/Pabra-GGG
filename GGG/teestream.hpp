#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;

typedef tee_device<ofstream, ofstream> TeeDevice;
typedef stream<TeeDevice> TeeStream;

/* EXAMPLE */
/*
class bio_teed_stream : bio::filtering_ostream
{

  bio::filtering_ostream teeOut;
  teeOut.push( io::tee(io::file_sink("out1.txt")) );
  teeOut.push( io::tee(io::file_sink("out2.txt")) );
}
*/
