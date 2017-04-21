#ifndef MACROS_HPP
#define MACROS_HPP

#include <stdio.h>
#include <stdexcept>

#define SUB_LOCATION(msg) msg=__FILE__;\
                          msg+=":";\
			  char line[10];\
			  sprintf(line, "%d", __LINE__);\
                          msg+=line;\
                          msg+=": error : in ";\
                          msg+=__FUNCTION__;\
			  msg+="\n";

#define THROW_ERROR(msg) std::string loc;\
  SUB_LOCATION(loc);			 \
  loc = loc + msg + "\n";\
  throw std::runtime_error(loc);

/*
#define SUB_LOCATION(msg) msg=__FILE__;\
                          msg+=":";\
                          msg+=__FUNCTION__;\
                          msg+=":";\
			  char line[10];\
			  sprintf(line, "%d", __LINE__);\
                          msg+=line;
*/
/*
#define SUB_LOCATION(msg) msg="\n";\
                          msg+="  File \""+__FILE__+"\", line ";\
			  char line[10];\
			  sprintf(line, "%d", __LINE__);\
                          msg+=line;\
                          msg+=", in ";\
                          msg+=__FUNCTION__;\
                          msg+="\n";
*/

#endif


