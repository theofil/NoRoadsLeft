#include <string>
#include <iomanip>
#include <sstream>

#ifndef any2string_h
#define any2string_h
template<class T> std::string any2string(T i){std::ostringstream buffer;    buffer << i;    return buffer.str();}
template<class T> std::string any2string2(T i){std::ostringstream buffer;    buffer << setprecision(3) << i;    return buffer.str();}
#endif
